//  MatrixSpotDistortionTests.swift
//  DiffractoKit#96: the detector's distortion fitted on the matrix's spots.
//
//  Planted: aluminium on [001], turned 12 degrees in the plane, its spots within 0.75 1/Å carried by a known stretch (and,
//  in one case, a wrong scale), on a small scan. The orientation the fit is given is the true one or one 1.5 degrees off in
//  the plane (the template match gives the angle to half an azimuthal bin). The numbers of the 5 % case are those measured
//  on a public SPED scan (Zenodo 10.5281/zenodo.2652906): 1.051 at 166 degrees.

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct MatrixSpotDistortionTests {
    static let q = 0.012
    static let origin: Float = 63.5
    static let reach = 0.75

    /// The symmetric stretch: `ratio` along the direction at `axisDeg`, 1 across it, times `scale`.
    static func stretch(ratio: Double, axisDeg: Double, scale: Double = 1) -> simd_double2x2 {
        let a = axisDeg * .pi / 180, c = cos(a), s = sin(a)
        let u = simd_double2x2(columns: (SIMD2(c, s), SIMD2(-s, c)))
        return u * simd_double2x2(diagonal: SIMD2(scale * ratio, scale)) * u.transpose
    }

    static func fit(distortion: simd_double2x2, offDeg: Double = 0, keep: (SIMD2<Double>) -> Bool = { _ in true },
                    positions: Int = 16) throws -> MatrixSpotDistortion? {
        let al = try PhaseOrientationCandidate(name: "matrix", crystal: CrystalPresets.aluminum, kMax: reach)
        func orientation(_ deg: Double) -> simd_double3x3 {
            ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(deg * Double.pi / 180))
        }
        let truth = PhaseVectorMatcher.matrixSpots(reflections: al.reflections, orientation: orientation(12), excitationSlabInvAngstrom: 0.02)
            .filter { simd_length($0) > 0.15 && simd_length($0) < reach && keep($0) }
        let peaks = ([SIMD2<Double>(0, 0)] + truth.map { distortion * $0 }).map { SIMD2(Double(origin) + $0.x / q, Double(origin) + $0.y / q) }
        let scan = [[SIMD2<Double>]](repeating: peaks, count: positions)
        let field = MatrixOrientationField(width: positions, height: 1, results: [OrientationResult](repeating: .empty, count: positions),
                                           orientations: [simd_double3x3?](repeating: orientation(12 + offDeg), count: positions))
        return PhaseVectorMatcher.matrixSpotDistortion(scan: scan, field: field, reflections: al.reflections,
                                                       originX: origin, originY: origin, invAngstromPerPixel: q,
                                                       maximumVectorInvAngstrom: 0.85)
    }

    /// A stretch of 5 % at 166 degrees comes back, from the true orientation and from one 1.5 degrees off in the plane.
    @Test(arguments: [0.0, 1.5, -1.5]) func aPlantedStretchComesBack(_ offDeg: Double) throws {
        let d = try #require(try Self.fit(distortion: Self.stretch(ratio: 1.0506, axisDeg: 166.3), offDeg: offDeg))
        print("MATRIX SPOT DISTORTION | off \(offDeg) deg | ratio \(d.stretchRatio), axis \(d.stretchAxisRad * 180 / .pi) deg, mean scale \(d.meanScale), pairs \(d.pairCount), residual \(d.residualInvAngstrom), without stretch \(d.residualWithoutStretchInvAngstrom)")
        #expect(abs(d.stretchRatio - 1.0506) < 0.0005)
        #expect(abs(d.stretchAxisRad * 180 / .pi - 166.3) < 0.5)
        #expect(abs(d.meanScale - 1.0506.squareRoot()) < 0.0005)
        #expect(d.residualInvAngstrom < 1e-6)
        #expect(d.residualWithoutStretchInvAngstrom > 0.005, "the stretch must show in the residual of rotation and scale alone")
        #expect(d.pairCount >= 16 * 8)
        // Taking the stretch out lays the measured spots back on a square lattice of the mean scale.
        let back = d.removingStretch * Self.stretch(ratio: 1.0506, axisDeg: 166.3)
        let s = MatrixSpotDistortionProbe.singularValues(back)
        #expect(abs(s.0 / s.1 - 1) < 0.0005)
        #expect(abs((s.0 * s.1).squareRoot() - 1.0506.squareRoot()) < 0.0005)
    }

    /// No distortion: ratio 1 and scale 1. A wrong scale alone: ratio 1 and that scale.
    @Test func anUndistortedPatternAndAWrongScale() throws {
        let none = try #require(try Self.fit(distortion: matrix_identity_double2x2, offDeg: 1.5))
        #expect(abs(none.stretchRatio - 1) < 1e-6)
        #expect(abs(none.meanScale - 1) < 1e-6)
        #expect(none.residualInvAngstrom < 1e-6 && none.residualWithoutStretchInvAngstrom < 1e-6)
        let scaled = try #require(try Self.fit(distortion: Self.stretch(ratio: 1, axisDeg: 0, scale: 1.03)))
        #expect(abs(scaled.stretchRatio - 1) < 1e-6)
        #expect(abs(scaled.meanScale - 1.03) < 1e-6)
    }

    /// Refusals: spots along one line fix no stretch; too few pairs; lengths that disagree.
    @Test func whatCannotBeFittedIsRefused() throws {
        // Only the spots along one row of the pattern (and its opposite): one direction.
        let al = try PhaseOrientationCandidate(name: "matrix", crystal: CrystalPresets.aluminum, kMax: Self.reach)
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(12 * Double.pi / 180))
        let all = PhaseVectorMatcher.matrixSpots(reflections: al.reflections, orientation: m, excitationSlabInvAngstrom: 0.02)
            .filter { simd_length($0) > 0.15 && simd_length($0) < Self.reach }
        let first = try #require(all.min { simd_length($0) < simd_length($1) })
        let direction = simd_normalize(first)
        #expect(try Self.fit(distortion: Self.stretch(ratio: 1.05, axisDeg: 30),
                             keep: { abs($0.x * direction.y - $0.y * direction.x) < 1e-6 }) == nil)
        #expect(try Self.fit(distortion: Self.stretch(ratio: 1.05, axisDeg: 30), positions: 2) == nil, "16 pairs are fewer than the minimum of 50")
        let field = MatrixOrientationField(width: 1, height: 1, results: [.empty], orientations: [m])
        #expect(PhaseVectorMatcher.matrixSpotDistortion(scan: [], field: field, reflections: al.reflections,
                                                        originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q) == nil)
    }
}

/// The singular values of a 2 × 2 matrix for the test's own check (larger first), by the closed form of M Mᵀ's eigenvalues.
enum MatrixSpotDistortionProbe {
    static func singularValues(_ m: simd_double2x2) -> (Double, Double) {
        let a = m * m.transpose
        let half = (a[0][0] + a[1][1]) / 2, diff = (a[0][0] - a[1][1]) / 2
        let root = (diff * diff + a[1][0] * a[1][0]).squareRoot()
        return ((half + root).squareRoot(), max(half - root, 0).squareRoot())
    }
}
