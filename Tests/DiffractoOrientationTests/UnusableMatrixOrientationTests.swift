//  UnusableMatrixOrientationTests.swift
//  DiffractoKit#94: a position where the matrix removal could lay the matrix orientation on none of the position's vectors
//  has no usable matrix orientation. Its residual is the matrix pattern itself, and a precipitate indexed on it is the
//  matrix's spots read as a phase.
//
//  Seen on a public SPED scan of an Al-Mg-Si alloy (Zenodo 10.5281/zenodo.2652906): the peak list handed over carried a
//  per-position origin error of about one pixel (0.026 1/Å) in two corners of the scan. The removal at 0.02 1/Å took no
//  aluminium spot there, all eight stayed in the residual, and about 25,000 matrix positions were called beta''.
//
//  Planted here: Al on [001] with beta'' at Andersen's relationship (the pair of `TwoPassOtherMaterialsTests`), enough
//  rich positions to learn the relationship, and matrix-only positions whose every peak is shifted by 0.03 1/Å, with
//  nothing removed, as the removal leaves them.

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct UnusableMatrixOrientationTests {
    static let q = 0.012
    static let origin: Float = 63.5
    static let tolerance = 0.02
    static let reach = 0.75

    struct Outcome {
        var learned = 0
        var rich = 0, richCalled = 0
        var shifted = 0, shiftedCalled = 0
        var matrix = 0, matrixIsMatrix = 0
    }

    static func run(shiftInvAngstrom: SIMD2<Double>) throws -> Outcome {
        let al = CrystalPresets.aluminum, beta = CrystalPresets.betaDoublePrime
        let relation = try OrientationRelationship.rotation(
            OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
            OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(3, 1, 0))),
            precipitate: beta, matrix: al)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(al.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(beta.laueClass))
        let planted = LearnedOrientationRelationship(relation: relation, positionCount: 0, spreadRad: 0, matrixSymmetry: ms, precipitateSymmetry: ps)
        let candidate = try PhaseOrientationCandidate(name: "precipitate", crystal: beta, kMax: reach)
        let matrixCandidate = try PhaseOrientationCandidate(name: "matrix", crystal: al, kMax: reach)
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(12 * Double.pi / 180))
        func inReach(_ spots: [SIMD2<Double>]) -> [SIMD2<Double>] {
            spots.filter { simd_length($0) > 0.15 && simd_length($0) < reach }.sorted { simd_length($0) < simd_length($1) }
        }
        let matrixSpots = inReach(PhaseVectorMatcher.matrixSpots(reflections: matrixCandidate.reflections, orientation: m, excitationSlabInvAngstrom: 0.02))
        func spots(_ variant: simd_double3x3) -> [SIMD2<Double>] {
            inReach(PhaseVectorMatcher.matrixSpots(reflections: candidate.reflections, orientation: variant.transpose * m, excitationSlabInvAngstrom: 0.02))
        }
        func own(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] {
            list.filter { p in !matrixSpots.contains { simd_length($0 - p) <= tolerance } }
        }
        func pixels(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] { list.map { SIMD2(Double(origin) + $0.x / q, Double(origin) + $0.y / q) } }
        let byRichness = planted.variants.sorted { own(spots($0)).count > own(spots($1)).count }
        let most = own(spots(byRichness[0])).count
        let rich = byRichness.filter { own(spots($0)).count == most }

        var scan: [[SIMD2<Double>]] = [], residual: [[SIMD2<Double>]] = [], removed: [Int] = [], kind: [String] = []
        for k in 0..<24 {
            let p = spots(rich[k % rich.count])
            scan.append(pixels(matrixSpots + p)); residual.append(pixels(own(p))); removed.append(matrixSpots.count); kind.append("rich")
        }
        // The matrix alone, every peak shifted: the removal takes none of them.
        let shifted = pixels(matrixSpots.map { $0 + shiftInvAngstrom })
        for _ in 0..<12 { scan.append(shifted); residual.append(shifted); removed.append(0); kind.append("shifted") }
        for _ in 0..<12 { scan.append(pixels(matrixSpots)); residual.append([]); removed.append(matrixSpots.count); kind.append("matrix") }

        let result = try #require(PhaseVectorMatcher.classifyInTwoPasses(
            residual: MatrixResidualScan(scan: residual, removedCounts: removed, inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count)),
            matrixOrientations: [simd_double3x3?](repeating: m, count: scan.count),
            phases: [TwoPassPhase(candidate: candidate, symmetry: ps)], matrixSymmetry: ms,
            settings: PhaseOrientationSearchSettings(reachInvAngstrom: reach),
            originX: origin, originY: origin, invAngstromPerPixel: q,
            matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scan.count, minimumPositions: 10))
        var out = Outcome()
        out.learned = result.relationships[0].count
        for (i, k) in kind.enumerated() {
            let v = result.verdicts[i]
            switch k {
            case "rich": out.rich += 1; if v.verdict == .indexed { out.richCalled += 1 }
            case "shifted": out.shifted += 1; if v.verdict == .indexed { out.shiftedCalled += 1 }
            default: out.matrix += 1; if v.verdict == .matrix { out.matrixIsMatrix += 1 }
            }
        }
        return out
    }

    /// The matrix pattern with every peak 0.03 1/Å off, along x, along y and along the diagonal: never a precipitate, while
    /// the rich positions are still the phase and the true matrix positions still matrix.
    @Test(arguments: [SIMD2(0.03, 0.0), SIMD2(0.0, 0.03), SIMD2(0.0212, -0.0212), SIMD2(-0.015, 0.026)])
    func aShiftedMatrixPatternIsNeverAPrecipitate(_ shift: SIMD2<Double>) throws {
        let o = try Self.run(shiftInvAngstrom: shift)
        print("UNUSABLE MATRIX | shift \(shift) | \(o)")
        #expect(o.learned >= 1, "\(o)")
        #expect(o.richCalled == o.rich, "\(o)")
        #expect(o.matrixIsMatrix == o.matrix, "\(o)")
        #expect(o.shiftedCalled == 0, "the matrix's own spots, shifted, were called a precipitate. \(o)")
    }
}
