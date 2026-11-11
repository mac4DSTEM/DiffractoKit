//
//  TwoPassClassificationTests.swift
//  Step 4 of "find phase orientations from the structure alone" (DiffractoKit#81): the two passes as one call. The free search
//  teaches the relationships; they are refined on the strong positions and every position is classified against the variants.
//  The free-search verdict is kept beside the final one, with the angle of its orientation to the learned relationship, so a
//  caller can say "this is the phase, but not at the relationship it takes elsewhere". Here on the residual alone; the whole
//  pattern, the firm rule and the neighbour rule are in TwoPassOtherMaterialsTests.
//  Planted truth: beta'' in aluminium at Andersen's relationship, as in OrientationRelationshipLearningTests.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct TwoPassClassificationTests {
    typealias Learning = OrientationRelationshipLearningTests

    /// The matrix orientation with the removal's in-plane correction applied predicts the spots the removal fitted: a field
    /// 3 degrees off and a correction of -3 degrees give the planted matrix orientation back.
    /// Mutation: the correction applied with the other sign.
    @Test func theCorrectedMatrixOrientationIsWhereTheRemovalPutIt() throws {
        let s = try Learning.setup()
        let off = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(15 * Double.pi / 180))
        let field = MatrixOrientationField(width: 2, height: 1, results: [.empty, .empty], orientations: [off, nil])
        let corrected = field.correctedOrientations(inPlaneCorrectionsRad: [-3 * .pi / 180, 0])
        #expect(corrected.count == 2 && corrected[1] == nil)
        let m = try #require(corrected[0])
        for c in 0..<3 { for r in 0..<3 { #expect(abs(m[c][r] - s.matrixOrientation[c][r]) < 1e-6) } }
        #expect(field.correctedOrientations(inPlaneCorrectionsRad: [0]).isEmpty, "a correction list of another length is refused")
    }

    /// The whole scan in one call. Twenty-four rich positions of the rich variants, twelve weak positions (three spots each) of
    /// other variants, six rich positions of the same crystal turned 45 degrees (a stand-in for a precipitate that does not
    /// keep the relationship), and plain matrix positions.
    /// - The relationship comes back within 0.5 degrees of the planted one, with twelve variants.
    /// - Every weak position is beta''; the free search alone had refused it.
    /// - The turned positions: the free search finds the crystal, and its orientation is reported far (more than 20 degrees)
    ///   from the learned relationship; at the rich and weak positions that angle is small (under 5 degrees) or absent.
    /// - A plain matrix position is matrix.
    @Test func theScanIsClassifiedInTwoPassesAndSaysWhereAPhaseLeavesItsRelationship() throws {
        let s = try Learning.setup()
        let rich = Learning.richVariants(s)
        func same(_ a: simd_double3x3, _ b: simd_double3x3) -> Bool {
            (0..<3).allSatisfy { c in (0..<3).allSatisfy { r in abs(a[c][r] - b[c][r]) < 1e-9 } }
        }
        let others = s.variants.filter { v in !same(v, rich[0]) && Learning.peaks(s, variant: v).count >= 3 }
        #expect(!others.isEmpty)
        var scan: [[SIMD2<Double>]] = []
        var kind: [String] = []
        for k in 0..<24 { scan.append(Learning.peaks(s, variant: rich[k % rich.count])); kind.append("rich") }
        for k in 0..<12 { scan.append(Learning.peaks(s, variant: others[k % others.count], count: 3)); kind.append("weak") }
        for _ in 0..<6 { scan.append(Learning.peaks(s, variant: rich[0], turnedDeg: 45)); kind.append("turned") }
        for _ in 0..<6 { scan.append([]); kind.append("matrix") }
        let residual = MatrixResidualScan(scan: scan, removedCounts: [Int](repeating: 8, count: scan.count),
                                          inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count))
        let matrix = [simd_double3x3?](repeating: s.matrixOrientation, count: scan.count)
        let result = try #require(PhaseVectorMatcher.classifyInTwoPasses(
            residual: residual, matrixOrientations: matrix,
            phases: [TwoPassPhase(candidate: s.candidate, symmetry: s.precipitateSymmetry)],
            matrixSymmetry: s.matrixSymmetry, settings: PhaseOrientationSearchSettings(reachInvAngstrom: Learning.reach),
            originX: Learning.origin, originY: Learning.origin, invAngstromPerPixel: Learning.q,
            minimumPositions: 10))
        #expect(result.verdicts.count == scan.count && result.freeSearch.count == scan.count && result.relationAnglesRad.count == scan.count)
        #expect(result.neighbourSupported.count == scan.count && result.relationships.count == 1)
        let learned = try #require(result.relationships.first?.first)
        let distance = PhaseVectorMatcher.relationDistance(learned.relation, s.relation, matrixSymmetry: s.matrixSymmetry,
                                                           precipitateSymmetry: s.precipitateSymmetry) * 180 / .pi
        #expect(distance < 0.5, "the relationship is \(distance) deg from the planted one")
        #expect(learned.variants.count == 12)
        for (i, k) in kind.enumerated() {
            let v = result.verdicts[i], free = result.freeSearch[i], angle = result.relationAnglesRad[i] * 180 / .pi
            switch k {
            case "rich":
                #expect(v.verdict == .indexed && v.phaseIndex == 0, "rich \(i): \(v.verdict)")
                #expect(angle < 5, "rich \(i): free-search orientation \(angle) deg from the relationship")
            case "weak":
                #expect(free.verdict == .notIndexed, "weak \(i): the free search must refuse three spots")
                #expect(v.verdict == .indexed && v.phaseIndex == 0, "weak \(i): \(v.verdict)")
                #expect(angle.isNaN, "weak \(i): no free-search orientation, so no angle")
            case "turned":
                #expect(free.verdict == .indexed && free.phaseIndex == 0, "turned \(i): the free search finds the crystal")
                #expect(angle > 20, "turned \(i): \(angle) deg from the relationship")
            default:
                #expect(v.verdict == .matrix)
                #expect(angle.isNaN)
            }
        }
        // Lengths that do not agree are refused.
        #expect(PhaseVectorMatcher.classifyInTwoPasses(
            residual: residual, matrixOrientations: Array(matrix.dropLast()),
            phases: [TwoPassPhase(candidate: s.candidate, symmetry: s.precipitateSymmetry)],
            matrixSymmetry: s.matrixSymmetry, settings: PhaseOrientationSearchSettings(reachInvAngstrom: Learning.reach),
            originX: Learning.origin, originY: Learning.origin, invAngstromPerPixel: Learning.q) == nil)
    }
}
