//
//  MatrixToleranceInSecondPassTests.swift
//  DiffractoKit#106: the second pass uses the matrix tolerance the removal was run with. A caller who widens the tolerance
//  (a coarse detector, a scale a little off) gets one value in the whole run: the removal, the own-pattern helper and the
//  classification against the variants.
//  Planted truth: aluminium on [001] with beta'' at Andersen's relationship, as NeighbourGrainSpotsTests; here the matrix's
//  measured spots sit 0.03 1/Å from their predicted places and the removal (run by the caller at 0.04) took them.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct MatrixToleranceInSecondPassTests {
    static let q = 0.012, reach = 0.75
    static let origin: Float = 63.5

    /// With the tolerance of the removal (0.04) handed on, a position whose matrix spots are 0.03 1/Å off is the matrix and a
    /// rich precipitate position beside such spots is the phase. The default (0.02) gives the verdicts of a call that names
    /// 0.02, so nothing changes for a caller who does not pass it.
    /// Mutation: the tolerance handed to the own-pattern helper and not to the classification (the displaced matrix spots
    /// would be kept in the pattern and then counted as the phase's own vectors).
    @Test func theRemovalsToleranceIsTheSecondPasssTolerance() throws {
        let al = CrystalPresets.aluminum, beta = CrystalPresets.betaDoublePrime
        let relation = try OrientationRelationship.rotation(
            OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
            OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(3, 1, 0))),
            precipitate: beta, matrix: al)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(al.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(beta.laueClass))
        let variants = LearnedOrientationRelationship(relation: relation, positionCount: 0, spreadRad: 0, matrixSymmetry: ms, precipitateSymmetry: ps).variants
        let candidate = try PhaseOrientationCandidate(name: "beta''", crystal: beta, kMax: Self.reach)
        let matrixCandidate = try PhaseOrientationCandidate(name: "Al", crystal: al, kMax: Self.reach)
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(12 * Double.pi / 180))
        func inReach(_ spots: [SIMD2<Double>]) -> [SIMD2<Double>] {
            spots.filter { simd_length($0) > 0.15 && simd_length($0) < Self.reach }.sorted { simd_length($0) < simd_length($1) }
        }
        let predicted = inReach(PhaseVectorMatcher.matrixSpots(reflections: matrixCandidate.reflections, orientation: m, excitationSlabInvAngstrom: 0.02))
        // The matrix spots as measured: each 0.03 1/Å farther out than predicted.
        let displaced = predicted.map { $0 + 0.03 * simd_normalize($0) }
        func betaSpots(_ variant: simd_double3x3) -> [SIMD2<Double>] {
            inReach(PhaseVectorMatcher.matrixSpots(reflections: candidate.reflections, orientation: variant.transpose * m, excitationSlabInvAngstrom: 0.02))
        }
        // A precipitate spot within the removal's tolerance of a measured matrix spot is taken by it.
        func own(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] {
            list.filter { p in !predicted.contains { simd_length($0 - p) <= 0.04 } && !displaced.contains { simd_length($0 - p) <= 0.04 } }
        }
        func pixels(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] { list.map { SIMD2(Double(Self.origin) + $0.x / Self.q, Double(Self.origin) + $0.y / Self.q) } }
        let byRichness = variants.sorted { own(betaSpots($0)).count > own(betaSpots($1)).count }
        let most = own(betaSpots(byRichness[0])).count
        let rich = byRichness.filter { own(betaSpots($0)).count == most }

        var scan: [[SIMD2<Double>]] = [], residual: [[SIMD2<Double>]] = [], removed: [Int] = [], kind: [String] = []
        for k in 0..<24 {
            let spots = own(betaSpots(rich[k % rich.count]))
            scan.append(pixels(displaced + spots)); residual.append(pixels(spots)); removed.append(displaced.count); kind.append("rich")
        }
        for _ in 0..<12 { scan.append(pixels(displaced)); residual.append([]); removed.append(displaced.count); kind.append("matrix") }

        func run(_ tolerance: Double?) throws -> TwoPassClassification {
            let residualScan = MatrixResidualScan(scan: residual, removedCounts: removed, inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count))
            let orientations = [simd_double3x3?](repeating: m, count: scan.count)
            let phases = [TwoPassPhase(candidate: candidate, symmetry: ps)]
            let settings = PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach)
            if let tolerance {
                return try #require(PhaseVectorMatcher.classifyInTwoPasses(
                    residual: residualScan, matrixOrientations: orientations, phases: phases, matrixSymmetry: ms, settings: settings,
                    originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                    matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scan.count, minimumPositions: 10,
                    matrixToleranceInvAngstrom: tolerance))
            }
            return try #require(PhaseVectorMatcher.classifyInTwoPasses(
                residual: residualScan, matrixOrientations: orientations, phases: phases, matrixSymmetry: ms, settings: settings,
                originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scan.count, minimumPositions: 10))
        }
        let wide = try run(0.04)
        #expect(!wide.relationships[0].isEmpty, "fixture: the rich positions must teach the relationship")
        for (i, k) in kind.enumerated() {
            let v = wide.verdicts[i]
            if k == "rich" {
                #expect(v.verdict == .indexed, "rich \(i): \(v.verdict), matched \(v.matchedCount) of \(v.survivingCount)")
            } else {
                #expect(v.verdict == .matrix, "matrix \(i): \(v.verdict), matched \(v.matchedCount) of \(v.survivingCount)")
            }
        }
        let byDefault = try run(nil), named = try run(0.02)
        #expect(byDefault.verdicts == named.verdicts, "the default must be 0.02")
    }
}
