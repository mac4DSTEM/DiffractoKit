//
//  NeighbourGrainSpotsTests.swift
//  DiffractoKit#103: at a grain boundary a pattern holds the spots of two matrix grains. The removal takes both (the
//  neighbouring-grain rule), and the free search says "matrix". The second pass scores the whole pattern, and must not
//  offer the neighbouring grain's spots to a phase as vectors of its own: on the simulated three-grain demo a row of 29
//  boundary positions was called beta'' on two of the neighbour's six spots.
//  Planted truth: aluminium on [001] with beta'' at Andersen's relationship (as TwoPassOtherMaterialsTests), and
//  positions that hold the [001] grain's spots together with those of a second aluminium grain on another zone.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct NeighbourGrainSpotsTests {
    static let q = 0.012, reach = 0.75, tolerance = 0.02
    static let origin: Float = 63.5

    /// Rich beta'' positions teach the relationship; then, for a second grain on [011] and on [111] at twelve in-plane
    /// angles each, a position holding both grains' spots whose removal took all of them. Every such position is the
    /// matrix: no vector of it is a phase's.
    /// Mutation: the second pass given every peak of the position (the neighbour's spots counted as the phase's own).
    @Test func aSecondGrainsSpotsAreNotAPhasesOwnVectors() throws {
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
        func alSpots(_ orientation: simd_double3x3) -> [SIMD2<Double>] {
            inReach(PhaseVectorMatcher.matrixSpots(reflections: matrixCandidate.reflections, orientation: orientation, excitationSlabInvAngstrom: 0.02))
        }
        let matrixSpots = alSpots(m)
        func betaSpots(_ variant: simd_double3x3) -> [SIMD2<Double>] {
            inReach(PhaseVectorMatcher.matrixSpots(reflections: candidate.reflections, orientation: variant.transpose * m, excitationSlabInvAngstrom: 0.02))
        }
        func own(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] { list.filter { p in !matrixSpots.contains { simd_length($0 - p) <= Self.tolerance } } }
        func pixels(_ list: [SIMD2<Double>]) -> [SIMD2<Double>] { list.map { SIMD2(Double(Self.origin) + $0.x / Self.q, Double(Self.origin) + $0.y / Self.q) } }
        let byRichness = variants.sorted { own(betaSpots($0)).count > own(betaSpots($1)).count }
        let most = own(betaSpots(byRichness[0])).count
        let rich = byRichness.filter { own(betaSpots($0)).count == most }

        var scan: [[SIMD2<Double>]] = [], residual: [[SIMD2<Double>]] = [], removed: [Int] = [], kind: [String] = []
        for k in 0..<24 {
            let spots = betaSpots(rich[k % rich.count])
            scan.append(pixels(matrixSpots + spots)); residual.append(pixels(own(spots))); removed.append(matrixSpots.count); kind.append("rich")
        }
        // The boundary positions: the removal took the spots of both grains, nothing is left.
        var neighbourSpotCounts: [Int] = []
        for zone in [SIMD3<Double>(0, 1, 1), SIMD3<Double>(1, 1, 1)] {
            for step in 0..<12 {
                let second = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: zone), inPlaneAngle: Float(Double(step) * 15 * Double.pi / 180))
                let theirs = own(alSpots(second))
                neighbourSpotCounts.append(theirs.count)
                scan.append(pixels(matrixSpots + theirs)); residual.append([]); removed.append(matrixSpots.count + theirs.count)
                kind.append("boundary \(Int(zone.x))\(Int(zone.y))\(Int(zone.z)) at \(step * 15) deg, \(theirs.count) spots of the second grain")
            }
        }
        #expect(neighbourSpotCounts.allSatisfy { $0 >= 2 }, "fixture: every boundary position must hold at least two spots of the second grain")
        for _ in 0..<6 { scan.append(pixels(matrixSpots)); residual.append([]); removed.append(matrixSpots.count); kind.append("matrix") }

        let result = try #require(PhaseVectorMatcher.classifyInTwoPasses(
            residual: MatrixResidualScan(scan: residual, removedCounts: removed, inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count)),
            matrixOrientations: [simd_double3x3?](repeating: m, count: scan.count),
            phases: [TwoPassPhase(candidate: candidate, symmetry: ps)], matrixSymmetry: ms,
            settings: PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach),
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
            matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scan.count, minimumPositions: 10))
        #expect(!result.relationships[0].isEmpty, "fixture: the rich positions must teach the relationship")
        var called: [String] = []
        for (i, k) in kind.enumerated() {
            let v = result.verdicts[i]
            if k == "rich" {
                #expect(v.verdict == .indexed, "rich \(i): \(v.verdict)")
            } else {
                if v.verdict != .matrix { called.append("\(k): \(v.verdict), matched \(v.matchedCount) of \(v.survivingCount)") }
            }
        }
        #expect(called.isEmpty, "positions of two matrix grains not called matrix: \(called)")
    }
}
