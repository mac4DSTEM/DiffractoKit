//
//  IndistinguishableVariantsTests.swift
//  DiffractoKit#107: which variants a spot pattern cannot tell apart. Planted truth: beta'' in aluminium at Andersen's
//  relationship (needles along the matrix's cube axes), the matrix seen down [001].
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct IndistinguishableVariantsTests {
    static let reach = 0.75

    struct Setup {
        let variants: [simd_double3x3]
        let candidate: PhaseOrientationCandidate
        let matrix: simd_double3x3
        /// Per variant: true when its needle axis (b of beta'') lies along the beam.
        let endOn: [Bool]
    }

    static func setup(zone: SIMD3<Double> = SIMD3(0, 0, 1)) throws -> Setup {
        let al = CrystalPresets.aluminum, beta = CrystalPresets.betaDoublePrime
        let relation = try OrientationRelationship.rotation(
            OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
            OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(3, 1, 0))),
            precipitate: beta, matrix: al)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(al.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(beta.laueClass))
        let variants = LearnedOrientationRelationship(relation: relation, positionCount: 0, spreadRad: 0, matrixSymmetry: ms, precipitateSymmetry: ps).variants
        let candidate = try PhaseOrientationCandidate(name: "beta''", crystal: beta, kMax: reach)
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: zone), inPlaneAngle: Float(12 * Double.pi / 180))
        // b of the monoclinic cell in the crystal's Cartesian frame, carried into the matrix frame by the variant, against the beam.
        let b = simd_normalize(beta.cell.realLattice[1])
        let beam = simd_normalize(m.columns.2)
        let endOn = variants.map { abs(simd_dot(simd_normalize($0 * b), beam)) > 0.999 }
        return Setup(variants: variants, candidate: candidate, matrix: m, endOn: endOn)
    }

    /// Aluminium down [001]: the four variants whose needle is end-on are each their own group, and the eight whose needle
    /// lies in the plane fall into four pairs, each pair two needles of one direction. Every variant is in exactly one group.
    /// Mutation: the comparison made one way only (a variant whose spots are a subset of another's would join it).
    @Test func inPlaneNeedlesComeInPairsAndEndOnNeedlesAlone() throws {
        let s = try Self.setup()
        #expect(s.variants.count == 12 && s.endOn.filter { $0 }.count == 4)
        let groups = PhaseVectorMatcher.indistinguishableVariants(s.variants, matrixOrientation: s.matrix, candidate: s.candidate,
                                                                  reachInvAngstrom: Self.reach)
        #expect(groups.flatMap { $0 }.sorted() == Array(0..<12), "every variant once: \(groups)")
        #expect(groups == groups.sorted { $0[0] < $1[0] } && groups.allSatisfy { $0 == $0.sorted() })
        for g in groups {
            if g.contains(where: { s.endOn[$0] }) {
                #expect(g.count == 1, "an end-on variant shares a group: \(g)")
            } else {
                #expect(g.count == 2, "in-plane variants not in a pair: \(g)")
                // The two of a pair are needles of one direction in the matrix.
                let b = simd_normalize(CrystalPresets.betaDoublePrime.cell.realLattice[1])
                let d0 = simd_normalize(s.variants[g[0]] * b), d1 = simd_normalize(s.variants[g[1]] * b)
                #expect(abs(simd_dot(d0, d1)) > 0.999, "a pair of two needle directions: \(g)")
            }
        }
        #expect(groups.filter { $0.count == 1 }.count == 4 && groups.filter { $0.count == 2 }.count == 4, "\(groups)")
    }

    /// A tolerance of nothing finds no two variants the same only if their spots differ at all: the pairs of the test above
    /// predict the same spots exactly, so they stay pairs; and no reach means no spots, so nothing can be told apart.
    @Test func theGroupsFollowTheGeometryNotTheTolerance() throws {
        let s = try Self.setup()
        let tight = PhaseVectorMatcher.indistinguishableVariants(s.variants, matrixOrientation: s.matrix, candidate: s.candidate,
                                                                 reachInvAngstrom: Self.reach, toleranceInvAngstrom: 1e-6)
        #expect(tight.filter { $0.count == 2 }.count == 4 && tight.filter { $0.count == 1 }.count == 4, "\(tight)")
        let blind = PhaseVectorMatcher.indistinguishableVariants(s.variants, matrixOrientation: s.matrix, candidate: s.candidate,
                                                                 reachInvAngstrom: 0.1)
        #expect(blind == [Array(0..<12)], "with no spot in reach every variant is the same: \(blind)")
        #expect(PhaseVectorMatcher.indistinguishableVariants([], matrixOrientation: s.matrix, candidate: s.candidate, reachInvAngstrom: Self.reach).isEmpty)
    }
}
