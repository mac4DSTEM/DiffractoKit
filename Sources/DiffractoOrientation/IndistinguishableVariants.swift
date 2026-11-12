//
//  IndistinguishableVariants.swift
//  Role: which variants of an orientation relationship a spot pattern cannot tell apart (DiffractoKit#107).
//
//  A variant is an orientation of the precipitate in the matrix. Seen along the beam, two variants can predict the same
//  spots within the detector's reach: for a beta'' needle lying in the foil plane of aluminium on [001], the two variants
//  that share the needle's direction do (measured on the simulated three-grain demo: each explains 96 of 96 vectors of
//  the needle H2, and no spot belongs to only one of them). A classification by spot positions then returns one of the
//  two by the order of a tie, and a caller that shows that one number shows a coin toss. This file says which variants
//  are one choice. Geometry only: the predicted spots and the tolerance within which two spots are one spot; nothing
//  measured enters, and no threshold of its own.
//
//  Checked against planted truth (`IndistinguishableVariantsTests`).
//

import Foundation
import simd
import DiffractoCore

extension PhaseVectorMatcher {

    /// The variants grouped by the spots they predict: two variants are in one group when every spot either predicts
    /// within the reach has a spot of the other within `toleranceInvAngstrom` (and so on through a chain of such pairs).
    /// A group of one is a variant the pattern can tell from every other.
    ///
    /// - Parameters:
    ///   - variants: the relations (precipitate to matrix) of the variants, as `LearnedOrientationRelationship.variants`
    ///     and `OrientationVariantSet.variants` hold them.
    ///   - matrixOrientation: the matrix orientation at the position (columns lab x, y, z in crystal coordinates). The
    ///     groups depend on it: which needles lie in the plane is a matter of the beam.
    ///   - reflections: the phase's reflections (``PhaseOrientationCandidate/reflections(of:kMax:)``).
    ///   - excitationSlabInvAngstrom: the slab the phase was matched with; with the reflections it predicts the spots.
    ///   - reachInvAngstrom: spots farther from the centre are not recorded and tell nothing apart.
    ///   - toleranceInvAngstrom: two predicted spots nearer than this are one spot: the pair radius of the matching.
    ///   - directBeamRadiusInvAngstrom: spots inside it are the direct beam.
    /// - Returns: the groups, each ascending, ordered by their first variant. Every index 0..<variants.count is in one.
    public static func indistinguishableVariants(_ variants: [simd_double3x3], matrixOrientation: simd_double3x3,
                                                 reflections: [ACOMReflection], excitationSlabInvAngstrom: Double,
                                                 reachInvAngstrom: Double, toleranceInvAngstrom: Double = 0.02,
                                                 directBeamRadiusInvAngstrom: Double = 0.15) -> [[Int]] {
        let spots: [[SIMD2<Double>]] = variants.map { relation in
            matrixSpots(reflections: reflections, orientation: relation.transpose * matrixOrientation,
                        excitationSlabInvAngstrom: excitationSlabInvAngstrom)
                .filter { let r = simd_length($0); return r > directBeamRadiusInvAngstrom && r <= reachInvAngstrom }
        }
        func covers(_ a: [SIMD2<Double>], _ b: [SIMD2<Double>]) -> Bool {
            a.allSatisfy { s in b.contains { simd_length($0 - s) <= toleranceInvAngstrom } }
        }
        // Union-find over the pairs that predict the same spots.
        var parent = Array(variants.indices)
        func root(_ i: Int) -> Int { var r = i; while parent[r] != r { r = parent[r] }; return r }
        for i in variants.indices {
            for j in variants.indices where j > i {
                if covers(spots[i], spots[j]) && covers(spots[j], spots[i]) {
                    let a = root(i), b = root(j)
                    if a != b { parent[max(a, b)] = min(a, b) }
                }
            }
        }
        var groups: [Int: [Int]] = [:]
        for i in variants.indices { groups[root(i), default: []].append(i) }
        return groups.keys.sorted().map { groups[$0]!.sorted() }
    }
}
