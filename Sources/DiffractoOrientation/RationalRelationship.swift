//
//  RationalRelationship.swift
//  Role: the rational orientation relationships near a learned one, reported, never chosen (DiffractoKit#81).
//
//  A published orientation relationship is two statements of the kind "a low-index plane of the precipitate is parallel to a
//  low-index plane of the matrix" or the same for directions (Burgers, Bagaryatsky, Andersen and the like). A relationship
//  learned from a scan (OrientationRelationshipLearning.swift) is a measured rotation, a few degrees off such a statement.
//  This file lists the statements that lie within a stated angle of the measured rotation.
//
//  It reports and does not choose. Between many pairs of structures several rational relationships lie within a few degrees
//  of each other (hcp and bcc: Burgers, Pitsch-Schrader, Potter; ferrite and cementite: Bagaryatsky, Isaichev), and a
//  relationship measured less precisely than they lie apart cannot tell them apart: the search then says so
//  (`isAmbiguous`). Nothing in the classification uses it.
//

import Foundation
import simd
import DiffractoCore

/// A rational orientation relationship: two statements and the rotation they fix.
public struct RationalRelationship: Sendable {
    /// The statement that fixes an axis exactly.
    public let first: OrientationRelationship
    /// The statement that fixes the rotation about it.
    public let second: OrientationRelationship
    /// Precipitate Cartesian vectors to matrix Cartesian vectors (``OrientationRelationship/rotation(_:_:precipitate:matrix:)``).
    public let relation: simd_double3x3
    /// The angle between this relationship and the measured rotation it was found for, radians: the smallest over the
    /// statements of it that were found (`relation` itself may be another variant, farther off).
    public let angleRad: Double
    /// The sum of the absolute indices of the four lattice vectors.
    public let indexSum: Int
    /// How far the second statement is from holding under `relation`, radians: the two statements of a relationship are
    /// consistent only when the angle between the two precipitate vectors equals the angle between the two matrix vectors.
    public let secondMisfitRad: Double
}

/// What ``PhaseVectorMatcher/rationalRelationships(near:precipitate:matrix:matrixSymmetry:precipitateSymmetry:maximumIndex:maximumIndexSum:limitRad:)`` found.
public struct RationalRelationshipSearch: Sendable {
    /// The distinct rational relationships within the limit, the nearest to the measured relation first. Two statements of
    /// one rotation (up to both symmetries) are one relationship, kept in its statement of the smallest index sum.
    public let candidates: [RationalRelationship]
    /// More than one relationship lies within the limit: the measurement does not say which it is.
    public var isAmbiguous: Bool { candidates.count > 1 }
}

extension PhaseVectorMatcher {

    /// The rational relationships within `limitRad` of a measured one.
    ///
    /// Every plane and direction with coprime indices up to `maximumIndex` of the precipitate is carried by `relation` into the
    /// matrix and paired with the matrix planes (a plane) or directions (a direction) it lands on within `limitRad`: a plane
    /// with a plane and a direction with a direction, as relationships are stated. Every two such pairs that are not parallel
    /// fix a rotation; those within `limitRad` of `relation` whose four lattice vectors have an index sum of at most
    /// `maximumIndexSum` are returned, each rotation once.
    ///
    /// - Parameters:
    ///   - relation: precipitate Cartesian vectors to matrix Cartesian vectors (a learned relationship's `relation`).
    ///   - maximumIndexSum: the largest sum of absolute indices of the four lattice vectors. Without a bound there is a
    ///     rational relationship within any angle of any rotation. The published ones are simple: cube on cube 4, Burgers,
    ///     Bagaryatsky and Andersen's beta'' relationship 7. Default 8.
    ///   - limitRad: how far a rational relationship may lie from the measured one: the measurement's uncertainty. Below a
    ///     right angle.
    ///   - statementToleranceRad: two statements are a relationship only when both can hold at once: the second may be at most
    ///     this far from parallel under the rotation the first fixes (default 0.25°). Two relationships closer than this are
    ///     listed as one. Without it, two statements that merely both lie near the measured relation would be listed in each
    ///     of their two orders as two relationships.
    /// - Returns: no candidates when `relation` is not a rotation, or the limit or the bounds are not usable.
    public static func rationalRelationships(near relation: simd_double3x3,
                                             precipitate: CrystalStructure, matrix: CrystalStructure,
                                             matrixSymmetry: CrystalOrientationSymmetry, precipitateSymmetry: CrystalOrientationSymmetry,
                                             maximumIndex: Int = 3, maximumIndexSum: Int = 8,
                                             limitRad: Double,
                                             statementToleranceRad: Double = 0.25 * Double.pi / 180) -> RationalRelationshipSearch {
        guard maximumIndex >= 1, maximumIndexSum >= 4, limitRad > 0, limitRad < Double.pi / 2,
              statementToleranceRad >= 0, statementToleranceRad.isFinite else { return RationalRelationshipSearch(candidates: []) }
        // A rotation: orthonormal columns, right-handed.
        let gram = relation.transpose * relation
        let identity = matrix_identity_double3x3
        let orthonormal = (0..<3).allSatisfy { c in (0..<3).allSatisfy { r in abs(gram[c][r] - identity[c][r]) < 1e-6 } }
        guard orthonormal, relation.determinant > 0 else { return RationalRelationshipSearch(candidates: []) }
        func gcd(_ a: Int, _ b: Int) -> Int { b == 0 ? abs(a) : gcd(b, a % b) }
        // One of every +/- pair of coprime index triples.
        var triples: [SIMD3<Int>] = []
        // A lattice vector of a pair weighs at most the index sum less the three other vectors' least.
        let n = min(maximumIndex, maximumIndexSum - 3)
        for h in -n...n { for k in -n...n { for l in -n...n {
            guard (h, k, l) != (0, 0, 0), gcd(gcd(h, k), l) == 1, abs(h) + abs(k) + abs(l) <= maximumIndexSum - 3 else { continue }
            let firstNonZero = h != 0 ? h : (k != 0 ? k : l)
            if firstNonZero > 0 { triples.append(SIMD3(h, k, l)) }
        } } }
        struct Element { let vector: LatticeVector; let unit: SIMD3<Double>; let weight: Int }
        func elements(of crystal: CrystalStructure) -> [Element] {
            var out: [Element] = []
            for t in triples {
                let weight = abs(t.x) + abs(t.y) + abs(t.z)
                for vector in [LatticeVector.direction(t), LatticeVector.plane(t)] {
                    let c = vector.cartesian(in: crystal)
                    let length = simd_length(c)
                    if length > 0, length.isFinite { out.append(Element(vector: vector, unit: c / length, weight: weight)) }
                }
            }
            return out
        }
        let ofPrecipitate = elements(of: precipitate), ofMatrix = elements(of: matrix)
        let cosLimit = cos(limitRad)
        struct Pair { let candidate: LatticeVector; let matrix: LatticeVector; let unitP: SIMD3<Double>; let unitM: SIMD3<Double>; let weight: Int }
        func negated(_ v: LatticeVector) -> LatticeVector {
            switch v {
            case .direction(let d): return .direction(SIMD3(-d.x, -d.y, -d.z))
            case .plane(let p): return .plane(SIMD3(-p.x, -p.y, -p.z))
            }
        }
        var pairs: [Pair] = []
        for p in ofPrecipitate {
            let carried = relation * p.unit
            for m in ofMatrix where p.weight + m.weight <= maximumIndexSum - 2 {
                switch (p.vector, m.vector) {
                case (.plane, .plane), (.direction, .direction): break
                default: continue
                }
                let d = simd_dot(carried, m.unit)
                guard abs(d) >= cosLimit else { continue }
                pairs.append(Pair(candidate: p.vector, matrix: d >= 0 ? m.vector : negated(m.vector), unitP: p.unit,
                                  unitM: d >= 0 ? m.unit : -m.unit, weight: p.weight + m.weight))
            }
        }
        var found: [RationalRelationship] = []
        for (i, a) in pairs.enumerated() {
            for (j, b) in pairs.enumerated() where j != i && a.weight + b.weight <= maximumIndexSum {
                let cosP = simd_dot(a.unitP, b.unitP), cosM = simd_dot(a.unitM, b.unitM)
                guard abs(cosP) < 0.999, abs(cosM) < 0.999 else { continue }
                // Both statements at once: the same angle between the two vectors on either side.
                let misfit = abs(acos(min(1, max(-1, cosP))) - acos(min(1, max(-1, cosM))))
                guard misfit <= statementToleranceRad else { continue }
                let first = OrientationRelationship(candidate: a.candidate, matrix: a.matrix)
                let second = OrientationRelationship(candidate: b.candidate, matrix: b.matrix)
                guard let r = try? OrientationRelationship.rotation(first, second, precipitate: precipitate, matrix: matrix) else { continue }
                let t = relation.transpose * r
                let angle = acos(min(1, max(-1, (t.columns.0.x + t.columns.1.y + t.columns.2.z - 1) / 2)))
                guard angle <= limitRad else { continue }
                found.append(RationalRelationship(first: first, second: second, relation: r, angleRad: angle, indexSum: a.weight + b.weight,
                                                  secondMisfitRad: misfit))
            }
        }
        // Each rotation once (up to both symmetries), in its simplest statement; then the nearest first.
        found.sort { ($0.indexSum, $0.angleRad) < ($1.indexSum, $1.angleRad) }
        var distinct: [RationalRelationship] = []
        let same = statementToleranceRad + 0.01 * Double.pi / 180
        for f in found {
            if let k = distinct.firstIndex(where: {
                relationDistance($0.relation, f.relation, matrixSymmetry: matrixSymmetry, precipitateSymmetry: precipitateSymmetry) < same
            }) {
                // The same relationship in another statement (or another of its variants): its simplest statement stays, with
                // the smallest angle any statement of it has to the measured relation.
                let kept = distinct[k]
                if f.angleRad < kept.angleRad {
                    distinct[k] = RationalRelationship(first: kept.first, second: kept.second, relation: kept.relation, angleRad: f.angleRad,
                                                       indexSum: kept.indexSum, secondMisfitRad: kept.secondMisfitRad)
                }
            } else {
                distinct.append(f)
            }
        }
        distinct.sort { $0.angleRad < $1.angleRad }
        return RationalRelationshipSearch(candidates: distinct)
    }
}
