//
//  ReflectionList.swift
//  Role: every (h, k, l) with 0 < |g| ≤ k_max, its d, whether the space-group operations allow it, and the rings
//        (families) the list falls into. The enumeration is mac4DSTEM 7448401b `Crystal.reflections` without the
//        structure factor (the per-axis index bound, DEVIATION of 2026-09-14, complete for oblique cells); the absence
//        rule is registration 3.8; the family rule is macTEM a58c84c `DSpacingIndexing.families` (ADR 022), 3.9.
//

import Foundation
import simd

public enum ReflectionListError: Error, Equatable, Sendable, CustomStringConvertible {
    /// k_max must be finite, positive and at most `ReflectionList.maximumKMax` Å⁻¹.
    case invalidKMax(Double)

    public var description: String {
        switch self {
        case .invalidKMax(let k):
            return "k_max must be a number of Å⁻¹ above 0 and at most \(ReflectionList.maximumKMax); \(k) is not."
        }
    }
}

/// One (h, k, l) of the list.
public struct Reflection: Sendable, Equatable {
    public let hkl: MillerIndex
    /// g = h b₁ + k b₂ + l b₃ (Å⁻¹, Cartesian).
    public let g: SIMD3<Double>
    /// |g| (Å⁻¹).
    public let gLength: Double
    /// 1/|g| (Å).
    public var d: Double { 1 / gLength }
    /// False when an operation (R, t) with (h k l) R = (h k l) has (h k l) · t not an integer: a systematic absence.
    public let isAllowed: Bool
}

/// A ring: the reflections of one |g|.
public struct ReflectionFamily: Sendable, Equatable {
    /// The smallest member's |g| (Å⁻¹).
    public let g: Double
    /// Every member (Friedel mates included), descending lexicographic.
    public let members: [MillerIndex]
    /// The canonical forms of the members' orbits, descending lexicographic.
    public let orbits: [MillerIndex]
    /// The orbits' names joined by " + " ("{1 1 1}", "{1 0 -1 1}", "{3 3 3} + {5 1 1}").
    public let name: String

    public var d: Double { 1 / g }
    public var multiplicity: Int { members.count }
}

public struct ReflectionList: Sendable, Equatable {
    /// The largest k_max a list is built to (Å⁻¹): the index loops run for minutes far above it.
    public static let maximumKMax: Double = 4
    /// Two reflections are one family while their |g| differ by at most this, relative to the family's first.
    public static let familyTolerance: Double = 1e-9
    /// Two |g| are tied, for the list's order, within this relative difference.
    public static let orderTolerance: Double = 1e-12

    public let kMax: Double
    /// Increasing |g|; equal |g| by descending lexicographic (h, k, l).
    public let reflections: [Reflection]
    /// The Laue group the families' orbits are named by.
    public let laueGroup: PointSymmetry
    public let laueClass: LaueClass?

    /// The list for `cell` with the absences of `operations` and the orbits of `laueGroup`.
    public init(cell: UnitCell, operations: [SymmetryOperation], laueGroup: PointSymmetry, laueClass: LaueClass?,
                kMax: Double) throws(ReflectionListError) {
        guard kMax.isFinite, kMax > 0, kMax <= ReflectionList.maximumKMax else { throw .invalidKMax(kMax) }
        self.kMax = kMax
        self.laueGroup = laueGroup
        self.laueClass = laueClass

        let tile = cell.realLattice.map { max(0, Int((kMax * length($0)).rounded(.up))) }
        var out: [Reflection] = []
        for h in -tile[0]...tile[0] {
            for k in -tile[1]...tile[1] {
                for l in -tile[2]...tile[2] {
                    if h == 0 && k == 0 && l == 0 { continue }
                    let hkl = MillerIndex(h, k, l)
                    let g = cell.gVector(hkl)
                    let gLen = length(g)
                    if gLen > kMax { continue }
                    out.append(Reflection(hkl: hkl, g: g, gLength: gLen,
                                          isAllowed: ReflectionList.isAllowed(hkl, by: operations)))
                }
            }
        }
        self.reflections = ReflectionList.ordered(out)
    }

    /// The general reflection condition (registration 3.8).
    public static func isAllowed(_ hkl: MillerIndex, by operations: [SymmetryOperation]) -> Bool {
        for op in operations where hkl.applying(op.rotation) == hkl {
            let phase = hkl.dot(op.translation)
            if abs(phase - phase.rounded()) > 1e-9 { return false }
        }
        return true
    }

    /// Increasing |g|; within a run of |g| tied to `orderTolerance`, descending lexicographic (h, k, l).
    static func ordered(_ reflections: [Reflection]) -> [Reflection] {
        var sorted = reflections.sorted { lhs, rhs in
            if lhs.gLength != rhs.gLength { return lhs.gLength < rhs.gLength }
            return lhs.hkl > rhs.hkl
        }
        var start = 0
        while start < sorted.count {
            var end = start + 1
            let g0 = sorted[start].gLength
            while end < sorted.count, sorted[end].gLength - g0 <= orderTolerance * max(g0, sorted[end].gLength) {
                end += 1
            }
            if end - start > 1 { sorted[start..<end].sort { $0.hkl > $1.hkl } }
            start = end
        }
        return sorted
    }

    public var allowed: [Reflection] { reflections.filter(\.isAllowed) }

    /// The rings (registration 3.9), on the allowed reflections by default; `allowedOnly: false` groups every
    /// reflection of the lattice.
    public func families(allowedOnly: Bool = true) -> [ReflectionFamily] {
        let source = allowedOnly ? allowed : reflections
        let millerBravais = laueClass?.usesMillerBravais ?? false
        var families: [ReflectionFamily] = []
        var group: [Reflection] = []
        func close() {
            guard let first = group.first else { return }
            let members = group.map(\.hkl).sorted(by: >)
            let orbits = Set(members.map { laueGroup.canonical($0) }).sorted(by: >)
            families.append(ReflectionFamily(
                g: group.map(\.gLength).min() ?? first.gLength,
                members: members,
                orbits: orbits,
                name: orbits.map { PointSymmetry.orbitName(ofCanonical: $0, millerBravais: millerBravais) }
                    .joined(separator: " + ")))
            group = []
        }
        for r in source {
            if let first = group.first, r.gLength > first.gLength * (1 + ReflectionList.familyTolerance) { close() }
            group.append(r)
        }
        close()
        return families
    }
}
