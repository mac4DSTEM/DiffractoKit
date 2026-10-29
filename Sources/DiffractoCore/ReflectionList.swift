//
//  ReflectionList.swift
//  Role: every (h, k, l) with 0 < |g| ≤ k_max, its d, whether the space-group operations allow it, and the rings
//        (families) the list falls into. The enumeration is reference commit 7448401b `Crystal.reflections` without the
//        structure factor (the per-axis index bound, DEVIATION of 2026-09-14, complete for oblique cells); the absence
//        rule is the general reflection condition; the family rule is macTEM before v1.0.0 (history squashed) `DSpacingIndexing.families`.
//

import Foundation
import simd

public enum ReflectionListError: Error, Equatable, Sendable, CustomStringConvertible {
    /// k_max must be finite, positive and at most `ReflectionList.maximumKMax` Å⁻¹. (A reflection with |g| equal to k_max
    /// is in or out by the last bit of |g|, which can differ by one ulp between two settings of the same lattice: the
    /// lists of a file's cell and of its conventional cell may differ by such a reflection.)
    case invalidKMax(Double)
    /// The cell is not a valid lattice (`UnitCell.issues`), or so large that the list would exceed
    /// `ReflectionList.maximumReflections` (about (4/3)π k_max³ V reflections; 255 bytes each).
    case cellTooLargeOrInvalid(estimatedReflections: Double)

    public var description: String {
        switch self {
        case .invalidKMax(let k):
            return "k_max must be a number of Å⁻¹ above 0 and at most \(ReflectionList.maximumKMax); \(k) is not."
        case .cellTooLargeOrInvalid(let n):
            return n.isFinite
                ? String(format: "This cell at this k_max holds about %.3g reflections, or needs a wider range of indices, than a "
                         + "list is built to by default (1e7 reflections, %d indices); lower k_max, or build a larger "
                         + "list by setting ReflectionList.Policy.maximumReflections (and maximumEnumeration).",
                         n, ReflectionList.maximumEnumeration)
                : "The cell is not a valid finite lattice, so no reflection list can be built."
        }
    }
}

/// `localizedDescription` (what an app shows) is the same text as `description`.
extension ReflectionListError: LocalizedError {
    public var errorDescription: String? { description }
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

    package init(hkl: MillerIndex, g: SIMD3<Double>, gLength: Double, isAllowed: Bool) {
        self.hkl = hkl
        self.g = g
        self.gLength = gLength
        self.isAllowed = isAllowed
    }
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

    /// Public memberwise initialiser (0.10.0): a consumer can build a family, for example to test a view of one. `members` and `orbits` are
    /// kept as given (the list's own families are descending lexicographic).
    public init(g: Double, members: [MillerIndex], orbits: [MillerIndex], name: String) {
        self.g = g; self.members = members; self.orbits = orbits; self.name = name
    }
}

public struct ReflectionList: Sendable, Equatable {

    /// The refusal bounds and tolerances of the list (DiffractoKit 0.2.0): the suite's choices, with `Policy.default` the
    /// set every release before used. A cap that is not a positive number (NaN, zero, negative) refuses every list
    /// (`invalidKMax` for `maximumKMax`, `cellTooLargeOrInvalid` for the counts); a tolerance that is NaN or negative is
    /// taken as 0 (exact comparison: with `familyTolerance` 0 only equal |g| share a ring, with `absenceTolerance` 0 a
    /// phase must be an integer to the last bit).
    ///
    /// **The caps are the protection.** `maximumKMax`, `maximumReflections` and `maximumEnumeration` are what keep a list
    /// within the memory and time a caller expects. Removing them (an infinite `maximumKMax`, the two counts at `Int.max`)
    /// removes that protection: a list at k_max 100 was measured to reach 10 GB resident in 85 s before it was stopped.
    /// Raise a cap to what the machine can hold, not to infinity; the kinematic list's `retainingAll: false` holds far less
    /// than the default for the same reflections.
    public struct Policy: Sendable, Equatable {
        /// The largest k_max a list is built to (Å⁻¹). Default 4.
        public var maximumKMax: Double
        /// The most reflections a list is built to hold (estimated by (4/3)π k_max³ V). Default 10 million.
        public var maximumReflections: Int
        /// The most (h, k, l) the enumeration box may hold. Default 100 million.
        public var maximumEnumeration: Int
        /// Two reflections are one family (the ring rule) while their |g| differ by at most this, relative to the
        /// family's first. Default 1e-9.
        public var familyTolerance: Double
        /// A reflection is absent when an operation's phase h·t lies further than this from an integer. Default 1e-9.
        public var absenceTolerance: Double

        public init(maximumKMax: Double = ReflectionList.maximumKMax,
                    maximumReflections: Int = ReflectionList.maximumReflections,
                    maximumEnumeration: Int = ReflectionList.maximumEnumeration,
                    familyTolerance: Double = ReflectionList.familyTolerance,
                    absenceTolerance: Double = 1e-9) {
            self.maximumKMax = maximumKMax
            self.maximumReflections = maximumReflections
            self.maximumEnumeration = maximumEnumeration
            self.familyTolerance = familyTolerance
            self.absenceTolerance = absenceTolerance
        }

        /// The suite's values.
        public static let `default` = Policy()

        /// The tolerance used: `familyTolerance` when it is a number of at least 0, else 0.
        var effectiveFamilyTolerance: Double { familyTolerance.isFinite && familyTolerance >= 0 ? familyTolerance : 0 }
        package var effectiveAbsenceTolerance: Double { absenceTolerance.isFinite && absenceTolerance >= 0 ? absenceTolerance : 0 }
    }

    /// The largest k_max a list is built to (Å⁻¹): the index loops run for minutes far above it.
    public static let maximumKMax: Double = 4
    /// Two reflections are one family while their |g| differ by at most this, relative to the family's first.
    public static let familyTolerance: Double = 1e-9
    /// Two |g| are tied, for the list's order, within this relative difference.
    public static let orderTolerance: Double = 1e-12
    /// The most reflections a list is built to hold (measured about 800 B each with their structure factors: about 8 GB);
    /// above it the list is refused, never trapped.
    public static let maximumReflections = 10_000_000
    /// The most (h, k, l) the enumeration box may hold: a needle-shaped cell has few reflections by volume but a long box
    /// (a = 1e7, b = c = 1e-3 Å passed the volume estimate and enumerated 24 million).
    public static let maximumEnumeration = 100_000_000

    public let kMax: Double
    /// The policy the list was built under (its ring rule is the policy's `familyTolerance`).
    public let policy: Policy
    /// Increasing |g|; equal |g| by descending lexicographic (h, k, l).
    public let reflections: [Reflection]
    /// The Laue group the families' orbits are named by.
    public let laueGroup: PointSymmetry
    public let laueClass: LaueClass?

    /// Whether orbit names print four Miller–Bravais indices: a hexagonal-axes Laue class (-3, -3m, 6/m, 6/mmm) in a
    /// cell whose metric is hexagonal (a = b, γ = 120°). A trigonal structure in rhombohedral axes keeps three
    /// indices.
    public let usesMillerBravais: Bool

    /// The list for `cell` with the absences of `operations` and the orbits of `laueGroup`.
    public init(cell: UnitCell, operations: [SymmetryOperation], laueGroup: PointSymmetry, laueClass: LaueClass?,
                kMax: Double, policy: Policy = .default) throws(ReflectionListError) {
        try ReflectionList.checkBounds(cell: cell, kMax: kMax, policy: policy)
        self.kMax = kMax
        self.policy = policy
        self.laueGroup = laueGroup
        self.laueClass = laueClass
        self.usesMillerBravais = ReflectionList.usesMillerBravais(laueClass: laueClass, cell: cell)

        var out: [Reflection] = []
        let absenceTolerance = policy.effectiveAbsenceTolerance
        ReflectionList.enumerate(cell: cell, kMax: kMax) { hkl, g, gLen in
            out.append(Reflection(hkl: hkl, g: g, gLength: gLen,
                                  isAllowed: ReflectionList.isAllowed(hkl, by: operations, tolerance: absenceTolerance)))
        }
        self.reflections = ReflectionList.ordered(out)
    }

    // MARK: Shared with the kinematic list (DiffractoScattering) so that its streaming path enumerates the same set

    /// The refusals of `init`: k_max within the policy's cap, a valid cell, the estimated count and the index range within
    /// the policy's bounds. Nothing is enumerated.
    package static func checkBounds(cell: UnitCell, kMax: Double, policy: Policy) throws(ReflectionListError) {
        guard kMax.isFinite, kMax > 0, kMax <= policy.maximumKMax else { throw .invalidKMax(kMax) }
        guard cell.isValid else { throw .cellTooLargeOrInvalid(estimatedReflections: .infinity) }
        let estimate = 4.0 / 3.0 * Double.pi * kMax * kMax * kMax * cell.volume
        let box = cell.realLattice.map { 2 * (kMax * length($0)).rounded(.up) + 1 }.reduce(1, *)
        guard estimate.isFinite, box.isFinite, estimate <= Double(policy.maximumReflections),
              box <= Double(policy.maximumEnumeration)
        else { throw .cellTooLargeOrInvalid(estimatedReflections: estimate) }
    }

    /// Whether orbit names print four Miller–Bravais indices (see `usesMillerBravais`).
    package static func usesMillerBravais(laueClass: LaueClass?, cell: UnitCell) -> Bool {
        (laueClass?.usesMillerBravais ?? false)
            && cell.metricFamily(lengthTolerance: SymmetryDetection.defaultMetricLengthTolerance,
                                 angleToleranceDeg: SymmetryDetection.defaultMetricAngleToleranceDeg) == .hexagonal
    }

    /// Calls `body` with every (h, k, l) with 0 < |g| <= `kMax`, its g and |g|, in the enumeration's order (h, then k, then l
    /// ascending), which is not the list's order (`ordered`). Holds nothing: the cell must have passed `checkBounds`.
    package static func enumerate(cell: UnitCell, kMax: Double,
                                  _ body: (MillerIndex, SIMD3<Double>, Double) throws -> Void) rethrows {
        let tile = cell.realLattice.map { max(0, Int((kMax * length($0)).rounded(.up))) }
        for h in -tile[0]...tile[0] {
            for k in -tile[1]...tile[1] {
                for l in -tile[2]...tile[2] {
                    if h == 0 && k == 0 && l == 0 { continue }
                    let hkl = MillerIndex(h, k, l)
                    let g = cell.gVector(hkl)
                    let gLen = length(g)
                    if gLen > kMax { continue }
                    try body(hkl, g, gLen)
                }
            }
        }
    }

    /// The general reflection condition: absent when an operation (R, t) with (h k l) R = (h k l)
    /// has (h k l) · t further than 1e-9 from an integer. For such an operation h · t does not depend on the origin,
    /// and a detected translation is the least-squares one, whose rounding error (ē − R ē, ē the sites' mean rounding)
    /// has h · (ē − R ē) = 0 exactly when h R = h: so a shifted or coarsely written cell gives the exact answer.
    public static func isAllowed(_ hkl: MillerIndex, by operations: [SymmetryOperation],
                                 tolerance: Double = 1e-9) -> Bool {
        let bound = tolerance.isFinite && tolerance >= 0 ? tolerance : 0
        for op in operations where hkl.applying(op.rotation) == hkl {
            let phase = hkl.dot(op.translation)
            if abs(phase - phase.rounded()) > bound { return false }
        }
        return true
    }

    /// Increasing |g|; within a run of |g| tied to `orderTolerance`, descending lexicographic (h, k, l).
    static func ordered(_ reflections: [Reflection]) -> [Reflection] {
        ordered(reflections, by: { $0 })
    }

    /// The same order over any elements that carry a reflection (the kinematic list orders its kept reflections with it).
    package static func ordered<T>(_ items: [T], by reflection: (T) -> Reflection) -> [T] {
        var sorted = items.sorted { lhs, rhs in
            let l = reflection(lhs), r = reflection(rhs)
            if l.gLength != r.gLength { return l.gLength < r.gLength }
            return l.hkl > r.hkl
        }
        var start = 0
        while start < sorted.count {
            var end = start + 1
            let g0 = reflection(sorted[start]).gLength
            while end < sorted.count,
                  reflection(sorted[end]).gLength - g0 <= orderTolerance * max(g0, reflection(sorted[end]).gLength) {
                end += 1
            }
            if end - start > 1 { sorted[start..<end].sort { reflection($0).hkl > reflection($1).hkl } }
            start = end
        }
        return sorted
    }

    public var allowed: [Reflection] { reflections.filter(\.isAllowed) }

    /// The rings, on the allowed reflections by default; `allowedOnly: false` groups every
    /// reflection of the lattice.
    public func families(allowedOnly: Bool = true) -> [ReflectionFamily] {
        ReflectionList.families(of: allowedOnly ? allowed : reflections, laueGroup: laueGroup,
                                millerBravais: usesMillerBravais, familyTolerance: policy.effectiveFamilyTolerance)
    }

    /// The ring rule on any reflections in the list's order (increasing |g|): a reflection joins the current ring while
    /// |g| ≤ g₀(1 + `familyTolerance`); orbits and names by `laueGroup` (DiffractoScattering applies it after the |F| cut).
    public static func families(of reflections: [Reflection], laueGroup: PointSymmetry,
                                millerBravais: Bool, familyTolerance: Double = ReflectionList.familyTolerance)
        -> [ReflectionFamily] {
        let ringTolerance = familyTolerance.isFinite && familyTolerance >= 0 ? familyTolerance : 0
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
        for r in reflections {
            if let first = group.first, r.gLength > first.gLength * (1 + ringTolerance) { close() }
            group.append(r)
        }
        close()
        return families
    }
}
