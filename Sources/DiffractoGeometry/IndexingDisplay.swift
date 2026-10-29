//
//  IndexingDisplay.swift
//  Ported from macTEM TEMCore/Diffraction/IndexingDisplay.swift (the app's source is the oracle; pins in P23).
//  Role: a ring's indexing AS THE ROOM DISPLAYS IT: the kit's `RingIndexing` decision over a kinematic list, with
//        the numbers a reader sees derived once, here, from the kit's own fields. The decision is the kit's, untouched; this file
//        adds the percentages, the odds, σ_q from the typed tolerance and the verdict's name.
//
//  The tolerance is the user's stated total 1-sigma uncertainty of |q| in PERCENT as typed: σ_q = q · tolerance / 100. δ is the
//  kit's signed deviation g/q − 1 (d_measured − d_listed over d_listed), shown in percent; |Δd|/d is its magnitude; z is
//  |q − g| / σ_q; odds = exp(margin), the likelihood ratio of the best family over the next, shown "> 10⁶ : 1" beyond 10⁶.
//  The product self-reports `validation: "none"` (RingIndexing.validation): never called validated here.
//

import Foundation
import DiffractoCore
import DiffractoScattering

/// A ring's indexing as a room shows it: the kit's decision plus the numbers derived from its fields.
public struct DisplayedIndexing: Sendable, Equatable {
    /// The measured |q| (Å⁻¹) and the standard uncertainty derived from the typed tolerance (Å⁻¹; NaN when q or the tolerance is not finite).
    public let q: Double
    public let sigmaQ: Double
    /// "accepted", or the kit's refusal name: "invalidInput", "ambiguous", "beyondKMax", "noneWithinTolerance".
    public let verdict: String
    /// The accepted assignment, nil when refused.
    public let assignment: RingIndexing.Assignment?
    /// Why the indexing was refused, nil when accepted.
    public let refusal: RingIndexing.Refusal?

    /// The kit's own wording of the refusal, nil when accepted.
    public var refusalText: String? { refusal?.reason }
    /// True when the indexing was accepted.
    public var isAccepted: Bool { assignment != nil }
    /// The label the kit gives the product: "none: ...".
    public var validation: String { RingIndexing.validation }

    /// Of the accepted family; nil otherwise.
    public var family: ReflectionFamily? { assignment?.best.family }
    /// d of the listed family and d = 1/q measured (Å).
    public var dListed: Double? { assignment?.best.family.d }
    public var dMeasured: Double? { assignment.map { 1 / $0.q } }
    /// δ = g/q − 1 in percent (signed) and |Δd|/d (a fraction, unsigned).
    public var deltaPercent: Double? { assignment.map { 100 * $0.best.deviation } }
    public var relativeDeltaD: Double? { assignment?.best.relativeDeviation }
    public var z: Double? { assignment?.best.z }
    /// The log-likelihood ratio over the next family; nil with no next family.
    public var margin: Double? { assignment?.margin }
    /// exp(margin): the odds the best family has over the next; +infinity where exp overflows; nil with no margin.
    public var odds: Double? { margin.map { exp($0) } }
    /// True when the odds exceed 10⁶ : 1, where the room shows `oddsCapLabel` and not a number.
    public var oddsBeyondCap: Bool { odds.map { $0 > DisplayedIndexing.oddsCap } ?? false }

    /// Odds above this show as `oddsCapLabel`.
    public static let oddsCap = 1e6
    /// The label shown beyond the cap.
    public static let oddsCapLabel = "> 10⁶ : 1"

    /// "n : 1" with three significant digits, or `oddsCapLabel` beyond the cap; "—" where there is no next family to compare with.
    public var oddsText: String {
        guard let odds else { return "—" }
        if odds > DisplayedIndexing.oddsCap { return DisplayedIndexing.oddsCapLabel }
        return String(format: "%.3g : 1", odds)
    }
}

/// Indexes a measured |q| from a typed percent tolerance and derives the displayed numbers.
public enum IndexingDisplay {

    /// Indexes |q| (Å⁻¹), with the tolerance in percent of q, against a kinematic list; the list's k_max is the kit's guard.
    public static func index(q: Double, tolerancePercent: Double, in list: KinematicReflectionList,
                              coverage k: Double = RingIndexing.defaultCoverage) -> DisplayedIndexing {
        index(q: q, tolerancePercent: tolerancePercent, families: list.families.map(\.family), kMax: list.kMax, coverage: k)
    }

    /// The same over a list of families and the k_max it was built to.
    public static func index(q: Double, tolerancePercent: Double, families: [ReflectionFamily], kMax: Double,
                              coverage k: Double = RingIndexing.defaultCoverage) -> DisplayedIndexing {
        let sigma = q * tolerancePercent / 100
        do {
            let assignment = try RingIndexing.index(q: q, sigma: sigma, in: families, kMax: kMax, coverage: k)
            return DisplayedIndexing(q: q, sigmaQ: sigma, verdict: "accepted", assignment: assignment, refusal: nil)
        } catch {
            return DisplayedIndexing(q: q, sigmaQ: sigma, verdict: error.name, assignment: nil, refusal: error)
        }
    }
}
