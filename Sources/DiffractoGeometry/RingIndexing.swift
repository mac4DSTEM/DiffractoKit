//
//  RingIndexing.swift
//  Role: the (hkl) family a measured |q| belongs to, against a ring list (DiffractoScattering's `KinematicFamily`, or DiffractoCore's
//        `ReflectionFamily`), and the decision of when it is NOT to be named. The decision is reference commit 9482831
//        `DSpacingIndexing.index` (macTEM ADR 022 and ADR 024), generalised: the list is DiffractoCore's family rule over any
//        crystal system ("RingIndexing"), so the lattice, orbit and
//        family-construction code of the seed is not ported and `unsupportedLattice` / `unsupportedElement` cannot occur
//        here (the list builder refuses an element without a form factor, and every system has a Laue group).
//
//  THE INPUT. q = |q| in Å⁻¹ (no 2π; d = 1/q) and σ_q, the STANDARD uncertainty (1σ, Å⁻¹) of q's total error: in
//  quadrature the calibration's (q·u_cal), the measurement's, the lattice constant's (q·u_cell) and, when q came from a
//  linear Δq·R scale, the flat-detector nonlinearity ≈ 0.375 λ² |q² − q₀²| q. σ_q is the caller's; an understated one
//  gives a confident wrong index. k is a Gaussian coverage factor (default 3). `kMax` is the k_max the list was built to
//  (Å⁻¹): the list is complete only below it, which is what `beyondKMax` guards.
//
//  THE DEFINITION (tools/slice-3/index-pins.py is written from the same text):
//    z, δ       z = |q − g| / σ_q; the signed deviation δ = (d_meas − d)/d = g/q − 1; |Δd|/d = |δ|. Candidates order by
//               (z, g).
//    margin     between two candidates, their log-likelihood ratio (z₂² − z₁²)/2; accepting needs ≥ ln 19 (the next at
//               most 1/19 as likely), derived from the Gaussian model, not measured.
//    reach      r = √(k² + 2 ln 19)·σ_q: no family farther than r can change a verdict.
//    verdict    the first that applies: invalidInput; ambiguous (two or more families with z ≤ k, or one whose margin to
//               the next listed family is below ln 19); beyondKMax (q + R > k_max, R = max(r, |q − g| of every family the
//               verdict names), so every number a verdict carries is exact for the whole list, not only for the part of it
//               below k_max); noneWithinTolerance; accepted.
//  Intensity plays no part: every listed family counts, however weak.
//
//  WHAT `accepted` MEANS, and what it does not. Of the LISTED families, this one is at least 19 times as likely as the
//  next (a likelihood ratio at equal priors, Gaussian error of the stated σ_q), and no second family is within kσ_q. It
//  does not mean the ring IS that family: a ring of no listed family (a support film, another phase) is accepted whenever
//  it falls inside one window; an understated σ_q gives confident wrong names in dense cells; a lone spot is not a ring.
//  The decision is `validation: "none"`: parity with a NumPy reading of the rule only, chance acceptance and real σ_q
//  unmeasured (the reference implementation's T3 and T5 numbers stand).
//

import Foundation
import DiffractoCore
import DiffractoScattering

public enum RingIndexing {

    /// Self-reported label: chance acceptance on foreign phases and the real σ_q are unmeasured until a pattern of known
    /// phase is indexed.
    public static let validation =
        "none: parity with a NumPy reading of the rule only; chance acceptance and real sigma_q unmeasured"

    /// k's default: a Gaussian coverage factor of 3 (99.73 % for a correct standard uncertainty).
    public static let defaultCoverage: Double = 3

    /// The least margin an accepted family has over the next: ln 19, the default of every `minimumMargin:` parameter
    /// below. A policy of the suite, derived from the Gaussian model (the next family at most 1/19 as likely), not
    /// measured: a larger margin names fewer rings, 0 accepts a family whatever the next one's likelihood (the
    /// ambiguity within kσ_q still refuses).
    public static let minimumMargin: Double = log(19.0)

    /// The reach in units of σ_q, √(k² + 2 ln 19): 3.86 at k = 3, with the margin that was asked for. NaN for a coverage
    /// that is not a finite positive number or a margin that is negative or not finite: the arguments `index` refuses with
    /// `invalidInput`.
    public static func reach(coverage k: Double, minimumMargin margin: Double = minimumMargin) -> Double {
        guard k.isFinite, k > 0, margin.isFinite, margin >= 0 else { return .nan }
        return (k * k + 2 * margin).squareRoot()
    }

    /// A family judged against one measurement.
    public struct Candidate: Sendable, Equatable {
        public let family: ReflectionFamily
        /// |q − g| / σ_q.
        public let z: Double
        /// The signed deviation (d_meas − d)/d = g/q − 1.
        public let deviation: Double

        public init(family: ReflectionFamily, z: Double, deviation: Double) {
            self.family = family; self.z = z; self.deviation = deviation
        }

        /// |Δd|/d.
        public var relativeDeviation: Double { abs(deviation) }
    }

    /// A measurement named: its family, and the next candidate with the margin between them. `accepted` is "this family
    /// is at least 19 times as likely as the next listed one and no other is within kσ_q", not "the ring is this family".
    public struct Assignment: Sendable, Equatable {
        public let q: Double
        public let sigma: Double
        public let coverage: Double
        public let best: Candidate
        /// The family of second-smallest z; nil only when the list holds no other family.
        public let next: Candidate?
        /// (z_next² − z²)/2, the log-likelihood ratio of the best over the next; nil with `next`.
        public let margin: Double?

        /// `RingIndexing.validation`: "none: ...".
        public var validation: String { RingIndexing.validation }

        /// d_meas = 1/q (Å).
        public var dMeasured: Double { 1 / q }
    }

    /// Why a measurement is not named, by name.
    public enum Refusal: Error, Equatable, Sendable, CustomStringConvertible, LocalizedError {
        /// A number that is not finite and positive (q, σ_q, k, k_max): the caller's mistake, not data.
        case invalidInput(String)
        /// Two or more families within kσ_q, or one whose margin to the next is below ln 19: the candidates, by (z, g).
        case ambiguous([Candidate])
        /// The list, built to k_max, cannot settle the verdict: q + R (Å⁻¹) passes it.
        case beyondKMax(reach: Double, kMax: Double)
        /// No family within kσ_q; the nearest (nil for an empty list).
        case noneWithinTolerance(nearest: Candidate?)

        /// The case's name, without its numbers: a stable key for logs and tools.
        public var name: String {
            switch self {
            case .invalidInput: return "invalidInput"
            case .ambiguous: return "ambiguous"
            case .beyondKMax: return "beyondKMax"
            case .noneWithinTolerance: return "noneWithinTolerance"
            }
        }

        /// One sentence, for the person reading the result.
        public var reason: String {
            switch self {
            case .invalidInput(let detail):
                return "Cannot index: \(detail)."
            case .ambiguous(let candidates):
                let named = candidates.map { "\($0.family.name) (\(String(format: "%.2f", $0.z))σ)" }.joined(separator: ", ")
                return "More than one family fits within the stated uncertainty: \(named)."
            case .beyondKMax(let reach, let kMax):
                return String(format: "The families that could decide this lie up to %.4f Å⁻¹, past the list's limit of %.4f Å⁻¹.",
                              reach, kMax)
            case .noneWithinTolerance(let nearest):
                let causes = "possible causes, unranked: a scale error beyond the stated uncertainty, another phase or a "
                    + "support film, lens distortion, or for a kinematically forbidden position double diffraction or "
                    + "bonding charge"
                guard let nearest else { return "The list holds no family; \(causes)." }
                let off = String(format: "%.1fσ (d off by %+.3f %%)", nearest.z, 100 * nearest.deviation)
                return "No family lies within the stated uncertainty; the nearest is \(nearest.family.name) at \(off); \(causes)."
            }
        }

        public var description: String { reason }
        public var errorDescription: String? { reason }
    }

    /// The family `q` (Å⁻¹) belongs to among the kinematic rings (DiffractoScattering's list after the |F| cut), given its standard
    /// uncertainty `sigma` (Å⁻¹), the k_max the list was built to and the coverage factor k; or the first refusal that
    /// applies (the header's order).
    /// The decision against a kinematic list, with the k_max it was built to (a free k_max larger
    /// than the list's turns the beyondKMax guard off; this overload cannot be given one).
    public static func index(q: Double, sigma: Double, in list: KinematicReflectionList,
                             coverage k: Double = defaultCoverage,
                             minimumMargin margin: Double = minimumMargin) throws(Refusal) -> Assignment {
        try index(q: q, sigma: sigma, in: list.families, kMax: list.kMax, coverage: k, minimumMargin: margin)
    }

    public static func index(q: Double, sigma: Double, in families: [KinematicFamily], kMax: Double,
                             coverage k: Double = defaultCoverage,
                             minimumMargin margin: Double = minimumMargin) throws(Refusal) -> Assignment {
        try index(q: q, sigma: sigma, in: families.map(\.family), kMax: kMax, coverage: k, minimumMargin: margin)
    }

    /// The same over DiffractoCore's rings (any list of families).
    public static func index(q: Double, sigma: Double, in families: [ReflectionFamily], kMax: Double,
                             coverage k: Double = defaultCoverage,
                             minimumMargin margin: Double = minimumMargin) throws(Refusal) -> Assignment {
        guard q.isFinite, q > 0 else { throw .invalidInput("|q| must be a positive number of Å⁻¹") }
        guard sigma.isFinite, sigma > 0 else { throw .invalidInput("σ_q must be a positive number of Å⁻¹") }
        guard k.isFinite, k > 0 else { throw .invalidInput("the coverage factor must be a positive number") }
        guard kMax.isFinite, kMax > 0 else { throw .invalidInput("k_max must be a positive number of Å⁻¹") }
        guard margin.isFinite, margin >= 0 else { throw .invalidInput("the minimum margin must be a number of at least 0") }

        let candidates = families
            .map { Candidate(family: $0, z: abs(q - $0.g) / sigma, deviation: $0.g / q - 1) }
            .sorted { ($0.z, $0.family.g) < ($1.z, $1.family.g) }
        let within = candidates.prefix { $0.z <= k }
        if within.count >= 2 { throw .ambiguous(Array(within)) }

        let r = reach(coverage: k, minimumMargin: margin) * sigma
        guard let best = within.first else {
            let nearest = candidates.first
            let reachEnd = q + max(r, nearest.map { abs(q - $0.family.g) } ?? 0)
            if reachEnd > kMax { throw .beyondKMax(reach: reachEnd, kMax: kMax) }
            throw .noneWithinTolerance(nearest: nearest)
        }
        let next = candidates.count > 1 ? candidates[1] : nil
        let measured = next.map { ($0.z * $0.z - best.z * best.z) / 2 }
        if let next, let ratio = measured, ratio < margin { throw .ambiguous([best, next]) }
        let reachEnd = q + max(r, abs(q - best.family.g), next.map { abs(q - $0.family.g) } ?? 0)
        if reachEnd > kMax { throw .beyondKMax(reach: reachEnd, kMax: kMax) }
        return Assignment(q: q, sigma: sigma, coverage: k, best: best, next: next, margin: measured)
    }
}
