//
//  SymmetryDetection.swift
//  Role: the operations of a structure given as a site list, found by testing every element of the metric's
//        holohedry on the sites (registration section 3.4b). The permutation test is mac4DSTEM 7448401b's
//        `CIFImport.isSymmetry` (candidate translations from site 0 onto every compatible site, targets distinct,
//        positions within a FIXED 5e-3 fractional, occupancies within 0.02), generalised from two generators of two
//        families to the whole holohedry of seven, and extended to keep every translation that works, so a centred
//        cell yields its centring translations and a screw axis its own.
//

import Foundation
import simd

public enum SymmetryDetection {
    /// Fractional-coordinate agreement required to call two sites the same (the apps' reasoning: ≈10δ for 3-decimal
    /// data, an order of magnitude below real symmetry breaking; fixed on purpose so a file cannot widen its own gate).
    public static let positionTolerance = 5e-3
    /// Occupancies are dimensionless and refined ones differ by thousandths between equivalent sites.
    public static let occupancyTolerance = 0.02
    /// The metric family is proposed at the apps' CIF tolerances (lengths 1e-3 relative, angles 0.05°): the atoms,
    /// not the metric, decide what is taken.
    public static let defaultMetricLengthTolerance = 1e-3
    public static let defaultMetricAngleToleranceDeg = 0.05

    public struct Result: Sendable, Equatable {
        /// The holohedry that was tested.
        public let candidates: PointSymmetry
        /// Every (R, t) found, the identity first; translations mod 1.
        public let operations: [SymmetryOperation]
        /// Candidates that fail at `positionTolerance` but pass at 10 × the written half-step (reported, not taken).
        public let unresolvedRotations: [IntegerMatrix]
    }

    /// The operations of `sites` in `cell`. `coarsestHalfStep` is the coarsest written rounding half-step of the
    /// coordinates when known (0 when exact); it only decides whether a failed candidate is reported as unresolved.
    public static func detect(
        sites: [AtomSite], cell: UnitCell, coarsestHalfStep: Double = 0,
        metricLengthTolerance: Double = defaultMetricLengthTolerance,
        metricAngleToleranceDeg: Double = defaultMetricAngleToleranceDeg
    ) -> Result {
        let family = cell.metricFamily(lengthTolerance: metricLengthTolerance, angleToleranceDeg: metricAngleToleranceDeg)
        let candidates = PointSymmetry.holohedry(of: family)
        guard !sites.isEmpty else {
            return Result(candidates: candidates, operations: [.identity], unresolvedRotations: [])
        }
        var operations: [SymmetryOperation] = []
        var unresolved: [IntegerMatrix] = []
        let resolutionLimited = 10 * coarsestHalfStep
        for rotation in candidates.rotations {
            let translations = symmetryTranslations(of: rotation, sites: sites, positionTolerance: positionTolerance)
            if translations.isEmpty {
                if resolutionLimited > positionTolerance,
                   !symmetryTranslations(of: rotation, sites: sites, positionTolerance: resolutionLimited).isEmpty {
                    unresolved.append(rotation)
                }
                continue
            }
            for t in translations { operations.append(SymmetryOperation(rotation: rotation, translation: t)) }
        }
        // The identity first, then by rotation order as the group lists them; translations as found.
        operations.sort { lhs, rhs in
            if lhs.rotation.isIdentity != rhs.rotation.isIdentity { return lhs.rotation.isIdentity }
            let l = lhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }, r = rhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }
            if l != r { return l.lexicographicallyPrecedes(r) }
            return (lhs.translation.x, lhs.translation.y, lhs.translation.z)
                < (rhs.translation.x, rhs.translation.y, rhs.translation.z)
        }
        return Result(candidates: candidates, operations: operations, unresolvedRotations: unresolved)
    }

    /// Every translation t (mod 1, deduplicated within 1e-6) for which p ↦ R p + t is a permutation of the site set
    /// preserving Z and occupancy. Candidates come from mapping site 0 onto every compatible site, which is exhaustive:
    /// if the operation is a symmetry then R x₀ + t lands on some site.
    static func symmetryTranslations(of rotation: IntegerMatrix, sites: [AtomSite],
                                     positionTolerance: Double) -> [SIMD3<Double>] {
        let reference = sites[0]
        let rotatedReference = rotation.apply(reference.fractional)
        var found: [SIMD3<Double>] = []
        for candidate in sites where candidate.z == reference.z
            && abs(candidate.occupancy - reference.occupancy) <= occupancyTolerance {
            let translation = wrapFractional(candidate.fractional - rotatedReference)
            guard !found.contains(where: { translationsEqual($0, translation, 1e-6) }) else { continue }
            if isPermutation(rotation, translation, of: sites, positionTolerance: positionTolerance) {
                found.append(translation)
            }
        }
        return found
    }

    /// Surjectivity is checked, not argued: finite-tolerance matching is not an equivalence relation, so two distinct
    /// sites can both land within tolerance of the same target (an FCC cell with a split atom is at best tetragonal).
    static func isPermutation(_ rotation: IntegerMatrix, _ translation: SIMD3<Double>, of sites: [AtomSite],
                              positionTolerance: Double) -> Bool {
        var claimed = Set<Int>()
        for site in sites {
            let image = wrapFractional(rotation.apply(site.fractional) + translation)
            let target = sites.indices.first { index in
                !claimed.contains(index)
                    && sites[index].z == site.z
                    && abs(sites[index].occupancy - site.occupancy) <= occupancyTolerance
                    && periodicDistance(sites[index].fractional, image) < positionTolerance
            }
            guard let target else { return false }
            claimed.insert(target)
        }
        return true
    }

    static func periodicDistance(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double {
        let dx = fractionalAxisDelta(a.x, b.x), dy = fractionalAxisDelta(a.y, b.y), dz = fractionalAxisDelta(a.z, b.z)
        return (dx * dx + dy * dy + dz * dz).squareRoot()
    }

    static func translationsEqual(_ a: SIMD3<Double>, _ b: SIMD3<Double>, _ tolerance: Double) -> Bool {
        fractionalAxisDelta(a.x, b.x) <= tolerance && fractionalAxisDelta(a.y, b.y) <= tolerance
            && fractionalAxisDelta(a.z, b.z) <= tolerance
    }
}
