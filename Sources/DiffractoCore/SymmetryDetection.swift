//
//  SymmetryDetection.swift
//  Role: the operations of a structure from its sites, found by testing every element of the metric's holohedry
//        (and any rotation a file's loop names) on the sites (registration section 3.4b, amended 2026-10-14: always
//        run, also after a file's expansion). The permutation test is mac4DSTEM 7448401b's `CIFImport.isSymmetry`
//        (candidate translations from site 0 onto every compatible site, targets distinct, positions within a FIXED
//        5e-3 fractional, occupancies within 0.02), generalised from two generators of two families to the whole
//        holohedry of seven, and extended to keep every translation that works, so a centred cell yields its centring
//        translations and a screw axis its own.
//
//        Each translation is REFINED by least squares over the matched sites (the mean of the residuals), then taken
//        from the file's loop when the loop has the same operation, else moved to the nearest k/24 only when it lies
//        within three times its own uncertainty of it (Gate B amendment, 2026-10-15: a fixed 5e-3 snap moved genuine
//        off-grid translations of a shifted origin and lost or merged operations).
//

import Foundation
import simd

public enum SymmetryDetection {
    /// Fractional-coordinate agreement required to call two sites the same (the apps' reasoning: ≈10δ for 3-decimal
    /// data, an order of magnitude below real symmetry breaking; fixed on purpose so a file cannot widen its own gate).
    public static let positionTolerance = 5e-3
    /// Occupancies are dimensionless and refined ones differ by thousandths between equivalent sites.
    public static let occupancyTolerance = 0.02
    /// The metric family is proposed at the ring rule's tolerance (1e-9 relative on lengths, 1e-9° on angles; the
    /// registration's 3.2 default): a rotation that is a symmetry then preserves every |g| to the 1e-9 the families are
    /// grouped by, so an orbit never splits across two rings (Gate B amendment, 2026-10-15; the apps' 1e-3 / 0.05°
    /// called a cell with c/a = 1.0005 cubic and split its rings).
    public static let defaultMetricLengthTolerance = 1e-9
    public static let defaultMetricAngleToleranceDeg = 1e-9
    /// The conventional grid of space-group translations: halves, thirds, quarters, sixths, eighths.
    public static let translationGrid = 1.0 / 24

    public struct Result: Sendable, Equatable {
        /// The rotations that were tested: the holohedry of the metric plus any extra candidates.
        public let candidates: PointSymmetry
        /// Every (R, t) found, the identity first; translations mod 1.
        public let operations: [SymmetryOperation]
        /// Candidates that fail at `positionTolerance` but pass at 10 × the written half-step (reported, not taken).
        public let unresolvedRotations: [IntegerMatrix]
        /// The largest uncertainty of a translation that was neither exact from the file nor on the k/24 grid: its
        /// least-squares residual plus twice the coordinates' written half-step, at least 1e-9. Reported.
        public let translationUncertainty: Double
    }

    /// The operations of `sites` in `cell`. `coarsestHalfStep` is the coarsest written rounding half-step of the
    /// coordinates when known (0 when exact). `extraCandidates` are rotations to test beyond the metric's holohedry;
    /// `preferredOperations` (a file's loop) supply exact translations for the operations they name.
    public static func detect(
        sites: [AtomSite], cell: UnitCell, coarsestHalfStep: Double = 0,
        extraCandidates: [IntegerMatrix] = [], preferredOperations: [SymmetryOperation] = [],
        metricLengthTolerance: Double = defaultMetricLengthTolerance,
        metricAngleToleranceDeg: Double = defaultMetricAngleToleranceDeg
    ) -> Result {
        let family = cell.metricFamily(lengthTolerance: metricLengthTolerance, angleToleranceDeg: metricAngleToleranceDeg)
        let holohedry = PointSymmetry.holohedry(of: family)
        let candidates = extraCandidates.isEmpty
            ? holohedry : PointSymmetry(rotations: holohedry.rotations + extraCandidates)
        guard !sites.isEmpty else {
            return Result(candidates: candidates, operations: [.identity], unresolvedRotations: [],
                          translationUncertainty: 1e-9)
        }
        var operations: [SymmetryOperation] = []
        var unresolved: [IntegerMatrix] = []
        var uncertainty = 1e-9
        let resolutionLimited = 10 * coarsestHalfStep
        for rotation in candidates.rotations {
            let found = symmetryTranslations(of: rotation, sites: sites, positionTolerance: positionTolerance,
                                             coarsestHalfStep: coarsestHalfStep,
                                             preferred: preferredOperations.filter { $0.rotation == rotation })
            if found.isEmpty {
                if resolutionLimited > positionTolerance,
                   !symmetryTranslations(of: rotation, sites: sites, positionTolerance: resolutionLimited,
                                         coarsestHalfStep: coarsestHalfStep, preferred: []).isEmpty {
                    unresolved.append(rotation)
                }
                continue
            }
            for (t, u) in found {
                operations.append(SymmetryOperation(rotation: rotation, translation: t))
                uncertainty = max(uncertainty, u)
            }
        }
        operations.sort { lhs, rhs in
            if lhs.rotation.isIdentity != rhs.rotation.isIdentity { return lhs.rotation.isIdentity }
            let l = lhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }, r = rhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }
            if l != r { return l.lexicographicallyPrecedes(r) }
            return (lhs.translation.x, lhs.translation.y, lhs.translation.z)
                < (rhs.translation.x, rhs.translation.y, rhs.translation.z)
        }
        return Result(candidates: candidates, operations: closedSubset(operations), unresolvedRotations: unresolved,
                      translationUncertainty: uncertainty)
    }

    /// The operations a space group can be made of: the set is closed under composition (mod lattice translations,
    /// within `positionTolerance`). Pseudo-symmetry within the fractional tolerance can pass single operations whose
    /// products do not pass (the apps' split-atom gold: 40 operations, rotations of m-3m, not a group); an operation
    /// that composes out of the set is dropped and the test repeated until the set is closed (2026-10-15). Sets above
    /// `closureCheckLimit` operations are returned as found (a 500-atom symmetric supercell has 24000).
    public static let closureCheckLimit = 2048

    static func closedSubset(_ operations: [SymmetryOperation]) -> [SymmetryOperation] {
        guard operations.count <= closureCheckLimit else { return operations }
        var set = operations
        while true {
            var byRotation: [IntegerMatrix: [SIMD3<Double>]] = [:]
            for op in set { byRotation[op.rotation, default: []].append(op.translation) }
            func contains(_ op: SymmetryOperation) -> Bool {
                byRotation[op.rotation]?.contains { translationsEqual($0, op.translation, positionTolerance) } ?? false
            }
            let kept = set.filter { g in g.rotation.isIdentity && SymmetryOperation.isWholeNumber(g.translation)
                || set.allSatisfy { h in contains(g.composed(with: h)) } }
            if kept.count == set.count { return set }
            set = kept
        }
    }

    /// Every translation t (mod 1) for which p ↦ R p + t is a permutation of the site set preserving Z and occupancy,
    /// each with its uncertainty (0 when exact: from the file, or on the grid). Candidates come from mapping site 0
    /// onto every compatible site, which is exhaustive: if the operation is a symmetry then R x₀ + t lands on some site.
    static func symmetryTranslations(of rotation: IntegerMatrix, sites: [AtomSite], positionTolerance: Double,
                                     coarsestHalfStep: Double,
                                     preferred: [SymmetryOperation]) -> [(SIMD3<Double>, Double)] {
        let reference = sites[0]
        let rotatedReference = rotation.apply(reference.fractional)
        var found: [(SIMD3<Double>, Double)] = []
        for candidate in sites where candidate.z == reference.z
            && abs(candidate.occupancy - reference.occupancy) <= occupancyTolerance {
            let trial = wrapFractional(candidate.fractional - rotatedReference)
            guard let assignment = permutation(rotation, trial, of: sites, positionTolerance: positionTolerance)
            else { continue }
            // Least squares: the translation that best maps every site onto its target is the trial plus the mean
            // of the matched residuals; the spread of the residuals is the fit's own uncertainty.
            var residuals: [SIMD3<Double>] = []
            for (i, j) in assignment.enumerated() {
                let image = rotation.apply(sites[i].fractional) + trial
                residuals.append(minimumImage(sites[j].fractional - image))
            }
            let mean = residuals.reduce(SIMD3<Double>(0, 0, 0), +) / Double(residuals.count)
            let spread = residuals.map { abs($0 - mean).max() }.max() ?? 0
            let u = spread + 2 * coarsestHalfStep + 1e-12
            let refined = wrapFractional(trial + mean)
            // The candidates for the final translation, best first: the file's exact one, the k/24 grid point, the
            // least-squares one, the trial. Each is RE-VERIFIED as a permutation before it is taken: moving a
            // translation after the test is how a non-symmetry got in (found 2026-10-15 while fixing Gate B's F1).
            var options: [(SIMD3<Double>, Double)] = []
            if let exact = preferred.first(where: { SymmetryDetection.translationsEqual($0.translation, refined, 3 * u + 1e-9) }) {
                options.append((exact.translationMod1, 0))
            }
            let grid = snapped(refined, within: 3 * u)
            if grid.onGrid { options.append((grid.value, 0)) }
            options.append((refined, u))
            options.append((trial, max(u, abs(mean).max())))
            guard let (t, tu) = options.first(where: {
                permutation(rotation, $0.0, of: sites, positionTolerance: positionTolerance) != nil
            }) else { continue }
            guard !found.contains(where: { translationsEqual($0.0, t, max(1e-6, 2 * u)) }) else { continue }
            found.append((t, tu))
        }
        return found
    }

    /// Each component moved to the nearest multiple of 1/24 when ALL three are within `within` of it.
    static func snapped(_ t: SIMD3<Double>, within: Double) -> (value: SIMD3<Double>, onGrid: Bool) {
        let grid = SIMD3((t.x / translationGrid).rounded(), (t.y / translationGrid).rounded(),
                         (t.z / translationGrid).rounded()) * translationGrid
        let onGrid = abs(t - grid).max() <= within
        return (onGrid ? wrapFractional(grid) : t, onGrid)
    }

    /// Each component in [−½, ½).
    static func minimumImage(_ d: SIMD3<Double>) -> SIMD3<Double> {
        d - (d + 0.5).rounded(.down)
    }

    /// The target of each site when p ↦ R p + t permutes the set, else nil. Surjectivity is checked, not argued:
    /// finite-tolerance matching is not an equivalence relation, so two distinct sites can both land within tolerance
    /// of the same target (the apps' split-atom gold is at best tetragonal).
    static func permutation(_ rotation: IntegerMatrix, _ translation: SIMD3<Double>, of sites: [AtomSite],
                            positionTolerance: Double) -> [Int]? {
        var claimed = Set<Int>()
        var assignment: [Int] = []
        assignment.reserveCapacity(sites.count)
        for site in sites {
            let image = wrapFractional(rotation.apply(site.fractional) + translation)
            let target = sites.indices.first { index in
                !claimed.contains(index)
                    && sites[index].z == site.z
                    && abs(sites[index].occupancy - site.occupancy) <= occupancyTolerance
                    && periodicDistance(sites[index].fractional, image) < positionTolerance
            }
            guard let target else { return nil }
            claimed.insert(target)
            assignment.append(target)
        }
        return assignment
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
