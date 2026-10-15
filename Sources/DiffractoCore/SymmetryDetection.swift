//
//  SymmetryDetection.swift
//  Role: the operations of a structure from its sites, found by testing every element of the metric's holohedry
//        (and any rotation a file's loop names) on the sites (registration section 3.4b, amended 2026-10-14: always
//        run, also after a file's expansion). The permutation test is reference commit 7448401b's `CIFImport.isSymmetry`
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
    /// Fractional-coordinate agreement required to call two sites the same (the reference implementations' reasoning: ≈10δ for 3-decimal
    /// data, an order of magnitude below real symmetry breaking; fixed on purpose so a file cannot widen its own gate).
    public static let positionTolerance = 5e-3
    /// Occupancies are dimensionless and refined ones differ by thousandths between equivalent sites.
    public static let occupancyTolerance = 0.02
    /// Through `CrystalStructure.detectingSymmetry(conventional: true)` (the default) only the LENGTH tolerance reaches the
    /// lattice-symmetry test that finds the cell and its Bravais lattice; the angle tolerance has no effect there (Gate B
    /// on 4876b18, F3), and spglib's symprec corresponds to a length tolerance of about 1e-5, not 1e-9: a bcc primitive
    /// cell written with 109.4712° is hR at the default and Im-3m from 1e-5.
    /// The metric family is proposed at the ring rule's tolerance (1e-9 relative on lengths, 1e-9° on angles; the
    /// registration's 3.2 default): a rotation that is a symmetry then preserves every |g| to the 1e-9 the families are
    /// grouped by, so an orbit never splits across two rings (Gate B amendment, 2026-10-15; the reference implementations' 1e-3 / 0.05°
    /// called a cell with c/a = 1.0005 cubic and split its rings).
    public static let defaultMetricLengthTolerance = 1e-9
    public static let defaultMetricAngleToleranceDeg = 1e-9
    /// The conventional grid of space-group translations: halves, thirds, quarters, sixths, eighths.
    public static let translationGrid = 1.0 / 24

    /// How two positions are called the same site (slice 1b). `.fractional(5e-3)` is slice 1's rule and the default:
    /// the Euclidean distance in fractional coordinates, whatever the cell's size. `.cartesian(Å)` is spglib's sense:
    /// the minimum-image distance in Å, so a pseudo-symmetry of 0.013 Å reads as symmetry in a 50 Å cell at the
    /// fractional default and not at `.cartesian(1e-3)`. The default stays fractional until a measured comparison is
    /// registered.
    ///
    /// The Cartesian value is spglib's `symprec` in the sense that it is a distance in Å, but it is NOT the same
    /// threshold: measured against spglib (Gate B on 4876b18, F4), an operation is kept here only while the displacement
    /// that breaks it is below about HALF of `t`, and is lost at displacements of 0.5 to 0.9 `t`, where spglib with
    /// `symprec = t` still finds it up to 1.0 to 1.1 `t`. The set of operations is also closed under composition at this
    /// bound (`closedSubset`). Roughly, `.cartesian(2 * s)` here corresponds to spglib's symprec `s`, and the
    /// proposed `1e-3` to a symprec of about 5e-4 Å. The default is unchanged.
    public enum Tolerance: Sendable, Equatable {
        case fractional(Double)
        case cartesian(Double)

        public static let `default`: Tolerance = .fractional(SymmetryDetection.positionTolerance)
        /// 10 × spglib's default symprec of 1e-4 Å: the value proposed for a Cartesian default (its effective bound on a
        /// displacement is about half of it, see `Tolerance`).
        public static let proposedCartesian: Tolerance = .cartesian(1e-3)
    }

    /// How the permutation test finds the target of an image: a cell list (bins of the position tolerance, O(1) per
    /// image) or slice 1's scan of every site (the reference path, O(N)). They return identical operations bit for bit.
    public enum Search: Sendable, Equatable {
        case cellList
        case bruteForce
    }

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
        metricAngleToleranceDeg: Double = defaultMetricAngleToleranceDeg,
        tolerance: Tolerance = .default, search: Search = .cellList
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
        let matcher = PositionMatcher(cell: cell, tolerance: tolerance)
        let index = search == .cellList ? SiteIndex(sites: sites, window: matcher.window) : nil
        for rotation in candidates.rotations {
            let found = symmetryTranslations(of: rotation, sites: sites, matcher: matcher, index: index,
                                             coarsestHalfStep: coarsestHalfStep,
                                             preferred: preferredOperations.filter { $0.rotation == rotation })
            if found.isEmpty {
                if resolutionLimited > matcher.fractionalBound {
                    let loose = PositionMatcher(cell: cell, tolerance: .fractional(resolutionLimited))
                    if !symmetryTranslations(of: rotation, sites: sites, matcher: loose, index: nil,
                                             coarsestHalfStep: coarsestHalfStep, preferred: []).isEmpty {
                        unresolved.append(rotation)
                    }
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
        return Result(candidates: candidates, operations: closedSubset(operations, tolerance: matcher.fractionalBound),
                      unresolvedRotations: unresolved,
                      translationUncertainty: uncertainty)
    }

    /// The operations a space group can be made of: the set is closed under composition (mod lattice translations,
    /// within `positionTolerance`). Pseudo-symmetry within the fractional tolerance can pass single operations whose
    /// products do not pass (the reference implementations' split-atom gold: 40 operations, rotations of m-3m, not a group); an operation
    /// that composes out of the set is dropped and the test repeated until the set is closed (2026-10-15). Sets above
    /// `closureCheckLimit` operations are returned as found (a 500-atom symmetric supercell has 24000).
    public static let closureCheckLimit = 2048

    static func closedSubset(_ operations: [SymmetryOperation], tolerance: Double = positionTolerance) -> [SymmetryOperation] {
        guard operations.count <= closureCheckLimit else { return operations }
        var set = operations
        while true {
            var byRotation: [IntegerMatrix: [SIMD3<Double>]] = [:]
            for op in set { byRotation[op.rotation, default: []].append(op.translation) }
            func contains(_ op: SymmetryOperation) -> Bool {
                byRotation[op.rotation]?.contains { translationsEqual($0, op.translation, tolerance) } ?? false
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
    static func symmetryTranslations(of rotation: IntegerMatrix, sites: [AtomSite], matcher: PositionMatcher,
                                     index: SiteIndex?, coarsestHalfStep: Double,
                                     preferred: [SymmetryOperation]) -> [(SIMD3<Double>, Double)] {
        let reference = sites[0]
        let rotatedReference = rotation.apply(reference.fractional)
        var found: [(SIMD3<Double>, Double)] = []
        for candidate in sites where candidate.z == reference.z
            && abs(candidate.occupancy - reference.occupancy) <= occupancyTolerance {
            let trial = wrapFractional(candidate.fractional - rotatedReference)
            guard let assignment = permutation(rotation, trial, of: sites, matcher: matcher, index: index)
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
                permutation(rotation, $0.0, of: sites, matcher: matcher, index: index) != nil
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
    /// of the same target (the reference implementations' split-atom gold is at best tetragonal).
    static func permutation(_ rotation: IntegerMatrix, _ translation: SIMD3<Double>, of sites: [AtomSite],
                            matcher: PositionMatcher, index: SiteIndex?) -> [Int]? {
        guard let index else { return permutationReference(rotation, translation, of: sites, matcher: matcher) }
        var claimed = [Bool](repeating: false, count: sites.count)
        var assignment: [Int] = []
        assignment.reserveCapacity(sites.count)
        for site in sites {
            let image = wrapFractional(rotation.apply(site.fractional) + translation)
            // The same target slice 1 picks: the LOWEST index among the unclaimed sites that match.
            var target = Int.max
            index.forEachCandidate(near: image) { j in
                if j < target, !claimed[j], sites[j].z == site.z,
                   abs(sites[j].occupancy - site.occupancy) <= occupancyTolerance,
                   matcher.matches(sites[j].fractional, image) { target = j }
            }
            guard target != Int.max else { return nil }
            claimed[target] = true
            assignment.append(target)
        }
        return assignment
    }

    /// Slice 1's scan of every site for every image (the CPU reference path of the cell list).
    static func permutationReference(_ rotation: IntegerMatrix, _ translation: SIMD3<Double>, of sites: [AtomSite],
                                     matcher: PositionMatcher) -> [Int]? {
        var claimed = Set<Int>()
        var assignment: [Int] = []
        assignment.reserveCapacity(sites.count)
        for site in sites {
            let image = wrapFractional(rotation.apply(site.fractional) + translation)
            let target = sites.indices.first { index in
                !claimed.contains(index)
                    && sites[index].z == site.z
                    && abs(sites[index].occupancy - site.occupancy) <= occupancyTolerance
                    && matcher.matches(sites[index].fractional, image)
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


// MARK: - Position matching and the cell list (slice 1b)

/// The rule that calls two fractional positions the same site.
struct PositionMatcher: Sendable {
    let tolerance: SymmetryDetection.Tolerance
    let lattice: [SIMD3<Double>]
    /// The largest fractional offset along each axis a match can have (the tolerance itself, or the Cartesian one
    /// times |bᵢ|), at most ½: the half-width of the cell list's search window.
    let window: SIMD3<Double>
    /// The largest component of `window`: the fractional tolerance that closure and equality tests use.
    let fractionalBound: Double
    private let cartesianSquared: Double
    /// Whether a single wrapped difference can be taken as the minimum image (the tolerance is far below half the
    /// shortest lattice vector).
    private let singleImage: Bool

    init(cell: UnitCell, tolerance: SymmetryDetection.Tolerance) {
        self.tolerance = tolerance
        self.lattice = cell.realLattice
        switch tolerance {
        case .fractional(let t):
            window = SIMD3(repeating: min(0.5, t))
            cartesianSquared = 0
            singleImage = true
        case .cartesian(let t):
            let b = cell.reciprocalLattice.map { length($0) }
            window = SIMD3(min(0.5, t * b[0]), min(0.5, t * b[1]), min(0.5, t * b[2]))
            cartesianSquared = t * t
            let shortest = cell.realLattice.map { length($0) }.min() ?? 1
            singleImage = t < 0.25 * shortest
        }
        fractionalBound = max(window.x, window.y, window.z)
    }

    /// Whether `a` and `b` (fractional, any image) are the same position within the tolerance. The argument order is
    /// slice 1's: a site, then the image looked for.
    func matches(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Bool {
        switch tolerance {
        case .fractional(let t):
            return SymmetryDetection.periodicDistance(a, b) < t
        case .cartesian:
            var d = SIMD3(a.x - b.x, a.y - b.y, a.z - b.z)
            d = d - d.rounded(.toNearestOrAwayFromZero)
            if singleImage {
                if abs(d.x) > window.x || abs(d.y) > window.y || abs(d.z) > window.z { return false }
                return cartesianLengthSquared(d) < cartesianSquared
            }
            var best = Double.infinity
            for sx in -1...1 { for sy in -1...1 { for sz in -1...1 {
                best = min(best, cartesianLengthSquared(d + SIMD3(Double(sx), Double(sy), Double(sz))))
            } } }
            return best < cartesianSquared
        }
    }

    private func cartesianLengthSquared(_ f: SIMD3<Double>) -> Double {
        let v = f.x * lattice[0] + f.y * lattice[1] + f.z * lattice[2]
        return (v * v).sum()
    }
}

/// A cell list over the sites: bins at least twice the search window wide, so a lookup touches at most two bins per
/// axis; the candidates of a bin are in increasing site index.
struct SiteIndex: Sendable {
    let bins: SIMD3<Int>
    let window: SIMD3<Double>
    private let start: [Int32]
    private let items: [Int32]

    init(sites: [AtomSite], window: SIMD3<Double>) {
        self.window = window
        // about one site per bin, never a bin narrower than twice the window
        let wanted = max(1, Int(cbrt(Double(sites.count)).rounded(.up)))
        var nb = SIMD3<Int>(repeating: 1)
        for axis in 0..<3 {
            let limit = window[axis] > 0 ? Int((1 / (2 * window[axis])).rounded(.down)) : 1_000_000
            nb[axis] = max(1, min(limit, wanted))
        }
        bins = nb
        let total = nb.x * nb.y * nb.z
        var counts = [Int32](repeating: 0, count: total + 1)
        var binOf = [Int](repeating: 0, count: sites.count)
        for (i, s) in sites.enumerated() {
            let b = SiteIndex.binCoordinates(wrapFractional(s.fractional), nb)
            let flat = (b.x * nb.y + b.y) * nb.z + b.z
            binOf[i] = flat
            counts[flat + 1] += 1
        }
        for k in 0..<total { counts[k + 1] += counts[k] }
        var cursor = counts
        var out = [Int32](repeating: 0, count: sites.count)
        for i in 0..<sites.count {                   // increasing index within every bin
            out[Int(cursor[binOf[i]])] = Int32(i)
            cursor[binOf[i]] += 1
        }
        start = counts
        items = out
    }

    private static func binCoordinates(_ p: SIMD3<Double>, _ nb: SIMD3<Int>) -> SIMD3<Int> {
        SIMD3(min(nb.x - 1, Int(p.x * Double(nb.x))), min(nb.y - 1, Int(p.y * Double(nb.y))),
              min(nb.z - 1, Int(p.z * Double(nb.z))))
    }

    /// Calls `body` with every site whose bin the window around `p` (wrapped) touches.
    func forEachCandidate(near p: SIMD3<Double>, _ body: (Int) -> Void) {
        func range(_ axis: Int) -> (lo: Int, count: Int) {
            let n = bins[axis]
            let lo = Int(((p[axis] - window[axis]) * Double(n)).rounded(.down))
            let hi = Int(((p[axis] + window[axis]) * Double(n)).rounded(.down))
            return hi - lo + 1 >= n ? (0, n) : (lo, hi - lo + 1)
        }
        let rx = range(0), ry = range(1), rz = range(2)
        for ix in 0..<rx.count {
            let bx = ((rx.lo + ix) % bins.x + bins.x) % bins.x
            for iy in 0..<ry.count {
                let by = ((ry.lo + iy) % bins.y + bins.y) % bins.y
                for iz in 0..<rz.count {
                    let bz = ((rz.lo + iz) % bins.z + bins.z) % bins.z
                    let flat = (bx * bins.y + by) * bins.z + bz
                    var k = Int(start[flat])
                    let end = Int(start[flat + 1])
                    while k < end { body(Int(items[k])); k += 1 }
                }
            }
        }
    }
}
