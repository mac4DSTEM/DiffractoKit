//
//  PhaseVectorMatching.swift
//  Role: the MATCHER half of vector-matched phase mapping: for each scan position, remove the matrix's reflections, score
//        what survives against every candidate entry of a `PhaseReferenceLibrary`, and return a phase label, a distance, and
//        an explicit "not indexed".
//
//  Method: Thronsen et al., Ultramicroscopy 255 (2024) 113861, CC BY 4.0 (attribution in NOTICE). The decision rule of
//  `.knownVariants` and its two constants (direct-matrix count 1, residual cutoff 0.07 Å⁻¹) are as the paper describes them;
//  no code was copied from the authors' notebook (their repository carries no licence). Three deviations from what the paper describes are deliberate and are named where
//  they occur: the matrix mask is applied in VECTOR space rather than to the image, the matrix orientation is fitted ONCE
//  for the scan rather than per pattern, and the verdict carries a completeness requirement the paper's mean-distance
//  score does not have.
//
//  WHY THIS EXISTS AT ALL. On 2026-09-11 a template-matched route was refuted on a measurement: a bare best-score argmax
//  over phases returned GOLD for a pattern containing only ALUMINIUM, because one radial bin is 0.05 Å⁻¹ and Al-Au (111)
//  differ by 6 % of one bin. Vector matching works on the peak positions themselves, at the precision the sub-pixel
//  refinement already delivers, so the discriminating quantity is a distance in Å⁻¹ that a reader can compare against the
//  spacing they care about, not a correlation score that hides its own resolution.
//
//  WHAT IT CANNOT DO, stated here rather than discovered later: it is zone-axis only; it cannot separate phases
//  overlapping along the beam; and strain shifts Bragg positions, which Thronsen et al. name as the main cause of their
//  interface errors. One matrix orientation is fitted for the whole scan, so a second matrix grain would be mis-removed:
//  that is a limit of this implementation, not of the method.
//
//  NOT SCORED AGAINST THE PUBLISHED GROUND TRUTH. Scoring this against the published ground truth of Thronsen et al. needs
//  a ~7.4 GB download and has not been run; the verdict rules have no reference implementation (py4DSTEM's phase tool is
//  a different method), so every number this file produces is an output of the definitions written in this file, not a
//  measurement against an independent oracle.
//
//  THE INPUT. A scan is `[[SIMD2<Double>]]`, one array of peak positions per scan position, row-major (`scanWidth` ×
//  `scanHeight`), in one of two forms (`PeakCoordinates`): detector pixels with an origin and a scale (converted with a
//  `Float` subtraction first, so peaks held as `Float` positions give bit-identical vectors), or calibrated Å⁻¹ vectors
//  about the origin.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

// MARK: - The matcher

public enum PhaseVectorMatcher {

    /// One pattern's peaks as calibrated in-plane vectors (Å⁻¹), with the
    /// direct beam removed.
    ///
    /// The frame is `dx = x − originX`, `dy = y − originY`, times
    /// `invAngstromPerPixel`, the frame of template matching. Ellipse
    /// correction and per-position origin collapse are the caller's job: they
    /// belong to calibration, not to matching. The R-Q rotation is deliberately
    /// NOT applied: the in-plane rotation is an axis of this search, fitted from
    /// the data, so applying a stored rotation as well would double-count it.
    ///
    /// Detector-pixel form: `peaks` are pixel positions (an application holding
    /// `Float` positions passes `SIMD2<Double>(Double(x), Double(y))` and gets
    /// bit-identical vectors); the subtraction is done in `Float`, then `Double`,
    /// then the scale, as in the apps' original.
    public static func experimentalVectors(peaks: [SIMD2<Double>],
                                            originX: Float, originY: Float,
                                            invAngstromPerPixel: Double,
                                            directBeamRadiusInvAngstrom: Double,
                                            maximumVectorInvAngstrom: Double = 0)
        -> [SIMD2<Double>] {
        experimentalVectors(
            peaks: peaks,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
            maximumVectorInvAngstrom: maximumVectorInvAngstrom)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func experimentalVectors(peaks: [SIMD2<Double>],
                                            coordinates: PeakCoordinates,
                                            directBeamRadiusInvAngstrom: Double,
                                            maximumVectorInvAngstrom: Double = 0)
        -> [SIMD2<Double>] {
        var out: [SIMD2<Double>] = []
        out.reserveCapacity(peaks.count)
        for p in peaks {
            if let q = experimentalVector(
                peak: p, coordinates: coordinates,
                directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: maximumVectorInvAngstrom) {
                out.append(q)
            }
        }
        return out
    }

    /// One peak as a calibrated vector, or nil when the matcher leaves it
    /// out (the direct beam, past the maximum vector, not finite). The single
    /// definition of "usable", shared by `experimentalVectors` and
    /// `claims(...)`, so a caller's idea of which disks were matched cannot
    /// drift from the map's.
    public static func experimentalVector(peak p: SIMD2<Double>,
                                           originX: Float, originY: Float,
                                           invAngstromPerPixel: Double,
                                           directBeamRadiusInvAngstrom: Double,
                                           maximumVectorInvAngstrom: Double = 0)
        -> SIMD2<Double>? {
        experimentalVector(
            peak: p,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
            maximumVectorInvAngstrom: maximumVectorInvAngstrom)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func experimentalVector(peak p: SIMD2<Double>,
                                           coordinates: PeakCoordinates,
                                           directBeamRadiusInvAngstrom: Double,
                                           maximumVectorInvAngstrom: Double = 0)
        -> SIMD2<Double>? {
        let q = coordinates.vector(of: p)
        let len = simd_length(q)
        guard len.isFinite, len > directBeamRadiusInvAngstrom,
              maximumVectorInvAngstrom <= 0 || len < maximumVectorInvAngstrom else { return nil }
        return q
    }

    // MARK: Nearest reference, with a length-band prune

    /// Index of and distance to the nearest vector of `refs` within `radius`.
    ///
    /// `refs` must be sorted by `length`, which `PhaseReferenceLibrary`
    /// guarantees: the reverse triangle inequality gives
    /// `| |u| − |v| | ≤ |u − v|`, so only the band `|u| ± radius` can hold a
    /// match and the rest is skipped without a distance computation.
    public static func nearest(_ u: SIMD2<Double>, in refs: [ReferenceVector],
                                radius: Double) -> (index: Int, distance: Double)? {
        guard !refs.isEmpty else { return nil }
        let uLen = simd_length(u)
        var lo = 0, hi = refs.count
        let lowerLength = uLen - radius
        while lo < hi {
            let mid = (lo + hi) / 2
            if refs[mid].length < lowerLength { lo = mid + 1 } else { hi = mid }
        }
        let upperLength = uLen + radius
        var best = -1
        var bestD = Double.infinity
        var i = lo
        while i < refs.count, refs[i].length <= upperLength {
            let d = simd_distance(u, refs[i].q)
            if d < bestD { bestD = d; best = i }
            i += 1
        }
        return (best >= 0 && bestD <= radius) ? (best, bestD) : nil
    }

    /// Per-worker scratch: one slot per reference vector, so scoring an entry
    /// allocates nothing.
    public final class Scratch: @unchecked Sendable {
        var bestPerRef: [Double]
        var touched: [Int]
        /// Room for `capacity` reference vectors (it grows on demand).
        public init(capacity: Int) {
            bestPerRef = [Double](repeating: .infinity, count: max(1, capacity))
            touched = []
            touched.reserveCapacity(max(1, capacity))
        }
    }

    /// Score `vectors` against one entry.
    ///
    /// For each experimental vector keep only its closest reference vector
    /// (their per-pattern reference subset, so a thin precipitate showing two
    /// spots is not penalised against a phase with forty), then average over
    /// the UNIQUE references hit — a reference claimed twice contributes its
    /// better distance once. Returns nil when nothing matched.
    public static func score(vectors: [SIMD2<Double>],
                              against entry: PhaseOrientationReference,
                              pairRadius: Double,
                              scratch: Scratch)
        -> (score: Double, matched: Int, uniqueReferences: Int)? {
        if scratch.bestPerRef.count < entry.vectors.count {
            scratch.bestPerRef = [Double](repeating: .infinity, count: entry.vectors.count)
        }
        scratch.touched.removeAll(keepingCapacity: true)
        var matched = 0
        for u in vectors {
            guard let hit = nearest(u, in: entry.vectors, radius: pairRadius) else { continue }
            matched += 1
            if scratch.bestPerRef[hit.index].isInfinite { scratch.touched.append(hit.index) }
            if hit.distance < scratch.bestPerRef[hit.index] {
                scratch.bestPerRef[hit.index] = hit.distance
            }
        }
        defer { for i in scratch.touched { scratch.bestPerRef[i] = .infinity } }
        guard !scratch.touched.isEmpty else { return nil }
        var sum = 0.0
        for i in scratch.touched { sum += scratch.bestPerRef[i] }
        return (sum / Double(scratch.touched.count), matched, scratch.touched.count)
    }

    /// Do `vectors` hold a Friedel pair — some u and v with |u + v| within
    /// `radius`? A ±g pair is one reciprocal-lattice row, the smallest thing
    /// that is evidence of a lattice rather than of two coincidences.
    public static func containsFriedelPair(_ vectors: [SIMD2<Double>], radius: Double) -> Bool {
        guard vectors.count >= 2 else { return false }
        for i in 0..<(vectors.count - 1) {
            for j in (i + 1)..<vectors.count
            where simd_length(vectors[i] + vectors[j]) <= radius && simd_length(vectors[i]) > radius {
                return true
            }
        }
        return false
    }

    // MARK: The whole scan

    /// The phase map of a scan: one `PhaseVectorResult` per position, row-major, or nil when the run was refused or
    /// cancelled (never a partial map).
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`), one array of peak positions per position; a malformed scan
    /// (width or height not positive, or a count that is not width · height) is refused (nil), as is a
    /// `matrixEntryIndex` outside the library, and a library whose candidate set is empty after the orientation
    /// relationships. `maximumWorkers` bounds the threads (nil: the active processor count); each position's result
    /// depends only on its own vectors, so the worker count moves no number. `progress` receives 0...1 from several
    /// threads; `isCancelled` is polled between positions.
    public static func map(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                            library: PhaseReferenceLibrary,
                            settings: PhaseVectorSettings,
                            originX: Float, originY: Float,
                            invAngstromPerPixel: Double,
                            matrixEntryIndex: Int? = nil,
                            maximumWorkers: Int? = nil,
                            progress: (@Sendable (Double) -> Void)? = nil,
                            isCancelled: (@Sendable () -> Bool)? = nil) -> PhaseMap? {
        map(scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            matrixEntryIndex: matrixEntryIndex, maximumWorkers: maximumWorkers,
            progress: progress, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func map(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                            library: PhaseReferenceLibrary,
                            settings: PhaseVectorSettings,
                            coordinates: PeakCoordinates,
                            matrixEntryIndex: Int? = nil,
                            maximumWorkers: Int? = nil,
                            progress: (@Sendable (Double) -> Void)? = nil,
                            isCancelled: (@Sendable () -> Bool)? = nil) -> PhaseMap? {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight) else { return nil }
        guard library.isConsistent else { return nil }
        if let matrixEntryIndex, !library.entries.indices.contains(matrixEntryIndex) { return nil }

        let fitted = matrixEntryIndex ?? fitMatrixOrientation(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: coordinates, isCancelled: isCancelled
        )?.entryIndex
        guard isCancelled?() != true else { return nil }

        let matrixEntry = fitted.map { library.entries[$0] }
        // The orientation relationship (2026-09-15): once the matrix is
        // fitted, a candidate entry whose phase lists an orientation
        // relationship is scored only when the relationship's derived
        // azimuth agrees with the fitted matrix entry's — this is what keeps
        // a θ′ edge-on entry, rotated 45° so its (002) at 0.345 Å⁻¹ sits on
        // face-on's (110) at 0.350, from being offered as face-on. A phase
        // with an empty list stays free, so the default (no matrix, or no
        // list anywhere) moves nothing.
        let candidates: [Int]
        if let matrixEntry {
            let matrixCrystal = library.phases[library.matrixPhaseIndex].crystal
            candidates = library.candidateEntryIndices.filter { index in
                let entry = library.entries[index]
                let relationships = library.phases[entry.phaseIndex].orientationRelationships
                // An empty list means free: it would otherwise empty the
                // candidate set and refuse the whole map.
                guard !relationships.isEmpty else { return true }
                return Self.orientationConsistent(
                    candidate: entry, candidateCrystal: library.phases[entry.phaseIndex].crystal,
                    matrix: matrixEntry, matrixCrystal: matrixCrystal,
                    relationships: relationships, toleranceDeg: settings.orientationRelationshipToleranceDeg)
            }
        } else {
            candidates = library.candidateEntryIndices
        }
        guard !candidates.isEmpty else { return nil }

        // Built once for the scan, scored only at positions that would
        // otherwise be labelled — see `classify` step 5 for what it is for.
        let challenge = matrixChallengeBases(library: library)

        var map = PhaseMap(width: scanWidth, height: scanHeight,
                           matrixEntryIndex: fitted ?? -1,
                           phaseNames: library.phases.map(\.displayName),
                           matrixPhaseIndex: library.matrixPhaseIndex)

        let capacity = library.entries.map(\.vectors.count).max() ?? 1
        let total = scan.count
        let counter = ProgressCounter()

        map.results.withUnsafeMutableBufferPointer { out in
            let ptr = SendableBox(out.baseAddress!)
            let workers = max(1, min(maximumWorkers ?? ProcessInfo.processInfo.activeProcessorCount, total))
            let progressBatch = max(1, total / (workers * 100))
            DispatchQueue.concurrentPerform(iterations: workers) { worker in
                guard isCancelled?() != true else { return }
                let scratch = Scratch(capacity: capacity)
                var pending = 0
                for position in Swift.stride(from: worker, to: total, by: workers) {
                    if isCancelled?() == true { break }
                    let vectors = experimentalVectors(
                        peaks: scan[position], coordinates: coordinates,
                        directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                        maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom
                    )
                    ptr.value[position] = classify(
                        vectors: vectors, library: library, settings: settings,
                        matrixEntry: matrixEntry, candidateEntryIndices: candidates,
                        scratch: scratch, matrixChallenge: challenge
                    )
                    pending += 1
                    if progress != nil && pending >= progressBatch {
                        let completed = counter.add(pending)
                        pending = 0
                        progress?(Double(completed) / Double(total))
                    }
                }
                if progress != nil && pending > 0 {
                    let completed = counter.add(pending)
                    progress?(Double(completed) / Double(total))
                }
            }
        }
        guard isCancelled?() != true else { return nil }
        return map
    }

    /// Minimal wrapper to pass a mutable pointer into the concurrent closure.
    private struct SendableBox: @unchecked Sendable {
        let value: UnsafeMutablePointer<PhaseVectorResult>
        init(_ v: UnsafeMutablePointer<PhaseVectorResult>) { value = v }
    }

    /// The positions done so far, shared by the workers: a lock around a captured counter, as in the apps' original.
    private final class ProgressCounter: @unchecked Sendable {
        private let lock = NSLock()
        private var done = 0
        func add(_ n: Int) -> Int {
            lock.lock(); done += n; let completed = done; lock.unlock()
            return completed
        }
    }
}

// MARK: - Which detected disk did which phase claim?

/// What the phase map made of one detected Bragg disk at one scan position.
public enum PhaseDiskClaim: Sendable, Equatable {
    /// Never offered to the matcher: inside the direct-beam radius, past the
    /// maximum vector, or not a finite position.
    case notMatched
    /// Removed as the fitted matrix orientation's (step 1 of `classify`).
    case matrix
    /// Paired, within the pair radius, with a reference vector of the winning
    /// entry of this phase (its index in the library's `phases`).
    case phase(Int)
    /// Offered to the matcher, and claimed by nothing: it survived matrix
    /// removal and the winning entry has no reference vector within the pair
    /// radius — or the position has no winning entry (matrix by exclusion, not
    /// indexed), so nothing beyond the matrix could claim it.
    case unexplained
}

extension PhaseVectorMatcher {

    /// The claim of every peak of one scan position, in `peaks` order.
    ///
    /// A REPLAY OF THE MATCHER'S OWN PAIRING, not a second opinion: usable
    /// vectors from `experimentalVector`, matrix removal by the same
    /// `nearest(_:in:radius:)` at `matrixToleranceInvAngstrom` that
    /// `classify` step 1 (and `classifyKnownVariants` step a) uses, and the
    /// phase claim by the same `nearest` at `pairRadiusInvAngstrom` against
    /// the entry the map recorded for this position
    /// (`result.entryIndex`) — the pairing `score(vectors:against:…)` counts
    /// as `matchedCount`. So, for an `.indexed` position, the number of
    /// `.phase` claims equals `result.matchedCount`, the `.matrix` claims
    /// `removedCount`, and matrix + phase + unexplained = the usable vectors.
    ///
    /// Nothing is re-decided here: the verdict, the winning phase and entry
    /// are read from `result`, so a claim can name only the phase the map
    /// drew. A position the map did not index gives no `.phase` claim.
    ///
    /// `peaks` must be the CALIBRATED peaks the run matched, and `matrixEntry`
    /// the run's fitted matrix entry, nil when none was fitted.
    ///
    /// LIMIT: where the verdict is `.matrix` because the matrix challenge
    /// (`classify` step 5) found the matrix on another zone axis, that axis
    /// is not recorded (`entryIndex` is cleared), so the disks it explains
    /// are reported as unexplained rather than matrix.
    public static func claims(peaks: [SIMD2<Double>],
                               originX: Float, originY: Float,
                               invAngstromPerPixel: Double,
                               settings: PhaseVectorSettings,
                               matrixEntry: PhaseOrientationReference?,
                               result: PhaseVectorResult,
                               library: PhaseReferenceLibrary) -> [PhaseDiskClaim] {
        claims(peaks: peaks,
               coordinates: .detectorPixels(originX: originX, originY: originY,
                                            invAngstromPerPixel: invAngstromPerPixel),
               settings: settings, matrixEntry: matrixEntry, result: result, library: library)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func claims(peaks: [SIMD2<Double>],
                               coordinates: PeakCoordinates,
                               settings: PhaseVectorSettings,
                               matrixEntry: PhaseOrientationReference?,
                               result: PhaseVectorResult,
                               library: PhaseReferenceLibrary) -> [PhaseDiskClaim] {
        var winning: PhaseOrientationReference?
        if result.verdict == .indexed, result.entryIndex >= 0,
           Int(result.entryIndex) < library.entries.count {
            winning = library.entries[Int(result.entryIndex)]
        }
        return peaks.map { peak in
            guard let u = experimentalVector(
                peak: peak, coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            else { return .notMatched }
            if let matrixEntry,
               nearest(u, in: matrixEntry.vectors,
                       radius: settings.matrixToleranceInvAngstrom) != nil {
                return .matrix
            }
            if let winning,
               nearest(u, in: winning.vectors, radius: settings.pairRadiusInvAngstrom) != nil {
                return .phase(Int(result.phaseIndex))
            }
            return .unexplained
        }
    }
}
