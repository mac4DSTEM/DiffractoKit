//
//  PhaseVectorMatcher+ZoneAxisFit.swift
//  Role: which zone axis the specimen is on, and the matrix orientation fitted once for the scan.
//  Split out of PhaseVectorMatching.swift, which carries the method, its sources and its limits.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

extension PhaseVectorMatcher {

    // MARK: Which zone axis is the specimen on?

    /// One candidate beam direction, scored against the scan's own peaks.
    public struct ZoneAxisFit: Sendable, Equatable {
        public init(zoneAxis: SIMD3<Int>, inPlaneRotationRad: Double,
                     matchedVectors: Int, totalVectors: Int, meanDistance: Double,
                     chanceMatchedVectors: Double, sweepMedianFraction: Double) {
            self.zoneAxis = zoneAxis
            self.inPlaneRotationRad = inPlaneRotationRad
            self.matchedVectors = matchedVectors
            self.totalVectors = totalVectors
            self.meanDistance = meanDistance
            self.chanceMatchedVectors = chanceMatchedVectors
            self.sweepMedianFraction = sweepMedianFraction
        }

        public let zoneAxis: SIMD3<Int>
        public let inPlaneRotationRad: Double
        /// Experimental vectors this axis accounts for, across the sample.
        public let matchedVectors: Int
        /// Experimental vectors in the sample, total.
        public let totalVectors: Int
        /// Mean |u − v| over the pairs it made, Å⁻¹.
        public let meanDistance: Double
        /// How many of `totalVectors` an entry this dense would match by
        /// accident, summed over the sampled patterns at each pattern's own
        /// reach. The number that says whether `matchedVectors` is evidence.
        ///
        /// WHY IT IS HERE (2026-09-14). Without it the sweep presented a
        /// ⟨112⟩ family at 8 % exactly the way it presents a ⟨110⟩ family at
        /// 38 %, and a reader took a coin toss for a fit. It is computed from
        /// the same `chanceMatchFraction` the matcher's own guard uses, so
        /// there is one definition of chance in this file rather than two.
        ///
        /// **IT DOES NOT MARK THAT MOTIVATING CASE, and saying so is the
        /// point** (a review, 2026-09-14). At the default reference
        /// settings aluminium's ⟨112⟩ entries carry 12 to 16 vectors, not the
        /// 48 the cap allows, so the expectation is 0.37–1.3 % and 8 % clears
        /// five times it. The mark fires only while each pattern's
        /// second-largest |u| stays under about 0.55–0.63 Å⁻¹, and Al {220}
        /// alone is at 0.699. What this buys is the number itself, reported
        /// instead of absent. That case is caught by the SECOND null
        /// below, `sweepMedianFraction`, not by a threshold on this one.
        public let chanceMatchedVectors: Double
        /// The sweep's own null (a diagnostic run, 2026-09-15): the explained
        /// fraction of the MEDIAN axis in the sweep this fit came from.
        ///
        /// WHY A SECOND NULL. On a real crystal every wrong axis explains a
        /// large share of the vectors by SHARED reflections, not by chance:
        /// measured on planted aluminium at 2° steps with 30 % of spots
        /// missing, 0.004 Å⁻¹ of jitter and three spurious peaks per
        /// pattern, the median wrong axis explains 11–16 % (⟨110⟩, ⟨112⟩
        /// and ⟨001⟩ plants) and the worst wrong family up to 25 %, while
        /// every one of the 49 axes clears five times the disc-chance
        /// expectation. That is why a ⟨112⟩ at 8 % could never be
        /// marked by `chanceMatchedVectors`: it was not chance, it was the
        /// true axis's reflections seen through a wrong projection. Against
        /// the sweep's median those true families sat at 4.8–6.8× and every
        /// wrong family at 0.7–1.6×; the bar is placed at 2×. On vectors
        /// pointing nowhere the sweep median is ~0.5 % and the disc rule is
        /// what refuses; the two rules cover different failures and both
        /// must pass.
        ///
        /// WHAT THE BAR IS NOT (a review, 2026-09-15): a clean separator
        /// for every truth. ⟨111⟩, ⟨012⟩ and ⟨210⟩ plants put the true family
        /// at 2.7–2.8× and the worst wrong at 1.0–2.0×; a ⟨122⟩ plant puts
        /// the true family at 2.8× and the ⟨100⟩ family at 2.7–2.9× — a
        /// wrong family clearing the bar, sometimes above the truth. A
        /// high-index truth with few reflections shares too much with a
        /// low-index axis for a ratio to tell them apart; both rows then
        /// read as informative and the tie caption is the honest thing on
        /// screen. Heavier degradation (to 92 % missing, 15 spurious) never
        /// made the sweep rule bind before the disc rule; a two-grain scan
        /// left both true families at 5×. DEVIATION: py4DSTEM's zone-axis
        /// tools have no sweep null and no informativeness verdict.
        ///
        /// Fewer than five axes in the sweep is no null at all — with two,
        /// the "median" is the top axis's own fraction and the ratio
        /// saturates at 1 — so the median is reported as 0 (ratio infinite,
        /// the sweep rule inert) below that count.
        public let sweepMedianFraction: Double
        /// Explained fraction relative to the sweep's median axis.
        public var sweepRatio: Double {
            sweepMedianFraction > 0 ? explainedFraction / sweepMedianFraction : .infinity
        }
        /// Below this ratio an axis explains no more than a wrong axis does
        /// on the same data (measured: wrong ≤ 1.6×, true ≥ 4.8×).
        public static let informativeSweepRatio = 2.0
        public var isAboveSweep: Bool { sweepRatio >= Self.informativeSweepRatio }
        /// Both nulls: chance vectors (the disc model) and wrong axes on a
        /// real crystal (the sweep median).
        public func isInformative(multiple: Double) -> Bool {
            isAboveChance(multiple: multiple) && isAboveSweep
        }

        public var explainedFraction: Double {
            totalVectors > 0 ? Double(matchedVectors) / Double(totalVectors) : 0
        }
        /// The same fraction chance alone would explain.
        public var chanceFraction: Double {
            totalVectors > 0 ? chanceMatchedVectors / Double(totalVectors) : 0
        }
        /// Does this fit carry information, by the matcher's own standard?
        /// `multiple` is `PhaseVectorSettings.chanceMatchMultiple`.
        public func isAboveChance(multiple: Double) -> Bool {
            Double(matchedVectors) >= multiple * chanceMatchedVectors
        }

        /// Mean distances closer than this, in Å⁻¹, rank as equal.
        ///
        /// WHY (2026-09-30): symmetry-equivalent axes tie on matched count
        /// by construction, and their mean distances then differ only in the
        /// last bits of the arithmetic that rotated their references — so a
        /// tied family's order was set by rounding: deterministic for one build
        /// and peak set, but not a property of the data. Hardening, not the
        /// fix of a recorded "run-dependent tie order": that came from a
        /// caller ranking through a Dictionary (hash order per process).
        /// The bar is far under
        /// anything a detector resolves: float32 peaks carry ~1e-7 Å⁻¹.
        public static let meanDistanceTieResolution = 1e-9

        /// The ranking `fitZoneAxis` returns, as a TOTAL order: more matched
        /// vectors first; then the lower mean distance (compared on the
        /// `meanDistanceTieResolution` grid); then the axis indices,
        /// lexicographically ascending — the order `lowIndexZoneAxes` already
        /// lists one symmetry family in, so a tied family reads as it did
        /// when it tied exactly (an exact tie ACROSS families may not). The
        /// last key only ever decides between fits the first two cannot
        /// tell apart; it never outranks a better count or a closer fit.
        public static func ranksBefore(_ a: ZoneAxisFit, _ b: ZoneAxisFit) -> Bool {
            if a.matchedVectors != b.matchedVectors { return a.matchedVectors > b.matchedVectors }
            let da = rankedDistance(a.meanDistance), db = rankedDistance(b.meanDistance)
            if da != db { return da < db }
            if a.zoneAxis.x != b.zoneAxis.x { return a.zoneAxis.x < b.zoneAxis.x }
            if a.zoneAxis.y != b.zoneAxis.y { return a.zoneAxis.y < b.zoneAxis.y }
            return a.zoneAxis.z < b.zoneAxis.z
        }

        /// A grid index, so equality is transitive (a tolerance comparison is
        /// not, and `sort` needs a strict weak order). Non-finite ranks last.
        private static func rankedDistance(_ d: Double) -> Double {
            d.isFinite ? (d / meanDistanceTieResolution).rounded() : .infinity
        }
    }

    /// Rank beam directions by how much of the scan's measured peaks each one
    /// explains — the question "what orientation is this specimen on?", asked
    /// of the data instead of assumed.
    ///
    /// WHY IT EXISTS. A matrix viewed down an axis it is not on presents no
    /// reflections to remove, nothing is removed, and every position comes
    /// back "not indexed" — which on screen is indistinguishable from a method
    /// that does not work. That happened on 2026-09-12: Al set to
    /// [001], matrix verdict count **zero** across 108 900 positions. Measured
    /// afterwards, the specimen was on **⟨110⟩**, and
    /// all five sampled ⟨110⟩ equivalents tied at 39.0 % — an exact tie across
    /// the family being what cubic symmetry requires, and therefore a check
    /// that the sweep is behaving rather than merely returning something.
    ///
    /// Ranked by explained fraction, tie-broken on the lower mean distance,
    /// then on the axis indices (`ZoneAxisFit.ranksBefore`, a total order).
    /// NOT by mean distance alone, for the same reason `fitMatrixOrientation`
    /// is not: an axis explaining one vector at 0.001 Å⁻¹ would beat one
    /// explaining nine at 0.01.
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`); a count that is not
    /// width · height gives no fits, as no peaks do. Throws when the kit refuses
    /// the structure (`PhaseReferenceLibrary.Failure.structureRefused`,
    /// `.unsupportedElements`), and `.libraryTooLarge(requested:limit:)` when `inPlaneStepDeg` asks for more than
    /// `PhaseReferenceLibrary.maximumInPlaneStepCount` rotations (a step finer than 0.001 degrees), with the count asked
    /// for and the cap.
    public static func fitZoneAxis(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                    crystal: CrystalStructure,
                                    referenceSettings: PhaseReferenceSettings,
                                    settings: PhaseVectorSettings,
                                    originX: Float, originY: Float,
                                    invAngstromPerPixel: Double,
                                    candidateAxes: [SIMD3<Int>]
                                        = PhaseReferenceLibrary.lowIndexZoneAxes,
                                    inPlaneStepDeg: Double = 5,
                                    sampleLimit: Int = 400,
                                    isCancelled: (@Sendable () -> Bool)? = nil)
        throws(PhaseReferenceLibrary.Failure) -> [ZoneAxisFit] {
        try fitZoneAxis(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, crystal: crystal,
            referenceSettings: referenceSettings, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            candidateAxes: candidateAxes, inPlaneStepDeg: inPlaneStepDeg,
            sampleLimit: sampleLimit, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func fitZoneAxis(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                    crystal: CrystalStructure,
                                    referenceSettings: PhaseReferenceSettings,
                                    settings: PhaseVectorSettings,
                                    coordinates: PeakCoordinates,
                                    candidateAxes: [SIMD3<Int>]
                                        = PhaseReferenceLibrary.lowIndexZoneAxes,
                                    inPlaneStepDeg: Double = 5,
                                    sampleLimit: Int = 400,
                                    isCancelled: (@Sendable () -> Bool)? = nil)
        throws(PhaseReferenceLibrary.Failure) -> [ZoneAxisFit] {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight),
              !scan.isEmpty, !candidateAxes.isEmpty else { return [] }

        // A subsample, spread across the scan: this answers a question about
        // the GRAIN, and a grain does not change between neighbouring probe
        // positions. Sampling is what keeps a 62-axis sweep interactive.
        let step = max(1, scan.count / max(1, sampleLimit))
        var sample: [[SIMD2<Double>]] = []
        var index = 0
        while index < scan.count {
            let v = experimentalVectors(
                peaks: scan[index], coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            if !v.isEmpty { sample.append(v) }
            index += step
        }
        guard !sample.isEmpty else { return [] }
        let totalVectors = sample.reduce(0) { $0 + $1.count }

        // Each pattern's own reach, for the chance expectation: the
        // SECOND-largest |u|, the same rule `classify` uses and for the same
        // measured reason — with `max`, one spurious maximum far out inflates
        // the area and weakens the estimate for every other vector.
        let reach: [Double] = sample.map { vectors in
            var largest = 0.0, second = 0.0
            for u in vectors {
                let r = simd_length(u)
                if r > largest { second = largest; largest = r }
                else if r > second { second = r }
            }
            return second > 0 ? second : largest
        }

        var reference = referenceSettings
        reference.inPlaneStepDeg = inPlaneStepDeg
        // A step finer than 360 / maximumInPlaneStepCount degrees is refused, with the count asked for (saturated at
        // Int.max) and the cap, not held at the cap: a sweep at a different step than the one asked for would be a
        // different answer.
        let stepCount = PhaseReferenceLibrary.inPlaneStepCount(reference)
        guard stepCount <= Double(PhaseReferenceLibrary.maximumInPlaneStepCount) else {
            throw PhaseReferenceLibrary.Failure.libraryTooLarge(
                requested: stepCount < 9.0e18 ? Int(stepCount) : Int.max,
                limit: PhaseReferenceLibrary.maximumInPlaneStepCount)
        }
        let reflections = try PhaseReferenceLibrary.reflections(of: crystal, kMaxInvAngstrom: reference.kMaxInvAngstrom)
        let rotations = PhaseReferenceLibrary.inPlaneSteps(reference)
        let scratch = Scratch(capacity: max(1, reference.maximumVectorsPerEntry))

        // Every axis's best rotation first; the sweep median is known only
        // once all of them are, and each fit carries it.
        var winners: [(axis: SIMD3<Int>, theta: Double, matched: Int, meanDistance: Double, chance: Double)] = []
        for axis in candidateAxes {
            if isCancelled?() == true { break }
            let base = PhaseReferenceLibrary.projectedVectors(
                reflections: reflections, crystal: crystal,
                zoneAxis: axis, settings: reference)
            guard !base.isEmpty else { continue }
            // Constant across this axis's rotations: rotating an entry moves
            // its references, it does not change how many there are.
            let unrotated = PhaseOrientationReference(
                phaseIndex: 0, zoneAxis: axis, inPlaneRotationRad: 0, vectors: base)
            var chanceMatched = 0.0
            for (index, vectors) in sample.enumerated() {
                chanceMatched += unrotated.chanceMatchFraction(
                    pairRadius: settings.matrixToleranceInvAngstrom,
                    accessibleRadius: reach[index]) * Double(vectors.count)
            }
            var best: (theta: Double, matched: Int, meanDistance: Double)?
            for theta in rotations {
                let entry = PhaseOrientationReference(
                    phaseIndex: 0, zoneAxis: axis, inPlaneRotationRad: theta,
                    vectors: PhaseReferenceLibrary.rotate(base, by: theta))
                var matched = 0, pairs = 0
                var total = 0.0
                for vectors in sample {
                    guard let sc = score(vectors: vectors, against: entry,
                                         pairRadius: settings.matrixToleranceInvAngstrom,
                                         scratch: scratch) else { continue }
                    matched += sc.matched
                    total += sc.score * Double(sc.uniqueReferences)
                    pairs += sc.uniqueReferences
                }
                guard pairs > 0 else { continue }
                let meanDistance = total / Double(pairs)
                if best == nil || matched > best!.matched
                    || (matched == best!.matched && meanDistance < best!.meanDistance) {
                    best = (theta, matched, meanDistance)
                }
            }
            if let best {
                winners.append((axis, best.theta, best.matched, best.meanDistance, chanceMatched))
            }
        }
        let fractions = winners.map { Double($0.matched) / Double(totalVectors) }.sorted()
        let sweepMedian = fractions.count >= 5 ? fractions[fractions.count / 2] : 0
        var out = winners.map {
            ZoneAxisFit(zoneAxis: $0.axis, inPlaneRotationRad: $0.theta,
                        matchedVectors: $0.matched, totalVectors: totalVectors,
                        meanDistance: $0.meanDistance, chanceMatchedVectors: $0.chance,
                        sweepMedianFraction: sweepMedian)
        }
        out.sort(by: ZoneAxisFit.ranksBefore)
        return out
    }

    // MARK: The matrix orientation, fitted once

    /// The single matrix entry that best explains the whole scan.
    ///
    /// DEVIATION, deliberate: the paper masks the matrix per pattern. Here the
    /// matrix is one grain across the scan, so its in-plane rotation is one
    /// number — fitting it once is both faster and steadier than re-deciding
    /// it at every position from a handful of peaks. The cost is stated in the
    /// file header: a second matrix grain would be mis-removed.
    ///
    /// Chosen by total matched vectors, tie-broken by lower mean distance. Not
    /// by mean distance alone: an entry matching one vector at 0.001 Å⁻¹ would
    /// win over one matching nine at 0.01.
    ///
    /// The candidate side is guarded by `minimumMatchedVectors`
    /// and `chanceMatchMultiple`, which are floors on COUNT, not on the
    /// fraction of an entry's accessible vectors that matched. So the trap IS
    /// open on the candidate side: the cross-phase winner below is still
    /// chosen by mean distance alone, and a 2-vector Friedel pair at 0.004
    /// beats a 10-vector match at 0.012. That is measured, not hypothetical
    /// (measured, and recorded as an open item).
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`); a count that is not
    /// width · height is refused (nil), as is a scan with no peaks.
    public static func fitMatrixOrientation(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                             library: PhaseReferenceLibrary,
                                             settings: PhaseVectorSettings,
                                             originX: Float, originY: Float,
                                             invAngstromPerPixel: Double,
                                             sampleLimit: Int = 2000,
                                             isCancelled: (@Sendable () -> Bool)? = nil)
        -> (entryIndex: Int, matched: Int, meanDistance: Double)? {
        fitMatrixOrientation(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            sampleLimit: sampleLimit, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func fitMatrixOrientation(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                             library: PhaseReferenceLibrary,
                                             settings: PhaseVectorSettings,
                                             coordinates: PeakCoordinates,
                                             sampleLimit: Int = 2000,
                                             isCancelled: (@Sendable () -> Bool)? = nil)
        -> (entryIndex: Int, matched: Int, meanDistance: Double)? {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight) else { return nil }
        let matrixEntries = library.matrixEntryIndices
        guard !matrixEntries.isEmpty, !scan.isEmpty else { return nil }

        let stride = max(1, scan.count / max(1, sampleLimit))
        var sample: [[SIMD2<Double>]] = []
        sample.reserveCapacity(min(scan.count, sampleLimit) + 1)
        var index = 0
        while index < scan.count {
            let v = experimentalVectors(
                peaks: scan[index], coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom
            )
            if !v.isEmpty { sample.append(v) }
            index += stride
        }
        guard !sample.isEmpty else { return nil }

        let capacity = matrixEntries.map { library.entries[$0].vectors.count }.max() ?? 1
        let scratch = Scratch(capacity: capacity)
        var best: (entryIndex: Int, matched: Int, meanDistance: Double)?
        for entryIndex in matrixEntries {
            if isCancelled?() == true { return best }
            let entry = library.entries[entryIndex]
            var matched = 0
            var sum = 0.0
            var pairs = 0
            for vectors in sample {
                guard let s = score(vectors: vectors, against: entry,
                                    pairRadius: settings.matrixToleranceInvAngstrom,
                                    scratch: scratch) else { continue }
                matched += s.matched
                sum += s.score * Double(s.uniqueReferences)
                pairs += s.uniqueReferences
            }
            guard matched > 0, pairs > 0 else { continue }
            let mean = sum / Double(pairs)
            if let b = best {
                if matched > b.matched || (matched == b.matched && mean < b.meanDistance) {
                    best = (entryIndex, matched, mean)
                }
            } else {
                best = (entryIndex, matched, mean)
            }
        }
        return best
    }

    /// Is `scan` a scan: width and height positive, and exactly width · height positions? A malformed scan is refused by
    /// `map` and `fitMatrixOrientation` (nil) and gives no zone-axis fits, as the apps' original refused one. A width · height that
    /// overflows is not a scan.
    static func scanIsWellFormed(_ scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int) -> Bool {
        guard scanWidth > 0, scanHeight > 0 else { return false }
        let (positions, overflow) = scanWidth.multipliedReportingOverflow(by: scanHeight)
        return !overflow && scan.count == positions
    }

    /// True when `a` outranks `b` as the cross-phase winner. Shipped
    /// (`completenessAware: false`): mean distance alone, ascending —
    /// `a.score < b.score`. Candidate (`completenessAware: true`): matched
    /// count first, descending, mean distance only as the tiebreak — the
    /// same rule `fitMatrixOrientation` already uses above for candidate
    /// zone axes. Extracted from `classify`'s sort so it is testable without
    /// a crystal library: `PhaseVectorSettings.completenessAwareCrossPhaseRanking`'s
    /// own doc comment has the measurement and the open item.
    public static func crossPhaseWinsOver(
        _ a: (matched: Int, score: Double), _ b: (matched: Int, score: Double),
        completenessAware: Bool
    ) -> Bool {
        guard completenessAware else { return a.score < b.score }
        return a.matched != b.matched ? a.matched > b.matched : a.score < b.score
    }
}
