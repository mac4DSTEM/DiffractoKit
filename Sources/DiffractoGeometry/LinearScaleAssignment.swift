//
//  LinearScaleAssignment.swift
//  Role: which of a crystal's shells is the first observed ring? `LinearScale.estimate` takes its reference shell from the caller; on Al [001]
//        the {111} reflections are not in the zone, the first visible ring is {200}, and assuming the first shell read Δq 13.4 % low.
//        `LinearScale.assign` tests the first few shells against the SEQUENCE of observed ring radii and refuses when the sequence cannot tell.
//        One case is assigned although the sequence cannot tell: two hypotheses that predict the SAME ratios (Al [001] with two rings is
//        {200}, {220} and equally {220}, {400}). No measurement of those rings separates them, so the lower shell is taken and the other is
//        named in `equallyFitting` (the single-ring rule's reasoning: the innermost ring is the lowest shell that fits).
//        The rings are read twice (position by position, and as the modes of all radii together) and nothing is assigned unless both readings
//        give the same scale: read one way only, a polycrystal's broad {111} ring was two rings that fit {311}, {222} (mac4DSTEM#110).
//        Ported from mac4DSTEM's `ShellAssignment` (calibration design 2026-10-22, section 6), which this replaces.
//
//  NO Python reference exists (py4DSTEM does not do this): checked against planted truth only, `validation: "none"`.
//

import Foundation

/// The shell found for the first observed ring.
public struct LinearScaleAssignment: Sendable, Equatable {
    /// Index into the crystal's distinct shells (0 = the shortest |g|).
    public let shellIndex: Int
    public let label: String
    /// The assigned shell's |g|, Å⁻¹: the reference radius of `LinearScale.estimate`.
    public let referenceRadiusInvAngstrom: Double
    /// The shell matched to the second observed ring, or the next listed shell when there is no second ring.
    public let nextShellInvAngstrom: Double?
    /// Median radii (pixels) of the observed rings, innermost first.
    public let observedRadiiPixels: [Double]
    /// Observed rings matched by the assignment (the first included), of `observedRings`.
    public let matchedRings: Int
    public let observedRings: Int
    /// RMS relative misfit over the matched rings after the first (0 when there are none).
    public let rmsMisfit: Double
    /// Labels of the higher starting shells that predict the same ratios as the assigned one and so fit these rings exactly as well
    /// (empty for a decided assignment). The lower shell was taken; a reader shows these beside the label.
    public let equallyFitting: [String]
    /// False with one observed ring: nothing could be tested and the first ring is assumed to be the first shell.
    public var isChecked: Bool { observedRings >= 2 }
    /// Always "none": no standard has validated the assignment.
    public var validation: String { "none" }
}

/// Why no shell was assigned. The sentences name what to do.
public enum LinearScaleAssignmentError: Error, LocalizedError, Equatable, Sendable {
    /// No starting shell explains the second ring: the observed ratios and what the crystal predicts for the first two hypotheses.
    case noMatch(crystal: String, observedRatios: [Double], predicted: [(label: String, ratios: [Double])])
    /// Two assignments fit equally well.
    case ambiguous(first: String, second: String)
    /// The two readings of the rings (position by position, and pooled over the scan) give different scales: the first ring each one found, with its radius.
    case readingsDisagree(perPosition: String, pooled: String)

    public static func == (a: Self, b: Self) -> Bool {
        switch (a, b) {
        case let (.noMatch(c1, o1, p1), .noMatch(c2, o2, p2)):
            return c1 == c2 && o1 == o2 && p1.count == p2.count && zip(p1, p2).allSatisfy { $0.label == $1.label && $0.ratios == $1.ratios }
        case let (.ambiguous(a1, b1), .ambiguous(a2, b2)): return a1 == a2 && b1 == b2
        case let (.readingsDisagree(a1, b1), .readingsDisagree(a2, b2)): return a1 == a2 && b1 == b2
        default: return false
        }
    }

    public var errorDescription: String? {
        switch self {
        case .noMatch(let crystal, let observed, let predicted):
            let seen = observed.map { String(format: "%.2f", $0) }.joined(separator: ", ")
            let sentences = predicted.map { "\($0.ratios.map { String(format: "%.2f", $0) }.joined(separator: ", ")) (\($0.label) first)" }
            return "The rings do not match any starting shell of \(crystal): the first ring gives ratios \(seen) but the crystal predicts \(sentences.joined(separator: " or "))."
        case .ambiguous(let first, let second):
            return "Two assignments fit equally well (\(first) and \(second) first). Pick the first ring in the crystal's list, or measure the pixel size from a distance."
        case .readingsDisagree(let perPosition, let pooled):
            return "The rings do not decide the first shell: read position by position the first ring is \(perPosition), read over the whole scan it is \(pooled), and the two scales differ. Broad or split rings do that. Correct the elliptical distortion first, pick the first ring in the crystal's list, or measure the pixel size from a distance."
        }
    }
}

extension LinearScale {
    /// Cluster width (relative) of the equivalents of one ring at one position, measured from the ring's first peak; also the half width of
    /// the window the pooled radii are smoothed with.
    static let assignmentClusterBand = 0.03
    /// The pooled radii are counted in bins this wide (relative) and smoothed with a triangular window of the cluster band on each side.
    static let assignmentBinWidth = 0.01
    /// Two neighbouring modes of the pooled radii are separate rings when the valley between them is at most this share of the lower mode.
    static let assignmentValley = 0.5
    /// The two readings agree when their scales are within this (relative): the tolerance a ring is matched to a shell with.
    static let assignmentReadingAgreement = 0.03
    /// Two hypotheses predict the same ratios when every matched shell ratio agrees to this (relative): floating point, not noise.
    static let assignmentSameRatios = 1e-9
    /// The first K observed rings are read; a ring needs `assignmentMinimumPositionShare` of the positions that hold any peak.
    static let assignmentMaximumRings = 3
    static let assignmentMinimumPositionShare = 0.2
    /// How many starting shells are tried.
    static let assignmentHypotheses = 4
    /// How many of the crystal's shells a hypothesis may match rings to: the starting shells plus the rings read, plus one. Without the bound a
    /// cubic crystal explains any ring set 1 : √2 : 2 twice ({200}, {220}, {400} and {220}, {400}, {440}); the bound keeps the matches to the
    /// low-order shells a real pattern shows (a pick, pinned by the planted Al [001] case).
    static let assignmentShellWindow = assignmentHypotheses + assignmentMaximumRings + 1

    /// Which of the crystal's first shells is the first observed ring.
    ///
    /// The rings are read twice, position by position (`ringsPerPosition`) and over the whole scan (`ringsPooled`), and each reading is scored
    /// on its own. An assignment is returned only when both readings give one and their scales agree within 3 %: the per-position one, whose
    /// first ring is the radius `estimate` measures. A broad ring, a split ring or a missing ring misleads one reading and not the other in the
    /// same way, so a disagreement is refused (`readingsDisagree`) instead of answered. What follows describes one reading.
    ///
    /// Per position the peaks are clustered (the equivalents of one ring: within 3 % of its first peak), the first three rings are
    /// read, and each ring's radius is the median over the positions that have it (at least 20 % of those with any peak). The hypotheses H_j "the first ring is shell j" (j below 4) predict
    /// the later rings; each observed ring is matched to the nearest later shell within max(3 % of its radius, 1 px). The best hypothesis wins
    /// if it matches strictly more rings than the runner-up, or the same number with a misfit under half and 0.5 % lower. Hypotheses that
    /// predict the same ratios as the best one cannot be told apart by these rings: the lowest of them is assigned and the others are named
    /// in `equallyFitting`. One observed ring
    /// cannot be tested: shell 0 is assumed, unless `priorScale` agrees with exactly one hypothesis within 10 % (it never overrides a sequence).
    ///
    /// - Parameters:
    ///   - radiiPerPosition: peak radii in pixels from the origin, one array per scan position.
    ///   - shellLengths: the crystal's distinct |g| (Å⁻¹), ascending.
    ///   - labels: one family name per shell for the sentences ("{200}"); "shell k" when absent.
    ///   - crystalName: for the "no match" sentence.
    ///   - priorScale: a Δq (Å⁻¹ per pixel) from the file or typed.
    /// - Returns: nil when no position holds a radius above `minimumRadiusPixels` (or the shell list is empty or has a non-positive length).
    /// - Throws: `LinearScaleAssignmentError` when the sequence matches no starting shell, two assignments with different predictions fit
    ///   equally well, or the two readings of the rings give different scales.
    public static func assign(
        radiiPerPosition: [[Double]], shellLengths: [Double], labels: [String]? = nil, crystalName: String = "the crystal",
        priorScale: Double? = nil, minimumRadiusPixels: Double = 2
    ) throws -> LinearScaleAssignment? {
        guard shellLengths.count > 0, shellLengths.allSatisfy({ $0.isFinite && $0 > 0 }) else { return nil }
        func name(_ k: Int) -> String { labels.flatMap { k < $0.count ? $0[k] : nil } ?? "shell \(k + 1)" }

        // Two readings of the same peaks, with different weaknesses, must give the same scale; a reading that refuses refuses for both.
        let perPosition = ringsPerPosition(radiiPerPosition, minimumRadiusPixels: minimumRadiusPixels)
        let pooled = ringsPooled(radiiPerPosition, minimumRadiusPixels: minimumRadiusPixels)
        guard let a = try assign(rings: perPosition, shellLengths: shellLengths, name: name, crystalName: crystalName, priorScale: priorScale),
              let b = try assign(rings: pooled, shellLengths: shellLengths, name: name, crystalName: crystalName, priorScale: priorScale) else { return nil }
        let scaleA = a.referenceRadiusInvAngstrom / a.observedRadiiPixels[0], scaleB = b.referenceRadiusInvAngstrom / b.observedRadiiPixels[0]
        guard abs(scaleA / scaleB - 1) <= assignmentReadingAgreement else {
            func at(_ x: LinearScaleAssignment) -> String { "\(x.label) at \(String(format: "%.1f", x.observedRadiiPixels[0])) px" }
            throw LinearScaleAssignmentError.readingsDisagree(perPosition: at(a), pooled: at(b))
        }
        return a
    }

    /// The rings read position by position: at each position the peaks are clustered (the equivalents of one ring: within 3 % of its first
    /// peak) and the first three clusters are its rings; a ring's radius is the median over the positions that have it (at least 20 % of
    /// those with any peak). Its weakness: a ring that spreads more than the band is split, and a position that lacks a ring shifts the
    /// count of the later ones.
    static func ringsPerPosition(_ radiiPerPosition: [[Double]], minimumRadiusPixels: Double) -> [Double] {
        var perRing = [[Double]](repeating: [], count: assignmentMaximumRings)
        var positionsWithPeaks = 0
        for peaks in radiiPerPosition {
            let radii = peaks.filter { $0.isFinite && $0 > minimumRadiusPixels }.sorted()
            guard !radii.isEmpty else { continue }
            positionsWithPeaks += 1
            var i = 0, ring = 0
            while i < radii.count, ring < assignmentMaximumRings {
                let limit = radii[i] * (1 + assignmentClusterBand)
                var j = i, sum = 0.0
                while j < radii.count, radii[j] <= limit { sum += radii[j]; j += 1 }
                perRing[ring].append(sum / Double(j - i))
                i = j; ring += 1
            }
        }
        guard positionsWithPeaks > 0 else { return [] }
        let needed = max(1, Int((assignmentMinimumPositionShare * Double(positionsWithPeaks)).rounded(.up)))
        var observed: [Double] = []
        for ring in 0 ..< assignmentMaximumRings {
            guard perRing[ring].count >= needed else { break }
            observed.append(median(perRing[ring]))
        }
        return observed
    }

    /// The rings read over the whole scan: the modes of the radii of all positions together (1 % bins on a logarithmic axis, smoothed over the
    /// 3 % band; two modes are separate when the valley between them is at most half the lower one). A ring needs at least 20 % of the
    /// positions that have any peak; its radius is the median of its radii. Its weakness: grains of several orientations add their rings
    /// into one sequence no single position shows, and an elliptical distortion splits a ring into two modes.
    static func ringsPooled(_ radiiPerPosition: [[Double]], minimumRadiusPixels: Double) -> [Double] {
        let positions = radiiPerPosition.map { $0.filter { $0.isFinite && $0 > minimumRadiusPixels && $0 > 0 } }.filter { !$0.isEmpty }
        guard let smallest = positions.joined().min(), let largest = positions.joined().max() else { return [] }
        let step = log(1 + assignmentBinWidth)
        let first = (log(smallest) / step).rounded(.down)
        func bin(_ r: Double) -> Int { Int((log(r) / step).rounded(.down) - first) }
        let half = Int((assignmentClusterBand / assignmentBinWidth).rounded())
        var counts = [Double](repeating: 0, count: bin(largest) + 1)
        for radii in positions { for r in radii { counts[bin(r)] += 1 } }
        // A triangular window of the cluster band on each side: one ring's equivalents within 3 % are one mode.
        var smooth = [Double](repeating: 0, count: counts.count)
        for i in counts.indices where counts[i] > 0 {
            for d in -half ... half where smooth.indices.contains(i + d) { smooth[i + d] += counts[i] * Double(half + 1 - abs(d)) }
        }
        // Local maxima; two neighbours are two modes only when the valley between them is at most half the lower one.
        var peaks: [Int] = [], cuts: [Int] = []
        for i in smooth.indices where smooth[i] > 0 && (i == 0 || smooth[i] > smooth[i - 1]) && (i == smooth.count - 1 || smooth[i] >= smooth[i + 1]) {
            guard let last = peaks.last else { peaks.append(i); continue }
            let valley = (last ... i).min { smooth[$0] < smooth[$1] }!
            if smooth[valley] <= assignmentValley * min(smooth[last], smooth[i]) {
                peaks.append(i); cuts.append(valley)
            } else if smooth[i] > smooth[last] {
                peaks[peaks.count - 1] = i
            }
        }
        // Each mode's radii and the positions that have it; a ring needs its share of the positions.
        func mode(of r: Double) -> Int { let b = bin(r); return cuts.firstIndex { b < $0 } ?? cuts.count }
        var members = [[Double]](repeating: [], count: peaks.count), holders = [Int](repeating: 0, count: peaks.count)
        for radii in positions {
            var seen = Set<Int>()
            for r in radii { let m = mode(of: r); members[m].append(r); seen.insert(m) }
            for m in seen { holders[m] += 1 }
        }
        let needed = max(1, Int((assignmentMinimumPositionShare * Double(positions.count)).rounded(.up)))
        return members.indices.filter { holders[$0] >= needed }.prefix(assignmentMaximumRings).map { median(members[$0]) }
    }

    /// The hypotheses scored against one reading of the rings (`assign`'s documentation). Nil without a first ring.
    private static func assign(rings observed: [Double], shellLengths: [Double], name: (Int) -> String, crystalName: String,
                               priorScale: Double?) throws -> LinearScaleAssignment? {
        guard let r1 = observed.first, r1 > 0 else { return nil }

        let shells = Array(shellLengths.prefix(assignmentShellWindow))
        let tried = min(assignmentHypotheses, shells.count)
        func assigned(_ j: Int, matched: Int, rms: Double, secondShell: Int?, equallyFitting: [String] = []) -> LinearScaleAssignment {
            LinearScaleAssignment(
                shellIndex: j, label: name(j), referenceRadiusInvAngstrom: shells[j],
                nextShellInvAngstrom: secondShell.map { shells[$0] } ?? (j + 1 < shells.count ? shells[j + 1] : nil),
                observedRadiiPixels: observed, matchedRings: matched, observedRings: observed.count, rmsMisfit: rms,
                equallyFitting: equallyFitting)
        }

        if observed.count == 1 || shells.count == 1 {
            var j = 0
            if let prior = priorScale, prior.isFinite, prior > 0 {
                let agreeing = (0 ..< tried).filter { abs(shells[$0] / r1 / prior - 1) <= 0.10 }
                if agreeing.count == 1 { j = agreeing[0] }
            }
            return assigned(j, matched: 1, rms: 0, secondShell: nil)
        }

        /// `shellPerRing[i]` is the shell matched to observed ring i + 1, nil when none was.
        struct Score { let j: Int; let matched: Int; let rms: Double; let secondShell: Int?; let shellPerRing: [Int?] }
        var scores: [Score] = []
        for j in 0 ..< tried {
            let scale = shells[j] / r1
            var last = j, matched = 1, squares = 0.0, second: Int?
            var shellPerRing = [Int?](repeating: nil, count: observed.count - 1)
            for ring in 1 ..< observed.count {
                let r = observed[ring]
                var best: (k: Int, error: Double)?
                for k in (last + 1) ..< max(shells.count, last + 1) {
                    let error = abs(shells[k] / scale - r) / r
                    if best == nil || error < best!.error { best = (k, error) }
                }
                guard let hit = best, abs(shells[hit.k] / scale - r) <= max(0.03 * r, 1) else { continue }
                matched += 1; squares += hit.error * hit.error; last = hit.k
                shellPerRing[ring - 1] = hit.k
                if ring == 1 { second = hit.k }
            }
            scores.append(Score(j: j, matched: matched, rms: matched > 1 ? (squares / Double(matched - 1)).squareRoot() : 0, secondShell: second,
                                shellPerRing: shellPerRing))
        }
        // A misfit below 0.2 % is rounding, not a difference between assignments.
        func effective(_ s: Score) -> Double { max(s.rms, 0.002) }
        let ranked = scores.sorted { $0.matched != $1.matched ? $0.matched > $1.matched : effective($0) < effective($1) }
        let best = ranked[0]
        if best.matched == 1 {
            let seen = (1 ..< observed.count).prefix(2).map { observed[$0] / r1 }
            let predicted = scores.sorted { $0.j < $1.j }.prefix(2).map { s in
                (label: name(s.j), ratios: (1 ... 2).compactMap { k in s.j + k < shells.count ? shells[s.j + k] / shells[s.j] : nil })
            }
            throw LinearScaleAssignmentError.noMatch(crystal: crystalName, observedRatios: Array(seen), predicted: predicted)
        }
        // Equal match counts are told apart only by a misfit that is both half and 0.5 % lower: ellipticity or an origin offset moves every
        // ring by about that much, so a smaller difference is noise.
        func beaten(_ s: Score) -> Bool {
            best.matched > s.matched
                || (best.matched == s.matched && effective(best) < 0.5 * effective(s) && effective(s) - effective(best) > 0.005)
        }
        // The same rings matched to shells at the same ratios: the two hypotheses are one prediction at two scales, and no radius of these
        // rings separates them.
        func sameRatios(_ s: Score) -> Bool {
            zip(best.shellPerRing, s.shellPerRing).allSatisfy { mine, theirs in
                switch (mine, theirs) {
                case (nil, nil): return true
                case let (k?, m?):
                    let p = shells[k] / shells[best.j], q = shells[m] / shells[s.j]
                    return abs(p - q) <= assignmentSameRatios * p
                default: return false
                }
            }
        }
        let rivals = ranked.dropFirst().filter { !beaten($0) }
        guard let other = rivals.first(where: { !sameRatios($0) }) else {
            // No rival, or only the same prediction at a higher scale: the lowest shell that fits is the first ring.
            let tied = ([best] + rivals).sorted { $0.j < $1.j }
            let lowest = tied[0]
            return assigned(lowest.j, matched: lowest.matched, rms: lowest.rms, secondShell: lowest.secondShell,
                            equallyFitting: tied.dropFirst().map { name($0.j) })
        }
        let (a, b) = best.j < other.j ? (best.j, other.j) : (other.j, best.j)
        throw LinearScaleAssignmentError.ambiguous(first: name(a), second: name(b))
    }

    /// `assign`, then `estimate` against the assigned shell, the shell the second ring matched as the check.
    public static func estimateAssigned(
        radiiPerPosition: [[Double]], shellLengths: [Double], labels: [String]? = nil, crystalName: String = "the crystal",
        priorScale: Double? = nil, minimumRadiusPixels: Double = 2
    ) throws -> (estimate: LinearScaleEstimate, assignment: LinearScaleAssignment)? {
        guard let a = try assign(radiiPerPosition: radiiPerPosition, shellLengths: shellLengths, labels: labels, crystalName: crystalName,
                                 priorScale: priorScale, minimumRadiusPixels: minimumRadiusPixels),
              let e = estimate(radiiPerPosition: radiiPerPosition, referenceRadiusInvAngstrom: a.referenceRadiusInvAngstrom,
                               secondShellRadiusInvAngstrom: a.nextShellInvAngstrom, minimumRadiusPixels: minimumRadiusPixels) else { return nil }
        return (e, a)
    }

    /// "{2 0 0}" is written "{200}" when every index is a single digit; anything else is kept.
    package static func compactShellLabel(_ name: String) -> String {
        let inner = name.dropFirst().dropLast()
        let indices = inner.split(separator: " ")
        guard name.hasPrefix("{"), name.hasSuffix("}"), !name.contains("+"), indices.allSatisfy({ $0.count == 1 && $0.first!.isNumber }) else {
            return name
        }
        return "{" + indices.joined() + "}"
    }
}
