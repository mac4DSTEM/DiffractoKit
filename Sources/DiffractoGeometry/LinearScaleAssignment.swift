//
//  LinearScaleAssignment.swift
//  Role: which of a crystal's shells is the first observed ring? `LinearScale.estimate` takes its reference shell from the caller; on Al [001]
//        the {111} reflections are not in the zone, the first visible ring is {200}, and assuming the first shell read Δq 13.4 % low.
//        `LinearScale.assign` tests the first few shells against the SEQUENCE of observed ring radii and refuses when the sequence cannot tell.
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

    public static func == (a: Self, b: Self) -> Bool {
        switch (a, b) {
        case let (.noMatch(c1, o1, p1), .noMatch(c2, o2, p2)):
            return c1 == c2 && o1 == o2 && p1.count == p2.count && zip(p1, p2).allSatisfy { $0.label == $1.label && $0.ratios == $1.ratios }
        case let (.ambiguous(a1, b1), .ambiguous(a2, b2)): return a1 == a2 && b1 == b2
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
        }
    }
}

extension LinearScale {
    /// Cluster width (relative) of the equivalents of one ring at one position.
    static let assignmentClusterBand = 0.03
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
    /// Per position the peaks are clustered (the equivalents of one ring), the first three rings are read, and each ring's radius is the median
    /// over the positions that have it (at least 20 % of those with any peak). The hypotheses H_j "the first ring is shell j" (j below 4) predict
    /// the later rings; each observed ring is matched to the nearest later shell within max(3 % of its radius, 1 px). The best hypothesis wins
    /// if it matches strictly more rings than the runner-up, or the same number with a misfit under half and 0.5 % lower. One observed ring
    /// cannot be tested: shell 0 is assumed, unless `priorScale` agrees with exactly one hypothesis within 10 % (it never overrides a sequence).
    ///
    /// - Parameters:
    ///   - radiiPerPosition: peak radii in pixels from the origin, one array per scan position.
    ///   - shellLengths: the crystal's distinct |g| (Å⁻¹), ascending.
    ///   - labels: one family name per shell for the sentences ("{200}"); "shell k" when absent.
    ///   - crystalName: for the "no match" sentence.
    ///   - priorScale: a Δq (Å⁻¹ per pixel) from the file or typed.
    /// - Returns: nil when no position holds a radius above `minimumRadiusPixels` (or the shell list is empty or has a non-positive length).
    /// - Throws: `LinearScaleAssignmentError` when the sequence matches no starting shell or two assignments fit equally well.
    public static func assign(
        radiiPerPosition: [[Double]], shellLengths: [Double], labels: [String]? = nil, crystalName: String = "the crystal",
        priorScale: Double? = nil, minimumRadiusPixels: Double = 2
    ) throws -> LinearScaleAssignment? {
        guard shellLengths.count > 0, shellLengths.allSatisfy({ $0.isFinite && $0 > 0 }) else { return nil }
        func name(_ k: Int) -> String { labels.flatMap { k < $0.count ? $0[k] : nil } ?? "shell \(k + 1)" }

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
        guard positionsWithPeaks > 0 else { return nil }
        let needed = max(1, Int((assignmentMinimumPositionShare * Double(positionsWithPeaks)).rounded(.up)))
        var observed: [Double] = []
        for ring in 0 ..< assignmentMaximumRings {
            guard perRing[ring].count >= needed else { break }
            observed.append(median(perRing[ring]))
        }
        guard let r1 = observed.first, r1 > 0 else { return nil }

        let shells = Array(shellLengths.prefix(assignmentShellWindow))
        let tried = min(assignmentHypotheses, shells.count)
        func assigned(_ j: Int, matched: Int, rms: Double, secondShell: Int?) -> LinearScaleAssignment {
            LinearScaleAssignment(
                shellIndex: j, label: name(j), referenceRadiusInvAngstrom: shells[j],
                nextShellInvAngstrom: secondShell.map { shells[$0] } ?? (j + 1 < shells.count ? shells[j + 1] : nil),
                observedRadiiPixels: observed, matchedRings: matched, observedRings: observed.count, rmsMisfit: rms)
        }

        if observed.count == 1 || shells.count == 1 {
            var j = 0
            if let prior = priorScale, prior.isFinite, prior > 0 {
                let agreeing = (0 ..< tried).filter { abs(shells[$0] / r1 / prior - 1) <= 0.10 }
                if agreeing.count == 1 { j = agreeing[0] }
            }
            return assigned(j, matched: 1, rms: 0, secondShell: nil)
        }

        struct Score { let j: Int; let matched: Int; let rms: Double; let secondShell: Int? }
        var scores: [Score] = []
        for j in 0 ..< tried {
            let scale = shells[j] / r1
            var last = j, matched = 1, squares = 0.0, second: Int?
            for ring in 1 ..< observed.count {
                let r = observed[ring]
                var best: (k: Int, error: Double)?
                for k in (last + 1) ..< max(shells.count, last + 1) {
                    let error = abs(shells[k] / scale - r) / r
                    if best == nil || error < best!.error { best = (k, error) }
                }
                guard let hit = best, abs(shells[hit.k] / scale - r) <= max(0.03 * r, 1) else { continue }
                matched += 1; squares += hit.error * hit.error; last = hit.k
                if ring == 1 { second = hit.k }
            }
            scores.append(Score(j: j, matched: matched, rms: matched > 1 ? (squares / Double(matched - 1)).squareRoot() : 0, secondShell: second))
        }
        // A misfit below 0.2 % is rounding, not a difference between assignments.
        func effective(_ s: Score) -> Double { max(s.rms, 0.002) }
        let ranked = scores.sorted { $0.matched != $1.matched ? $0.matched > $1.matched : effective($0) < effective($1) }
        let best = ranked[0], runnerUp = ranked[1]
        if best.matched == 1 {
            let seen = (1 ..< observed.count).prefix(2).map { observed[$0] / r1 }
            let predicted = scores.sorted { $0.j < $1.j }.prefix(2).map { s in
                (label: name(s.j), ratios: (1 ... 2).compactMap { k in s.j + k < shells.count ? shells[s.j + k] / shells[s.j] : nil })
            }
            throw LinearScaleAssignmentError.noMatch(crystal: crystalName, observedRatios: Array(seen), predicted: predicted)
        }
        // Equal match counts are told apart only by a misfit that is both half and 0.5 % lower: ellipticity or an origin offset moves every
        // ring by about that much, so a smaller difference is noise.
        let decisive = best.matched > runnerUp.matched
            || (best.matched == runnerUp.matched && effective(best) < 0.5 * effective(runnerUp)
                && effective(runnerUp) - effective(best) > 0.005)
        guard decisive else {
            let (a, b) = best.j < runnerUp.j ? (best.j, runnerUp.j) : (runnerUp.j, best.j)
            throw LinearScaleAssignmentError.ambiguous(first: name(a), second: name(b))
        }
        return assigned(best.j, matched: best.matched, rms: best.rms, secondShell: best.secondShell)
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
