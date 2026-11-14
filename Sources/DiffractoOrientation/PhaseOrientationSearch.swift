//
//  PhaseOrientationSearch.swift
//  Role: which phase explains the spots the matrix did not, at which orientation, found from the structures alone. Step 2 of
//        "find phase orientations from the structure alone" (DiffractoKit#81).
//
//  How: for the residual spots of one scan position, each candidate phase's templates (an `ACOMMatcher` over the phase's
//  fundamental zone) PROPOSE a few distinct orientations. The phase is then DECIDED by the vector score at those orientations,
//  the known-variants score of `PhaseVectorMatcher.classifyKnownVariants`: every surviving vector's distance to the nearest
//  predicted spot, summed and divided by the number of distinct spots hit, in Å⁻¹. The template score never decides the
//  phase (its radial bins are wider than the difference between neighbouring phases; it once called aluminium gold). A CHANCE
//  GUARD then asks whether the matches could be an accident among all the orientations that were tried.
//
//  NOT A PHASE DECIDER ON ITS OWN. A free search over the whole orientation space is right where a position carries much
//  evidence and cannot be right where it carries little: with two to five matched vectors (most precipitate positions of a
//  real scan) either the guard refuses a true precipitate, or, without the guard, the phase with the largest cell explains
//  everything. Measured 2026-11-11 on Thronsen et al. dataset A (stride 3, θ′ and T1 from their structures alone, θ′ edge-on
//  and face-on merged): 81.6 % of positions equal the truth with the guard and 94.5 % without, against 98.7 % for the typed
//  zone axes of `PhaseVectorMatcher.map`. On a planted scan with rich patterns it finds every variant and refuses a phase
//  that is not in the list. Use it as the first pass of a two-pass classification (strong positions teach the orientation
//  relationship, the rest is classified against its variants), not as the map.
//
//  This is how template-matching phase mapping works in ASTAR, pyxem and py4DSTEM in outline (templates over the fundamental
//  zone, then a second number); no Python call does these steps as one, and the decision rule is this kit's own. Checked
//  against planted truth (`PhaseOrientationSearchTests`).
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// A phase to look for: its templates and its reflections.
public struct PhaseOrientationCandidate: Sendable {
    public let name: String
    /// A matcher on this phase's template library (its fundamental zone).
    public let matcher: ACOMMatcher
    /// The phase's reflections; −g must be there for every g (a kinematic list has both).
    public let reflections: [ACOMReflection]
    /// Half the thickness of the slab of reciprocal space that an orientation brings into the pattern, Å⁻¹.
    public var excitationSlabInvAngstrom: Double

    public init(name: String, matcher: ACOMMatcher, reflections: [ACOMReflection], excitationSlabInvAngstrom: Double = 0.05) {
        self.name = name
        self.matcher = matcher
        self.reflections = reflections
        self.excitationSlabInvAngstrom = excitationSlabInvAngstrom
    }

    /// The candidate of a crystal: its kinematic reflections up to `kMax` and a template library over its fundamental zone.
    ///
    /// - Throws: ``ACOMError`` (the structure refused, or its Laue class undetermined).
    public init(name: String, crystal: CrystalStructure, kMax: Double,
                templateSettings: ACOMTemplateSettings = ACOMTemplateSettings(),
                excitationSlabInvAngstrom: Double = 0.05) throws(ACOMError) {
        let reflections = try Self.reflections(of: crystal, kMax: kMax)
        let library = try ACOMTemplateLibrary.generate(crystal: crystal, kMax: kMax, settings: templateSettings)
        guard let matcher = ACOMMatcher(library: library) else { throw .invalidSettings("no FFT plan for the library's azimuthal length") }
        self.init(name: name, matcher: matcher, reflections: reflections, excitationSlabInvAngstrom: excitationSlabInvAngstrom)
    }

    /// The reflections a candidate of `crystal` holds: its kinematic reflections up to `kMax`, each with its intensity. For a
    /// caller that predicts spots and needs no template library
    /// (``PhaseVectorMatcher/indistinguishableVariants(_:matrixOrientation:reflections:excitationSlabInvAngstrom:reachInvAngstrom:toleranceInvAngstrom:directBeamRadiusInvAngstrom:)``).
    public static func reflections(of crystal: CrystalStructure, kMax: Double) throws(ACOMError) -> [ACOMReflection] {
        let list: KinematicReflectionList
        do { list = try crystal.kinematicReflections(kMax: kMax, retainingAll: false) }
        catch { throw .kinematics(error.description) }
        return list.kept.map { ACOMReflection(g: $0.reflection.g, intensity: $0.structureFactor.intensity) }
    }
}

/// The settings of ``PhaseVectorMatcher/searchOrientations(peaks:removedCount:candidates:settings:originX:originY:invAngstromPerPixel:)``.
/// The pair radius, the cut-off, the direct-matrix count, the direct-beam radius and the maximum vector mean what they mean
/// in `PhaseVectorSettings`, with the same defaults.
public struct PhaseOrientationSearchSettings: Sendable, Equatable {
    /// How many distinct orientations each phase's templates propose per position.
    public var proposalsPerPhase: Int
    /// A vector within this distance of a predicted spot is matched, Å⁻¹ (`PhaseVectorSettings.pairRadiusInvAngstrom`).
    public var pairRadiusInvAngstrom: Double
    /// A winner whose score is above this is not indexed, Å⁻¹ (`PhaseVectorSettings.residualCutoffInvAngstrom`).
    public var residualCutoffInvAngstrom: Double
    /// With at most this many vectors left there is nothing to search (`PhaseVectorSettings.directMatrixMaximumVectors`).
    public var directMatrixMaximumVectors: Int
    /// Vectors inside it are the direct beam and are not searched, Å⁻¹.
    public var directBeamRadiusInvAngstrom: Double
    /// 0 uses every vector.
    public var maximumVectorInvAngstrom: Double
    /// The expectation of orientations, among all tried, that match as many vectors by chance must be BELOW this. `.infinity`
    /// switches the guard off; a value that is not a positive number refuses everything.
    public var chanceLevel: Double
    /// The radius of the disc of the pattern in which a spot can be recorded (the detector's reach from the origin), Å⁻¹:
    /// the area a stray vector could have fallen anywhere in.
    public var reachInvAngstrom: Double

    public init(reachInvAngstrom: Double, proposalsPerPhase: Int = 3, pairRadiusInvAngstrom: Double = 0.02,
                residualCutoffInvAngstrom: Double = 0.07, directMatrixMaximumVectors: Int = 1,
                directBeamRadiusInvAngstrom: Double = 0.15, maximumVectorInvAngstrom: Double = 0, chanceLevel: Double = 0.05) {
        self.reachInvAngstrom = reachInvAngstrom
        self.proposalsPerPhase = proposalsPerPhase
        self.pairRadiusInvAngstrom = pairRadiusInvAngstrom
        self.residualCutoffInvAngstrom = residualCutoffInvAngstrom
        self.directMatrixMaximumVectors = directMatrixMaximumVectors
        self.directBeamRadiusInvAngstrom = directBeamRadiusInvAngstrom
        self.maximumVectorInvAngstrom = maximumVectorInvAngstrom
        self.chanceLevel = chanceLevel
    }
}

/// What the search made of one position.
public struct PhaseOrientationVerdict: Sendable, Equatable {
    /// `.indexed`: a candidate won within the cut-off and the chance guard. `.notIndexed`: vectors are left and nothing
    /// explains them well enough. `.matrix`: too few vectors are left to search and the matrix took spots here.
    /// `.noData`: no usable vector and nothing removed.
    public var verdict: PhaseVerdict = .noData
    /// The winner's index in the candidates, also for a refused winner (its numbers are kept for diagnosis); −1 when no
    /// candidate was scored.
    public var phaseIndex: Int = -1
    /// The winner's orientation with the in-plane fit applied (columns lab x, y, z in crystal coordinates).
    public var orientation: simd_double3x3?
    /// The winner's vector score, Å⁻¹; NaN when none was scored.
    public var score: Double = .nan
    /// Vectors within the pair radius of a predicted spot of the winner.
    public var matchedCount: Int = 0
    /// Usable vectors searched.
    public var survivingCount: Int = 0
    /// The winner's template correlation (what proposed it; it decides nothing).
    public var templateScore: Float = 0
    /// Orientations the templates stand for, over all candidates: templates × azimuthal bins × 2 (the mirror pass).
    public var orientationsTried: Int = 0
    /// The fewest matched vectors the chance guard accepts here (0 with the guard off).
    public var requiredMatches: Int = 0
    /// The best other phase, and its score; −1 and NaN when there is none.
    public var runnerUpPhaseIndex: Int = -1
    public var runnerUpScore: Double = .nan

    public init() {}

    // Every stored property, with NaN equal to NaN (an unscored position equals itself). A new property must be added here.
    public static func == (a: Self, b: Self) -> Bool {
        func same(_ x: Double, _ y: Double) -> Bool { x == y || (x.isNaN && y.isNaN) }
        return a.verdict == b.verdict && a.phaseIndex == b.phaseIndex && a.orientation == b.orientation && same(a.score, b.score)
            && a.matchedCount == b.matchedCount && a.survivingCount == b.survivingCount && a.templateScore == b.templateScore
            && a.orientationsTried == b.orientationsTried && a.requiredMatches == b.requiredMatches
            && a.runnerUpPhaseIndex == b.runnerUpPhaseIndex && same(a.runnerUpScore, b.runnerUpScore)
    }
}

extension PhaseVectorMatcher {

    /// The in-plane rotation (2D least squares, twice, over the pairs within three tolerances) that lays predicted spots on
    /// measured vectors, at most ``matrixInPlaneCorrectionLimitRad``: the fit of
    /// ``removingMatrixSpots(scan:field:reflections:originX:originY:invAngstromPerPixel:toleranceInvAngstrom:excitationSlabInvAngstrom:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:)``.
    static func fittingInPlane(spots input: [SIMD2<Double>], vectors: [SIMD2<Double>], toleranceInvAngstrom: Double)
        -> (spots: [SIMD2<Double>], rotationRad: Double) {
        var spots = input
        var total = 0.0
        for _ in 0..<2 {
            var cross = 0.0, dot = 0.0, pairs = 0
            for u in vectors {
                var best: SIMD2<Double>?, bestDistance = Double.infinity
                for s in spots { let d = simd_length(u - s); if d < bestDistance { bestDistance = d; best = s } }
                guard let n = best, bestDistance <= 3 * toleranceInvAngstrom else { continue }
                cross += n.x * u.y - n.y * u.x
                dot += n.x * u.x + n.y * u.y
                pairs += 1
            }
            guard pairs >= 2, dot > 0 else { break }
            let delta = atan2(cross, dot)
            guard abs(total + delta) <= matrixInPlaneCorrectionLimitRad else { break }
            total += delta
            let c = cos(delta), s = sin(delta)
            spots = spots.map { SIMD2(c * $0.x - s * $0.y, s * $0.x + c * $0.y) }
        }
        return (spots, total)
    }

    /// The fewest matches that are unlikely to be chance among `orientationsTried` orientations: the smallest k for which
    /// orientationsTried × P(at least k of `vectors` stray vectors land on a spot) is below `chanceLevel`, where one stray
    /// vector lands on one of `spots` spots with probability spots × (pairRadius / reach)² (the area of the pair discs over
    /// the area of the pattern, at most 1). A union bound: the orientations are not independent, so it errs on the side of
    /// asking for more. `vectors + 1` when no count is enough (also for a level that is not a positive number); 0 with the
    /// guard off (`chanceLevel` = +infinity). The spot count is the winning orientation's, the orientations tried are all
    /// candidates'.
    public static func requiredMatches(vectors m: Int, spots n: Int, orientationsTried: Int,
                                       settings: PhaseOrientationSearchSettings) -> Int {
        if settings.chanceLevel == .infinity { return 0 }
        guard m > 0, n > 0, orientationsTried > 0, settings.reachInvAngstrom > 0, settings.chanceLevel > 0 else { return m + 1 }
        let ratio = settings.pairRadiusInvAngstrom / settings.reachInvAngstrom
        let p = min(1, Double(n) * ratio * ratio)
        guard p > 0 else { return 1 }       // nothing can match by chance: one match is one more than chance gives
        // P(Binomial(m, p) >= k), from the top down.
        var tail = 0.0
        var required = m + 1
        for k in stride(from: m, through: 0, by: -1) {
            var logChoose = 0.0
            if k > 0 { for j in 1...k { logChoose += log(Double(m - k + j)) - log(Double(j)) } }
            let term = p >= 1 ? (k == m ? 1.0 : 0.0)
                : exp(logChoose + Double(k) * log(p) + Double(m - k) * log(1 - p))
            tail += term
            if Double(orientationsTried) * tail < settings.chanceLevel { required = k } else { break }
        }
        return required
    }

    /// Which candidate explains the residual spots of one position, and at which orientation.
    ///
    /// - Parameters:
    ///   - peaks: the position's peaks the matrix did not explain (``MatrixResidualScan/scan``), detector pixels.
    ///   - removedCount: how many spots the matrix took at this position (``MatrixResidualScan/removedCounts``).
    ///   - candidates: the phases to look for; the result's `phaseIndex` is an index into this array.
    public static func searchOrientations(peaks: [SIMD2<Double>], removedCount: Int,
                                          candidates: [PhaseOrientationCandidate],
                                          settings: PhaseOrientationSearchSettings,
                                          originX: Float, originY: Float, invAngstromPerPixel: Double) -> PhaseOrientationVerdict {
        var out = PhaseOrientationVerdict()
        var usablePeaks: [SIMD2<Double>] = [], vectors: [SIMD2<Double>] = []
        for p in peaks {
            if let u = experimentalVector(peak: p, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                          directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                                          maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom) {
                usablePeaks.append(p); vectors.append(u)
            }
        }
        out.survivingCount = vectors.count
        out.orientationsTried = candidates.reduce(0) { $0 + $1.matcher.library.count * $1.matcher.library.geometry.nAzimuthal * 2 }
        if vectors.count <= settings.directMatrixMaximumVectors {
            out.verdict = removedCount > 0 ? .matrix : (vectors.isEmpty ? .noData : .notIndexed)
            return out
        }
        let acomPeaks = usablePeaks.map { ACOMPeak(x: Float($0.x), y: Float($0.y), intensity: 1) }

        struct Scored { var phase: Int; var orientation: simd_double3x3; var score: Double; var matched: Int; var spots: Int; var template: Float }
        var bestPerPhase: [Scored] = []
        for (phase, candidate) in candidates.enumerated() {
            let library = candidate.matcher.library
            var best: Scored?
            for proposal in candidate.matcher.matches(peaks: acomPeaks, originX: originX, originY: originY,
                                                      invAngstromPerPixel: invAngstromPerPixel, count: settings.proposalsPerPhase) {
                let m = library.orientationMatrix(templateIndex: proposal.templateIndex, inPlaneAngle: proposal.inPlaneAngle,
                                                  mirrored: proposal.mirrored)
                let predicted = matrixSpots(reflections: candidate.reflections, orientation: m,
                                            excitationSlabInvAngstrom: candidate.excitationSlabInvAngstrom)
                guard !predicted.isEmpty else { continue }
                let fitted = fittingInPlane(spots: predicted, vectors: vectors, toleranceInvAngstrom: settings.pairRadiusInvAngstrom)
                // The known-variants score: every vector's distance to its true nearest spot, over the distinct spots hit.
                var sum = 0.0, matched = 0
                var hits = Set<Int>()
                for u in vectors {
                    var nearest = 0, distance = Double.infinity
                    for (k, s) in fitted.spots.enumerated() { let d = simd_length(u - s); if d < distance { distance = d; nearest = k } }
                    sum += distance
                    hits.insert(nearest)
                    if distance <= settings.pairRadiusInvAngstrom { matched += 1 }
                }
                let score = sum / Double(hits.count)
                if best == nil || score < best!.score {
                    // The orientation with the fit applied: the spots turned by +δ are g·(cos δ x − sin δ y), g·(sin δ x + cos δ y).
                    let c = cos(fitted.rotationRad), s = sin(fitted.rotationRad)
                    let turned = simd_double3x3(columns: (c * m.columns.0 - s * m.columns.1, s * m.columns.0 + c * m.columns.1, m.columns.2))
                    best = Scored(phase: phase, orientation: turned, score: score, matched: matched, spots: fitted.spots.count,
                                  template: proposal.score)
                }
            }
            if let best { bestPerPhase.append(best) }
        }
        guard var winner = bestPerPhase.first else { out.verdict = .notIndexed; return out }
        for s in bestPerPhase.dropFirst() where s.score < winner.score { winner = s }
        var runnerUp: Scored?
        for s in bestPerPhase where s.phase != winner.phase { if runnerUp == nil || s.score < runnerUp!.score { runnerUp = s } }

        out.phaseIndex = winner.phase
        out.orientation = winner.orientation
        out.score = winner.score
        out.matchedCount = winner.matched
        out.templateScore = winner.template
        if let runnerUp { out.runnerUpPhaseIndex = runnerUp.phase; out.runnerUpScore = runnerUp.score }
        out.requiredMatches = requiredMatches(vectors: vectors.count, spots: winner.spots, orientationsTried: out.orientationsTried,
                                              settings: settings)
        out.verdict = winner.score <= settings.residualCutoffInvAngstrom && winner.matched >= out.requiredMatches ? .indexed : .notIndexed
        return out
    }

    /// ``searchOrientations(peaks:removedCount:candidates:settings:originX:originY:invAngstromPerPixel:)`` at every position of
    /// a residual scan, on worker threads. Each position depends only on its own peaks, so the worker count moves no number.
    ///
    /// - Returns: nil when the residual scan's two arrays do not have the same length, or the run was cancelled.
    public static func searchOrientations(residual: MatrixResidualScan, candidates: [PhaseOrientationCandidate],
                                          settings: PhaseOrientationSearchSettings,
                                          originX: Float, originY: Float, invAngstromPerPixel: Double,
                                          maximumWorkers: Int? = nil,
                                          isCancelled: (@Sendable () -> Bool)? = nil,
                                          progress: (@Sendable (Double) -> Void)? = nil) -> [PhaseOrientationVerdict]? {
        let total = residual.scan.count
        guard residual.removedCounts.count == total else { return nil }
        guard total > 0 else { return [] }
        var results = [PhaseOrientationVerdict](repeating: PhaseOrientationVerdict(), count: total)
        results.withUnsafeMutableBufferPointer { buffer in
            let box = VerdictBox(buffer.baseAddress!)
            let workers = max(1, min(maximumWorkers ?? ProcessInfo.processInfo.activeProcessorCount, total))
            // The positions done so far, reported about every 1 % of the scan (mac4DSTEM#333: this step ran for minutes unreported).
            let counter = PhaseVectorMatcher.ProgressCounter()
            let batch = max(1, total / (workers * 100))
            DispatchQueue.concurrentPerform(iterations: workers) { worker in
                var pending = 0
                for position in Swift.stride(from: worker, to: total, by: workers) {
                    if isCancelled?() == true { break }
                    box.value[position] = searchOrientations(
                        peaks: residual.scan[position], removedCount: residual.removedCounts[position], candidates: candidates,
                        settings: settings, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel)
                    pending += 1
                    if let progress, pending >= batch { counter.add(pending, of: total, reporting: progress); pending = 0 }
                }
                if let progress, pending > 0 { counter.add(pending, of: total, reporting: progress) }
            }
        }
        guard isCancelled?() != true else { return nil }
        return results
    }

    private struct VerdictBox: @unchecked Sendable {
        let value: UnsafeMutablePointer<PhaseOrientationVerdict>
        init(_ v: UnsafeMutablePointer<PhaseOrientationVerdict>) { value = v }
    }
}
