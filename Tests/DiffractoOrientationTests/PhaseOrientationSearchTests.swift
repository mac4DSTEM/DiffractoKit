//
//  PhaseOrientationSearchTests.swift
//  Step 2 of "find phase orientations from the structure alone" (DiffractoKit#81): for the spots the matrix did not explain, each
//  candidate phase's templates propose a few orientations, the vector score decides the phase at those orientations, and a
//  chance guard counts the orientations tried. Planted truth: the patterns are the reference vectors of `PhaseReferenceLibrary`
//  entries (another projection code than the templates'), at an in-plane angle no template bin holds exactly.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct PhaseOrientationSearchTests {
    static let q = 0.012
    static let origin: Float = 63.5
    static let reach = 0.75

    /// The planted pattern of a crystal down a lattice zone axis at 34 degrees in-plane: the phase library's own entry.
    static func plantedVectors(_ crystal: CrystalStructure, zone: SIMD3<Int>) throws -> [SIMD2<Double>] {
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = reach
        reference.minimumIntensityFraction = 0
        let library = try PhaseReferenceLibrary.build(phases: [
            PhaseDefinition(id: "m", displayName: "matrix", crystal: CrystalPresets.aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "c", displayName: "candidate", crystal: crystal, role: .candidate, zoneAxes: [zone]),
        ], settings: reference)
        let entry = try #require(library.entries.first { $0.phaseIndex == 1 && abs($0.inPlaneRotationRad - 34 * .pi / 180) < 1e-9 })
        return entry.vectors.map(\.q).filter { simd_length($0) > 0.15 && simd_length($0) < reach }
    }

    static func peaks(_ vectors: [SIMD2<Double>]) -> [SIMD2<Double>] {
        vectors.map { SIMD2(Double(origin) + $0.x / q, Double(origin) + $0.y / q) }
    }

    static func candidates() throws -> [PhaseOrientationCandidate] {
        [try PhaseOrientationCandidate(name: "beta''", crystal: CrystalPresets.betaDoublePrime, kMax: reach),
         try PhaseOrientationCandidate(name: "Mg", crystal: CrystalPresets.magnesium, kMax: reach)]
    }

    static func settings() -> PhaseOrientationSearchSettings {
        PhaseOrientationSearchSettings(reachInvAngstrom: reach)
    }

    static func search(_ vectors: [SIMD2<Double>], removed: Int = 0, _ candidates: [PhaseOrientationCandidate],
                       _ settings: PhaseOrientationSearchSettings) -> PhaseOrientationVerdict {
        PhaseVectorMatcher.searchOrientations(peaks: peaks(vectors), removedCount: removed, candidates: candidates, settings: settings,
                                              originX: origin, originY: origin, invAngstromPerPixel: q)
    }

    /// A phase is found from its structure alone, at the orientation it has: the end-on and the in-plane pattern of beta'' are
    /// both indexed as beta'', every planted vector is matched, and the orientation returned predicts the planted spots.
    /// Mutations: only the best template tried without the in-plane fit (the score rises above the pair radius at 34 degrees);
    /// the orientation returned without the in-plane correction (the last expectation).
    @Test(arguments: [SIMD3<Int>(0, 1, 0), SIMD3<Int>(0, 0, 1)])
    func aPhaseIsFoundAtItsOwnOrientation(zone: SIMD3<Int>) throws {
        let candidates = try Self.candidates()
        // Turned 1.2 degrees off the 34 degree entry: 35.2 degrees is half an azimuthal bin (2.8 degrees at 128 bins) from the
        // nearest bin, so without the in-plane fit a spot at 0.7 Å⁻¹ is 0.017 Å⁻¹ off.
        let turn = 1.2 * Double.pi / 180
        let planted = try Self.plantedVectors(CrystalPresets.betaDoublePrime, zone: zone).map {
            SIMD2(cos(turn) * $0.x - sin(turn) * $0.y, sin(turn) * $0.x + cos(turn) * $0.y)
        }
        #expect(planted.count >= 6)
        let v = Self.search(planted, candidates, Self.settings())
        #expect(v.verdict == .indexed, "\(v)")
        #expect(v.phaseIndex == 0)
        #expect(v.survivingCount == planted.count && v.matchedCount == planted.count, "matched \(v.matchedCount) of \(planted.count)")
        #expect(v.score < 0.003, "score \(v.score)")
        #expect(v.matchedCount >= v.requiredMatches)
        let orientation = try #require(v.orientation)
        let spots = PhaseVectorMatcher.matrixSpots(reflections: candidates[0].reflections, orientation: orientation,
                                                   excitationSlabInvAngstrom: candidates[0].excitationSlabInvAngstrom)
        for u in planted {
            #expect(spots.contains { simd_length($0 - u) < 0.003 }, "the orientation returned does not predict the planted \(u)")
        }
    }

    /// The vector score decides between phases: a magnesium [0001] pattern is magnesium, and beta'' is the runner-up with the
    /// worse score. Mutation: the template score decides (not pinned by this case alone; see the kit's note on why it must not).
    @Test func theVectorScoreDecidesThePhase() throws {
        let candidates = try Self.candidates()
        let planted = try Self.plantedVectors(CrystalPresets.magnesium, zone: SIMD3(0, 0, 1))
        let v = Self.search(planted, candidates, Self.settings())
        #expect(v.verdict == .indexed && v.phaseIndex == 1, "\(v)")
        #expect(v.runnerUpPhaseIndex == 0)
        #expect(v.runnerUpScore > v.score)
    }

    /// The chance guard counts the orientations tried. Three true beta'' spots fit a beta'' orientation perfectly, and three
    /// vectors are too few to tell that from chance among the orientations searched: not indexed, with the count the guard
    /// asked for. With the guard switched off the same three are indexed, which is what the guard is there to prevent.
    /// Mutation: the guard computed for one orientation instead of all tried (three matches then pass).
    @Test func threeMatchingVectorsAreNotEnoughAmongAllTheOrientationsTried() throws {
        let candidates = try Self.candidates()
        let planted = Array(try Self.plantedVectors(CrystalPresets.betaDoublePrime, zone: SIMD3(0, 1, 0)).prefix(3))
        let v = Self.search(planted, candidates, Self.settings())
        let tried = candidates.reduce(0) { $0 + $1.matcher.library.count * $1.matcher.library.geometry.nAzimuthal * 2 }
        #expect(v.orientationsTried == tried)
        #expect(v.verdict == .notIndexed, "\(v)")
        #expect(v.requiredMatches > 3)
        var off = Self.settings()
        off.chanceLevel = .infinity
        let unguarded = Self.search(planted, candidates, off)
        #expect(unguarded.verdict == .indexed && unguarded.matchedCount == 3 && unguarded.requiredMatches == 0, "\(unguarded)")
    }

    /// The guard's count against the binomial tail written out. One stray vector lands on a spot with p = 2 x (0.1 / 1)^2 =
    /// 0.02; of 10 vectors, P(at least 5) = 7.3e-7 and P(at least 6) = 1.2e-8. Among 100 000 orientations the first is 0.073
    /// (not below 0.05) and the second 0.0012: six are asked for. With one orientation P(at least 2) = 0.016 is enough.
    @Test func theChanceGuardIsTheBinomialTailTimesTheOrientationsTried() {
        var s = PhaseOrientationSearchSettings(reachInvAngstrom: 1, pairRadiusInvAngstrom: 0.1)
        func tail(_ m: Int, _ p: Double, _ k: Int) -> Double {
            (k...m).reduce(0.0) { sum, j in
                var choose = 1.0
                for i in 0..<j { choose *= Double(m - i) / Double(i + 1) }
                return sum + choose * pow(p, Double(j)) * pow(1 - p, Double(m - j))
            }
        }
        #expect(100_000 * tail(10, 0.02, 5) > 0.05 && 100_000 * tail(10, 0.02, 6) < 0.05)
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 10, spots: 2, orientationsTried: 100_000, settings: s) == 6)
        #expect(tail(10, 0.02, 1) > 0.05 && tail(10, 0.02, 2) < 0.05)
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 10, spots: 2, orientationsTried: 1, settings: s) == 2)
        // No count is enough: more is asked for than there are vectors.
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 2, spots: 2, orientationsTried: 100_000_000, settings: s) == 3)
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 0, spots: 2, orientationsTried: 10, settings: s) == 1)
        // The guard is off for +infinity only; a level that is not a positive number refuses.
        s.chanceLevel = .infinity
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 10, spots: 2, orientationsTried: 100_000, settings: s) == 0)
        for level in [Double.nan, 0.0, -1.0, -Double.infinity] {
            s.chanceLevel = level
            #expect(PhaseVectorMatcher.requiredMatches(vectors: 10, spots: 2, orientationsTried: 100_000, settings: s) == 11)
        }
        // A pair radius of 0: nothing matches by chance, one match is asked for.
        s = PhaseOrientationSearchSettings(reachInvAngstrom: 1, pairRadiusInvAngstrom: 0)
        #expect(PhaseVectorMatcher.requiredMatches(vectors: 10, spots: 2, orientationsTried: 100_000, settings: s) == 1)
    }

    /// The several distinct matches and the match near a zone, directly: best first and equal to `match`, pairwise farther
    /// apart than the library's distinct angle, an over-long request ends, too few peaks give nothing; a zone with no
    /// template near it gives the empty result.
    @Test func severalDistinctMatchesAndTheMatchNearAZone() throws {
        let matcher = try Self.candidates()[0].matcher
        let library = matcher.library
        let peaks = Self.peaks(try Self.plantedVectors(CrystalPresets.betaDoublePrime, zone: SIMD3(0, 1, 0)))
            .map { ACOMPeak(x: Float($0.x), y: Float($0.y), intensity: 1) }
        let single = matcher.match(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        let five = matcher.matches(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, count: 5)
        #expect(five.count == 5)
        #expect(five[0].templateIndex == single.templateIndex && five[0].score == single.score && five[0].inPlaneAngle == single.inPlaneAngle)
        for i in 1..<five.count { #expect(five[i].score <= five[i - 1].score) }
        let limit = cos(library.distinctOrientationRadians)
        for i in 0..<five.count { for j in (i + 1)..<five.count {
            #expect(simd_dot(library.zoneAxes[five[i].templateIndex], library.zoneAxes[five[j].templateIndex]) <= limit)
        } }
        let all = matcher.matches(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, count: 100_000)
        #expect(all.count >= 5 && all.count < library.count)
        #expect(matcher.matches(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, count: 0).isEmpty)
        #expect(matcher.matches(peaks: Array(peaks.prefix(1)), originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, count: 3).isEmpty)

        let zone = library.zoneAxes[single.templateIndex]
        let near = matcher.match(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, near: zone, withinRadians: 1e-3)
        #expect(near.templateIndex == single.templateIndex && near.score == single.score)
        let nowhere = matcher.match(peaks: peaks, originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                                    near: -zone, withinRadians: 1e-3)
        #expect(nowhere == .empty)
    }

    /// Too few vectors to search: the matrix when the matrix took spots here, no data when there was nothing, and not
    /// indexed for a lone vector the matrix did not explain.
    @Test func tooFewVectors() throws {
        let candidates = try Self.candidates()
        let lone = [SIMD2(0.31, 0.05)]
        #expect(Self.search([], removed: 8, candidates, Self.settings()).verdict == .matrix)
        #expect(Self.search(lone, removed: 8, candidates, Self.settings()).verdict == .matrix)
        #expect(Self.search([], removed: 0, candidates, Self.settings()).verdict == .noData)
        #expect(Self.search(lone, removed: 0, candidates, Self.settings()).verdict == .notIndexed)
        #expect(Self.search(try Self.plantedVectors(CrystalPresets.betaDoublePrime, zone: SIMD3(0, 1, 0)), [], Self.settings()).verdict == .notIndexed)
    }

    /// The scan form is the single-position form at every position, whatever the worker count.
    @Test func theScanFormEqualsThePositionForm() throws {
        let candidates = try Self.candidates()
        let a = Self.peaks(try Self.plantedVectors(CrystalPresets.betaDoublePrime, zone: SIMD3(0, 1, 0)))
        let b = Self.peaks(try Self.plantedVectors(CrystalPresets.magnesium, zone: SIMD3(0, 0, 1)))
        let residual = MatrixResidualScan(scan: [a, b, [], a], removedCounts: [0, 8, 8, 2], inPlaneCorrectionsRad: [0, 0, 0, 0])
        let one = try #require(PhaseVectorMatcher.searchOrientations(
            residual: residual, candidates: candidates, settings: Self.settings(),
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, maximumWorkers: 1))
        let four = try #require(PhaseVectorMatcher.searchOrientations(
            residual: residual, candidates: candidates, settings: Self.settings(),
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, maximumWorkers: 4))
        #expect(one == four)
        #expect(one.map(\.verdict) == [.indexed, .indexed, .matrix, .indexed])
        #expect(one.map(\.phaseIndex) == [0, 1, -1, 0])
        for (i, peaks) in residual.scan.enumerated() {
            #expect(one[i] == PhaseVectorMatcher.searchOrientations(
                peaks: peaks, removedCount: residual.removedCounts[i], candidates: candidates, settings: Self.settings(),
                originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        }
    }
}
