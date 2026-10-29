//
//  PhaseVectorMatcher+Classify.swift
//  Role: one pattern classified: the default rule, the known-variants rule (Thronsen et al.) and the orientation relationship check.
//  Split out of PhaseVectorMatching.swift, which carries the method, its sources and its limits.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

extension PhaseVectorMatcher {

    // MARK: One pattern

    /// Classify one pattern's already-calibrated vectors.
    public static func classify(vectors: [SIMD2<Double>],
                                 library: PhaseReferenceLibrary,
                                 settings: PhaseVectorSettings,
                                 matrixEntry: PhaseOrientationReference?,
                                 candidateEntryIndices: [Int],
                                 scratch: Scratch,
                                 matrixChallenge: [PhaseOrientationReference] = []) -> PhaseVectorResult {
        // Dispatch only, and only here: everything below this line is
        // `.search`, untouched by `.knownVariants` existing. The second rule
        // is a second function, not a branch threaded through this one, so
        // `.search`'s numbers cannot drift by way of a shared code path —
        // pinned by `testKnownVariantsRuleDoesNotChangeSearchsOwnResult`.
        if settings.classificationRule == .knownVariants {
            return classifyKnownVariants(vectors: vectors, library: library, settings: settings,
                                         matrixEntry: matrixEntry,
                                         candidateEntryIndices: candidateEntryIndices)
        }
        var result = PhaseVectorResult()
        guard !vectors.isEmpty else { return result }
        // Caller-supplied indices outside the library give no verdict (`map` checks the library's own consistency once).
        guard candidateEntryIndices.allSatisfy({ library.entries.indices.contains($0) }) else { return result }

        // 1 — matrix removal in vector space.
        var surviving: [SIMD2<Double>] = []
        surviving.reserveCapacity(vectors.count)
        if let matrixEntry {
            for u in vectors {
                if nearest(u, in: matrixEntry.vectors,
                           radius: settings.matrixToleranceInvAngstrom) == nil {
                    surviving.append(u)
                }
            }
        } else {
            surviving = vectors
        }
        result.survivingCount = Int32(surviving.count)
        result.removedCount = Int32(vectors.count - surviving.count)

        // 2 — the matrix, by exclusion.
        if surviving.count < settings.minimumVectors {
            guard Int(result.removedCount) >= settings.minimumMatrixEvidence else {
                result.verdict = .notIndexed
                return result
            }
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            return result
        }

        // 3 — score every candidate entry; keep the best per phase.
        //
        // The area chance is measured against is the one the DATA occupies,
        // not the library's kMax: if every surviving vector sits inside
        // 0.7 Å⁻¹, references beyond it cannot be hit by accident and must not
        // dilute the estimate.
        // The radius is the SECOND largest |u|, not the largest. Measured
        // 2026-09-12 (a review finding): with `max`, one spurious maximum far
        // out inflates the area the chance expectation is computed over,
        // weakens the guard for every other vector in the pattern, and took a
        // 1024-pattern false-positive rate from 0 to 10. Second-largest costs
        // nothing, is exact whenever the outermost peak is real, and cannot be
        // moved by a single outlier.
        var largest = 0.0, accessibleRadius = 0.0
        for u in surviving {
            let r = simd_length(u)
            if r > largest { accessibleRadius = largest; largest = r }
            else if r > accessibleRadius { accessibleRadius = r }
        }
        if accessibleRadius <= 0 { accessibleRadius = largest }
        var bestPerPhase: [Int: (entryIndex: Int, score: Double, matched: Int)] = [:]
        // The floor a candidate must clear: `minimumMatchedVectors`, or the
        // Friedel-pair floor when the survivors hold u and −u (see the
        // setting's note). Computed once per pattern, not per entry.
        let matchedFloor = containsFriedelPair(surviving, radius: settings.pairRadiusInvAngstrom)
            ? min(settings.minimumMatchedVectors, settings.friedelPairMinimumMatchedVectors)
            : settings.minimumMatchedVectors
        for entryIndex in candidateEntryIndices {
            let entry = library.entries[entryIndex]
            guard let s = score(vectors: surviving, against: entry,
                                pairRadius: settings.pairRadiusInvAngstrom,
                                scratch: scratch) else { continue }
            let chance = entry.chanceMatchFraction(
                pairRadius: settings.pairRadiusInvAngstrom,
                accessibleRadius: accessibleRadius) * Double(surviving.count)
            guard s.matched >= matchedFloor,
                  Double(s.matched) >= settings.chanceMatchMultiple * chance
            else { continue }
            let current = bestPerPhase[entry.phaseIndex]
            if current == nil || s.score < current!.score {
                bestPerPhase[entry.phaseIndex] = (entryIndex, s.score, s.matched)
            }
        }
        guard !bestPerPhase.isEmpty else {
            if settings.matrixFallbackExplainedFraction > 0, !vectors.isEmpty,
               Double(result.removedCount) / Double(vectors.count)
                   >= settings.matrixFallbackExplainedFraction {
                result.verdict = .matrix
                result.phaseIndex = Int32(library.matrixPhaseIndex)
                result.entryIndex = -1
                return result
            }
            result.verdict = .notIndexed
            return result
        }

        // An exact tie goes to the lower phase index (the order the library lists the phases), so the label never depends on
        // the Dictionary's per-process hash order; identical to the apps' original
        // wherever it is deterministic.
        let ranked = bestPerPhase.sorted {
            let a = (matched: $0.value.matched, score: $0.value.score), b = (matched: $1.value.matched, score: $1.value.score)
            let aWins = crossPhaseWinsOver(a, b, completenessAware: settings.completenessAwareCrossPhaseRanking)
            let bWins = crossPhaseWinsOver(b, a, completenessAware: settings.completenessAwareCrossPhaseRanking)
            return aWins != bWins ? aWins : $0.key < $1.key
        }
        let winner = ranked[0]
        result.score = Float(winner.value.score)
        result.matchedCount = Int32(winner.value.matched)
        if ranked.count > 1 {
            result.runnerUpPhaseIndex = Int32(ranked[1].key)
            result.runnerUpScore = Float(ranked[1].value.score)
        }

        // 4 — the two refusals. Both are recorded with the winner's numbers
        // still attached, so "not indexed" can be inspected rather than
        // guessed at.
        if winner.value.score > settings.notIndexedAboveInvAngstrom {
            result.verdict = .notIndexed
            return result
        }
        if settings.minimumPhaseContrastInvAngstrom > 0, ranked.count > 1 {
            let margin = ranked[1].value.score - winner.value.score
            if margin < settings.minimumPhaseContrastInvAngstrom {
                result.verdict = .notIndexed
                return result
            }
        }

        // 5 — the matrix gets the last word.
        //
        // WHY, and what it fixes (2026-09-14, a diagnostic run). Steps 1-4 ask only
        // which CANDIDATE explains the surviving vectors best. They never ask
        // whether those vectors are unexplained at all. At a position where
        // the matrix crystal is present on a zone axis other than the fitted
        // one, step 1 removes almost nothing and the whole pattern is offered
        // to the candidates — and a candidate coherent with the matrix wins by
        // default, because nothing else is in the competition. Measured on the
        // demonstration scan: 2 250 positions of pure aluminium on [011] came back
        // 100 % "β″ [001]". β″ [001] covers 10 of the 16 [011]Al reflections
        // at 0.0059 Å⁻¹; Al [011] covers all 16 at 0.0000 Å⁻¹, and was never
        // asked.
        //
        // So before a position may be called a candidate phase, the matrix
        // crystal is offered every low-index zone axis
        // (`matrixChallengeBases`, the fixed 49-direction list — NOT a list
        // the scan selects, and it ignores any zone axis the user restricted
        // the matrix phase to), scored by the same rule and held to the same
        // eligibility as any candidate.
        //
        // THE RULE: the challenger must explain STRICTLY MORE vectors, at a
        // smaller mean distance. See `challengeByMatrix` for why "strictly
        // more" rather than "at least as many" — it is what makes a fully
        // explained precipitate impossible to erase.
        //
        // THE IN-PLANE ROTATION IS DERIVED HERE, NOT FITTED FOR THE SCAN. An
        // earlier implementation took one rotation per axis from
        // `fitZoneAxis` — a whole-scan fit — and the fixture refuted it: at a
        // grain-B position the pool carried [-1 1 0] at 100°, matching NOTHING,
        // while the same axis at 130° matches 8 of 8 at 0.0001 Å⁻¹. A
        // whole-scan fit ranks a rotation by matches summed over every grain,
        // so a rotation picking up scattered matches everywhere outranks the
        // one that fits the grain exactly. A grain's rotation is a property of
        // the GRAIN, and must be asked at the position.
        //
        // Searching 180 rotations × 49 axes per position would cost more than
        // the whole map. Instead the rotations that can possibly win are
        // enumerated: if an orientation explains this pattern, it puts each
        // reference under an observed vector of the same length, so taking the
        // longest few observed vectors and every reference of matching length
        // gives every rotation worth testing — exactly, not on a 2° grid.
        if !matrixChallenge.isEmpty,
           let challenged = challengeByMatrix(
               surviving: surviving, bases: matrixChallenge, settings: settings,
               accessibleRadius: accessibleRadius,
               beating: (matched: winner.value.matched, score: winner.value.score),
               scratch: scratch) {
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            result.entryIndex = -1
            // The CHALLENGER's numbers, not the rejected candidate's: they are
            // the ones that describe the verdict given. A review caught the first
            // version discarding them while a comment claimed the verdict could
            // be inspected.
            result.score = Float(challenged.score)
            result.matchedCount = Int32(challenged.matched)
            return result
        }

        result.verdict = .indexed
        result.phaseIndex = Int32(winner.key)
        result.entryIndex = Int32(winner.value.entryIndex)
        return result
    }

    // MARK: The known-variants rule (Thronsen et al.)

    /// The per-position rule Thronsen et al. describe
    /// (written from the paper, no code copied): score every non-matrix library entry against the
    /// survivors and take the argmin, no floors and no guards. `classify`
    /// dispatches here for `settings.classificationRule == .knownVariants`;
    /// this function never runs for `.search`.
    ///
    /// DEVIATION from the paper's rule (which assigns Al directly when few vectors survive): the paper's
    /// five references are θ′ and T1 only — Al is never among them. Here Al
    /// (the matrix) is decided ENTIRELY by step b, by survivor count, and a
    /// matrix-phase entry is excluded from scoring in step c for that same
    /// reason: it is not that it would score badly, it is that the rule
    /// never asks.
    public static func classifyKnownVariants(vectors: [SIMD2<Double>],
                                              library: PhaseReferenceLibrary,
                                              settings: PhaseVectorSettings,
                                              matrixEntry: PhaseOrientationReference?,
                                              candidateEntryIndices: [Int]) -> PhaseVectorResult {
        var result = PhaseVectorResult()
        guard !vectors.isEmpty else { return result }
        guard candidateEntryIndices.allSatisfy({ library.entries.indices.contains($0) }) else { return result }

        // a — matrix removal, identical to `classify` step 1: same
        // tolerance, same survivors, same removedCount. DEVIATION: the
        // paper masks in IMAGE space at a fixed 8 detector pixels; this
        // removes in Å⁻¹ VECTOR space at `matrixToleranceInvAngstrom` — the
        // same vector-space departure the file header already states for
        // `.search`, inherited here because this step IS `.search`'s step
        // 1, not a reimplementation of it.
        var surviving: [SIMD2<Double>] = []
        surviving.reserveCapacity(vectors.count)
        if let matrixEntry {
            for u in vectors {
                if nearest(u, in: matrixEntry.vectors,
                           radius: settings.matrixToleranceInvAngstrom) == nil {
                    surviving.append(u)
                }
            }
        } else {
            surviving = vectors
        }
        result.survivingCount = Int32(surviving.count)
        result.removedCount = Int32(vectors.count - surviving.count)

        // b — the matrix, by exclusion (survivors at or below `directMatrixMaximumVectors`).
        if surviving.count <= settings.directMatrixMaximumVectors {
            guard Int(result.removedCount) >= settings.minimumMatrixEvidence else {
                result.verdict = .notIndexed
                return result
            }
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            return result
        }

        // c — score every entry whose phase is NOT the matrix's. The paper's
        // per-entry scoring: for every survivor, the TRUE nearest
        // reference vector of the entry (the argmin over the WHOLE set, no
        // radius — see `nearestReferenceVector` below for why `nearest(_:in:
        // radius:)` above is not reused), summed over every survivor, and
        // divided by the count of DISTINCT reference vectors hit
        // (the count of unique reference vectors) — a reference claimed
        // twice costs the sum without buying a second denominator credit,
        // which is what makes two survivors splitting one spot score worse
        // than two survivors hitting two spots at the same distances.
        //
        // DEVIATION (the paper's loop carries none of these — each omission
        // named against the paper's rule, which has no such term):
        //  - no `minimumMatchedVectors` / Friedel-pair floor (the paper: every
        //    survivor is scored, however few clear any particular entry);
        //  - no `chanceMatchMultiple` guard (the paper: no chance expectation
        //    is computed at all);
        //  - no matrix challenge, `classify` step 5 (the paper: the matrix is
        //    never re-offered once step b has passed it by);
        //  - no `minimumPhaseContrastInvAngstrom` margin (the paper:
        //    the phase with the lowest score is taken
        //    outright, no runner-up margin considered);
        //  - no completeness-aware cross-phase ranking (the paper: the argmin
        //    is by score alone, matched count never enters the ordering).
        //
        // DEVIATION: the paper's five entries are hand-built at fixed
        // orientations (T1 at two in-plane rotations under the stated
        // orientation relationship, θ′ edge-on at two, face-on at one);
        // every non-matrix entry THIS library holds is scored here instead —
        // for this app, whatever zone axes and in-plane steps the phase list
        // and `PhaseReferenceSettings` produced. T1 is a free in-plane sweep
        // until a later session derives its stated relationship
        // ((0001)T1 ∥ (111)Al, [1‑10]Al ∥ [10‑10]T1) into the constrained
        // form θ′ already has.
        var bestPerPhase: [Int: (entryIndex: Int, score: Double, matched: Int)] = [:]
        for entryIndex in candidateEntryIndices {
            let entry = library.entries[entryIndex]
            guard entry.phaseIndex != library.matrixPhaseIndex, !entry.vectors.isEmpty
            else { continue }
            var sum = 0.0
            var uniqueHits = Set<Int>()
            var matchedWithinRadius = 0
            for u in surviving {
                guard let hit = nearestReferenceVector(u, in: entry.vectors) else { continue }
                sum += hit.distance
                uniqueHits.insert(hit.index)
                if hit.distance <= settings.pairRadiusInvAngstrom { matchedWithinRadius += 1 }
            }
            guard !uniqueHits.isEmpty else { continue }
            let score = sum / Double(uniqueHits.count)
            if bestPerPhase[entry.phaseIndex] == nil || score < bestPerPhase[entry.phaseIndex]!.score {
                bestPerPhase[entry.phaseIndex] = (entryIndex, score, matchedWithinRadius)
            }
        }
        guard !bestPerPhase.isEmpty else {
            result.verdict = .notIndexed
            return result
        }

        // d — argmin over ALL entries, ties to the first in library order.
        // Reduced through `bestPerPhase` above (the global minimum is always
        // some phase's own minimum), then picked by entry index ascending so
        // a tie keeps the earliest entry, the way an argmin over a flat
        // array would.
        let ranked = bestPerPhase.map {
            (phaseIndex: $0.key, entryIndex: $0.value.entryIndex,
             score: $0.value.score, matched: $0.value.matched)
        }.sorted { $0.entryIndex < $1.entryIndex }
        var winner = ranked[0]
        for candidate in ranked.dropFirst() where candidate.score < winner.score {
            winner = candidate
        }
        var runnerUp: (phaseIndex: Int, score: Double)?
        for candidate in ranked where candidate.phaseIndex != winner.phaseIndex {
            if runnerUp == nil || candidate.score < runnerUp!.score {
                runnerUp = (candidate.phaseIndex, candidate.score)
            }
        }

        // Filled before the verdict is decided: a refused position must
        // still carry the would-be winner's numbers, so a caller can print
        // the distribution of scores that missed the cutoff.
        result.phaseIndex = Int32(winner.phaseIndex)
        result.entryIndex = Int32(winner.entryIndex)
        result.score = Float(winner.score)
        result.matchedCount = Int32(winner.matched)
        if let runnerUp {
            result.runnerUpPhaseIndex = Int32(runnerUp.phaseIndex)
            result.runnerUpScore = Float(runnerUp.score)
        }

        // The paper's residual cutoff (0.07 Å⁻¹): a winner above it is not indexed.
        result.verdict = winner.score > settings.residualCutoffInvAngstrom ? .notIndexed : .indexed

        // e — the evidence guard (DEVIATION, `knownVariantsMinimumSpecificReflections`):
        // an accepted winner with too few phase-specific reflections is the
        // matrix. Only `.indexed` is judged; a refusal stays a refusal. As in
        // the other matrix fallbacks, `entryIndex` is cleared; the rejected
        // winner's score, matched count and runner-up stay for diagnosis.
        let minimumSpecific = settings.knownVariantsMinimumSpecificReflections
        if result.verdict == .indexed, minimumSpecific > 0,
           specificReflectionCount(surviving: surviving, entry: library.entries[winner.entryIndex],
                                   matrixEntry: matrixEntry, settings: settings) < minimumSpecific {
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            result.entryIndex = -1
        }
        return result
    }

    /// Distinct reference vectors of `entry` hit by some survivor within
    /// `pairRadiusInvAngstrom`, whose own `q` is not within
    /// `matrixToleranceInvAngstrom` of any matrix reference. With no matrix
    /// entry every hit counts. Transcribes the inline guard of the measuring
    /// tool the rule was developed with, so this rule reproduces its measured
    /// numbers.
    public static func specificReflectionCount(surviving: [SIMD2<Double>],
                                                entry: PhaseOrientationReference,
                                                matrixEntry: PhaseOrientationReference?,
                                                settings: PhaseVectorSettings) -> Int {
        var hits = Set<Int>()
        for u in surviving {
            if let hit = nearest(u, in: entry.vectors, radius: settings.pairRadiusInvAngstrom) {
                hits.insert(hit.index)
            }
        }
        guard let matrixEntry else { return hits.count }
        return hits.filter { index in
            let q = entry.vectors[index].q
            return !matrixEntry.vectors.contains {
                simd_distance($0.q, q) <= settings.matrixToleranceInvAngstrom
            }
        }.count
    }

    /// The TRUE nearest reference vector of `refs` — brute force, no radius
    /// and no length-band prune.
    ///
    /// `nearest(_:in:radius:)` above is deliberately NOT reused: its
    /// length-band prune assumes a candidate outside `|u| ± radius` can be
    /// discarded without being computed, which is exactly the restriction
    /// this rule does not have — a survivor here must find its true nearest
    /// reference however far away that reference sits, because the residual
    /// itself is what the cutoff judges.
    private static func nearestReferenceVector(_ u: SIMD2<Double>, in refs: [ReferenceVector])
        -> (index: Int, distance: Double)? {
        guard !refs.isEmpty else { return nil }
        var best = 0
        var bestDistance = simd_distance(u, refs[0].q)
        for i in 1..<refs.count {
            let d = simd_distance(u, refs[i].q)
            if d < bestDistance { bestDistance = d; best = i }
        }
        return (best, bestDistance)
    }

    /// The MATRIX crystal projected down every low-index zone axis, unrotated —
    /// the pool `classify` challenges a candidate label with. One entry per
    /// axis (49 at the shipped list), built once per scan; the in-plane
    /// rotation is derived at each position by `challengeByMatrix`, because a
    /// grain's rotation is a property of the grain and not of the scan.
    public static func matrixChallengeBases(
        library: PhaseReferenceLibrary
    ) -> [PhaseOrientationReference] {
        guard library.phases.indices.contains(library.matrixPhaseIndex) else { return [] }
        let crystal = library.phases[library.matrixPhaseIndex].crystal
        // `build` made this very list for this phase and structure, so it cannot be refused here; a library assembled by
        // hand around a structure the kit refuses has no challenge pool.
        guard let reflections = try? PhaseReferenceLibrary.reflections(
            of: crystal, kMaxInvAngstrom: library.settings.kMaxInvAngstrom),
              !reflections.isEmpty else { return [] }
        var out: [PhaseOrientationReference] = []
        out.reserveCapacity(PhaseReferenceLibrary.lowIndexZoneAxes.count)
        for axis in PhaseReferenceLibrary.lowIndexZoneAxes {
            let base = PhaseReferenceLibrary.projectedVectors(
                reflections: reflections, crystal: crystal,
                zoneAxis: axis, settings: library.settings)
            guard !base.isEmpty else { continue }
            out.append(PhaseOrientationReference(
                phaseIndex: library.matrixPhaseIndex, zoneAxis: axis,
                inPlaneRotationRad: 0, vectors: base))
        }
        return out
    }

    /// Does any orientation of the matrix crystal explain `surviving` at least
    /// as well as the winning candidate does? Returns the orientation that
    /// does, or nil.
    ///
    /// THE RULE, and why it is what it is. A challenger must explain STRICTLY
    /// MORE of the pattern than the winning candidate, at a smaller mean
    /// distance, having cleared the same floors a candidate clears.
    ///
    /// "Strictly more" is not fussiness; it is what makes the guard unable to
    /// erase a precipitate. A candidate that already explains every surviving
    /// vector cannot be beaten on count by anything, so a genuine precipitate
    /// pattern is safe BY CONSTRUCTION. The first version required only "at
    /// least as many", and a review measured what that costs: with a sparse
    /// precipitate of 3 to 6 vectors and 0.004 Å⁻¹ of jitter — a third of a
    /// detector pixel on the demonstration scan, inside this file's own stated
    /// measurement budget — a 49-axis search found a coincidental matrix fit
    /// with the same count and a smaller distance, and took 26.5 % of
    /// three-vector precipitates. A second matrix grain is the opposite case:
    /// the candidate explains part of it and the matrix explains all of it.
    ///
    /// Rotations are enumerated, not searched on a grid: an orientation that
    /// explains the pattern must put some reference of matching length under
    /// each observed vector, so every rotation worth testing is
    /// `angle(u) − angle(v)` for an observed `u` and a reference `v` with
    /// `||u| − |v||` inside the pair radius, taken over every surviving
    /// vector. The enumerated angle is exact only for a zero-distance match —
    /// its error is about `asin(ρ/|u|)` — and that error biases toward keeping
    /// the candidate label, which is the safe direction.
    public static func challengeByMatrix(
        surviving: [SIMD2<Double>],
        bases: [PhaseOrientationReference],
        settings: PhaseVectorSettings,
        accessibleRadius: Double,
        beating winner: (matched: Int, score: Double),
        scratch: Scratch
    ) -> (entry: PhaseOrientationReference, score: Double, matched: Int)? {
        guard !surviving.isEmpty, !bases.isEmpty else { return nil }
        let radius = settings.pairRadiusInvAngstrom
        // EVERY surviving vector seeds a rotation, not a chosen few. The first
        // version took the three LONGEST, on the reasoning that an angle is
        // most precise far out — and a review measured that backwards: spurious
        // maxima are typically FARTHER out than any real reflection, so the
        // longest-three rule hands the seeds to exactly them. With three
        // spurious peaks beyond the outermost real one, the catch rate for a
        // second matrix grain went from 100 % to 0 % (2026-09-14, 200 trials
        // per scenario). Seeding on all of them removes the failure mode and
        // the arbitrary constant with it; the cost is linear in the vector
        // count, measured at ~1.6 µs per surviving vector per challenged
        // position.
        //
        // The BEST qualifying orientation is returned, not the first. The
        // numbers it carries reach the result and the user, so "whichever
        // zone axis came first in the list" is not good enough for them.
        var best: (entry: PhaseOrientationReference, score: Double, matched: Int)?
        var rotated = surviving
        for base in bases {
            let chance = base.chanceMatchFraction(
                pairRadius: radius, accessibleRadius: accessibleRadius)
                * Double(surviving.count)
            for u in surviving {
                let uLength = simd_length(u)
                guard uLength > 0 else { continue }
                let uAngle = atan2(u.y, u.x)
                for v in base.vectors where abs(v.length - uLength) <= radius {
                    // Rotating the references by θ puts v under u; scoring the
                    // observations rotated by −θ against the unrotated base is
                    // the same comparison and rotates 8 vectors instead of 48.
                    let theta = uAngle - atan2(v.q.y, v.q.x)
                    let c = cos(theta), s = sin(theta)
                    for i in surviving.indices {
                        rotated[i] = SIMD2(surviving[i].x * c + surviving[i].y * s,
                                           -surviving[i].x * s + surviving[i].y * c)
                    }
                    guard let candidate = score(vectors: rotated, against: base,
                                                pairRadius: radius, scratch: scratch)
                    else { continue }
                    guard candidate.matched > winner.matched,
                          candidate.score < winner.score,
                          candidate.matched >= settings.minimumMatchedVectors,
                          Double(candidate.matched) >= settings.chanceMatchMultiple * chance
                    else { continue }
                    if let current = best,
                       current.matched > candidate.matched
                        || (current.matched == candidate.matched
                            && current.score <= candidate.score) { continue }
                    // Rotated, because the entry claims `theta`: a
                    // `PhaseOrientationReference`'s vectors are defined to be
                    // the ones AFTER its in-plane rotation, and handing back
                    // the unrotated base with a non-zero angle would be a trap
                    // for the next caller.
                    best = (PhaseOrientationReference(
                        phaseIndex: base.phaseIndex, zoneAxis: base.zoneAxis,
                        inPlaneRotationRad: theta,
                        vectors: PhaseReferenceLibrary.rotate(base.vectors, by: theta)),
                        candidate.score, candidate.matched)
                }
            }
        }
        return best
    }

    // MARK: The orientation relationship

    /// The in-plane azimuth, radians, of a Cartesian vector projected into
    /// an entry's detector frame (`ZoneAxisFrame.detectorBasis` of its
    /// zone axis, then the entry's own in-plane rotation) — nil when the
    /// vector is not in the zone (|v·n̂| > 1e-6·|v|, `n̂` the unit zone axis)
    /// or is zero.
    public static func projectedAzimuth(of v: SIMD3<Double>, zoneAxis n: SIMD3<Double>,
                                         inPlaneRotationRad: Double) -> Double? {
        let vLen = simd_length(v)
        guard vLen.isFinite, vLen > 1e-12 else { return nil }
        let nLen = simd_length(n)
        guard nLen.isFinite, nLen > 1e-12 else { return nil }
        let nHat = n / nLen
        guard abs(simd_dot(v, nHat)) <= 1e-6 * vLen else { return nil }

        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: nHat)
        let x0 = simd_dot(v, basis.columns.0)
        let y0 = simd_dot(v, basis.columns.1)
        let c = cos(inPlaneRotationRad), s = sin(inPlaneRotationRad)
        // The library's own rotate convention (`PhaseReferenceLibrary.rotate`):
        // x' = c·x − s·y, y' = s·x + c·y, i.e. azimuth + θ.
        let x = c * x0 - s * y0
        let y = s * x0 + c * y0
        return atan2(y, x)
    }

    /// Is a candidate entry consistent with any listed relationship, given
    /// the fitted matrix entry?
    ///
    /// For each pair: the candidate's azimuth (its vector in the candidate
    /// crystal, projected through the candidate entry's own zone and
    /// rotation) minus the matrix's azimuth (likewise, through the matrix
    /// entry's own zone and rotation) must be within `toleranceDeg` of 0
    /// modulo 180° (see `PhaseDefinition.orientationRelationships` for why
    /// 180 and not 360). A pair whose vector is not in the respective zone
    /// does not apply to this (candidate, matrix) pair of entries; if no
    /// listed pair applies, the entry is free (true). An empty list is free.
    public static func orientationConsistent(candidate: PhaseOrientationReference,
                                               candidateCrystal: CrystalStructure,
                                               matrix: PhaseOrientationReference,
                                               matrixCrystal: CrystalStructure,
                                               relationships: [OrientationRelationship],
                                               toleranceDeg: Double) -> Bool {
        guard !relationships.isEmpty else { return true }
        guard let candidateZone = PhaseReferenceLibrary.cartesianZoneAxis(
                candidate.zoneAxis, crystal: candidateCrystal),
              let matrixZone = PhaseReferenceLibrary.cartesianZoneAxis(
                matrix.zoneAxis, crystal: matrixCrystal)
        else { return true } // a degenerate zone axis never reaches an entry; defensive.

        var anyApplicable = false
        for relationship in relationships {
            let candidateVector = relationship.candidate.cartesian(in: candidateCrystal)
            let matrixVector = relationship.matrix.cartesian(in: matrixCrystal)
            guard let candidateAz = projectedAzimuth(
                    of: candidateVector, zoneAxis: candidateZone,
                    inPlaneRotationRad: candidate.inPlaneRotationRad),
                  let matrixAz = projectedAzimuth(
                    of: matrixVector, zoneAxis: matrixZone,
                    inPlaneRotationRad: matrix.inPlaneRotationRad)
            else { continue }
            anyApplicable = true
            var d = (candidateAz - matrixAz) * 180 / .pi
            d = d.truncatingRemainder(dividingBy: 180)
            if d < 0 { d += 180 }
            let distance = min(d, 180 - d)
            if distance < toleranceDeg { return true }
        }
        return !anyApplicable
    }
}
