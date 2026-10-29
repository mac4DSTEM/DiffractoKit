//
//  MatrixChallengeTests.swift
//  Role: the matrix's last word before a position may be called a precipitate.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct MatrixChallengeTests {
    /// Mutation this names: deleting the step-5 block in `classify`, or
    /// passing an empty `matrixChallenge` from `map`. Both restore the defect
    /// exactly — the second assertion here is the defect, kept as the control.
    ///
    /// The planted rotation is 37.3°, deliberately NOT a multiple of the 2°
    /// library step. MEASURED, not assumed: snapping the derived rotation to
    /// that grid does NOT turn this test red — 0.7° at |q| ≤ 0.8 Å⁻¹ is
    /// 0.0098 Å⁻¹, still inside the 0.02 pair radius, so the verdict survives
    /// at a worse distance. The grid mutation is caught by
    /// `testAMatrixOrientationExplainingLessDoesNotWin`, which compares the
    /// distance rather than the verdict.
    @Test func testAMatrixGrainOnAnotherZoneAxisIsNotLabelledAsACandidate() throws {
        let library = try challengeLibrary()
        let planted = matrixVectorsOnAnotherAxis(library, degrees: 37.3)
        #expect(planted.count >= 4,
            "the fixture planted too few vectors to decide anything")
        let settings = PhaseVectorSettings()
        let bases = PhaseVectorMatcher.matrixChallengeBases(library: library)
        #expect(!(bases.isEmpty),
            "the challenge pool is empty; nothing was tested")

        // The control: without the challenge this IS the defect.
        let unchallenged = PhaseVectorMatcher.classify(
            vectors: planted, library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 256))
        #expect(unchallenged.verdict == .indexed,
            "the fixture no longer reproduces the defect it was built for, so the challenged case below proves nothing")

        let challenged = PhaseVectorMatcher.classify(
            vectors: planted, library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 256),
            matrixChallenge: bases)
        #expect(challenged.verdict == .matrix,
            "aluminium on a second zone axis was labelled a candidate phase")
        #expect(Int(challenged.phaseIndex) == library.matrixPhaseIndex)
    }

    /// Mutation this names: relaxing `candidate.matched > winner.matched` to
    /// `>=`. That single character is what makes a precipitate impossible to
    /// erase — a candidate explaining every surviving vector cannot be beaten
    /// on count — and a review measured what `>=` costs before it was tightened:
    /// **26.5 % of three-vector precipitates taken by the matrix** at
    /// 0.004 Å⁻¹ of jitter, a third of a detector pixel on the demo cube.
    ///
    /// THE JITTER IS THE POINT. The first version of this test planted a
    /// candidate entry's own vectors verbatim, so the winner's mean distance
    /// was exactly 0 and `candidate.score < winner.score` could never hold for
    /// any challenger — it passed under every loosening of the guard it
    /// claimed to hold. A real measurement is never exact.
    @Test func testTheChallengeDoesNotTakeATruePrecipitate() throws {
        let library = try challengeLibrary()
        let candidate = library.entries[library.candidateEntryIndices[0]]
        let bases = PhaseVectorMatcher.matrixChallengeBases(library: library)
        let settings = PhaseVectorSettings()

        // Deterministic jitter, well inside the pair radius: a real peak is
        // measured, not looked up.
        var seed: UInt64 = 0xD1B54A32D192ED03
        func jitter() -> Double {
            seed = seed &* 6364136223846793005 &+ 1442695040888963407
            return (Double(seed >> 11) / Double(UInt64(1) << 53) - 0.5) * 0.008
        }
        // Sparse, because sparse is the case that was being stolen: the file
        // header names "a thin precipitate showing two spots", and
        // `minimumMatchedVectors` admits three.
        for count in [3, 4, 5, 6] {
            let planted = candidate.vectors.prefix(count).map {
                $0.q + SIMD2(jitter(), jitter())
            }
            let result = PhaseVectorMatcher.classify(
                vectors: Array(planted), library: library, settings: settings,
                matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
                scratch: PhaseVectorMatcher.Scratch(capacity: 256),
                matrixChallenge: bases)
            #expect(result.verdict == .indexed,
            "a \(count)-vector precipitate was taken by the matrix")
            #expect(Int(result.phaseIndex) != library.matrixPhaseIndex,
            "a \(count)-vector precipitate was relabelled as the matrix")
        }
    }

    /// Mutation this names: dropping `candidate.matched >= winner.matched`, or
    /// turning `candidate.score < winner.score` into `<=`. The rule is that a
    /// challenger must explain AT LEAST AS MUCH of the pattern, not merely
    /// explain a corner of it more precisely.
    @Test func testAMatrixOrientationExplainingLessDoesNotWin() throws {
        let library = try challengeLibrary()
        let planted = matrixVectorsOnAnotherAxis(library, degrees: 37.3)
        let bases = PhaseVectorMatcher.matrixChallengeBases(library: library)
        let settings = PhaseVectorSettings()
        let scratch = PhaseVectorMatcher.Scratch(capacity: 256)

        // ONE named orientation, not the whole pool: the pool holds six
        // ⟨011⟩ variants and 49 axes in all, and a boundary asserted against
        // all of them would be decided by whichever one happens to win.
        // Restricting it to the planted axis makes "the challenger's numbers"
        // exact.
        // The planted axis itself, up to sign — NOT the ⟨011⟩ family, which
        // has six members in `lowIndexZoneAxes` and would reintroduce the same
        // first-versus-best ambiguity.
        let oneBase = bases.filter {
            $0.zoneAxis == SIMD3(0, 1, 1) || $0.zoneAxis == SIMD3(0, -1, -1)
        }
        #expect(oneBase.count == 1,
            "expected exactly one [011] base in the pool")

        // What that orientation can actually do here, measured rather than
        // assumed, so the comparisons below sit exactly on its boundary.
        guard let best = PhaseVectorMatcher.challengeByMatrix(
            surviving: planted, bases: oneBase, settings: settings,
            accessibleRadius: 1.0, beating: (matched: 0, score: .infinity),
            scratch: scratch) else {
            Issue.record("the matrix explains every planted vector exactly and still lost")
            return
        }
        #expect(best.matched == planted.count,
            "the planted grain was not fully explained by its own orientation")

        // ON the boundary: a winner the challenger merely TIES on count is not
        // beaten, because the rule is strictly more. This is the assertion that
        // makes a fully explained precipitate safe, and it is what turns
        // `>` into `>=` red.
        #expect(PhaseVectorMatcher.challengeByMatrix(
            surviving: planted, bases: oneBase, settings: settings,
            accessibleRadius: 1.0, beating: (matched: best.matched, score: 1.0),
            scratch: scratch) == nil,
            "a challenger that merely tied the winner's count took the position")

        // And a winner it explains more than, it takes.
        #expect((PhaseVectorMatcher.challengeByMatrix(
            surviving: planted, bases: oneBase, settings: settings,
            accessibleRadius: 1.0, beating: (matched: best.matched - 1, score: 1.0),
            scratch: scratch)) != nil,
            "a challenger explaining one more vector than the winner did not take it")

        // The distance half of the rule, at its own boundary.
        #expect((PhaseVectorMatcher.challengeByMatrix(
            surviving: planted, bases: oneBase, settings: settings,
            accessibleRadius: 1.0, beating: (matched: best.matched - 1, score: best.score),
            scratch: scratch)) == nil,
            "a challenger tied the winner's distance and took the position anyway")
    }

    /// THE FRIEDEL-PAIR FLOOR (2026-09-15, Thronsen step 3). A candidate that
    /// explains exactly the two survivors u and −u is indexed; two survivors
    /// that are not a pair are not enough. Mutations this names: the
    /// `containsFriedelPair` clause dropped (the pair case goes not indexed),
    /// or the floor applied to any two survivors (the control goes indexed
    /// — its two survivors both match the candidate, so it is the floor and
    /// nothing else that keeps it out).
    @Test func testAFriedelPairAloneIndexesACandidateAndTwoStraySurvivorsDoNot() throws {
        let library = try challengeLibrary()
        let settings = PhaseVectorSettings()
        let candidate = library.entries[library.candidateEntryIndices[0]]
        let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
        // A reference vector of the candidate that no matrix vector sits on,
        // and its Friedel partner.
        let matrix = library.entries[library.matrixEntryIndices[0]]
        let own = try #require(candidate.vectors.first { v in
            !matrix.vectors.contains { simd_distance($0.q, v.q) < 2 * settings.pairRadiusInvAngstrom }
                && candidate.vectors.contains { simd_distance($0.q, -v.q) < 1e-6 }
        },
            "the candidate has no Friedel pair clear of the matrix")
        let pair = [own.q, -own.q]
        let indexed = PhaseVectorMatcher.classify(
            vectors: pair, library: library, settings: settings, matrixEntry: matrix,
            candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        #expect(indexed.verdict == .indexed,
            "a Friedel pair of the candidate's own reflections was not indexed: \(indexed.verdict)")
        #expect(Int(indexed.phaseIndex) == candidate.phaseIndex)

        // Two survivors that are NOT a pair but are both the candidate's own
        // reflections, clear of the matrix, no shorter than `own` (so the
        // chance guard passes as it did for the pair): two matched of two,
        // which only the pair floor could admit. A stray that matches one
        // reference cannot tell the floor from the count (a review of 2026-09-15).
        let other = try #require(candidate.vectors.first { v in
            simd_length(v.q) >= simd_length(own.q) - 1e-9
                && simd_distance(v.q, own.q) > 2 * settings.pairRadiusInvAngstrom
                && simd_distance(v.q, -own.q) > 2 * settings.pairRadiusInvAngstrom
                && !matrix.vectors.contains { simd_distance($0.q, v.q) < 2 * settings.pairRadiusInvAngstrom }
        },
            "the candidate has no second reflection clear of the matrix")
        let stray = [own.q, other.q]
        #expect(!(PhaseVectorMatcher.containsFriedelPair(stray, radius: settings.pairRadiusInvAngstrom)))
        let control = PhaseVectorMatcher.classify(
            vectors: stray, library: library, settings: settings, matrixEntry: matrix,
            candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        #expect(control.verdict != .indexed,
            "two stray survivors were indexed as a candidate under the pair floor")
    }

    // MARK: Fixtures

    /// Al as the matrix on [001], β″ [001] as the candidate — the pair that
    /// produced the false label, because β″ is coherent with Al and its [001]
    /// net covers most of Al's [011] net.
    private func challengeLibrary() throws -> PhaseReferenceLibrary {
        var settings = PhaseReferenceSettings()
        settings.kMaxInvAngstrom = 0.8
        settings.inPlaneStepDeg = 2
        return try PhaseReferenceLibrary.build(phases: [
            PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum,
                            role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime,
                            role: .candidate, zoneAxes: [SIMD3(0, 0, 1)]),
        ], settings: settings)
    }

    /// The matrix crystal's [011] net — a second aluminium grain, the case the
    /// library has no entry for because the user named one zone axis.
    private func matrixVectorsOnAnotherAxis(_ library: PhaseReferenceLibrary,
                                            degrees: Double) -> [SIMD2<Double>] {
        let crystal = library.phases[library.matrixPhaseIndex].crystal
        let reflections = refl(crystal, library.settings.kMaxInvAngstrom)
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: reflections, crystal: crystal,
            zoneAxis: SIMD3(0, 1, 1), settings: library.settings)
        return PhaseReferenceLibrary.rotate(base, by: degrees * .pi / 180).map(\.q)
    }
}
