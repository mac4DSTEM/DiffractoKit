//
//  PhaseVectorMatchingTests.swift
//  Role: the planted-truth and unit-level tests of the reference library and the vector matcher, moved to Swift Testing from the
//        application tests of the ported types with every expectation kept. Each test's comment
//        names the mutation it catches.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct PhaseVectorMatchingTests {
    // MARK: The reciprocal lattice the library is built from

    /// Every lattice point inside kMax, enumerated over a tile far larger than
    /// any bound the code could pick. Ground truth is the lattice itself: a
    /// one-atom P1 cell has no extinction and `f_e > 0`, so every one of these
    /// must be returned; an fcc cell must return exactly the all-even /
    /// all-odd subset.
    private func latticePoints(_ crystal: CrystalStructure, kMax: Double, tile: Int = 48) -> Set<SIMD3<Int>> {
        var out = Set<SIMD3<Int>>()
        for h in -tile...tile { for k in -tile...tile { for l in -tile...tile {
            if h == 0 && k == 0 && l == 0 { continue }
            let g = Double(h) * crystal.cell.reciprocalLattice[0] + Double(k) * crystal.cell.reciprocalLattice[1]
                + Double(l) * crystal.cell.reciprocalLattice[2]
            if simd_length(g) <= kMax { out.insert(SIMD3(h, k, l)) }
        } } }
        return out
    }

    /// Runs on the kit's reflection list (`PhaseReferenceLibrary.reflections`), whose box is the per-axis bound.
    ///
    /// Mutation: the tile bound `ceil(kMax / kMin)` that stood until
    /// 2026-09-14 (and still stands in py4DSTEM). The exact bound on an index
    /// is |h| = |g·a₁| ≤ kMax·|a₁|, and for a b-unique monoclinic cell
    /// 1/kMin ≤ a·sin β < a, so the old tile fell short as β left 90° and
    /// reflections went missing with no signal. A review measured 6 lost at
    /// β = 110°, 48 at 115°, 198 at 125° on the β″-shaped cell
    /// (2026-09-12); this pins the 125° case against the lattice and the
    /// shipped cells against their own extinction rules, so a bound that
    /// over-tiles is not caught here — only one that under-tiles.
    @Test func testReflectionsCoverEveryLatticePointInsideKMaxOnObliqueCells() throws {
        let kMax = 1.6
        for betaDeg in [105.3, 115.0, 125.0] {
            let cell = makeCrystal(a: 15.16, b: 4.05, c: 6.74, betaDeg: betaDeg,
                               sites: [AtomSite(z: 13, fractional: SIMD3(0, 0, 0))])
            let expected = latticePoints(cell, kMax: kMax)
            let got = Set(refl(cell, kMax).map { SIMD3($0.reflection.hkl.h, $0.reflection.hkl.k, $0.reflection.hkl.l) })
            #expect(got.count == expected.count,
            "β = \(betaDeg)°: \(expected.subtracting(got).count) lattice points inside kMax are missing, e.g. \(expected.subtracting(got).prefix(3))")
            #expect(got.isSubset(of: expected),
            "β = \(betaDeg)°: a reflection outside kMax")
        }
        // The shipped cells, against their space groups: fcc keeps hkl of one
        // parity; β″ (C2/m) keeps h + k even. Both must be complete AND exact.
        let al = CrystalStructure.fcc(a: 4.05, z: 13)
        let alExpected = latticePoints(al, kMax: kMax).filter {
            ($0.x & 1) == ($0.y & 1) && ($0.y & 1) == ($0.z & 1)
        }
        let alGot = Set(refl(al, kMax).map { SIMD3($0.reflection.hkl.h, $0.reflection.hkl.k, $0.reflection.hkl.l) })
        #expect(alGot == alExpected)
        // Tolerance 1e-9, not the default 1e-4: 38 allowed β″ reflections are
        // accidentally weak enough to fall under the default (measured
        // 2026-09-14), and this asserts the LATTICE, not the intensity cut.
        let beta = CrystalStructure.betaDoublePrime
        let betaExpected = latticePoints(beta, kMax: kMax).filter { ($0.x + $0.y) & 1 == 0 }
        let betaGot = Set(try! beta.kinematicReflections(kMax: kMax, cut: 1e-9).kept.map { SIMD3($0.reflection.hkl.h, $0.reflection.hkl.k, $0.reflection.hkl.l) })
        #expect(betaGot == betaExpected)
    }

    // MARK: The library

    /// Mutation: `cartesianZoneAxis` returning `SIMD3(Double(uvw.x), ...)`
    /// normalised — the naive reading of [u v w]. Invisible on every cubic
    /// phase, wrong on every monoclinic one, which is the case the feature
    /// exists for.
    @Test func testZoneAxisIsARealSpaceDirectionNotACartesianTriple() throws {
        let beta = CrystalStructure.betaDoublePrime
        let n = PhaseReferenceLibrary.cartesianZoneAxis(SIMD3(0, 0, 1), crystal: beta)
        let axis = try #require(n)
        let betaRad = beta.cell.betaDeg * .pi / 180
        #expect(abs(axis.x - cos(betaRad)) <= 1e-12)
        #expect(abs(axis.y - 0) <= 1e-12)
        #expect(abs(axis.z - sin(betaRad)) <= 1e-12)
        // And it really differs from the naive triple, so the mutation is
        // detectable rather than merely wrong.
        #expect((acos(abs(axis.z)) * 180 / .pi) > 10)
    }

    /// Mutation: dropping the `guard n > 1e-12` in `cartesianZoneAxis`.
    @Test func testZeroZoneAxisYieldsNoDirection() throws {
        #expect(PhaseReferenceLibrary.cartesianZoneAxis(
            SIMD3(0, 0, 0), crystal: CrystalStructure.aluminum) == nil)
    }

    /// Mutation: `matrixIndices.count == 1` relaxed to `>= 1`, or the guard
    /// removed. Without exactly one matrix there is no "by exclusion" answer,
    /// and silently picking the first would make the verdict depend on the
    /// order the user added phases.
    @Test func testLibraryRefusesUnlessExactlyOneMatrixPhase() throws {
        func build(_ roles: [PhaseDefinition.Role]) throws -> PhaseReferenceLibrary {
            try PhaseReferenceLibrary.build(phases: roles.enumerated().map {
                PhaseDefinition(id: "p\($0.offset)", displayName: "P\($0.offset)",
                                crystal: .aluminum, role: $0.element,
                                zoneAxes: [SIMD3(0, 0, 1)])
            }, settings: coarseSettings())
        }
        #expect(throws: PhaseReferenceLibrary.Failure.matrixPhaseNotUnique(count: 0)) {
            try build([.candidate, .candidate])
        }
        #expect(throws: PhaseReferenceLibrary.Failure.matrixPhaseNotUnique(count: 2)) {
            try build([.matrix, .matrix, .candidate])
        }
        #expect(throws: PhaseReferenceLibrary.Failure.noCandidatePhases) {
            try build([.matrix])
        }
        #expect(throws: Never.self) { try build([.matrix, .candidate]) }
    }

    /// Mutation: the `requested <= maximumEntries` guard removed.
    ///
    /// It does NOT catch the guard being MOVED after the build loop, and two
    /// attempts to make it do so both failed (measured 2026-09-12): building
    /// 720 entries takes under a second, and so does building 44 640, because
    /// each entry is a rotation of an already-projected set of at most 48
    /// vectors. Timing cannot separate "refused before any work" from
    /// "refused after cheap work". The ordering is stated in the source and is
    /// what makes the refusal meaningful for a library that WOULD be expensive;
    /// this test pins the refusal and its numbers, which is what it can see.
    /// Saying so beats keeping an assertion that cannot fail — this repo has
    /// shipped one of those before (the L4
    /// control that everyone believed had teeth).
    @Test func testOversizedLibraryIsRefusedWithTheNumbersThatExplainIt() throws {
        var settings = coarseSettings()
        settings.inPlaneStepDeg = 1
        settings.maximumEntries = 10
        // Empty `zoneAxes` means every low-index axis, so this asks for tens of
        // thousands of entries of a 22-atom monoclinic cell. The size is the
        // point: at 720 entries the timing bound below had no teeth — building
        // them took under a second and the test passed against the mutation it
        // names (measured 2026-09-12).
        let axes = PhaseReferenceLibrary.lowIndexZoneAxes.count
        let expected = 2 * axes * 360
        #expect(expected > 20000)
        let started = Date()
        do {
            _ = try PhaseReferenceLibrary.build(phases: [
                PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum,
                                role: .matrix, zoneAxes: []),
                PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime,
                                role: .candidate, zoneAxes: []),
            ], settings: settings)
            Issue.record("expected libraryTooLarge, got a library")
        } catch let error as PhaseReferenceLibrary.Failure {
            guard case .libraryTooLarge(let requested, let limit) = error else {
                Issue.record("expected libraryTooLarge, got \(error)")
                return
            }
            #expect(requested == expected)
            #expect(limit == 10)
        }
        // Kept as a smoke bound only — it is NOT the ordering check its first
        // version claimed to be (see the comment above).
        #expect(Date().timeIntervalSince(started) < 30.0)
    }

    /// Mutation: the cap applied before the strongest-first sort, or the tie
    /// break dropped so the kept set depends on enumeration order.
    @Test func testDensityCapKeepsTheStrongestAndIsDeterministic() throws {
        var settings = coarseSettings()
        settings.minimumIntensityFraction = 0
        settings.maximumVectorsPerEntry = 0
        let beta = CrystalStructure.betaDoublePrime
        let uncapped = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(beta, settings.kMaxInvAngstrom), crystal: beta,
            zoneAxis: SIMD3(0, 1, 0), settings: settings)
        settings.maximumVectorsPerEntry = 12
        let capped = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(beta, settings.kMaxInvAngstrom), crystal: beta,
            zoneAxis: SIMD3(0, 1, 0), settings: settings)
        #expect(uncapped.count > 12)
        #expect(capped.count == 12)
        let weakestKept = capped.map(\.relativeIntensity).min() ?? 0
        let strongestDropped = uncapped
            .filter { kept in !capped.contains { $0.h == kept.h && $0.l == kept.l } }
            .map(\.relativeIntensity).max() ?? 0
        #expect(weakestKept >= strongestDropped,
            "the cap dropped a reflection stronger than one it kept")
        // Determinism: the same inputs give the same set, every time.
        let again = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(beta, settings.kMaxInvAngstrom), crystal: beta,
            zoneAxis: SIMD3(0, 1, 0), settings: settings)
        #expect(capped.map(\.length) == again.map(\.length))
    }

    /// Mutation: `rotate` using `[c, s; -s, c]` (the transpose), or recomputing
    /// `length` from the rotated components and letting rounding reorder the
    /// array. The matcher's binary search assumes the sort, so a reordering is
    /// a wrong answer, not a slow one.
    @Test func testRotationPreservesLengthsAndSortOrder() throws {
        let al = CrystalStructure.aluminum
        let settings = coarseSettings()
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, settings.kMaxInvAngstrom), crystal: al,
            zoneAxis: SIMD3(0, 0, 1), settings: settings)
        let turned = PhaseReferenceLibrary.rotate(base, by: 37.2 * .pi / 180)
        #expect(base.count == turned.count)
        for (b, t) in zip(base, turned) {
            #expect(abs(b.length - t.length) <= 1e-15)
            #expect(abs(simd_length(t.q) - t.length) <= 1e-12)
        }
        #expect(turned == (turned.sorted { $0.length < $1.length }.map { $0 }),
            "the rotated set is no longer sorted by |q|")
        // The sense of the rotation: +37.2° takes (x, 0) into the first
        // quadrant, not the fourth. A transposed matrix fails here.
        let probe = PhaseReferenceLibrary.rotate(
            [ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(1, 0), length: 1, relativeIntensity: 1)],
            by: 37.2 * .pi / 180)
        #expect(probe[0].q.y > 0)
    }

    /// Mutation: `inPlaneSteps` returning `stride(from: 0, through: 360)`,
    /// which duplicates 0 and 360, or dropping the guard so a step of 0 loops.
    @Test func testInPlaneStepsCoverTheCircleOnceAndIncludeZero() throws {
        var settings = coarseSettings()
        settings.inPlaneStepDeg = 90
        let steps = PhaseReferenceLibrary.inPlaneSteps(settings)
        #expect(steps.count == 4)
        #expect(steps.first == 0)
        #expect(steps.last! < (2 * .pi))
        settings.inPlaneStepDeg = 0
        #expect(PhaseReferenceLibrary.inPlaneSteps(settings) == [0])
        settings.inPlaneStepDeg = .nan
        #expect(PhaseReferenceLibrary.inPlaneSteps(settings) == [0])
    }

    /// Mutation: `chanceMatchFraction` using radius instead of radius², or
    /// dropping the `min(1, …)` clamp so it reports a probability above one.
    @Test func testChanceMatchFractionIsAnAreaRatioAndIsClamped() throws {
        let entry = PhaseOrientationReference(
            phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: (0..<10).map {
                ReferenceVector(h: $0, k: 0, l: 0, q: SIMD2(Double($0) * 0.1, 0),
                                length: Double($0) * 0.1, relativeIntensity: 1)
            })
        // 10 discs of radius 0.02 in a disc of radius 1: 10 × 0.0004 / 1.
        #expect(abs(entry.chanceMatchFraction(pairRadius: 0.02, accessibleRadius: 1) - 0.004) <= 1e-12)
        // Doubling the radius quadruples it — an area, not a length.
        #expect(abs(entry.chanceMatchFraction(pairRadius: 0.04, accessibleRadius: 1) - 0.016) <= 1e-12)
        #expect(entry.chanceMatchFraction(pairRadius: 1, accessibleRadius: 1) == 1)
        #expect(entry.chanceMatchFraction(pairRadius: 0.02, accessibleRadius: 0) == 0)
    }

    // MARK: The matcher

    /// Mutation: `len > directBeamRadius` weakened to `len >= 0`, or the
    /// origin subtraction dropped. The direct beam is the brightest peak in
    /// every pattern and would match the shortest reference vector of every
    /// phase at once.
    @Test func testExperimentalVectorsSubtractTheOriginAndDropTheDirectBeam() throws {
        let peaks = [
            pixel(128, 128),      // the direct beam
            pixel(148, 128),        // +20 px in x
            pixel(128, 108),        // −20 px in y
        ]
        let vectors = PhaseVectorMatcher.experimentalVectors(
            peaks: peaks, originX: 128, originY: 128, invAngstromPerPixel: 0.01,
            directBeamRadiusInvAngstrom: 0.05)
        #expect(vectors.count == 2)
        #expect(abs(vectors[0].x - 0.2) <= 1e-9)
        #expect(abs(vectors[0].y - 0) <= 1e-9)
        // y is NOT flipped: a peak at a smaller row gives a negative q_y, the
        // same sense as template matching.
        #expect(abs(vectors[1].y - (-0.2)) <= 1e-9)
    }

    /// The outer reach (`maximumVectorInvAngstrom`, 2026-09-15): a vector at
    /// or beyond it is dropped, 0 keeps everything. Mutations: `<` to `<=`
    /// (the vector AT the reach survives), the guard dropped (both survive).
    @Test func testTheOuterReachDropsVectorsAtAndBeyondIt() throws {
        let peaks = [
            pixel(128, 128),
            pixel(148, 128),        // 0.20 Å⁻¹
            pixel(168, 128),        // 0.40 Å⁻¹, at the reach
            pixel(128, 178),        // 0.50 Å⁻¹, beyond it
        ]
        func vectors(reach: Double) -> [SIMD2<Double>] {
            PhaseVectorMatcher.experimentalVectors(
                peaks: peaks, originX: 128, originY: 128, invAngstromPerPixel: 0.01,
                directBeamRadiusInvAngstrom: 0.05, maximumVectorInvAngstrom: reach)
        }
        #expect(vectors(reach: 0).count == 3,
            "0 is not a reach")
        let capped = vectors(reach: 0.4)
        #expect(capped.count == 1,
            "the vector at the reach or beyond survived: \(capped)")
        #expect(abs((capped.first?.x ?? 0) - 0.2) <= 1e-9)
    }

    /// Mutation: the `bestD <= radius` test dropped, so `nearest` always
    /// returns its closest vector however far away it is — which turns every
    /// pattern into a match for every phase.
    ///
    /// THE ANGULAR PROBE IS THE ONE THAT BITES, and the first version of this
    /// test did not have it. Every radial probe is caught by the length-band
    /// prune before the clamp is ever consulted — a reference outside
    /// `|u| ± radius` is never even visited, so `best` stays −1 and the
    /// function returns nil whether the clamp is there or not. Only a
    /// reference at the SAME |q| and a different azimuth reaches the clamp.
    /// Measured 2026-09-12: with only radial probes this test passed against
    /// the mutation it names.
    @Test func testNearestRespectsThePairRadius() throws {
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.5, 0), length: 0.5, relativeIntensity: 1),
            ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.9, 0), length: 0.9, relativeIntensity: 1),
        ]
        let hit = PhaseVectorMatcher.nearest(SIMD2(0.51, 0), in: refs, radius: 0.02)
        #expect(hit?.index == 0)
        #expect(abs((hit?.distance ?? .nan) - 0.01) <= 1e-12)
        #expect(PhaseVectorMatcher.nearest(SIMD2(0.7, 0), in: refs, radius: 0.02) == nil)
        #expect(PhaseVectorMatcher.nearest(SIMD2(0.5, 0), in: [], radius: 0.02) == nil)
        // Same |q| as the first reference, 90° away: inside the length band,
        // 0.707 Å⁻¹ from any reference. Only the clamp refuses it.
        #expect(PhaseVectorMatcher.nearest(SIMD2(0, 0.5), in: refs, radius: 0.02) == nil,
            "a reference at the same |q| but 90° away was accepted as a match")
        // And a near miss in angle alone, just outside the radius.
        let a = 0.05 / 0.5      // arc ≈ chord at this scale
        #expect((PhaseVectorMatcher.nearest(SIMD2(0.5 * cos(a), 0.5 * sin(a)),
                                                in: refs, radius: 0.02)) == nil)
    }

    /// Mutation: the binary search's lower bound computed from `uLen` rather
    /// than `uLen − radius`, which silently skips the reference just inside
    /// the band. Only vectors whose |q| is slightly ABOVE the reference's
    /// expose it, so both directions are probed.
    @Test func testTheLengthBandPruneFindsMatchesOnBothSidesOfTheBand() throws {
        let refs = (1...40).map {
            ReferenceVector(h: $0, k: 0, l: 0, q: SIMD2(Double($0) * 0.05, 0),
                            length: Double($0) * 0.05, relativeIntensity: 1)
        }
        for i in 1...40 {
            let exact = Double(i) * 0.05
            for offset in [-0.015, 0.015] {
                let hit = PhaseVectorMatcher.nearest(
                    SIMD2(exact + offset, 0), in: refs, radius: 0.02)
                #expect(hit?.index == (i - 1),
            "missed the reference at |q| = \(exact) from \(offset)")
            }
        }
    }

    /// Mutation: `minimumVectors` compared with `>` rather than `>=`, or the
    /// matrix branch returning `.notIndexed`. A pattern that is pure matrix is
    /// the commonest case in a real scan — it must not be an unknown.
    @Test func testTooLittleSurvivingMeansTheMatrixByExclusion() throws {
        let library = try smallLibrary()
        let matrixEntry = library.entries[library.matrixEntryIndices[0]]
        let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
        // Every vector taken from the matrix's own reference set.
        let vectors = matrixEntry.vectors.prefix(6).map(\.q)
        let result = PhaseVectorMatcher.classify(
            vectors: Array(vectors), library: library, settings: PhaseVectorSettings(),
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(result.verdict == .matrix)
        #expect(Int(result.phaseIndex) == library.matrixPhaseIndex)
        #expect(Int(result.removedCount) == 6)
        #expect(Int(result.survivingCount) == 0)
    }

    /// Mutation: the fall-back firing at its default of 0; the fraction
    /// compared with `>` so a position exactly at the bar is lost;
    /// `survivingCount` used in place of `removedCount`, which inverts it.
    ///
    /// The fall-back exists because a position is otherwise called matrix only
    /// by EXCLUSION, never because the matrix explains it, so improving
    /// detection turns matrix positions into refusals (measured 2026-09-16,
    /// a recorded open item). It ships OFF, and this test pins that first: the
    /// default must be a refusal, because a silent matrix verdict overstates
    /// the phase fraction this feature exists to report.
    @Test func testTheMatrixFallBackIsOffByDefaultAndFiresOnlyAboveItsFraction() throws {
        let library = try smallLibrary()
        let matrixEntry = library.entries[library.matrixEntryIndices[0]]
        let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
        // Six vectors the matrix explains, plus three on a ring no phase has a
        // reflection near: too many survive for the by-exclusion branch, and
        // nothing clears the candidate guards. Explained fraction 6 / 9 = 0.667.
        let explained = matrixEntry.vectors.prefix(6).map(\.q)
        let unexplainable = (0..<3).map { i -> SIMD2<Double> in
            let a = Double(i) * .pi / 4
            return SIMD2(0.313 * cos(a), 0.313 * sin(a))
        }
        let vectors = Array(explained) + unexplainable

        var settings = PhaseVectorSettings()
        #expect(settings.matrixFallbackExplainedFraction == 0,
            "the fall-back must ship off -- it converts an honest refusal into a positive matrix claim")
        let refused = PhaseVectorMatcher.classify(
            vectors: vectors, library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(refused.verdict == .notIndexed,
            "at the default the position must be refused, not called matrix")
        #expect(Int(refused.removedCount) == 6)
        #expect(Int(refused.survivingCount) == 3)

        settings.matrixFallbackExplainedFraction = 0.5
        let taken = PhaseVectorMatcher.classify(
            vectors: vectors, library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(taken.verdict == .matrix,
            "0.667 explained is above the 0.5 bar, so the matrix takes it")
        #expect(Int(taken.phaseIndex) == library.matrixPhaseIndex)
        #expect(Int(taken.entryIndex) == (-1))

        settings.matrixFallbackExplainedFraction = 0.8
        let stillRefused = PhaseVectorMatcher.classify(
            vectors: vectors, library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(stillRefused.verdict == .notIndexed,
            "0.667 explained is below the 0.8 bar -- the fraction must be compared, not merely tested for being armed")
    }

    /// Mutation: the `notIndexedAboveInvAngstrom` comparison inverted or
    /// removed, so an unexplained pattern is forced into whichever phase
    /// happened to score least badly. The refusal IS the feature. Also: the
    /// score's mean replaced by its max (case 1b).
    ///
    /// TWO CASES, and the first version of this test had only the second.
    /// Vectors pointing nowhere are refused by the ELIGIBILITY guards
    /// (`minimumMatchedVectors`, `chanceMatchMultiple`) long before the
    /// threshold is consulted, so that case passes against the mutation it
    /// names — measured 2026-09-12. What reaches the threshold is a pattern
    /// that matches a phase's vectors well enough to be eligible and not well
    /// enough to be called: an offset INSIDE the pair radius and ABOVE the
    /// not-indexed distance.
    @Test func testUnexplainedVectorsAreNotIndexedRatherThanLabelled() throws {
        let library = try smallLibrary()
        let settings = PhaseVectorSettings()
        let scratch = PhaseVectorMatcher.Scratch(capacity: 256)

        // Case 1 — eligible, but too far to call. Every vector is a real β″
        // reference displaced by 0.018 Å⁻¹: inside the 0.020 pair radius, above
        // the 0.015 not-indexed distance (0.75 of the radius since 2026-09-15).
        let entry = library.entries[library.candidateEntryIndices[0]]
        #expect(settings.pairRadiusInvAngstrom > 0.018)
        #expect(settings.notIndexedAboveInvAngstrom < 0.018)
        let displaced = entry.vectors.prefix(10).map { $0.q + SIMD2(0.018, 0) }
        let nearMiss = PhaseVectorMatcher.classify(
            vectors: Array(displaced), library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(nearMiss.verdict == .notIndexed,
            "a pattern 0.018 Å⁻¹ off every reference was called a phase")
        #expect(nearMiss.matchedCount > 0,
            "the near-miss case never reached the threshold")
        #expect(nearMiss.score.isFinite,
            "a refused position must still carry the numbers behind the refusal")
        #expect(nearMiss.phaseIndex == (-1))
        #expect(nearMiss.entryIndex == (-1))

        // Case 1b — NON-uniform residuals (a review of 2026-09-15: every case above
        // displaces every vector identically, so a `max` or a median in place
        // of the mean was invisible). Half the vectors 0.012 off, half 0.017
        // off: mean 0.0145 is under the 0.015 verdict distance, the largest is
        // over it. The mean is the rule, so this is called.
        let mixed = entry.vectors.prefix(10).enumerated().map { i, v in
            v.q + SIMD2(i < 5 ? 0.012 : 0.017, 0)
        }
        let called = PhaseVectorMatcher.classify(
            vectors: Array(mixed), library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(called.verdict == .indexed,
            "a pattern whose MEAN residual is under the verdict distance was refused -- mutation: the score is not the mean (max or median)")

        // Case 2 — vectors on a ring no phase has a reflection near. Refused
        // by the eligibility guards, which is a different mechanism.
        let elsewhere = (0..<8).map { i -> SIMD2<Double> in
            let a = Double(i) * .pi / 4
            return SIMD2(0.313 * cos(a), 0.313 * sin(a))
        }
        let unknown = PhaseVectorMatcher.classify(
            vectors: elsewhere, library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(unknown.verdict == .notIndexed)
        #expect(unknown.phaseIndex == (-1))
    }

    /// Mutation: `classify` returning `.notIndexed` for an empty pattern.
    /// "No peaks here" and "peaks I cannot explain" are different facts about
    /// the specimen and must not share a colour on the map.
    @Test func testAnEmptyPatternIsNoDataNotNotIndexed() throws {
        let library = try smallLibrary()
        let result = PhaseVectorMatcher.classify(
            vectors: [], library: library, settings: PhaseVectorSettings(),
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 64))
        #expect(result.verdict == .noData)
        #expect(result.survivingCount == 0)
    }

    /// Mutation: `score` averaging over MATCHED vectors instead of unique
    /// references. Two experimental vectors claiming one reference must count
    /// once, at the better distance — otherwise a split peak halves the score
    /// of whichever phase happens to sit between its two halves.
    @Test func testScoreAveragesOverUniqueReferencesAtTheirBestDistance() throws {
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.5, 0), length: 0.5, relativeIntensity: 1),
        ]
        let entry = PhaseOrientationReference(phaseIndex: 1, zoneAxis: SIMD3(0, 0, 1),
                                              inPlaneRotationRad: 0, vectors: refs)
        let scored = PhaseVectorMatcher.score(
            vectors: [SIMD2(0.504, 0), SIMD2(0.512, 0)], against: entry,
            pairRadius: 0.02, scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(scored?.matched == 2)
        #expect(scored?.uniqueReferences == 1)
        #expect(abs((scored?.score ?? .nan) - 0.004) <= 1e-12)
    }

    /// Mutation: `Scratch.bestPerRef` not reset after a call, so the second
    /// entry scored with the same scratch inherits the first's distances.
    /// Nothing else in the suite would notice — every other test scores once.
    @Test func testScratchIsReusableWithoutCarryingStateBetweenCalls() throws {
        let entry = PhaseOrientationReference(
            phaseIndex: 1, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: [ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.5, 0),
                                      length: 0.5, relativeIntensity: 1)])
        let scratch = PhaseVectorMatcher.Scratch(capacity: 8)
        let first = PhaseVectorMatcher.score(vectors: [SIMD2(0.5, 0)], against: entry,
                                             pairRadius: 0.02, scratch: scratch)
        #expect(abs((first?.score ?? .nan) - 0) <= 1e-15)
        let second = PhaseVectorMatcher.score(vectors: [SIMD2(0.515, 0)], against: entry,
                                              pairRadius: 0.02, scratch: scratch)
        #expect(abs((second?.score ?? .nan) - 0.015) <= 1e-12)
        // 90° away at the same |q|: inside the length band, so this reaches the
        // radius clamp rather than being pruned before it (see
        // `testNearestRespectsThePairRadius`).
        let third = PhaseVectorMatcher.score(vectors: [SIMD2(0, 0.5)], against: entry,
                                             pairRadius: 0.02, scratch: scratch)
        #expect(third == nil,
            "a reference outside the radius still counted as matched")
    }

    /// Mutation: `map` accepting a scan whose peak count does not
    /// match `scanWidth * scanHeight`, which would index out of bounds on the
    /// first short row rather than refuse.
    @Test func testMapRefusesAMalformedScan() throws {
        let library = try smallLibrary()
        let malformed = ScanFixture(scanWidth: 4, scanHeight: 4, peaks: [[], [], []])
        #expect(PhaseVectorMatcher.map(
            scan: malformed, library: library, settings: PhaseVectorSettings(),
            originX: 0, originY: 0, invAngstromPerPixel: 0.01) == nil)
    }

    /// Mutation: BOTH `guard cancellation?.isCancelled != true` lines in `map`
    /// removed, so a cancelled run returns a half-filled map presented as a
    /// finished one.
    ///
    /// Both, and stated that way because it was measured: removing either one
    /// alone leaves this test passing (2026-09-12), because the other still
    /// returns nil. They are redundant BY DESIGN — one saves the scan's work,
    /// one guarantees the answer — and a test with a pre-cancelled token can
    /// only see that at least one survives. What it does guarantee is the
    /// property that matters to a caller: a cancelled run never hands back a
    /// partial map wearing the look of a finished one.
    @Test func testACancelledRunReturnsNothingRatherThanAPartialMap() throws {
        let library = try smallLibrary()
        let scan = ScanFixture(scanWidth: 2, scanHeight: 2,
                                 peaks: Array(repeating: [], count: 4))
        #expect(PhaseVectorMatcher.map(
            scan: scan, library: library, settings: PhaseVectorSettings(),
            originX: 0, originY: 0, invAngstromPerPixel: 0.01, isCancelled: { true }) == nil)
    }

    /// Mutation: any coordinate of the published β″ structure, or the C2/m
    /// expansion collapsing a 4i site to 2. The cell content is what the paper
    /// states and is the only independent check on a hand-entered structure.
    @Test func testBetaDoublePrimeMatchesThePublishedCellContent() throws {
        let beta = CrystalStructure.betaDoublePrime
        #expect(beta.sites.count == 22)
        #expect((beta.sites.filter { $0.z == 12 }.count) == 10)   // Mg
        #expect((beta.sites.filter { $0.z == 14 }.count) == 12)   // Si
        #expect(abs(beta.cell.a - 15.16) <= 1e-12)
        #expect(abs(beta.cell.b - 4.05) <= 1e-12)
        #expect(abs(beta.cell.c - 6.74) <= 1e-12)
        #expect(abs(beta.cell.betaDeg - 105.3) <= 1e-12)
        // No two sites coincide — a collapsed expansion would double an atom
        // and multiply straight into every structure factor.
        for i in beta.sites.indices {
            for j in (i + 1)..<beta.sites.count {
                let d = beta.sites[i].fractional - beta.sites[j].fractional
                #expect(simd_length(d) > 1e-6)
            }
        }
        #expect(PhaseReferenceLibrary.unsupportedElements(of: beta).isEmpty)

        // THE CELL CONTENT DOES NOT PIN THE COORDINATES, and the first version
        // of this test stopped above. Measured 2026-09-12: changing Si3's z
        // from 0.617 to 0.671 — a transposed pair of digits, the likeliest
        // transcription error there is — left every assertion above passing,
        // and the gated measuring tool too, because |q(h0l)| depends on the CELL and
        // not on where the atoms sit inside it. So the six published sites are
        // pinned directly against Andersen et al. 1998, Table 3 set 3: a
        // reader can check these six rows against the paper, which is the only
        // independent check a hand-entered structure can have.
        let published: [(Int, Double, Double)] = [
            (12, 0.0,    0.0),      // Mg1, 2a
            (12, 0.3459, 0.089),    // Mg2
            (12, 0.430,  0.652),    // Mg3
            (14, 0.0565, 0.649),    // Si1
            (14, 0.1885, 0.224),    // Si2
            (14, 0.2171, 0.617),    // Si3
        ]
        for (z, x, zc) in published {
            let present = beta.sites.contains {
                $0.z == z && abs($0.fractional.x - x) < 1e-12
                    && abs($0.fractional.y) < 1e-12 && abs($0.fractional.z - zc) < 1e-12
            }
            #expect(present,
            "the published site Z=\(z) (\(x), 0, \(zc)) is not in the cell")
        }

        // And one physical consequence, so a coordinate that is wrong in a way
        // the list above would not catch (a whole site replaced, say) still
        // has somewhere to fail: no two atoms may sit closer than a real Mg-Si
        // or Si-Si contact.
        var shortest = Double.infinity
        for i in beta.sites.indices {
            for j in (i + 1)..<beta.sites.count {
                for dx in -1...1 {
                    for dz in -1...1 {
                        var d = beta.sites[i].fractional - beta.sites[j].fractional
                        d.x += Double(dx); d.z += Double(dz)
                        let cart = d.x * beta.cell.realLattice[0] + d.y * beta.cell.realLattice[1]
                            + d.z * beta.cell.realLattice[2]
                        shortest = min(shortest, simd_length(cart))
                    }
                }
            }
        }
        #expect(shortest > 2.0,
            "shortest interatomic contact \(shortest) Å is not chemistry")
        #expect(shortest < 3.2,
            "shortest interatomic contact \(shortest) Å — the cell is not bonded")
    }

    /// Mutation: `classify` reading `pairRadiusInvAngstrom` for the matrix
    /// removal, or `matrixToleranceInvAngstrom` for the candidate scoring.
    /// Both ship at 0.02, so the conflation was invisible to every check
    /// (a review, mutation E, 2026-09-12) — this fixture sets them 0.005 and
    /// 0.020 and displaces every vector by 0.012, between the two.
    @Test func testMatrixRemovalAndCandidateScoringReadTheirOwnRadii() throws {
        let library = try smallLibrary()
        var settings = PhaseVectorSettings()
        settings.matrixToleranceInvAngstrom = 0.005
        settings.pairRadiusInvAngstrom = 0.020
        let scratch = PhaseVectorMatcher.Scratch(capacity: 256)
        let matrixEntry = library.entries[library.matrixEntryIndices[0]]
        let candidate = library.entries[library.candidateEntryIndices[0]]
        let d = SIMD2(0.012, 0.0)
        let offMatrix = matrixEntry.vectors.prefix(6).map { $0.q + d }
        let offCandidate = candidate.vectors.prefix(10).map { $0.q + d }
        let result = PhaseVectorMatcher.classify(
            vectors: offMatrix + offCandidate, library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(result.removedCount == 0,
            "a vector 0.012 off the matrix, with a 0.005 tolerance, was removed: matrix removal read the pair radius")
        #expect(result.survivingCount == 16)
        #expect(result.matchedCount >= 3,
            "candidate vectors 0.012 off, within a 0.020 pair radius, were not matched: scoring read the matrix tolerance")
    }

    // MARK: Fixtures

    private func coarseSettings() -> PhaseReferenceSettings {
        var s = PhaseReferenceSettings()
        s.kMaxInvAngstrom = 1.0
        s.inPlaneStepDeg = 90
        return s
    }

    private func matrixAluminium() -> PhaseDefinition {
        PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum,
                        role: .matrix, zoneAxes: [SIMD3(0, 0, 1)])
    }

    private func candidateBeta() -> PhaseDefinition {
        PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime,
                        role: .candidate, zoneAxes: [SIMD3(0, 1, 0)])
    }

    /// An exact tie between phases is broken by the lower phase index, not by the process's hash seed.
    /// Eight identical candidate phases score identically on a pattern of one of them; the label
    /// must be the first, the runner-up the second. Under the source (a Dictionary sorted by score alone) both hold in about
    /// 1 run in 56. Mutation: the phase-index tie-break dropped.
    @Test func anExactTieBetweenPhasesGoesToTheLowerPhaseIndex() throws {
        let beta = candidateBeta()
        let copies = (0..<8).map {
            PhaseDefinition(id: "\(beta.id)-\($0)", displayName: "\(beta.displayName) \($0)", crystal: beta.crystal,
                            role: .candidate, zoneAxes: beta.zoneAxes,
                            orientationRelationships: beta.orientationRelationships,
                            excitationSlabInvAngstrom: beta.excitationSlabInvAngstrom)
        }
        let library = try PhaseReferenceLibrary.build(phases: [matrixAluminium()] + copies, settings: coarseSettings())
        let entry = library.entries[library.candidateEntryIndices[0]]
        let result = PhaseVectorMatcher.classify(
            vectors: entry.vectors.prefix(10).map(\.q), library: library, settings: PhaseVectorSettings(),
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 256))
        #expect(result.verdict == .indexed)
        #expect(result.phaseIndex == 1)
        #expect(result.runnerUpPhaseIndex == 2)
    }

    private func smallLibrary() throws -> PhaseReferenceLibrary {
        try PhaseReferenceLibrary.build(phases: [matrixAluminium(), candidateBeta()],
                                        settings: coarseSettings())
    }
}
