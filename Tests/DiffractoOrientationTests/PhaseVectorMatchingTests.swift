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

    /// The orientation relationship, stated as parallel lattice vectors
    /// (2026-09-15, replacing a library-frame degree list; the diagnostic record
    /// is kept with the open items):
    /// `PhaseDefinition.orientationRelationships` restricts a candidate
    /// entry to ones whose DERIVED azimuth (`PhaseVectorMatcher.
    /// projectedAzimuth`, from each entry's own zone and rotation) agrees
    /// with a listed pair's, relative to the fitted matrix entry. Mutations
    /// this names: the 180° fold dropped (a candidate 180° from the listed
    /// pair — equally valid, since a flat-Ewald ZOLZ excites g and −g alike
    /// — refused); the matrix azimuth ignored (as if `matrixRad` were
    /// always 0); the candidate zone's own `detectorBasis` replaced by the
    /// matrix's; a pair whose vector is not in an entry's own zone treated
    /// as applying to it anyway instead of being skipped.
    @Test func testTheOrientationRelationshipIsStatedAsParallelVectorsAndDerivedPerZone() throws {
        // --- 1. `LatticeVector.cartesian` ---------------------------------
        let al = CrystalStructure.aluminum
        let g200 = LatticeVector.plane(SIMD3(2, 0, 0)).cartesian(in: al)
        #expect(abs(simd_length(g200) - (2 / al.cell.a)) <= 1e-9,
            "(200)Al should have |g| = 2/a")
        let d100 = LatticeVector.direction(SIMD3(1, 0, 0)).cartesian(in: al)
        #expect(abs((simd_length(d100 - SIMD3(al.cell.a, 0, 0))) - 0) <= 1e-9,
            "[100]Al should be a·x̂")

        // --- 2. `projectedAzimuth` -----------------------------------------
        let zoneAxis001 = try #require(PhaseReferenceLibrary.cartesianZoneAxis(SIMD3(0, 0, 1), crystal: al),
            "Al [001] should not be degenerate")
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: zoneAxis001)
        let expectedAzimuth200 = atan2(simd_dot(g200, basis.columns.1),
                                       simd_dot(g200, basis.columns.0))
        let azimuth200 = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: g200, zoneAxis: zoneAxis001, inPlaneRotationRad: 0),
            "(200) is in the [001] zone and must project")
        #expect(abs(azimuth200 - expectedAzimuth200) <= 1e-9,
            "the (200) azimuth must come from the entry's own detectorBasis, not a hardcoded angle")

        let g002 = LatticeVector.plane(SIMD3(0, 0, 2)).cartesian(in: al)
        #expect(PhaseVectorMatcher.projectedAzimuth(
            of: g002, zoneAxis: zoneAxis001, inPlaneRotationRad: 0) == nil,
            "(002) lies along the [001] zone axis -- not in the zone, must not project")

        let theta = 30 * Double.pi / 180
        let rotatedAzimuth = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: g200, zoneAxis: zoneAxis001, inPlaneRotationRad: theta))
        var advance = (rotatedAzimuth - azimuth200).truncatingRemainder(dividingBy: 2 * .pi)
        if advance < 0 { advance += 2 * .pi }
        #expect(abs(advance - theta) <= 1e-9,
            "rotating the entry by θ should advance its azimuth by +θ")

        // --- 3. `orientationConsistent` on its own, with the expected
        // answer computed independently of the function under test --------
        //
        // 3a. Candidate and matrix share the SAME crystal, zone axis and
        // relationship vector, so the azimuth difference `orientationConsistent`
        // must compare is EXACTLY (candidate rotation − matrix rotation), no
        // crystal geometry to derive by hand.
        let sameGeometryRelationship = OrientationRelationship(
            candidate: .plane(SIMD3(2, 0, 0)), matrix: .plane(SIMD3(2, 0, 0)))
        func sameZoneEntry(rotationDeg: Double) -> PhaseOrientationReference {
            PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1),
                                      inPlaneRotationRad: rotationDeg * .pi / 180, vectors: [])
        }
        #expect(PhaseVectorMatcher.orientationConsistent(
            candidate: sameZoneEntry(rotationDeg: 8), candidateCrystal: al,
            matrix: sameZoneEntry(rotationDeg: 0), matrixCrystal: al,
            relationships: [sameGeometryRelationship], toleranceDeg: 10),
            "8° apart, within tolerance 10, was refused")
        #expect(!(PhaseVectorMatcher.orientationConsistent(
            candidate: sameZoneEntry(rotationDeg: 12), candidateCrystal: al,
            matrix: sameZoneEntry(rotationDeg: 0), matrixCrystal: al,
            relationships: [sameGeometryRelationship], toleranceDeg: 10)),
            "12° apart, outside tolerance 10, was allowed")
        #expect(PhaseVectorMatcher.orientationConsistent(
            candidate: sameZoneEntry(rotationDeg: 171), candidateCrystal: al,
            matrix: sameZoneEntry(rotationDeg: 0), matrixCrystal: al,
            relationships: [sameGeometryRelationship], toleranceDeg: 10),
            "171° apart folds to 9° from 180° -- mutation: the 180° fold dropped")
        #expect(!(PhaseVectorMatcher.orientationConsistent(
            candidate: sameZoneEntry(rotationDeg: 169), candidateCrystal: al,
            matrix: sameZoneEntry(rotationDeg: 0), matrixCrystal: al,
            relationships: [sameGeometryRelationship], toleranceDeg: 10)),
            "169° apart folds to 11° from 180°, outside tolerance")
        #expect(PhaseVectorMatcher.orientationConsistent(
            candidate: sameZoneEntry(rotationDeg: 95), candidateCrystal: al,
            matrix: sameZoneEntry(rotationDeg: 90), matrixCrystal: al,
            relationships: [sameGeometryRelationship], toleranceDeg: 10),
            "a candidate 5° past a matrix at 90° was refused -- mutation: the matrix azimuth ignored")

        // 3b. Candidate and matrix now use DIFFERENT zone axes (still the
        // same crystal, to keep this arithmetic rather than crystallographic):
        // the expected offset is computed by calling `projectedAzimuth`
        // directly -- already verified in part 2 above -- not by calling
        // `orientationConsistent` and trusting whatever it says.
        let zoneAxis100 = try #require(PhaseReferenceLibrary.cartesianZoneAxis(SIMD3(1, 0, 0), crystal: al))
        let candidateVectorB = LatticeVector.plane(SIMD3(0, 2, 0)).cartesian(in: al)
        let matrixVectorB = LatticeVector.plane(SIMD3(2, 0, 0)).cartesian(in: al)
        let az0cB = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: candidateVectorB, zoneAxis: zoneAxis100, inPlaneRotationRad: 0)) * 180 / .pi
        let az0mB = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: matrixVectorB, zoneAxis: zoneAxis001, inPlaneRotationRad: 0)) * 180 / .pi
        let differentZoneRelationship = OrientationRelationship(
            candidate: .plane(SIMD3(0, 2, 0)), matrix: .plane(SIMD3(2, 0, 0)))
        let thetaMB = 20.0
        // Solved so the TRUE azimuth difference is exactly 0: if the
        // candidate's own [100] basis were replaced by the matrix's [001]
        // basis, az0cB would be wrong and this would stop landing at 0.
        let thetaCB = (thetaMB + az0mB - az0cB).truncatingRemainder(dividingBy: 360)
        #expect(PhaseVectorMatcher.orientationConsistent(
            candidate: PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(1, 0, 0),
                                                 inPlaneRotationRad: thetaCB * .pi / 180, vectors: []),
            candidateCrystal: al,
            matrix: PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1),
                                              inPlaneRotationRad: thetaMB * .pi / 180, vectors: []),
            matrixCrystal: al,
            relationships: [differentZoneRelationship], toleranceDeg: 10),
            "the exact azimuth match (solved from projectedAzimuth's own frames) was refused -- mutation: the candidate zone's own detectorBasis replaced by the matrix's")
        #expect(!(PhaseVectorMatcher.orientationConsistent(
            candidate: PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(1, 0, 0),
                                                 inPlaneRotationRad: (thetaCB + 45) * .pi / 180, vectors: []),
            candidateCrystal: al,
            matrix: PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1),
                                              inPlaneRotationRad: thetaMB * .pi / 180, vectors: []),
            matrixCrystal: al,
            relationships: [differentZoneRelationship], toleranceDeg: 10)),
            "45° off the exact match was allowed")

        // --- 4. `map` end to end: θ′ edge-on [100] as candidate on Al
        // [001] matrix, inPlaneStepDeg = 45, the relationship
        // (002)θ′ ∥ (200)Al. The candidate entry AT the OR angle must be
        // indexed; a neighbour 45° off must not -- and the same holds with
        // the matrix fitted away from 0°, because a grain's rotation is not
        // always 0 (a review of 2026-09-15: a previous test fit the matrix at 0°
        // throughout and so never caught the matrix azimuth being ignored).
        // This is an INTEGRATION check that the pieces compose through
        // `map`; the mutations above are what part 3 independently catches.
        let thetaPrime = makeCrystal(a: 4.04, b: 4.04, c: 5.80, sites: [
            AtomSite(z: 13, fractional: [0.000000, 0.000000, 0.000000]),
            AtomSite(z: 13, fractional: [0.500000, 0.500000, 0.500000]),
            AtomSite(z: 13, fractional: [0.000000, 0.000000, 0.500000]),
            AtomSite(z: 13, fractional: [0.500000, 0.500000, 0.000000]),
            AtomSite(z: 29, fractional: [0.000000, 0.500000, 0.250000]),
            AtomSite(z: 29, fractional: [0.500000, 0.000000, 0.750000]),
        ])
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 0.8
        reference.inPlaneStepDeg = 45
        let relationship = OrientationRelationship(
            candidate: .plane(SIMD3(0, 0, 2)), matrix: .plane(SIMD3(2, 0, 0)))
        let matrixPhase = PhaseDefinition(id: "al", displayName: "Al", crystal: al,
                                          role: .matrix, zoneAxes: [SIMD3(0, 0, 1)])
        let candidatePhase = PhaseDefinition(
            id: "theta-edge", displayName: "θ′ edge-on", crystal: thetaPrime,
            role: .candidate, zoneAxes: [SIMD3(1, 0, 0)],
            orientationRelationships: [relationship])
        let library = try PhaseReferenceLibrary.build(
            phases: [matrixPhase, candidatePhase], settings: reference)

        var fixtureSettings = PhaseVectorSettings()
        fixtureSettings.directBeamRadiusInvAngstrom = 0
        fixtureSettings.orientationRelationshipToleranceDeg = 10
        let scale = 0.005

        func braggFor(entryIndex: Int, in aLibrary: PhaseReferenceLibrary) -> ScanFixture {
            let peaks: [SIMD2<Double>] = aLibrary.entries[entryIndex].vectors.map {
                pixel(Float($0.q.x / scale), Float($0.q.y / scale))
            }
            return ScanFixture(scanWidth: 1, scanHeight: 1, peaks: [peaks])
        }

        func check(matrixDeg: Double) throws {
            let matrixIndex = try #require(library.matrixEntryIndices.first {
                abs(library.entries[$0].inPlaneRotationRad * 180 / .pi - matrixDeg) < 1e-6
            },
            "no matrix entry at \(matrixDeg)° in the fixture")
            let matrixEntry = library.entries[matrixIndex]

            let passing = library.candidateEntryIndices.filter { index in
                PhaseVectorMatcher.orientationConsistent(
                    candidate: library.entries[index], candidateCrystal: thetaPrime,
                    matrix: matrixEntry, matrixCrystal: al,
                    relationships: [relationship], toleranceDeg: 10)
            }
            let hit = try #require(passing.first,
            "no candidate rotation of the 45° sweep satisfies the OR at matrix \(matrixDeg)°")
            let hitDeg = library.entries[hit].inPlaneRotationRad * 180 / .pi
            let miss = try #require(library.candidateEntryIndices.first { index in
                guard !passing.contains(index) else { return false }
                let deg = library.entries[index].inPlaneRotationRad * 180 / .pi
                let delta = abs(deg - hitDeg).truncatingRemainder(dividingBy: 360)
                let wrapped = min(delta, 360 - delta)
                return wrapped > 44 && wrapped < 46
            },
            "no 45°-off, non-passing neighbour of the OR entry at matrix \(matrixDeg)°")

            let hitMap = try #require(PhaseVectorMatcher.map(
                scan: braggFor(entryIndex: hit, in: library), library: library, settings: fixtureSettings,
                originX: 0, originY: 0, invAngstromPerPixel: scale, matrixEntryIndex: matrixIndex))
            #expect(hitMap.results[0].verdict == .indexed,
            "the candidate entry AT the OR angle (matrix \(matrixDeg)°) was not indexed")

            let missMap = try #require(PhaseVectorMatcher.map(
                scan: braggFor(entryIndex: miss, in: library), library: library, settings: fixtureSettings,
                originX: 0, originY: 0, invAngstromPerPixel: scale, matrixEntryIndex: matrixIndex))
            #expect(missMap.results[0].verdict != .indexed,
            "the candidate entry 45° off the OR angle (matrix \(matrixDeg)°) was indexed anyway")
        }
        try check(matrixDeg: 0)
        try check(matrixDeg: 45)

        // --- 5. An empty list is free -------------------------------------
        let freePhase = PhaseDefinition(
            id: "theta-edge", displayName: "θ′ edge-on", crystal: thetaPrime,
            role: .candidate, zoneAxes: [SIMD3(1, 0, 0)], orientationRelationships: [])
        let freeLibrary = try PhaseReferenceLibrary.build(
            phases: [matrixPhase, freePhase], settings: reference)
        let freeMatrixIndex = try #require(freeLibrary.matrixEntryIndices.first {
            freeLibrary.entries[$0].inPlaneRotationRad.magnitude < 1e-9
        })
        // Any grid rotation will do -- 45°, which the constrained sweep
        // above found to be off the OR angle at matrix 0°.
        let freeCandidateIndex = try #require(freeLibrary.candidateEntryIndices.first {
            abs(freeLibrary.entries[$0].inPlaneRotationRad * 180 / .pi - 45) < 1e-6
        })
        let freeMap = try #require(PhaseVectorMatcher.map(
            scan: braggFor(entryIndex: freeCandidateIndex, in: freeLibrary), library: freeLibrary,
            settings: fixtureSettings, originX: 0, originY: 0, invAngstromPerPixel: scale,
            matrixEntryIndex: freeMatrixIndex))
        #expect(freeMap.results[0].verdict == .indexed,
            "an empty orientation-relationship list should not constrain the map")

        // --- 6. A pair whose vector is not in the entry's zone does not
        // constrain it: (200)θ′ lies exactly along the [100] zone axis, so
        // it can never be projected, and the phase must stay free rather
        // than refusing every entry.
        let outOfZoneRelationship = OrientationRelationship(
            candidate: .plane(SIMD3(2, 0, 0)), matrix: .plane(SIMD3(2, 0, 0)))
        let outOfZonePhase = PhaseDefinition(
            id: "theta-edge", displayName: "θ′ edge-on", crystal: thetaPrime,
            role: .candidate, zoneAxes: [SIMD3(1, 0, 0)],
            orientationRelationships: [outOfZoneRelationship])
        let outOfZoneLibrary = try PhaseReferenceLibrary.build(
            phases: [matrixPhase, outOfZonePhase], settings: reference)
        let outOfZoneMatrixIndex = try #require(outOfZoneLibrary.matrixEntryIndices.first {
            outOfZoneLibrary.entries[$0].inPlaneRotationRad.magnitude < 1e-9
        })
        let outOfZoneCandidateIndex = try #require(outOfZoneLibrary.candidateEntryIndices.first {
            abs(outOfZoneLibrary.entries[$0].inPlaneRotationRad * 180 / .pi - 45) < 1e-6
        })
        let outOfZoneMap = try #require(PhaseVectorMatcher.map(
            scan: braggFor(entryIndex: outOfZoneCandidateIndex, in: outOfZoneLibrary), library: outOfZoneLibrary,
            settings: fixtureSettings, originX: 0, originY: 0, invAngstromPerPixel: scale,
            matrixEntryIndex: outOfZoneMatrixIndex))
        #expect(outOfZoneMap.results[0].verdict == .indexed,
            "a relationship pair not in either zone must not constrain the entry")

        // --- a review of 2026-09-15: every plane above has the same azimuth in
        // either crystal (orthogonal cells, (h00)-type), so the two crystals
        // SWAPPED inside `orientationConsistent` was invisible. A hexagonal
        // candidate's (110) sits 60° from a in the basal plane, a cubic
        // matrix's 45°, so the derived offset is 15° one way and −15° the
        // other; the swap flips its sign. (The exact-tolerance boundary is
        // real-valued after atan2 and is not pinned: `<` vs `<=` differ only
        // at equality, which no fixture reaches exactly.)
        let hexagonal = makeCrystal(a: 4.94775, b: 4.94775, c: 14.14499,
                                alphaDeg: 90, betaDeg: 90, gammaDeg: 120,
                                sites: [AtomSite(z: 13, fractional: SIMD3(0, 0, 0))])
        let cubic = CrystalStructure.aluminum
        let basal = SIMD3<Int>(0, 0, 1)
        let up = SIMD3<Double>(0, 0, 1)           // both cells carry c along z
        let azHex = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: LatticeVector.plane(SIMD3(1, 1, 0)).cartesian(in: hexagonal),
            zoneAxis: up, inPlaneRotationRad: 0))
        let azCubic = try #require(PhaseVectorMatcher.projectedAzimuth(
            of: LatticeVector.plane(SIMD3(1, 1, 0)).cartesian(in: cubic),
            zoneAxis: up, inPlaneRotationRad: 0))
        let derivedOffset = azCubic - azHex
        #expect((abs(derivedOffset) * 180 / .pi) > 12,
            "the two (110) azimuths must differ by more than the tolerance for this to bite")
        let oneOneZero = [OrientationRelationship(candidate: .plane(SIMD3(1, 1, 0)),
                                                  matrix: .plane(SIMD3(1, 1, 0)))]
        func hexCandidate(at rad: Double) -> PhaseOrientationReference {
            PhaseOrientationReference(phaseIndex: 1, zoneAxis: basal, inPlaneRotationRad: rad, vectors: [])
        }
        let cubicMatrix = PhaseOrientationReference(phaseIndex: 0, zoneAxis: basal,
                                                    inPlaneRotationRad: 0, vectors: [])
        #expect(PhaseVectorMatcher.orientationConsistent(
            candidate: hexCandidate(at: derivedOffset), candidateCrystal: hexagonal,
            matrix: cubicMatrix, matrixCrystal: cubic, relationships: oneOneZero, toleranceDeg: 10),
            "the offset derived from the two crystals was not accepted")
        #expect(!(PhaseVectorMatcher.orientationConsistent(
            candidate: hexCandidate(at: -derivedOffset), candidateCrystal: hexagonal,
            matrix: cubicMatrix, matrixCrystal: cubic, relationships: oneOneZero, toleranceDeg: 10)),
            "the offset with its sign flipped was accepted -- mutation: the two crystals swapped")
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

    // MARK: - Item K: cross-phase-completeness-guard (off by default)

    /// `PhaseVectorMatcher.crossPhaseWinsOver` in isolation -- the exact
    /// scenario measured (a review, 2026-09-15): a
    /// 2-vector match at 0.004 Å⁻¹ against a 10-vector match of another
    /// phase at 0.012 Å⁻¹. Shipped (mean distance alone) picks the sparse
    /// match; the completeness-aware candidate picks the dense one.
    ///
    /// Mutation this catches: the `completenessAware` guard deleted (the
    /// candidate rule would then always apply, changing shipped behavior),
    /// or `matched != matched` compared with `<` instead of `>` (the dense
    /// match would lose to the sparse one even with the guard on).
    @Test func testCrossPhaseWinsOverPicksTheSparseMatchByDefaultAndTheDenseMatchWithTheGuard() throws {
        let sparsePrecise = (matched: 2, score: 0.004)
        let denseImprecise = (matched: 10, score: 0.012)

        #expect(PhaseVectorMatcher.crossPhaseWinsOver(sparsePrecise, denseImprecise, completenessAware: false),
            "shipped default (mean distance alone): the sparse, precise match must still win")
        #expect(!(PhaseVectorMatcher.crossPhaseWinsOver(denseImprecise, sparsePrecise, completenessAware: false)),
            "the ordering must be antisymmetric: the dense match cannot ALSO outrank the sparse one")

        #expect(PhaseVectorMatcher.crossPhaseWinsOver(denseImprecise, sparsePrecise, completenessAware: true),
            "completeness-aware (off by default): the dense match must win when the flag is on")
        #expect(!(PhaseVectorMatcher.crossPhaseWinsOver(sparsePrecise, denseImprecise, completenessAware: true)),
            "and the sparse match must no longer outrank it")
    }

    /// Equal matched counts always fall back to mean distance, with the
    /// guard on or off -- the tiebreak, not a separate rule. Mutation: the
    /// tiebreak branch skipped when `completenessAware` is true (an equal
    /// matched count would then order arbitrarily by dictionary iteration).
    @Test func testCrossPhaseWinsOverTiesOnMatchedCountFallToMeanDistanceEitherWay() throws {
        let lowerDistance = (matched: 5, score: 0.006)
        let higherDistance = (matched: 5, score: 0.009)
        for completenessAware in [false, true] {
            #expect(PhaseVectorMatcher.crossPhaseWinsOver(
                    lowerDistance, higherDistance, completenessAware: completenessAware),
            "completenessAware=\(completenessAware): equal matched counts must still prefer the closer mean")
        }
    }

    /// `completenessAwareCrossPhaseRanking` ships off -- the item's own
    /// requirement ("moves a per-position verdict... measure and propose,
    /// never merge"). If this ever reads true, `classify` is silently using
    /// the unmeasured candidate ranking on every real dataset.
    @Test func testCompletenessAwareCrossPhaseRankingShipsOff() throws {
        #expect(!(PhaseVectorSettings().completenessAwareCrossPhaseRanking))
    }

    // MARK: - `.knownVariants` (Thronsen et al.'s own rule)

    /// A library built directly from hand-placed reference vectors, bypassing
    /// `PhaseReferenceLibrary.build`'s crystallography so every distance in a
    /// `.knownVariants` test below is exact rather than derived from a
    /// projected cell. Each of `candidateEntries` becomes one
    /// `PhaseOrientationReference` at its own phase index; phase 0 (index
    /// `matrixPhaseIndex`) is reserved for the matrix and carries no entry
    /// unless the test builds one itself and splices it in.
    private func knownVariantsLibrary(
        candidateEntries: [(phaseIndex: Int, vectors: [ReferenceVector])],
        matrixPhaseIndex: Int = 0
    ) -> PhaseReferenceLibrary {
        var phases: [PhaseDefinition] = [
            PhaseDefinition(id: "matrix", displayName: "Matrix", crystal: .aluminum,
                            role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
        ]
        var entries: [PhaseOrientationReference] = []
        for (offset, candidate) in candidateEntries.enumerated() {
            phases.append(PhaseDefinition(id: "c\(offset)", displayName: "C\(offset)",
                                          crystal: .aluminum, role: .candidate,
                                          zoneAxes: [SIMD3(0, 0, 1)]))
            entries.append(PhaseOrientationReference(
                phaseIndex: candidate.phaseIndex, zoneAxis: SIMD3(0, 0, 1),
                inPlaneRotationRad: 0, vectors: candidate.vectors))
        }
        return PhaseReferenceLibrary(phases: phases, settings: PhaseReferenceSettings(),
                                     entries: entries, matrixPhaseIndex: matrixPhaseIndex)
    }

    /// (1) A planted candidate pattern -- exactly a library entry's own
    /// vectors -- is labelled that phase, `.indexed`, at zero residual.
    /// Mutation: `classifyKnownVariants` never reached from `classify`, or
    /// the winner's phase/entry index mixed up.
    @Test func testKnownVariantsPlantedPatternIsLabelledThatPhase() throws {
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.30, 0), length: 0.30, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.40), length: 0.40, relativeIntensity: 1),
            ReferenceVector(h: 1, k: 1, l: 0, q: SIMD2(0.30, 0.40), length: 0.5, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        let result = PhaseVectorMatcher.classify(
            vectors: refs.map(\.q), library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(result.verdict == .indexed)
        #expect(Int(result.phaseIndex) == 1)
        #expect(Int(result.entryIndex) == library.candidateEntryIndices[0])
        #expect(abs(result.score - 0) <= 1e-12)
        #expect(Int(result.matchedCount) == 3)
    }

    /// (2) Survivor count alone decides the matrix branch: 0 and 1 survivors
    /// -> `.matrix`, 2 -> reaches scoring. Mutation: `<=` relaxed to `<`, or
    /// the matrix shortcut firing at 2 survivors too.
    @Test func testKnownVariantsSurvivorCountDecidesMatrixVsScored() throws {
        let matrixEntry = PhaseOrientationReference(
            phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: [ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.4, 0), length: 0.4, relativeIntensity: 1)])
        let candidateRefs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(1.0, 0), length: 1.0, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(1.2, 0), length: 1.2, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: candidateRefs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        let scratch = PhaseVectorMatcher.Scratch(capacity: 8)

        // 0 survivors: the only input vector is the matrix's own.
        let zero = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.4, 0)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        #expect(zero.verdict == .matrix)
        #expect(Int(zero.phaseIndex) == library.matrixPhaseIndex)
        #expect(Int(zero.survivingCount) == 0)
        #expect(Int(zero.removedCount) == 1)

        // 1 survivor: one matches the matrix and is removed, one does not.
        let one = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.4, 0), SIMD2(1.0, 0)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        #expect(one.verdict == .matrix)
        #expect(Int(one.survivingCount) == 1)
        #expect(Int(one.removedCount) == 1)

        // 2 survivors: neither matches the matrix -- must reach scoring,
        // not the by-exclusion shortcut.
        let two = PhaseVectorMatcher.classify(
            vectors: [SIMD2(1.0, 0), SIMD2(1.2, 0)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        #expect(two.verdict == .indexed,
            "2 survivors must reach scoring, not the matrix shortcut")
        #expect(Int(two.phaseIndex) == 1)
        #expect(abs(two.score - 0) <= 1e-12)
    }

    /// (3) The unique-hit denominator: two survivors that both land nearest
    /// the SAME reference spot must score WORSE than two survivors landing
    /// nearest two DISTINCT spots at the same two distances -- cell 16's
    /// `np.unique(ref_tmp[n], axis=0).shape[0]` denominator.
    /// Break-first M1: divide by `surviving.count` instead of the unique-hit
    /// count -- both cases then have denominator 2 and score identically,
    /// so `XCTAssertGreaterThan` below must go red.
    @Test func testKnownVariantsUniqueHitDenominatorPenalisesRepeatedHits() throws {
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.30, 0), length: 0.30, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0.60, 0), length: 0.60, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        let scratch = PhaseVectorMatcher.Scratch(capacity: 8)

        // Both survivors nearest the SAME reference (0.30): distances 0.01
        // and 0.01, unique hits 1 -> score 0.02 / 1 = 0.02.
        let sameSpot = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.29, 0), SIMD2(0.31, 0)], library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        // Survivors nearest two DIFFERENT references at the same two
        // distances: unique hits 2 -> score 0.02 / 2 = 0.01.
        let distinctSpots = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.29, 0), SIMD2(0.61, 0)], library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)

        #expect(sameSpot.verdict == .indexed)
        #expect(distinctSpots.verdict == .indexed)
        #expect(abs(sameSpot.score - 0.02) <= 1e-6)
        #expect(abs(distinctSpots.score - 0.01) <= 1e-6)
        #expect(sameSpot.score > distinctSpots.score,
            "two hits on one reference must score worse than two hits on distinct ones")
    }

    /// (4) The residual cutoff, right at the edge: mean residual just below
    /// `residualCutoffInvAngstrom` is `.indexed`, just above is
    /// `.notIndexed`. Break-first M2: flip the comparison -- both cases
    /// swap verdicts, so this must go red.
    @Test func testKnownVariantsCutoffBoundary() throws {
        let cutoff = PhaseVectorSettings().residualCutoffInvAngstrom
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.30, 0), length: 0.30, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.60), length: 0.60, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        // The cutoff alone: these survivors sit ~0.07 Å⁻¹ from their
        // references, outside the pair radius, so the evidence guard (tested
        // below) would send the under-cutoff case to the matrix.
        settings.knownVariantsMinimumSpecificReflections = 0
        let scratch = PhaseVectorMatcher.Scratch(capacity: 8)

        func result(offset: Double) -> PhaseVectorResult {
            PhaseVectorMatcher.classify(
                vectors: [SIMD2(0.30 + offset, 0), SIMD2(0, 0.60 + offset)],
                library: library, settings: settings, matrixEntry: nil,
                candidateEntryIndices: library.candidateEntryIndices, scratch: scratch)
        }

        let justBelow = result(offset: cutoff - 0.001)
        let justAbove = result(offset: cutoff + 0.001)
        #expect(justBelow.verdict == .indexed,
            "mean residual \(justBelow.score) is under the cutoff \(cutoff)")
        #expect(justAbove.verdict == .notIndexed,
            "mean residual \(justAbove.score) is over the cutoff \(cutoff)")
        // Still filled on a refusal, so the probe can print the distribution
        // of would-be winners among the not-indexed positions.
        #expect(Int(justAbove.phaseIndex) == 1)
        #expect(justAbove.score.isFinite)
    }

    // MARK: - `.knownVariants` evidence guard (a diagnostic run 2026-09-23)

    /// Matrix at (0.4, 0); candidate phase 1 holds a SHARED reflection at
    /// (0.415, 0) (0.015 from the matrix's, inside `matrixToleranceInvAngstrom`)
    /// and a SPECIFIC one at (0, 0.75). `second` places the second survivor.
    private func guardFixture(second: SIMD2<Double>, minimumSpecific: Int)
        -> PhaseVectorResult {
        let matrixEntry = PhaseOrientationReference(
            phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: [ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.4, 0), length: 0.4, relativeIntensity: 1)])
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.415, 0), length: 0.415, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.75), length: 0.75, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        settings.knownVariantsMinimumSpecificReflections = minimumSpecific
        // (0.43, 0) is 0.03 from the matrix's (0.4, 0): it survives removal,
        // and hits the SHARED reference 0.015 away, inside the pair radius.
        return PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.43, 0), second], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
    }

    /// (G1) A winner whose only in-radius hit is a reflection the matrix
    /// shares is the matrix under the guard, and `.indexed` without it. The
    /// second survivor (0, 0.70) is 0.05 from the specific reference, outside
    /// the pair radius: score (0.015 + 0.05) / 2 = 0.0325, under the cutoff.
    /// Break-first: (a) the guard removed, or (b) the shared-with-matrix filter
    /// removed from `specificReflectionCount` (the shared hit then counts), each
    /// leaves the k = 1 case `.indexed`, so this must go red.
    @Test func testKnownVariantsGuardSendsSharedOnlyWinnerToMatrix() throws {
        let guarded = guardFixture(second: SIMD2(0, 0.70), minimumSpecific: 1)
        #expect(guarded.verdict == .matrix)
        #expect(Int(guarded.phaseIndex) == 0,
            "the matrix's phase index")
        #expect(Int(guarded.entryIndex) == (-1),
            "cleared, as in every matrix fallback")
        #expect(abs(guarded.score - 0.0325) <= 1e-6,
            "the rejected winner's score is kept")

        let unguarded = guardFixture(second: SIMD2(0, 0.70), minimumSpecific: 0)
        #expect(unguarded.verdict == .indexed)
        #expect(Int(unguarded.phaseIndex) == 1)
        #expect(abs(unguarded.score - 0.0325) <= 1e-6)
    }

    /// (G2) The boundary: one specific hit ((0, 0.76), 0.01 from (0, 0.75))
    /// passes k = 1 and fails k = 2. Break-first: `<` changed to `<=` in the
    /// guard sends the k = 1 case to the matrix, so this must go red.
    @Test func testKnownVariantsGuardBoundaryAtOneSpecificReflection() throws {
        let atOne = guardFixture(second: SIMD2(0, 0.76), minimumSpecific: 1)
        #expect(atOne.verdict == .indexed,
            "exactly one specific reflection meets k = 1")
        #expect(Int(atOne.phaseIndex) == 1)
        let atTwo = guardFixture(second: SIMD2(0, 0.76), minimumSpecific: 2)
        #expect(atTwo.verdict == .matrix,
            "one specific reflection fails k = 2")
    }

    /// (G3) A refusal stays a refusal: a winner over the residual cutoff with
    /// no specific hit is `.notIndexed`, not the matrix. Second survivor
    /// (0, 0.90): 0.15 from (0, 0.75), score (0.015 + 0.15) / 2 = 0.0825 > 0.07.
    /// Break-first: the guard applied to every verdict (the `.indexed` test
    /// dropped) turns this into `.matrix`, so this must go red.
    @Test func testKnownVariantsGuardLeavesRefusalsAlone() throws {
        let refused = guardFixture(second: SIMD2(0, 0.90), minimumSpecific: 1)
        #expect(refused.verdict == .notIndexed)
        #expect(Int(refused.phaseIndex) == 1,
            "the would-be winner, as for every refusal")
    }

    /// (G5) "Shared with the matrix" is judged at the MATRIX tolerance, not
    /// the pair radius (a review of 2026-09-23: with both at 0.02 no test told them
    /// apart). Pair radius 0.03, matrix tolerance 0.01; the candidate's
    /// (0.42, 0) sits 0.02 from the matrix's (0.4, 0): outside the tolerance,
    /// so SPECIFIC. Survivor (0.43, 0) survives removal (0.03 > 0.01) and hits
    /// it at 0.01; (0, 0.70) is 0.05 from (0, 0.75), outside the radius.
    /// Score (0.01 + 0.05) / 2 = 0.03. Break-first: the shared filter at
    /// `pairRadiusInvAngstrom` makes (0.42, 0) shared, specific 0, matrix — so
    /// this must go red.
    @Test func testKnownVariantsGuardJudgesSharedAtTheMatrixTolerance() throws {
        let matrixEntry = PhaseOrientationReference(
            phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: [ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.4, 0), length: 0.4, relativeIntensity: 1)])
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.42, 0), length: 0.42, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.75), length: 0.75, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        settings.pairRadiusInvAngstrom = 0.03
        settings.matrixToleranceInvAngstrom = 0.01
        let result = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.43, 0), SIMD2(0, 0.70)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(result.verdict == .indexed,
            "one specific reflection at the matrix tolerance")
        #expect(Int(result.phaseIndex) == 1)
        #expect(abs(result.score - 0.03) <= 1e-6)
    }

    /// (G6) The guard judges the WINNER, not the first entry in library
    /// order (a review of 2026-09-23: every guard test had one candidate). Entry 0
    /// (phase 1) has a specific hit, (0.445, 0), but scores (0.015 + 0.10) / 2
    /// = 0.0575. Entry 1 (phase 2) wins at (0.015 + 0.05) / 2 = 0.0325 with
    /// only a shared hit, (0.415, 0), so the guard sends the position to the
    /// matrix. Break-first: the guard evaluated on `ranked[0]` (entry 0,
    /// specific 1) keeps it indexed as phase 2, so this must go red.
    @Test func testKnownVariantsGuardJudgesTheWinnerNotTheFirstEntry() throws {
        let matrixEntry = PhaseOrientationReference(
            phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1), inPlaneRotationRad: 0,
            vectors: [ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.4, 0), length: 0.4, relativeIntensity: 1)])
        let loser = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.445, 0), length: 0.445, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.60), length: 0.60, relativeIntensity: 1),
        ]
        let winner = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.415, 0), length: 0.415, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.75), length: 0.75, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [
            (phaseIndex: 1, vectors: loser), (phaseIndex: 2, vectors: winner),
        ])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        let scratch = PhaseVectorMatcher.Scratch(capacity: 8)
        let result = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.43, 0), SIMD2(0, 0.70)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(result.verdict == .matrix,
            "the winner (phase 2) has no specific reflection")
        #expect(abs(result.score - 0.0325) <= 1e-6,
            "phase 2's score: it was the winner")

        settings.knownVariantsMinimumSpecificReflections = 0
        let unguarded = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.43, 0), SIMD2(0, 0.70)], library: library, settings: settings,
            matrixEntry: matrixEntry, candidateEntryIndices: library.candidateEntryIndices,
            scratch: scratch)
        #expect(Int(unguarded.phaseIndex) == 2,
            "without the guard phase 2 wins")
    }

    /// (G4) Shipped on at 1 (decision of 2026-09-23). Break-first: the
    /// default changed to 0, so this must go red.
    @Test func testKnownVariantsGuardShipsOnAtOne() throws {
        #expect(PhaseVectorSettings().knownVariantsMinimumSpecificReflections == 1)
    }

    /// (5) Partial explanation loses: an entry explaining 2 of 12 survivors
    /// at 0.003 Å⁻¹ and leaving the other 10 far away must lose to an entry
    /// explaining all 12 at 0.008 Å⁻¹ -- the measured θ′ edge-on -> T1
    /// mechanism (a residual-detail record of 2026-09-21:
    /// a Friedel-pair floor lets a 2-of-12 partial match win under `.search`;
    /// `.knownVariants` has no such floor, and summing every survivor is what
    /// makes leaving ten unexplained costly instead of merely absent from the
    /// count). Break-first M3: sum only the matched-within-radius distances
    /// -- the ten far survivors drop out of entry X's sum, its score
    /// collapses to near 0.003, and it wins instead, so this must go red.
    @Test func testKnownVariantsPartialExplanationLoses() throws {
        let survivors = (0..<12).map { SIMD2<Double>(Double($0) * 1.0, 0) }
        // Entry X: two reference vectors. Survivors 0 and 1 land 0.003 from
        // one each; survivors 2...11 are all nearest the SECOND of the two,
        // "far away" as the residual record describes.
        let entryXRefs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.003, 0), length: 0.003, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(1.003, 0), length: 1.003, relativeIntensity: 1),
        ]
        // Entry Y: twelve reference vectors, each 0.008 from its own survivor.
        let entryYRefs = (0..<12).map { i in
            ReferenceVector(h: i, k: 0, l: 0, q: SIMD2(Double(i) * 1.0 + 0.008, 0),
                            length: Double(i) * 1.0 + 0.008, relativeIntensity: 1)
        }
        let library = knownVariantsLibrary(candidateEntries: [
            (phaseIndex: 1, vectors: entryXRefs),
            (phaseIndex: 2, vectors: entryYRefs),
        ])
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        let result = PhaseVectorMatcher.classify(
            vectors: survivors, library: library, settings: settings, matrixEntry: nil,
            candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 16))
        #expect(result.verdict == .indexed)
        #expect(Int(result.phaseIndex) == 2,
            "the entry explaining all 12 survivors must win over one explaining only 2 of them")
        #expect(abs(result.score - 0.008) <= 1e-6)
        #expect(Int(result.runnerUpPhaseIndex) == 1)
    }

    /// (6) A matrix-phase entry is never the winner under this rule, even
    /// when spliced into the candidate list with a perfect (zero-distance)
    /// score -- the matrix is decided only by step b, never scored.
    /// Mutation: the `entry.phaseIndex != library.matrixPhaseIndex` guard
    /// removed from the scoring loop.
    @Test func testKnownVariantsNeverPicksTheMatrixPhaseEvenWhenItWouldScoreBest() throws {
        let matrixRefs = [
            ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.50, 0), length: 0.50, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 2, l: 0, q: SIMD2(0, 0.70), length: 0.70, relativeIntensity: 1),
        ]
        let matrixEntry = PhaseOrientationReference(phaseIndex: 0, zoneAxis: SIMD3(0, 0, 1),
                                                    inPlaneRotationRad: 0, vectors: matrixRefs)
        // A real candidate, 0.01 Å⁻¹ off each spot -- the only legitimate winner.
        let candidateRefs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.51, 0), length: 0.51, relativeIntensity: 1),
            ReferenceVector(h: 0, k: 1, l: 0, q: SIMD2(0, 0.71), length: 0.71, relativeIntensity: 1),
        ]
        let base = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: candidateRefs)])
        // Splice the matrix's own entry into the SAME entries array, so the
        // guard under test is `classify`'s own, not merely a caller that
        // already knew to leave the matrix out.
        let library = PhaseReferenceLibrary(phases: base.phases, settings: base.settings,
                                            entries: [matrixEntry] + base.entries,
                                            matrixPhaseIndex: base.matrixPhaseIndex)
        var settings = PhaseVectorSettings()
        settings.classificationRule = .knownVariants
        // Deliberately pass EVERY entry as a candidate, the matrix's included.
        let result = PhaseVectorMatcher.classify(
            vectors: [SIMD2(0.50, 0), SIMD2(0, 0.70)], library: library, settings: settings,
            matrixEntry: nil, candidateEntryIndices: Array(library.entries.indices),
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(result.verdict == .indexed)
        #expect(Int(result.phaseIndex) == 1,
            "the matrix's own entry scored 0 and must still never win")
        #expect(abs(result.score - 0.01) <= 1e-6)
    }

    /// (7) `.search` is unchanged by `.knownVariants` existing: a hand-placed
    /// fixture's FULL `PhaseVectorResult` is pinned, and the default settings
    /// (`.search`) agree exactly with settings that name `.search` explicitly
    /// -- the dispatch adds a branch, it does not move which branch the
    /// default takes or any field `.search` computes. `PhaseVectorResult`'s
    /// synthesized `==` is not used directly: NaN != NaN would make a
    /// passing runner-up comparison look like a failure.
    @Test func testKnownVariantsRuleDoesNotChangeSearchsOwnResult() throws {
        let refs = [
            ReferenceVector(h: 1, k: 0, l: 0, q: SIMD2(0.10, 0), length: 0.10, relativeIntensity: 1),
            ReferenceVector(h: 2, k: 0, l: 0, q: SIMD2(0.20, 0), length: 0.20, relativeIntensity: 1),
            ReferenceVector(h: 3, k: 0, l: 0, q: SIMD2(0.30, 0), length: 0.30, relativeIntensity: 1),
            ReferenceVector(h: 5, k: 0, l: 0, q: SIMD2(0.50, 0), length: 0.50, relativeIntensity: 1),
        ]
        let library = knownVariantsLibrary(candidateEntries: [(phaseIndex: 1, vectors: refs)])
        let vectors = [SIMD2(0.101, 0.0), SIMD2(0.199, 0.0), SIMD2(0.302, 0.0), SIMD2(0.505, 0.0)]
        var settings = PhaseVectorSettings()
        #expect(settings.classificationRule == .search,
            "the shipped default must stay .search")

        let result = PhaseVectorMatcher.classify(
            vectors: vectors, library: library, settings: settings, matrixEntry: nil,
            candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(result.verdict == .indexed)
        #expect(Int(result.phaseIndex) == 1)
        #expect(Int(result.entryIndex) == library.candidateEntryIndices[0])
        #expect(abs(result.score - 0.00225) <= 1e-6)
        #expect(Int(result.matchedCount) == 4)
        #expect(Int(result.survivingCount) == 4)
        #expect(Int(result.removedCount) == 0)
        #expect(Int(result.runnerUpPhaseIndex) == (-1))
        #expect(result.runnerUpScore.isNaN)

        settings.classificationRule = .search
        let explicit = PhaseVectorMatcher.classify(
            vectors: vectors, library: library, settings: settings, matrixEntry: nil,
            candidateEntryIndices: library.candidateEntryIndices,
            scratch: PhaseVectorMatcher.Scratch(capacity: 8))
        #expect(explicit.verdict == result.verdict)
        #expect(explicit.phaseIndex == result.phaseIndex)
        #expect(explicit.entryIndex == result.entryIndex)
        #expect(explicit.score == result.score)
        #expect(explicit.matchedCount == result.matchedCount)
        #expect(explicit.survivingCount == result.survivingCount)
        #expect(explicit.removedCount == result.removedCount)
    }
}

// MARK: Diagnosis 2026-10-21: the cross-phase ranking and the matrix verdict by exclusion

struct PhaseDiagnosisTests {
    private func library() throws -> PhaseReferenceLibrary {
        var s = PhaseReferenceSettings(); s.kMaxInvAngstrom = 1.0; s.inPlaneStepDeg = 10
        let al = PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)])
        let beta = PhaseDefinition(id: "beta", displayName: "b", crystal: .betaDoublePrime, role: .candidate,
                                   zoneAxes: [SIMD3(0, 1, 0), SIMD3(0, 0, 1)])
        let mg = PhaseDefinition(id: "mg", displayName: "Mg", crystal: CrystalPresets.magnesium, role: .candidate,
                                 zoneAxes: [SIMD3(0, 0, 1), SIMD3(1, 0, 0)])
        return try PhaseReferenceLibrary.build(phases: [al, beta, mg], settings: s)
    }

    /// Wrong labels out of 400 planted patterns (two candidate phases, reference vectors of the true entry with Gaussian
    /// noise `sigma`, a fraction `dropout` removed, two random spurious vectors), per ranking rule. The measurement behind
    /// `completenessAwareCrossPhaseRanking` staying off (2026-10-21): the rule that counts matched vectors first loses at
    /// small noise (a denser library entry takes spurious vectors by chance) and wins only when the noise approaches the pair
    /// radius, so neither ranking dominates and the default stays the paper's mean distance.
    private func wrongLabels(sigma: Double, dropout: Double, aware: Bool) throws -> Int {
        let lib = try library()
        let cands = lib.candidateEntryIndices
        var seed: UInt64 = 777
        func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
        func gauss() -> Double { (-2 * log(max(rnd(), 1e-12))).squareRoot() * cos(2 * .pi * rnd()) }
        var settings = PhaseVectorSettings(); settings.completenessAwareCrossPhaseRanking = aware
        let scratch = PhaseVectorMatcher.Scratch(capacity: 512)
        var wrong = 0
        for trial in 0..<400 {
            let truth = trial % 2 == 0 ? 1 : 2
            let entries = cands.filter { lib.entries[$0].phaseIndex == truth }
            let e = lib.entries[entries[Int(rnd() * Double(entries.count)) % entries.count]]
            var v: [SIMD2<Double>] = []
            for r in e.vectors where r.length > 0.15 {
                if rnd() < dropout { continue }
                v.append(r.q + SIMD2(gauss(), gauss()) * sigma)
            }
            for _ in 0..<2 { let a = rnd() * 2 * .pi, r = 0.2 + rnd() * 0.7; v.append(SIMD2(cos(a), sin(a)) * r) }
            let res = PhaseVectorMatcher.classify(vectors: v, library: lib, settings: settings, matrixEntry: nil,
                                                  candidateEntryIndices: cands, scratch: scratch)
            if res.verdict == .indexed, Int(res.phaseIndex) != truth { wrong += 1 }
        }
        return wrong
    }

    @Test func theCompletenessRankingLosesAtSmallNoiseAndWinsNearThePairRadius() throws {
        // Measured with the seed above: default 0 / aware 23 wrong at sigma 0.003, dropout 0.4;
        // default 53 / aware 20 at sigma 0.010, dropout 0.4.
        #expect(try wrongLabels(sigma: 0.003, dropout: 0.4, aware: false) <= 2)
        #expect(try wrongLabels(sigma: 0.003, dropout: 0.4, aware: true) >= 10)
        #expect(try wrongLabels(sigma: 0.010, dropout: 0.4, aware: true) < (try wrongLabels(sigma: 0.010, dropout: 0.4, aware: false)))
    }

    /// One stray peak and none of the matrix's: the survivor count (1) calls it the matrix, the shipped rule. With
    /// `minimumMatrixEvidence` 1 the position is "not indexed", while a pattern the matrix explains stays matrix.
    /// RED UNDER: the setting absent (compile error), then its guard removed.
    @Test func matrixByExclusionCanRequireTheMatrixToHaveBeenSeen() throws {
        let lib = try library()
        let matrixEntry = lib.entries[lib.matrixEntryIndices[0]]
        let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
        let stray: [SIMD2<Double>] = [SIMD2(0.31, 0.12)]
        let explained = Array(matrixEntry.vectors.prefix(6).map(\.q))
        for rule in [PhaseVectorSettings.ClassificationRule.search, .knownVariants] {
            var off = PhaseVectorSettings(); off.classificationRule = rule
            var on = off; on.minimumMatrixEvidence = 1
            func verdict(_ v: [SIMD2<Double>], _ s: PhaseVectorSettings) -> PhaseVerdict {
                PhaseVectorMatcher.classify(vectors: v, library: lib, settings: s, matrixEntry: matrixEntry,
                                            candidateEntryIndices: lib.candidateEntryIndices, scratch: scratch).verdict
            }
            #expect(verdict(stray, off) == .matrix, "\(rule): shipped")
            #expect(verdict(stray, on) == .notIndexed, "\(rule): no matrix vector seen")
            #expect(verdict(explained, on) == .matrix, "\(rule): matrix explained everything")
        }
    }
}
