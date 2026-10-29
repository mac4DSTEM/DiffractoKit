//
//  RingIndexingTests.swift
//  The indexing decision ("RingIndexing", "Truth"), on the kit's own ring
//  lists (DiffractoScattering's `kinematicReflections(kMax: 1.2).families`). It ports every test of reference commit 9482831
//  `DSpacingIndexingTests` (28) and `MonoclinicIndexingTests` (13) that concerns the decision or the family lists,
//  re-expressed on the kit's rings; the ones that test what DiffractoCore and DiffractoScattering own are mapped to a comment at the end.
//  The dyadic lists are the kit's rings of an exactly-scaled orthorhombic cell (a 2, b 1, c 1/2 Å: g = h/2, k, 2 l, so
//  g in {1, 1.5, 2, 3, ...} is exact), so z, the k test, the margin and the reach are pinned to the last bit.
//  Each test names the source mutations that turn it red.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

struct RingIndexingTests {

    typealias Index = RingIndexing
    private static let kMax = 1.2

    // MARK: Helpers

    private func verdict(_ q: Double, _ sigma: Double, _ list: [ReflectionFamily], kMax: Double = 4,
                         k: Double = 3) -> Result<Index.Assignment, Index.Refusal> {
        do throws(Index.Refusal) { return .success(try Index.index(q: q, sigma: sigma, in: list, kMax: kMax, coverage: k)) }
        catch { return .failure(error) }
    }

    private func accepted(_ r: Result<Index.Assignment, Index.Refusal>, sourceLocation: SourceLocation = #_sourceLocation)
        -> Index.Assignment? {
        switch r {
        case .success(let a): return a
        case .failure(let e):
            Issue.record("expected an assignment, got \(e.name): \(e.reason)", sourceLocation: sourceLocation)
            return nil
        }
    }

    private func refusal(_ r: Result<Index.Assignment, Index.Refusal>, sourceLocation: SourceLocation = #_sourceLocation)
        -> Index.Refusal? {
        switch r {
        case .failure(let e): return e
        case .success:
            Issue.record("expected a refusal, got an assignment", sourceLocation: sourceLocation)
            return nil
        }
    }

    private func family(_ list: [ReflectionFamily], containing m: MillerIndex) throws -> ReflectionFamily {
        try #require(list.first { $0.members.contains(m) }, "no family holds \(m)")
    }

    private func hkl(_ h: Int, _ k: Int, _ l: Int) -> MillerIndex { MillerIndex(h, k, l) }

    // MARK: The lists, on the kit's rings (the reference implementation's list tests, re-expressed)

    /// Au has 8 rings and {111} holds 8 members; Al 8. (reference testGoldFamilies.) Mutants: absences kept (Au 8 -> 20),
    /// the multiplicity as the orbit count (8 -> 1).
    @Test func testGoldFamilies() throws {
        let a = 4.0782
        let au = try IndexPresets.rings(IndexPresets.gold)
        #expect(au.count == 8)
        let first = try #require(au.first)
        #expect(first.multiplicity == 8)
        #expect(first.orbits == [hkl(1, 1, 1)])
        #expect(first.name == "{1 1 1}")
        #expect(abs(first.g - 3.0.squareRoot() / a) <= 1e-12)
        #expect(abs(first.d - a / 3.0.squareRoot()) <= 1e-12)
        #expect(try IndexPresets.rings(IndexPresets.aluminium).count == 8)
        #expect(au.map(\.g) == au.map(\.g).sorted())
    }

    /// Mg (hcp): {10-10} 6 members, {0002} 2, {10-11} and {11-20} named in Miller-Bravais indices, {11-20} ONE orbit, 24
    /// rings, {10-15} and {21-32} stay apart (0.156 % in g). (reference testMagnesiumFamiliesAndOrbits.) Mutants: grouping at
    /// 1e-2 (24 -> 22), the multiplicity as the orbit count.
    @Test func testMagnesiumFamiliesAndOrbits() throws {
        let mg = try IndexPresets.rings(IndexPresets.magnesium)
        #expect(mg.count == 24)
        let prism = try family(mg, containing: hkl(1, 0, 0))
        #expect(prism.multiplicity == 6)
        #expect(prism.name == "{1 0 -1 0}")
        let basal = try family(mg, containing: hkl(0, 0, 2))
        #expect(basal.multiplicity == 2)
        #expect(basal.name == "{0 0 0 2}")
        let pyramid = try family(mg, containing: hkl(1, 0, 1))
        #expect(pyramid.name == "{1 0 -1 1}")
        #expect(pyramid.multiplicity == 12)
        let eleven = try family(mg, containing: hkl(1, 1, 0))
        #expect(eleven.name == "{1 1 -2 0}")
        #expect(eleven.orbits.count == 1, "(1 1 0) and (-1 -1 0) are one orbit")
        #expect(eleven.multiplicity == 6)
        let a = try family(mg, containing: hkl(1, 0, 5))
        let b = try family(mg, containing: hkl(2, 1, 2))
        #expect(a.g != b.g, "Mg {10-15} and {21-32} are two rings")
        #expect(abs(abs(b.g / a.g - 1) - 1.5e-3) <= 1e-4)
    }

    /// Si (diamond): {333} and {511} are one ring of 32 members and two orbits; 10 rings. (reference
    /// testSiliconCoincidentOrbitsAreOneFamily.) Mutants: multiplicity as the orbit count, absences kept (10 -> 36).
    @Test func testSiliconCoincidentOrbitsAreOneFamily() throws {
        let si = try IndexPresets.rings(IndexPresets.silicon)
        #expect(si.count == 10)
        let ring = try family(si, containing: hkl(3, 3, 3))
        #expect(ring.multiplicity == 32)
        #expect(ring.orbits == [hkl(5, 1, 1), hkl(3, 3, 3)])
        #expect(ring.name == "{5 1 1} + {3 3 3}")
        #expect(ring.members.contains(hkl(5, 1, 1)))
    }

    /// Ideal c/a hcp: {30-32} and {0006} coincide: one ring of 14 members (12 + 2), 23 rings. (reference
    /// testIdealHcpCoincidentOrbitsAreOneFamily.) Mutant: grouping at 0 splits them at the roundoff of |g|.
    @Test func testIdealHcpCoincidentOrbitsAreOneFamily() throws {
        let list = try IndexPresets.rings(IndexPresets.idealHCP)
        #expect(list.count == 23)
        let ring = try family(list, containing: hkl(3, 0, 2))
        #expect(ring.multiplicity == 14)
        #expect(ring.orbits.count == 2)
        #expect(ring.members.contains(hkl(0, 0, 6)))
    }

    /// The |F| cut keeps weak rings: the Weak-cut cell (a second Cu site at occupancy 0.99815) has 10. (reference
    /// testWeakFamiliesCount.) Mutant: the cut raised to 1e-2 (10 -> 6).
    @Test func testWeakFamiliesCount() throws {
        #expect(try IndexPresets.rings(IndexPresets.weakCut).count == 10)
    }

    /// Two reflections of one |g| to 1e-5, not 1e-9 (ideal hcp with c stretched by 1e-5): {30-32} and {0006} are two
    /// rings (24, not 23). (reference testNearCoincidentHcpOrbitsStaySeparate.) Mutant: the family tolerance 1e-4.
    @Test func testNearCoincidentHcpOrbitsStaySeparate() throws {
        let a = 3.2094
        let list = try IndexPresets.rings(IndexPresets.hcp(12, a: a, c: a * (8.0 / 3.0).squareRoot() * (1 + 1e-5)))
        #expect(list.count == 24)
        let x = try family(list, containing: hkl(3, 0, 2))
        let y = try family(list, containing: hkl(0, 0, 6))
        #expect(x.g != y.g)
        #expect(x.multiplicity == 12)
        #expect(y.multiplicity == 2)
    }

    // MARK: The deviation, z, margin and reach (exact)

    /// delta of q = 1.01 g is g/q - 1 = -0.00990099..., not 1 - q/g = -0.01. (reference testSignedDeviationIsOverTheMeasuredD.)
    /// Mutant: delta over the measured d reversed.
    @Test func testSignedDeviationIsOverTheMeasuredD() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        let g = au[0].g
        let q = 1.01 * g
        let a = try #require(accepted(verdict(q, 0.005 * q, au, kMax: Self.kMax)))
        #expect(abs(a.best.deviation - (1 / 1.01 - 1)) <= 1e-14)
        #expect(abs(a.best.deviation - -0.0099009900990099) <= 1e-14)
        #expect(abs(a.best.relativeDeviation - 0.0099009900990099) <= 1e-14)
        #expect(abs(a.dMeasured - 1 / q) <= 1e-15)
        #expect(abs(a.best.z - 0.01 * g / (0.005 * q)) <= 1e-12)
    }

    /// Exact z and margin: g in {1, 1.5, 2}, sigma 0.125, q = 1: z = 0, the next (1.5) at z = 4, margin (16 - 0)/2 = 8.
    /// (reference testDyadicZAndMargin.) Mutants: z as |q - g|/(sigma q), the margin as z2 - z1 (4), the next as the upper
    /// neighbour.
    @Test func testDyadicZAndMargin() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let a = try #require(accepted(verdict(1, 0.125, list)))
        #expect(a.best.z == 0)
        let next = try #require(a.next)
        #expect(next.z == 4)
        #expect(next.family.g == 1.5)
        #expect(next.deviation == 0.5)
        #expect(a.margin == 8)
        // z = |q - g|/sigma at a point that is not a family: q = 1.0625 -> z = 0.5
        let b = try #require(accepted(verdict(1.0625, 0.125, list)))
        #expect(b.best.z == 0.5)
        #expect(abs(b.best.deviation - (1 / 1.0625 - 1)) <= 1e-15)
    }

    /// z <= k, not z^2 <= k: q = 1.375, sigma = 0.125 on {1, 2, 3} is z = 3 and accepted at the default k = 3 (the next, 2,
    /// at z = 5: margin 8); at k = 2.9 it is `noneWithinTolerance`. Also pins the default k = 3 and z = 2.5 accepted.
    /// (reference testWithinIsZAtMostKAndDefaultKIsThree.) Mutants: z <= k to z < k (the z = 3 plant is refused), z^2 <= k,
    /// z over sigma q.
    @Test func testWithinIsZAtMostKAndDefaultKIsThree() throws {
        #expect(Index.defaultCoverage == 3)
        let list = try IndexPresets.dyadic([1, 2, 3])
        let at3 = try #require(accepted(verdict(1.375, 0.125, list)))
        #expect(at3.best.z == 3)
        #expect(at3.coverage == 3)
        #expect(try #require(at3.next).z == 5)
        let viaDefault = try Index.index(q: 1.375, sigma: 0.125, in: list, kMax: 4)
        #expect(viaDefault.coverage == 3)
        let none = try #require(refusal(verdict(1.375, 0.125, list, k: 2.9)))
        #expect(none.name == "noneWithinTolerance")
        let at25 = try #require(accepted(verdict(1.3125, 0.125, list)))
        #expect(at25.best.z == 2.5)
        let wantMargin: Double = (5.5 * 5.5 - 2.5 * 2.5) / 2
        #expect(at25.margin == wantMargin)
    }

    /// Two rings both at z = 2 <= k = 2 (g {1, 1.5}, q = 1.25, sigma 0.125) are `ambiguous`, payload in (z, g) order.
    /// (reference testBothWithinKIsAmbiguousInOrder.) Mutants: z <= k to z < k (none within), z^2 <= k.
    @Test func testBothWithinKIsAmbiguousInOrder() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let r = try #require(refusal(verdict(1.25, 0.125, list, k: 2)))
        guard case .ambiguous(let candidates) = r else { Issue.record("expected ambiguous, got \(r.name)"); return }
        #expect(candidates.map(\.family.g) == [1, 1.5])
        #expect(candidates.map(\.z) == [2, 2])
        // equal z in (z, g) order whichever way the list is given
        let rev = try #require(refusal(verdict(1.25, 0.125, Array(list.reversed()), k: 2)))
        guard case .ambiguous(let c2) = rev else { Issue.record("expected ambiguous, got \(rev.name)"); return }
        #expect(c2.map(\.family.g) == [1, 1.5])
    }

    /// Two rings within k with a margin far above ln 19 are still `ambiguous` (z = 0 and z = 4 at k = 4: margin 8).
    /// (reference testTwoWithinKIsAmbiguousEvenWithALargeMargin.) Mutant: the second-within check dropped.
    @Test func testTwoWithinKIsAmbiguousEvenWithALargeMargin() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let r = try #require(refusal(verdict(1, 0.125, list, k: 4)))
        guard case .ambiguous(let candidates) = r else { Issue.record("expected ambiguous, got \(r.name)"); return }
        #expect(candidates.map(\.family.g) == [1, 1.5])
    }

    /// One ring within k whose margin to the next is below ln 19 (q = 1, sigma 0.25, k = 1: z = 0 and 2, margin 2 < 2.944)
    /// is `ambiguous` with [best, next]. (reference testMarginBelowLn19IsAmbiguous.) Mutant: the margin clause dropped.
    @Test func testMarginBelowLn19IsAmbiguous() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let r = try #require(refusal(verdict(1, 0.25, list, k: 1)))
        guard case .ambiguous(let candidates) = r else { Issue.record("expected ambiguous, got \(r.name)"); return }
        #expect(candidates.map(\.family.g) == [1, 1.5])
        #expect(Index.minimumMargin == log(19.0))
    }

    /// The margin plant stated before the code: Au {331} -> {420}, sigma = gap/5.7, z = 2.5 (next 3.2, margin 1.995 < ln 19):
    /// `ambiguous`. (reference testAuMarginPlantIsAmbiguous.) Mutants: the margin clause dropped; the margin as z2 - z1.
    @Test func testAuMarginPlantIsAmbiguous() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        let f331 = try family(au, containing: hkl(3, 3, 1))
        let f420 = try family(au, containing: hkl(4, 2, 0))
        let sigma = (f420.g - f331.g) / 5.7
        let q = f331.g + 2.5 * sigma
        let r = try #require(refusal(verdict(q, sigma, au, kMax: Self.kMax)))
        guard case .ambiguous(let candidates) = r else { Issue.record("expected ambiguous, got \(r.name)"); return }
        #expect(candidates.count == 2)
        #expect(abs(candidates[0].z - 2.5) <= 1e-9)
        #expect(abs(candidates[1].z - 3.2) <= 1e-9)
    }

    /// The reach on {1, 1.5, 2}, q = 1, sigma 0.125: the next (1.5) is 0.5 away, so q + R = 1.5. At k_max 1.5 it is accepted
    /// (1.5 is not above k_max); at k_max 1.4375 it is `beyondKMax`, though q + k sigma = 1.375 is below the edge.
    /// (reference testReachIncludesTheDistanceToTheNextFamily.) Mutants: the reach at k sigma alone; reach without the next
    /// ring's distance.
    @Test func testReachIncludesTheDistanceToTheNextFamily() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let inside = try #require(accepted(verdict(1, 0.125, list, kMax: 1.5)))
        #expect(inside.best.family.g == 1)
        let r = try #require(refusal(verdict(1, 0.125, list, kMax: 1.4375)))
        guard case .beyondKMax(let reach, let kMax) = r else { Issue.record("expected beyondKMax, got \(r.name)"); return }
        #expect(reach == 1.5)
        #expect(kMax == 1.4375)
        #expect(abs(Index.reach(coverage: 3) - (9 + 2 * log(19.0)).squareRoot()) <= 1e-15)
    }

    /// q + R = k_max exactly is accepted (the clause is `q + R > k_max`): the same measurement with k_max one ulp below the
    /// reach is `beyondKMax`, one ulp above is accepted. (The "q + R = kMax exactly" boundary.) Mutant:
    /// `>=` for `>`.
    @Test func testReachEqualToKMaxIsAccepted() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        #expect(accepted(verdict(1, 0.125, list, kMax: 1.5)) != nil)
        #expect(accepted(verdict(1, 0.125, list, kMax: (1.5).nextUp)) != nil)
        let r = try #require(refusal(verdict(1, 0.125, list, kMax: (1.5).nextDown)))
        #expect(r.name == "beyondKMax")
        // the none path: nothing within, q + max(r sigma, nearest distance) = 1.5 + exactly the edge
        let none = Index.reach(coverage: 3) * 0.125
        #expect(none < 1)
        let n = try #require(refusal(verdict(1.25, 0.0625, list, kMax: 4)))
        #expect(n.name == "noneWithinTolerance" || n.name == "ambiguous")
    }

    // MARK: Which ring, which next

    /// Au {200} measured just above its g: the best is {200} and the next is {111}, the nearer of its two neighbours
    /// ({220} is farther). (reference testAuNextIsTheNearerNeighbour.) Mutants: the next as the upper neighbour.
    @Test func testAuNextIsTheNearerNeighbour() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        let f200 = try family(au, containing: hkl(2, 0, 0))
        let q = f200.g * 1.0001
        let a = try #require(accepted(verdict(q, 1e-3 * q, au, kMax: Self.kMax)))
        #expect(a.best.family.g == f200.g)
        #expect(try #require(a.next).family.g == (try family(au, containing: hkl(1, 1, 1))).g)
        #expect(try #require(a.next).z > 100)
    }

    /// The nearest ring of a `noneWithinTolerance`, with its z and signed deviation; a tie in z goes to the smaller g.
    /// (reference testNoneCarriesTheNearestFamily.) Mutant: the nearest as the first g >= q (g = 2 here).
    @Test func testNoneCarriesTheNearestFamily() throws {
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let r = try #require(refusal(verdict(1.75, 0.0625, list)))
        guard case .noneWithinTolerance(let nearest?) = r else { Issue.record("expected none with a nearest, got \(r.name)"); return }
        #expect(nearest.family.g == 1.5)
        #expect(nearest.z == 4)
        #expect(abs(nearest.deviation - (1.5 / 1.75 - 1)) <= 1e-15)
        // the order of the input list does not matter: a tie in z goes to the smaller g whichever comes first
        let rev = try #require(refusal(verdict(1.75, 0.0625, Array(list.reversed()))))
        guard case .noneWithinTolerance(let n2?) = rev else { Issue.record("expected none, got \(rev.name)"); return }
        #expect(n2.family.g == 1.5)
    }

    // MARK: The edges of the list

    /// A one-ring list (Au, k_max 0.45): accepted with no next and no margin; an empty list (k_max 0.3) is
    /// `noneWithinTolerance` with no nearest, and `beyondKMax` close to the edge. (reference testOneFamilyAndEmptyLists.)
    /// Mutants: the next as the upper neighbour (nil -> wrong), the refusals reordered.
    @Test func testOneFamilyAndEmptyLists() throws {
        let one = try IndexPresets.rings(IndexPresets.gold, kMax: 0.45)
        #expect(one.count == 1)
        let g = one[0].g
        let a = try #require(accepted(verdict(g, 1e-3 * g, one, kMax: 0.45)))
        #expect(a.next == nil)
        #expect(a.margin == nil)
        let empty = try IndexPresets.rings(IndexPresets.gold, kMax: 0.3)
        #expect(empty.isEmpty)
        let r = try #require(refusal(verdict(0.2, 2e-4, empty, kMax: 0.3)))
        guard case .noneWithinTolerance(let nearest) = r else { Issue.record("expected none, got \(r.name)"); return }
        #expect(nearest == nil)
        let beyond = try #require(refusal(verdict(0.299, 3e-4, empty, kMax: 0.3)))
        #expect(beyond.name == "beyondKMax")
    }

    /// Mg {20-25} (the list's last ring), Si {620} (its listed next, {531}, is nearer than the edge but the true next,
    /// {533}, is past it) at their own g, and q = 1.25 on Au, are `beyondKMax`. (reference testBeyondKMax.) Mutants: the reach
    /// at k sigma alone (Si {620} is accepted).
    @Test func testBeyondKMax() throws {
        let mg = try IndexPresets.rings(IndexPresets.magnesium)
        let last = try #require(mg.last)
        #expect(last.name == "{2 0 -2 5}")
        #expect(try #require(refusal(verdict(last.g, 1e-3 * last.g, mg, kMax: Self.kMax))).name == "beyondKMax")
        let si = try IndexPresets.rings(IndexPresets.silicon)
        let f620 = try family(si, containing: hkl(6, 2, 0))
        #expect(try #require(refusal(verdict(f620.g, 1e-3 * f620.g, si, kMax: Self.kMax))).name == "beyondKMax")
        let au = try IndexPresets.rings(IndexPresets.gold)
        #expect(try #require(refusal(verdict(1.25, 1.25e-3, au, kMax: Self.kMax))).name == "beyondKMax")
    }

    /// Kinematically forbidden positions are `noneWithinTolerance`: Si {200}, Si {222}, Mg {0001}, Au {100}, and Au {111}
    /// against Al's list. (reference testForbiddenPositionsAreNone.) Mutant: absences kept.
    @Test func testForbiddenPositionsAreNone() throws {
        let si = try IndexPresets.rings(IndexPresets.silicon)
        for g in [2 / 5.4309, 12.0.squareRoot() / 5.4309] {
            #expect(try #require(refusal(verdict(g, 1e-3 * g, si, kMax: Self.kMax))).name == "noneWithinTolerance")
        }
        let mg = try IndexPresets.rings(IndexPresets.magnesium)
        let q0001 = 1 / 5.2108
        #expect(try #require(refusal(verdict(q0001, 1e-3 * q0001, mg, kMax: Self.kMax))).name == "noneWithinTolerance")
        let au = try IndexPresets.rings(IndexPresets.gold)
        let q100 = 1 / 4.0782
        #expect(try #require(refusal(verdict(q100, 1e-3 * q100, au, kMax: Self.kMax))).name == "noneWithinTolerance")
        let al = try IndexPresets.rings(IndexPresets.aluminium)
        let g111 = 3.0.squareRoot() / 4.0782
        #expect(try #require(refusal(verdict(g111, 1e-3 * g111, al, kMax: Self.kMax))).name == "noneWithinTolerance")
    }

    /// Mg {10-15} and {21-32}, 0.156 % apart, with the measurement 40 % of the way: `ambiguous` at u = 1e-3. (reference
    /// testMgNearPairIsAmbiguous.) Mutants: grouping at 1e-2, the second-within check.
    @Test func testMgNearPairIsAmbiguous() throws {
        let mg = try IndexPresets.rings(IndexPresets.magnesium)
        let a = try family(mg, containing: hkl(1, 0, 5))
        let b = try family(mg, containing: hkl(2, 1, 2))
        let q = a.g + 0.4 * (b.g - a.g)
        #expect(try #require(refusal(verdict(q, 1e-3 * q, mg, kMax: Self.kMax))).name == "ambiguous")
    }

    /// A margin of exactly ln 19 is accepted (the clause is `margin < ln 19`, not `<=`). ln 19 is irrational, so the test
    /// scans sigma by ulps, at four q, for measurements on {1, 2} whose computed margin IS `minimumMargin`, and counts the
    /// ones that are accepted; it also counts the ones whose margin is one ulp short, which must be ambiguous. It fails if
    /// no exact measurement is found. (reference testMarginExactlyLn19IsAccepted.) Mutant: margin <= ln 19 for margin < ln 19
    /// (the exact ones are refused).
    @Test func testMarginExactlyLn19IsAccepted() throws {
        let list = try IndexPresets.dyadic([1, 2])
        let target = Index.minimumMargin
        var exact = 0, acceptedExact = 0, shortAmbiguous = 0
        for offset in [0.0, 0.125, 0.25, 0.375] {
            let q = 1 + offset
            // z_best = offset/sigma, z_next = (1 - offset)/sigma, margin = ((1 - offset)^2 - offset^2)/(2 sigma^2)
            var sigma = ((1 - 2 * offset) / (2 * target)).squareRoot()
            for _ in 0..<3000 { sigma = sigma.nextDown }
            for _ in 0..<6000 {
                sigma = sigma.nextUp
                let k = (offset / sigma + (1 - offset) / sigma) / 2         // z_best <= k < z_next
                let zb = abs(q - 1) / sigma, zn = abs(q - 2) / sigma
                let margin = (zn * zn - zb * zb) / 2
                if margin == target {
                    exact += 1
                    if case .success = verdict(q, sigma, list, k: k) { acceptedExact += 1 }
                } else if margin == target.nextDown {
                    if case .failure(let e) = verdict(q, sigma, list, k: k), e.name == "ambiguous" { shortAmbiguous += 1 }
                }
            }
        }
        #expect(exact > 0, "no measurement with a margin of exactly ln 19 was found")
        #expect(acceptedExact == exact, "a margin of exactly ln 19 must be accepted")
        #expect(shortAmbiguous > 0)
    }

    /// The label is carried: it begins "none", and an assignment reports it. (reference testValidationLabelIsNone.)
    @Test func testValidationLabelIsNone() throws {
        #expect(Index.validation.hasPrefix("none:"))
        let a = try #require(accepted(verdict(1, 0.125, try IndexPresets.dyadic([1, 1.5, 2]))))
        #expect(a.validation == Index.validation)
    }

    // MARK: The refusals, by name and in order

    /// `invalidInput` from `index`: q, sigma, k and k_max each at 0, -1, NaN, +/-infinity. (reference
    /// testInvalidInputIsRefusedByName; the list builder's own k_max refusal is DiffractoCore's `ReflectionListError.invalidKMax`,
    /// tested in ScatteringTruthTests.testRefusals.) Mutants: sigma = 0 accepted, k_max unchecked.
    @Test func testInvalidInputIsRefusedByName() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        let bad: [Double] = [0, -1, .nan, .infinity, -.infinity]
        for x in bad {
            #expect(try #require(refusal(verdict(x, 1e-3, au, kMax: Self.kMax))).name == "invalidInput", "q = \(x)")
            #expect(try #require(refusal(verdict(0.42, x, au, kMax: Self.kMax))).name == "invalidInput", "sigma = \(x)")
            #expect(try #require(refusal(verdict(0.42, 1e-3, au, kMax: Self.kMax, k: x))).name == "invalidInput", "k = \(x)")
            #expect(try #require(refusal(verdict(0.42, 1e-3, au, kMax: x))).name == "invalidInput", "k_max = \(x)")
        }
        #expect(accepted(verdict(0.42452, 1e-3, au, kMax: 4)) != nil || refusal(verdict(0.42452, 1e-3, au, kMax: 4)) != nil)
    }

    /// The first that applies: invalid input beats every other check, even where the answer would be `beyondKMax` or
    /// `ambiguous`; `ambiguous` beats `beyondKMax`; `beyondKMax` beats `noneWithinTolerance`. (reference
    /// testFirstThatAppliesPairs, re-expressed for the verdict order.) Mutants: any two checks swapped.
    @Test func testVerdictOrder() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        #expect(try #require(refusal(verdict(.nan, 0, au, kMax: Self.kMax))).name == "invalidInput")
        #expect(try #require(refusal(verdict(1.25, -1, au, kMax: Self.kMax))).name == "invalidInput")
        // two rings within k, and q + R past k_max: ambiguous wins
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        #expect(try #require(refusal(verdict(1, 0.125, list, kMax: 1.0, k: 4))).name == "ambiguous")
        // nothing within k, and q + R past k_max: beyondKMax wins over noneWithinTolerance
        #expect(try #require(refusal(verdict(1.75, 0.0625, list, kMax: 1.8))).name == "beyondKMax")
        #expect(try #require(refusal(verdict(1.75, 0.0625, list, kMax: 4))).name == "noneWithinTolerance")
    }

    /// The names are a stable key and the reason is a sentence naming no application. (reference testRefusalNamesAndReasons.)
    @Test func testRefusalNamesAndReasons() throws {
        let au = try IndexPresets.rings(IndexPresets.gold)
        let none = try #require(refusal(verdict(0.2452, 2.4e-4, au, kMax: Self.kMax)))
        #expect(none.name == "noneWithinTolerance")
        #expect(none.reason.contains("possible causes, unranked"))
        #expect(none.errorDescription == none.reason)
        #expect(none.description == none.reason)
        let list = try IndexPresets.dyadic([1, 1.5, 2])
        let all: [Index.Refusal] = [
            none,
            try #require(refusal(verdict(1.25, 0.125, list, k: 2))),
            try #require(refusal(verdict(1, 0.125, list, kMax: 1.4375))),
            try #require(refusal(verdict(0.42, 0, au, kMax: Self.kMax))),
            .noneWithinTolerance(nearest: nil),
        ]
        #expect(all.map(\.name) == ["noneWithinTolerance", "ambiguous", "beyondKMax", "invalidInput", "noneWithinTolerance"])
        for r in all {
            for app in ["mac4DSTEM", "macTEM", "ACOM"] { #expect(!r.reason.contains(app), "\(r.reason)") }
        }
        #expect(!Index.validation.contains("macTEM") && !Index.validation.contains("mac4DSTEM"))
    }

    // MARK: The overloads, and beta''

    /// The two overloads (DiffractoScattering's `KinematicFamily` and DiffractoCore's `ReflectionFamily`) give one answer, in the same
    /// payloads, for an accepted, an ambiguous and a refused plant.
    @Test func testKinematicAndReflectionOverloadsAgree() throws {
        let kin = try IndexPresets.mSq.kinematicReflections(kMax: 1.2).families
        let ref = kin.map(\.family)
        let g = ref[10].g
        for (q, u) in [(g, 1e-4), (g, 1e-2), (g * 1.5, 1e-3), (1.19, 1e-3), (0.1, 1e-3)] {
            do throws(Index.Refusal) {
                let a = try Index.index(q: q, sigma: u * q, in: kin, kMax: 1.2)
                let b = try Index.index(q: q, sigma: u * q, in: ref, kMax: 1.2)
                #expect(a == b)
            } catch {
                #expect(error == (verdict(q, u * q, ref, kMax: 1.2).failure))
            }
        }
    }

    /// On beta'': q on {2 0 0} at u = 1e-3 is accepted, next {0 0 1}; q midway between {3 1 5} and {7 3 0} (1.03e-5 apart)
    /// at u = 1e-3 is ambiguous, naming both. (reference testIndexOnTheBetaList.) Mutant: grouping at 1e-4 (the pair merges
    /// into one ring: accepted).
    @Test func testIndexOnTheBetaList() throws {
        let l = try IndexPresets.rings(IndexPresets.betaDoublePrime)
        let f200 = try family(l, containing: hkl(2, 0, 0))
        let a = try Index.index(q: f200.g, sigma: 1e-3 * f200.g, in: l, kMax: 1.2)
        #expect(a.best.family.name == "{2 0 0}")
        #expect(a.next?.family.name == "{0 0 1}")
        #expect(abs(a.best.z) <= 1e-12)
        let p = try family(l, containing: hkl(3, 1, 5)), r = try family(l, containing: hkl(7, 3, 0))
        #expect(abs(abs(p.g - r.g) / min(p.g, r.g) - 1.03e-5) <= 5e-8)
        let q = (p.g + r.g) / 2
        let e = try #require(refusal(verdict(q, 1e-3 * q, l, kMax: 1.2)))
        guard case .ambiguous(let c) = e else { Issue.record("a pair 1e-5 apart at sigma/q = 1e-3 must be ambiguous, got \(e.name)"); return }
        #expect(Set(c.map(\.family.name)).isSuperset(of: ["{3 1 5}", "{7 3 0}"]), "every ring within k sigma is named, these two among them")
    }

    /// reference ADR 024, stated before the code: at s = 0 (q = g of each ring) and sigma_q = u q, the number of beta'' rings that are
    /// accepted is 8 at u = 1e-2, 23 at 3e-3 and 4 at 3e-2; at 1e-2 they are {2 0 0}, {0 0 1}, {2 0 -1}, {2 0 1}, {1 1 0},
    /// {1 1 -1}, {3 1 0}, {3 1 -1}. Mutants: z <= k to z < k, the margin clause moved, reach without the next ring.
    @Test func testBetaDoublePrimeAcceptedCountsAtZeroOffset() throws {
        let l = try IndexPresets.rings(IndexPresets.betaDoublePrime)
        #expect(l.count == 417)
        func accepted(at u: Double) -> [String] {
            l.compactMap { f in
                if case .success(let a) = verdict(f.g, u * f.g, l, kMax: 1.2), a.best.family.name == f.name { return f.name }
                return nil
            }
        }
        #expect(accepted(at: 1e-2).count == 8)
        #expect(accepted(at: 3e-3).count == 23)
        #expect(accepted(at: 3e-2).count == 4)
        #expect(Set(accepted(at: 1e-2)) == ["{2 0 0}", "{0 0 1}", "{2 0 -1}", "{2 0 1}", "{1 1 0}", "{1 1 -1}", "{3 1 0}", "{3 1 -1}"])
        // the d of the eight is at least 3.01 A
        for name in accepted(at: 1e-2) { #expect(try #require(l.first { $0.name == name }).d >= 3.01) }
    }

    // MARK: Monoclinic cells, on the kit's rings

    /// Orbit sizes on M-P: {1 1 0} 4 members, {2 0 0} 2 (k = 0), {0 1 0} 2, {0 0 1} 2, {1 1 1} 4; each ring is one orbit
    /// and its members are the full orbit; 392 rings and 1384 reflections at k_max 1.2. (reference testOrbitSizesAndMultiplicity.)
    @Test func testOrbitSizesAndMultiplicityOnMP() throws {
        let l = try IndexPresets.rings(IndexPresets.mP)
        for (m, n) in [(hkl(1, 1, 0), 4), (hkl(2, 0, 0), 2), (hkl(0, 1, 0), 2), (hkl(0, 0, 1), 2), (hkl(1, 1, 1), 4)] {
            let f = try family(l, containing: m)
            #expect(f.multiplicity == n, "\(m)")
            #expect(f.orbits.count == 1)
            #expect(Set(f.members).count == n)
        }
        let f110 = try family(l, containing: hkl(1, 1, 0))
        #expect(Set(f110.members) == [hkl(1, 1, 0), hkl(-1, 1, 0), hkl(-1, -1, 0), hkl(1, -1, 0)])
        #expect(l.count == 392)
        let multiplicities1: [Int] = l.map(\.multiplicity)
        let multiplicitySum1: Int = multiplicities1.reduce(0, +)
        #expect(multiplicitySum1 == 1384)
    }

    /// The ring g is the list's own: the closed forms d(h00) = a sin(beta)/h, d(0k0) = b/k, d(00l) = c sin(beta)/l and
    /// 1/d^2 = (h^2/a^2 + l^2/c^2 - 2hl cos(beta)/(ac))/sin^2(beta) + k^2/b^2 for (101) and (10-1), at 1e-12 relative.
    /// (reference testClosedFormDSpacingsOnTheMonoclinicCell.)
    @Test func testClosedFormDSpacingsOnTheMonoclinicCell() throws {
        let (a, b, c) = (7.2, 4.9, 5.6), beta = 104.0 * Double.pi / 180
        let s = sin(beta), co = cos(beta)
        let l = try IndexPresets.rings(IndexPresets.mP)
        func d(_ h: Int, _ k: Int, _ ll: Int) throws -> Double { try family(l, containing: hkl(h, k, ll)).d }
        func closed(_ h: Double, _ k: Double, _ ll: Double) -> Double {
            1 / (((h * h / (a * a) + ll * ll / (c * c) - 2 * h * ll * co / (a * c)) / (s * s) + k * k / (b * b)).squareRoot())
        }
        #expect(abs(try d(1, 0, 0) - a * s) <= 1e-12 * a)
        #expect(abs(try d(0, 1, 0) - b) <= 1e-12 * b)
        #expect(abs(try d(0, 0, 1) - c * s) <= 1e-12 * c)
        #expect(abs(try d(1, 0, 1) - closed(1, 0, 1)) <= 1e-12 * 7)
        #expect(abs(try d(1, 0, -1) - closed(1, 0, -1)) <= 1e-12 * 7)
        #expect(abs(try d(1, 0, 1) - d(1, 0, -1)) > 1e-3)
        #expect(try family(l, containing: hkl(1, 0, 1)).members.contains(hkl(-1, 0, -1)))
        #expect(!(try family(l, containing: hkl(1, 0, 1)).members.contains(hkl(1, 0, -1))))
    }

    /// beta'': (2 0 1) and (2 0 -1) are two rings, 4.3246 and 5.6555 Å. (reference testPlusAndMinusLAreDifferentFamilies.)
    @Test func testPlusAndMinusLAreDifferentFamilies() throws {
        let l = try IndexPresets.rings(IndexPresets.betaDoublePrime)
        let p = try family(l, containing: hkl(2, 0, 1)), m = try family(l, containing: hkl(2, 0, -1))
        #expect(p.g != m.g)
        #expect(abs(p.d - 4.324567590028738) <= 1e-9)
        #expect(abs(m.d - 5.655535439744675) <= 1e-9)
        #expect(p.name == "{2 0 1}")
        #expect(m.name == "{2 0 -1}")
    }

    /// An acute beta: M-P with beta 76 deg is M-P with beta 104 deg and the sign of l mirrored, so the same rings carry the
    /// swapped names. (reference testAcuteBetaSwapsTheNamesOfPlusAndMinusL.)
    @Test func testAcuteBetaSwapsTheNamesOfPlusAndMinusL() throws {
        let a = try IndexPresets.rings(IndexPresets.mono(7.2, 4.9, 5.6, 76)), o = try IndexPresets.rings(IndexPresets.mP)
        #expect(a.count == o.count)
        let a101 = try family(a, containing: hkl(1, 0, 1)), o10m1 = try family(o, containing: hkl(1, 0, -1))
        let a10m1 = try family(a, containing: hkl(1, 0, -1)), o101 = try family(o, containing: hkl(1, 0, 1))
        #expect(abs(a101.g - o10m1.g) <= 1e-12 * o10m1.g)
        #expect(abs(a10m1.g - o101.g) <= 1e-12 * o101.g)
        #expect(abs(a101.g - a10m1.g) > 1e-3)
        #expect(a101.name == "{1 0 1}")
        #expect(a10m1.name == "{1 0 -1}")
        #expect(a101.members.contains(hkl(-1, 0, -1)))
        #expect(!a101.members.contains(hkl(1, 0, -1)))
        #expect(a101.d > a10m1.d, "at 76 deg (1 0 1) is the larger spacing; at 104 deg it is (1 0 -1)")
        #expect(o10m1.d > o101.d)
    }

    /// M-C: no member with h + k odd; beta'' (C-centred) the same, and its first ten rings are the predicted ones (d to 4
    /// decimals, multiplicity, name). (reference testCentredCellsAndTheFirstBetaFamilies.) Mutant: absences kept.
    @Test func testCentredCellsAndTheFirstBetaFamilies() throws {
        for s in [IndexPresets.mC, IndexPresets.betaDoublePrime] {
            for f in try IndexPresets.rings(s) {
                for m in f.members { #expect((m.h + m.k) % 2 == 0, "\(m)") }
            }
        }
        let first: [(String, Double, Int)] = [("{2 0 0}", 7.3113, 2), ("{0 0 1}", 6.5011, 2), ("{2 0 -1}", 5.6555, 2),
            ("{2 0 1}", 4.3246, 2), ("{1 1 0}", 3.9031, 4), ("{4 0 0}", 3.6557, 2), ("{4 0 -1}", 3.6207, 2),
            ("{1 1 -1}", 3.4554, 4), ("{2 0 -2}", 3.3124, 2), ("{0 0 2}", 3.2506, 2)]
        let fam = try IndexPresets.rings(IndexPresets.betaDoublePrime)
        for (i, e) in first.enumerated() {
            #expect(fam[i].name == e.0)
            #expect(abs(fam[i].d - e.1) <= 5e-5)
            #expect(fam[i].multiplicity == e.2)
        }
    }

    /// a = c, cos beta = -1/4: (h,k,l) and (l,k,h) have one |g|, so {1 0 0} + {0 0 1} is one ring of 4 and {1 1 0} +
    /// {0 1 1} one of 8; M-sq has 83 rings, 62 of them holding several orbits, 698 reflections. (reference
    /// testAccidentalCoincidenceIsOneFamilyNamedByItsOrbits; the kit reproduces its names.) Mutants: grouping at 1e-4 or
    /// 0, the name listing the first orbit only.
    @Test func testAccidentalCoincidenceIsOneFamilyNamedByItsOrbits() throws {
        let l = try IndexPresets.rings(IndexPresets.mSq)
        let a = try family(l, containing: hkl(1, 0, 0))
        #expect(a.name == "{1 0 0} + {0 0 1}")
        #expect(a.multiplicity == 4)
        // With its full symmetry the same cell's {1 0 0} and {0 0 1} are one orbit, not a coincidence.
        let full = try IndexPresets.rings(IndexPresets.mSqFull)
        let merged = try family(full, containing: hkl(1, 0, 0))
        #expect(merged.multiplicity == 4)
        #expect(merged.orbits.count == 1, "\(merged.name)")
        #expect(abs(a.d - 4.8412) <= 5e-5)
        let b = try family(l, containing: hkl(1, 1, 0))
        #expect(b.name == "{1 1 0} + {0 1 1}")
        #expect(b.multiplicity == 8)
        #expect(l.count == 83)
        let matching1: Int = l.filter { (x: ReflectionFamily) -> Bool in x.orbits.count > 1 }.count
        #expect(matching1 == 62)
        let multiplicities2: [Int] = l.map(\.multiplicity)
        let multiplicitySum2: Int = multiplicities2.reduce(0, +)
        #expect(multiplicitySum2 == 698)
    }

    /// beta = 90 deg + 1e-4 splits the (h,k,l)/(h,k,-l) pairs (402 rings, none with two orbits); + 1e-8 (past the 1e-9 deg
    /// tolerance of the metric family) merges them within the 1e-9 ring rule (239 rings, 163 with two orbits). (reference
    /// testNearOrthorhombicSplitsAndMerges.) Mutants: grouping at 1e-4 (the split pairs merge), at 0 (the merged pairs split).
    @Test func testNearOrthorhombicSplitsAndMerges() throws {
        let split = try IndexPresets.rings(IndexPresets.mono(7.2, 4.9, 5.6, 90 + 1e-4))
        #expect(split.count == 402)
        let matching2: Int = split.filter { (x: ReflectionFamily) -> Bool in x.orbits.count > 1 }.count
        #expect(matching2 == 0)
        let merged = try IndexPresets.rings(IndexPresets.mono(7.2, 4.9, 5.6, 90 + 1e-8))
        #expect(merged.count == 239)
        let matching3: Int = merged.filter { (x: ReflectionFamily) -> Bool in x.orbits.count > 1 }.count
        #expect(matching3 == 163)
    }

    /// beta'' through CIFImport: 22 sites, space group 12, 417 rings, 1442 reflections at k_max 1.2, 113 rings of
    /// multiplicity 2 and 304 of 4, none with several orbits. (reference testBetaDoublePrimeFromCIFText; the CIF reading
    /// itself is DiffractoCore's CIFImportTests.)
    @Test func testBetaDoublePrimeFromCIFText() throws {
        let s = IndexPresets.betaDoublePrime
        #expect(s.sites.count == 22)
        let matching4: Int = s.sites.filter { (x: AtomSite) -> Bool in x.z == 12 }.count
        #expect(matching4 == 10)
        let matching5: Int = s.sites.filter { (x: AtomSite) -> Bool in x.z == 14 }.count
        #expect(matching5 == 12)
        #expect(s.spaceGroupNumber == 12)
        #expect(s.cell.betaDeg == 105.3)
        #expect(s.cell.a == 15.16)
        let l = try IndexPresets.rings(s)
        #expect(l.count == 417)
        let multiplicities3: [Int] = l.map(\.multiplicity)
        let multiplicitySum3: Int = multiplicities3.reduce(0, +)
        #expect(multiplicitySum3 == 1442)
        let matching6: Int = l.filter { (x: ReflectionFamily) -> Bool in x.multiplicity == 2 }.count
        #expect(matching6 == 113)
        let matching7: Int = l.filter { (x: ReflectionFamily) -> Bool in x.multiplicity == 4 }.count
        #expect(matching7 == 304)
        let matching8: Int = l.filter { (x: ReflectionFamily) -> Bool in x.orbits.count > 1 }.count
        #expect(matching8 == 0)
    }

    // MARK: New in the kit: systems the reference implementation refused (no app truth: closed forms planted)

    /// A tetragonal cell (a 4, c 5 Å, one Cu) is indexed: g(h k l) = sqrt((h^2 + k^2)/a^2 + l^2/c^2) in closed form; the
    /// {1 0 0} ring (4 members) at 1/a, {0 0 1} (2) at 1/c, {1 1 0} (4) at sqrt(2)/a, {1 0 1} (8) at the Pythagorean sum;
    /// the nearest rings are named as the next. The reference implementation refused this cell (`unsupportedLattice`).
    @Test func testTetragonalCellIsIndexedByClosedForm() throws {
        let (a, c) = (4.0, 5.0)
        let cu = [AtomSite(z: 29, fractional: [0, 0, 0])]
        let l = try IndexPresets.rings(IndexPresets.structure(.tetragonal(a: a, c: c), cu))
        let closed: [(MillerIndex, Double, Int, String)] = [
            (hkl(1, 0, 0), 1 / a, 4, "{1 0 0}"), (hkl(0, 0, 1), 1 / c, 2, "{0 0 1}"),
            (hkl(1, 1, 0), 2.0.squareRoot() / a, 4, "{1 1 0}"),
            (hkl(1, 0, 1), (1 / (a * a) + 1 / (c * c)).squareRoot(), 8, "{1 0 1}"),
            (hkl(1, 1, 1), (2 / (a * a) + 1 / (c * c)).squareRoot(), 8, "{1 1 1}"),
        ]
        for (m, g, n, name) in closed {
            let f = try family(l, containing: m)
            #expect(abs(f.g - g) <= 1e-13, "\(name)")
            #expect(f.multiplicity == n, "\(name)")
            #expect(f.name == name)
        }
        // {0 0 1} (g = 0.2) is followed by {1 0 0} (0.25): a plant at {0 0 1} names it, with {1 0 0} as the next, z = 250 for
        // u = 1e-3 at q = 1/c: (1/a - 1/c)/(1e-3/c) = 250 (the margin 250^2/2 = 31250 is far above ln 19).
        let q = 1 / c
        let r = try Index.index(q: q, sigma: 1e-3 * q, in: l, kMax: 1.2)
        #expect(r.best.family.name == "{0 0 1}")
        #expect(r.next?.family.name == "{1 0 0}")
        #expect(abs(r.best.z) <= 1e-12)
        #expect(abs(try #require(r.next).z - 250) <= 1e-9)
        #expect(abs(try #require(r.margin) - 31250) <= 1e-6)
        // {0 0 1} and {1 0 0} are 0.05 apart: sigma_q/q = 0.1 at the midpoint puts both within 3 sigma (z = 1.1): ambiguous
        let mid = 0.225
        #expect(try #require(refusal(verdict(mid, 0.1 * mid, l, kMax: 1.2))).name == "ambiguous")
    }

    /// An orthorhombic cell (3, 4, 5 Å): the rings {1 0 0}, {0 1 0}, {0 0 1} are three single-orbit rings of 2 at 1/a, 1/b,
    /// 1/c; each plant names its ring.
    @Test func testOrthorhombicCellIsIndexedByClosedForm() throws {
        let cu = [AtomSite(z: 29, fractional: [0, 0, 0])]
        let l = try IndexPresets.rings(IndexPresets.structure(.orthorhombic(a: 3, b: 4, c: 5), cu))
        for (m, g, name) in [(hkl(1, 0, 0), 1.0 / 3, "{1 0 0}"), (hkl(0, 1, 0), 0.25, "{0 1 0}"), (hkl(0, 0, 1), 0.2, "{0 0 1}")] {
            let f = try family(l, containing: m)
            #expect(abs(f.g - g) <= 1e-13)
            #expect(f.multiplicity == 2)
            #expect(f.name == name)
            let a = try Index.index(q: g * (1 + 1e-4), sigma: 1e-3 * g, in: l, kMax: 1.2)
            #expect(a.best.family.name == name)
            #expect(abs(a.best.z - 0.1) <= 1e-9)
        }
    }

    // MARK: Mapped to the kit tests that own them (not re-run here)

    // reference testFamilyIntensityIsTheSumOverMembers  -> ScatteringPinsTests.testStructureMatchesNumPyPins (sumF2 per ring)
    // reference testLatticeRecognitionTolerance, MonoclinicIndexingTests.testLatticeRecognition
    //                                               -> UnitCellTests.testMetricFamilyOfEachConstructor and
    //                                                  RefuterRegressionTests.testNearlyCubicCellIsTetragonalAndNoRingIsNamedTwice
    //                                                  (the kit has no unsupportedLattice: every system is indexed)
    // reference testUnsupportedLatticeAndElementAreRefusedByName, testRefusalOrderWithMonoclinicCells,
    //   testReasonNamesMonoclinicAndStillRefusesOrthorhombic
    //                                               -> ScatteringTruthTests.testRefusals (unsupportedElements before any F,
    //                                                  invalidKMax); the lattice refusals do not exist in the kit
    // reference MonoclinicIndexingTests.testCanonicalFormAndNames
    //                                               -> PointSymmetryTests.testCanonicalFormsAndNamesMatchTheApps and
    //                                                  testCanonicalIsAnOrbitInvariant
}

private extension Result {
    var failure: Failure? {
        if case .failure(let e) = self { return e }
        return nil
    }
}
