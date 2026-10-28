//
//  ZoneAxisTests.swift
//  Zone-axis primitives ("ZoneAxis", "Truth"): the 49 low-index axes by
//  value, the Cartesian axis of a cell against closed forms and the zone law, the detector basis, and the projection with
//  the flat and the curved Ewald sphere. Source mutations (each broken once): the detector-basis reference threshold, the sign of
//  the curvature term, the slab test (<= to <), the Cartesian axis over the reciprocal rows.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

struct ZoneAxisTests {

    /// reference commit 7d6fdbd1 `PhaseReferenceLibrary.lowIndexZoneAxes`, compiled from its source and printed (49, in order).
    private static let appsLowIndex: [(Int, Int, Int)] = [
        (-1, -1, -1), (-2, -2, -1), (-1, -1, 0), (-2, -2, 1), (-1, -1, 1), (-2, -1, -2), (-2, -1, -1), (-2, -1, 0),
        (-2, -1, 1), (-2, -1, 2), (-1, 0, -1), (-2, 0, -1), (-1, 0, 0), (-2, 0, 1), (-1, 0, 1), (-2, 1, -2),
        (-2, 1, -1), (-2, 1, 0), (-2, 1, 1), (-2, 1, 2), (-1, 1, -1), (-2, 2, -1), (-1, 1, 0), (-2, 2, 1),
        (-1, 1, 1), (-1, -2, -2), (-1, -2, -1), (-1, -2, 0), (-1, -2, 1), (-1, -2, 2), (-1, -1, -2), (-1, -1, 2),
        (-1, 0, -2), (-1, 0, 2), (-1, 1, -2), (-1, 1, 2), (-1, 2, -2), (-1, 2, -1), (-1, 2, 0), (-1, 2, 1),
        (-1, 2, 2), (0, -1, -1), (0, -2, -1), (0, -1, 0), (0, -2, 1), (0, -1, 1), (0, -1, -2), (0, -1, 2),
        (0, 0, -1),
    ]

    private static func gcd(_ a: Int, _ b: Int) -> Int {
        var a = a, b = b
        while b != 0 { (a, b) = (b, a % b) }
        return a
    }

    /// The four cells whose axes are checked: cubic Au, the beta'' monoclinic cell, hexagonal Mg, a general triclinic cell.
    private static var cells: [(String, UnitCell)] {
        [("Au", .cubic(a: 4.0782)), ("beta''", .monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)),
         ("Mg", .hexagonal(a: 3.2094, c: 5.2108)),
         ("triclinic", .triclinic(a: 3.1, b: 4.2, c: 5.3, alphaDeg: 80, betaDeg: 95, gammaDeg: 105))]
    }

    /// |u a1 + v a2 + w a3| by the law of cosines, from the six cell parameters (no lattice matrix).
    private static func lengthOfDirection(_ z: ZoneAxis, _ c: UnitCell) -> Double {
        let r = Double.pi / 180
        let u = Double(z.u), v = Double(z.v), w = Double(z.w)
        let sq = u * u * c.a * c.a + v * v * c.b * c.b + w * w * c.c * c.c
            + 2 * u * v * c.a * c.b * cos(c.gammaDeg * r)
            + 2 * u * w * c.a * c.c * cos(c.betaDeg * r)
            + 2 * v * w * c.b * c.c * cos(c.alphaDeg * r)
        return sq.squareRoot()
    }

    // MARK: The 49 axes

    /// 49 axes, distinct, equal by value and in order to the reference implementation's; each coprime, |u|,|v|,|w| <= 2, one of each +- pair,
    /// none zero; the rule rebuilt here (u, v, w from -2 to 2, reduced by the true gcd with the floor at the end, the first of
    /// a +- pair kept) gives the same list. Mutants: the gcd floored inside (50 axes, [0 0 2] beside [0 0 1]), both of a +-
    /// pair kept (98), the range -2...2 changed.
    @Test func testLowIndexAxesAreTheAppsFortyNine() {
        let axes = ZoneAxis.lowIndex
        #expect(axes.count == 49)
        #expect(Set(axes).count == 49)
        #expect(axes.map { [$0.u, $0.v, $0.w] } == Self.appsLowIndex.map { [$0.0, $0.1, $0.2] }, "by value, in the reference implementations' order")
        for a in axes {
            let g = Self.gcd(Self.gcd(abs(a.u), abs(a.v)), abs(a.w))
            #expect(g == 1, "\(a) is not primitive")
            #expect(max(abs(a.u), abs(a.v), abs(a.w)) <= 2)
            #expect(!(a.u == 0 && a.v == 0 && a.w == 0))
            #expect(!axes.contains(ZoneAxis(-a.u, -a.v, -a.w)), "\(a) and its opposite are both listed")
        }
        // the rule, rebuilt
        var seen = Set<ZoneAxis>()
        var out: [ZoneAxis] = []
        for u in -2...2 { for v in -2...2 { for w in -2...2 {
            if u == 0 && v == 0 && w == 0 { continue }
            let g = max(1, Self.gcd(Self.gcd(abs(u), abs(v)), abs(w)))
            let r = ZoneAxis(u / g, v / g, w / g)
            if seen.contains(r) || seen.contains(ZoneAxis(-r.u, -r.v, -r.w)) { continue }
            seen.insert(r)
            out.append(r)
        } } }
        #expect(out == axes)
        #expect(axes.first == ZoneAxis(-1, -1, -1) && axes.last == ZoneAxis(0, 0, -1))
        #expect(ZoneAxis(2, -1, 0).description == "[2 -1 0]")
    }

    // MARK: The Cartesian axis

    /// beta'' [0 0 1] is (cos beta, 0, sin beta) in the kit frame (a along x, b in xy), to 1e-12; [1 0 0] is x, [0 1 0] is y.
    /// A zone axis is a real-space direction: [1 0 1] is a1 + a3 normalised, NOT (1, 0, 1)/sqrt 2. Mutant: the axis taken over
    /// the reciprocal rows (cos beta changes sign structure), the Cartesian triple.
    @Test func testBetaDoublePrimeAxesAreRealSpaceDirections() throws {
        let beta = 105.3 * Double.pi / 180
        let cell = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)
        let c001 = try #require(ZoneAxis(0, 0, 1).cartesian(in: cell))
        #expect(abs(c001.x - cos(beta)) <= 1e-12 && abs(c001.y) <= 1e-12 && abs(c001.z - sin(beta)) <= 1e-12, "\(c001)")
        let c100 = try #require(ZoneAxis(1, 0, 0).cartesian(in: cell))
        #expect(abs(c100.x - 1) <= 1e-12 && abs(c100.y) <= 1e-12 && abs(c100.z) <= 1e-12)
        let c010 = try #require(ZoneAxis(0, 1, 0).cartesian(in: cell))
        #expect(abs(c010.x) <= 1e-12 && abs(c010.y - 1) <= 1e-12 && abs(c010.z) <= 1e-12)
        let c101 = try #require(ZoneAxis(1, 0, 1).cartesian(in: cell))
        let a1 = SIMD3<Double>(15.16, 0, 0), a3 = SIMD3<Double>(6.74 * cos(beta), 0, 6.74 * sin(beta))
        let want = (a1 + a3) / simd_length(a1 + a3)
        #expect(simd_length(c101 - want) <= 1e-12)
        #expect(simd_length(c101 - SIMD3<Double>(1, 0, 1) / 2.0.squareRoot()) > 0.1)
        // negation negates
        let m = try #require(ZoneAxis(0, 0, -1).cartesian(in: cell))
        #expect(simd_length(m + c001) <= 1e-15)
    }

    /// Cubic: the Cartesian axis is (u, v, w)/|(u, v, w)| for every one of the 49. A zero axis has no direction (nil), as
    /// the reference implementation's `cartesianZoneAxis`.
    @Test func testCubicAxesAreTheCartesianTripleAndZeroHasNoDirection() throws {
        let cell = UnitCell.cubic(a: 4.0782)
        for z in ZoneAxis.lowIndex {
            let n = try #require(z.cartesian(in: cell))
            let t = SIMD3<Double>(Double(z.u), Double(z.v), Double(z.w))
            #expect(simd_length(n - t / simd_length(t)) <= 1e-12, "\(z)")
        }
        #expect(ZoneAxis(0, 0, 0).cartesian(in: cell) == nil)
    }

    /// The zone law, in four cells and for every one of the 49 axes: |n| = 1; g . n = (h u + k v + l w)/|r_uvw| for every
    /// (h k l) with |h|, |k|, |l| <= 3, |r_uvw| by the law of cosines from the six cell parameters. This ties the Cartesian
    /// axis to the reciprocal lattice, to 1e-12. Mutant: the axis over the reciprocal rows (the law fails off cubic cells).
    @Test func testZoneLawTiesTheAxisToTheReciprocalLattice() throws {
        for (name, cell) in Self.cells {
            for z in ZoneAxis.lowIndex {
                let n = try #require(z.cartesian(in: cell))
                #expect(abs(simd_length(n) - 1) <= 1e-13, "\(name) \(z)")
                let r = Self.lengthOfDirection(z, cell)
                for h in -3...3 { for k in -3...3 { for l in -3...3 {
                    let g = cell.gVector(MillerIndex(h, k, l))
                    let want = Double(h * z.u + k * z.v + l * z.w) / r
                    if abs(simd_dot(g, n) - want) > 1e-12 {
                        Issue.record("\(name) \(z) (\(h) \(k) \(l)): g.n \(simd_dot(g, n)), zone law \(want)")
                        return
                    }
                } } }
            }
        }
    }

    // MARK: The detector basis

    /// Exact cases: n = z gives e1 = y, e2 = -x; n = x (|n.x| >= 0.9, reference y) gives e1 = z, e2 = -y; n = y gives
    /// e1 = normalize(y x x) = -z..., computed by hand below.
    @Test func testDetectorBasisExactCases() {
        let z = ZoneAxis.detectorBasis([0, 0, 1])
        #expect(z.columns.0 == [0, 1, 0] && z.columns.1 == [-1, 0, 0] && z.columns.2 == [0, 0, 1])
        let x = ZoneAxis.detectorBasis([1, 0, 0])
        #expect(x.columns.0 == [0, 0, 1] && x.columns.1 == [0, -1, 0] && x.columns.2 == [1, 0, 0])
        // n = y: reference x; y x x = -z; e1 = -z; e2 = y x (-z) = -x
        let y = ZoneAxis.detectorBasis([0, 1, 0])
        #expect(y.columns.0 == [0, 0, -1] && y.columns.1 == [-1, 0, 0] && y.columns.2 == [0, 1, 0])
    }

    /// The reference switches at |n.x| = 0.9 exactly: below it the reference is x, at 0.9 and above it is y, and the two give
    /// opposite e1 for n = (0.9, sqrt(0.19), 0). Mutant: `<` to `<=` (the threshold), 0.9 changed, the reference swapped.
    @Test func testDetectorBasisReferenceThresholdIsExactlyPointNine() {
        let at = ZoneAxis.detectorBasis([0.9, 0.19.squareRoot(), 0])
        #expect(at.columns.0 == [0, 0, 1], "|n.x| = 0.9 uses the y reference: \(at.columns.0)")
        let below = ZoneAxis.detectorBasis([(0.9).nextDown, 0.19.squareRoot(), 0])
        #expect(below.columns.0 == [0, 0, -1], "|n.x| just below 0.9 uses the x reference: \(below.columns.0)")
        // the threshold is on |n.x|: n = (-0.9, -sqrt(0.19), 0) switches the same way
        let neg = ZoneAxis.detectorBasis([-0.9, -(0.19.squareRoot()), 0])
        #expect(neg.columns.0 == [0, 0, -1], "\(neg.columns.0)")
        let negBelow = ZoneAxis.detectorBasis([-(0.9).nextDown, -(0.19.squareRoot()), 0])
        #expect(negBelow.columns.0 == [0, 0, 1], "\(negBelow.columns.0)")
    }

    /// For every one of the 49 axes in four cells: the columns are orthonormal (1e-14), the determinant is +1, the third
    /// column is n, e1 x e2 = n, and e2 = n x e1. Mutants: e2 = e1 x n (det -1), the third column not n.
    @Test func testDetectorBasisIsOrthonormalRightHandedWithNLast() throws {
        for (name, cell) in Self.cells {
            for z in ZoneAxis.lowIndex {
                let n = try #require(z.cartesian(in: cell))
                let b = ZoneAxis.detectorBasis(n)
                let e1 = b.columns.0, e2 = b.columns.1, e3 = b.columns.2
                #expect(e3 == n, "\(name) \(z)")
                #expect(abs(simd_length(e1) - 1) <= 1e-14 && abs(simd_length(e2) - 1) <= 1e-14, "\(name) \(z)")
                #expect(abs(simd_dot(e1, e2)) <= 1e-14 && abs(simd_dot(e1, n)) <= 1e-14 && abs(simd_dot(e2, n)) <= 1e-14,
                        "\(name) \(z)")
                #expect(abs(simd_determinant(b) - 1) <= 1e-14, "\(name) \(z)")
                #expect(simd_length(simd_cross(e1, e2) - n) <= 1e-14, "\(name) \(z)")
                #expect(simd_length(simd_cross(n, e1) - e2) <= 1e-15, "\(name) \(z)")
            }
        }
    }

    // MARK: The projection

    private func kept(_ s: CrystalStructure, kMax: Double = 1.2) throws -> [KinematicReflection] {
        try s.kinematicReflections(kMax: kMax).kept
    }

    /// Flat sphere: a reflection with h u + k v + l w = 0 lies in the zero-order Laue zone and has s_g = 0 (to 1e-12); the
    /// slab |s_g| <= 1e-9 keeps exactly the ZOLZ reflections, for every one of the 49 axes on Au, beta'' and Mg; the others
    /// have |s_g| >= 1/|r_uvw|. Positions: x = g.e1, y = g.e2, and for a ZOLZ reflection |(x, y)| = |g|. Mutants: the slab
    /// test `<=` to `<` is checked in the boundary test; s_g without the dot (kept set changes).
    @Test func testFlatSphereKeepsExactlyTheZeroOrderLaueZone() throws {
        for (name, s) in [("Au", IndexPresets.gold), ("beta''", IndexPresets.betaDoublePrime), ("Mg", IndexPresets.magnesium)] {
            let list = try kept(s)
            #expect(!list.isEmpty)
            for z in ZoneAxis.lowIndex {
                let n = try #require(z.cartesian(in: s.cell))
                let projected = ZoneAxis.project(list, onto: n, sgMax: 1e-9)
                let zolz = list.filter { $0.reflection.hkl.h * z.u + $0.reflection.hkl.k * z.v + $0.reflection.hkl.l * z.w == 0 }
                #expect(projected.map { $0.reflection.reflection.hkl } == zolz.map { $0.reflection.hkl }, "\(name) \(z)")
                let basis = ZoneAxis.detectorBasis(n)
                for p in projected {
                    let g = p.reflection.reflection.g
                    #expect(abs(p.sg) <= 1e-12, "\(name) \(z) \(p.reflection.reflection.hkl): s_g \(p.sg)")
                    #expect(p.x == simd_dot(g, basis.columns.0) && p.y == simd_dot(g, basis.columns.1))
                    #expect(abs((p.x * p.x + p.y * p.y).squareRoot() - simd_length(g)) <= 1e-12 * simd_length(g))
                }
                // everything else is one Laue-zone spacing or more away
                let r = Self.lengthOfDirection(z, s.cell)
                let outside = ZoneAxis.project(list, onto: n, sgMax: 100).filter { p in
                    let h = p.reflection.reflection.hkl
                    return h.h * z.u + h.k * z.v + h.l * z.w != 0
                }
                for p in outside { #expect(abs(p.sg) >= 1 / r - 1e-12, "\(name) \(z) \(p.reflection.reflection.hkl)") }
            }
        }
    }

    /// Au [0 0 1] by hand: e1 = y, e2 = -x, so a (h k 0) reflection sits at (x, y) = (k, -h)/a; s_g = l/a flat. The (2 0 0)
    /// reflection is at (0, -2/a). Mutant: e1 and e2 swapped, the sign of e2.
    @Test func testAuOnTheCubeAxisByHand() throws {
        let a = 4.0782
        let list = try kept(IndexPresets.gold)
        let projected = ZoneAxis.project(list, onto: [0, 0, 1], sgMax: 1e-9)
        let inPlane: [KinematicReflection] = list.filter { (k: KinematicReflection) -> Bool in k.reflection.hkl.l == 0 }
        let projectedCount: Int = projected.count
        #expect(projectedCount == inPlane.count)
        let r200 = try #require(projected.first { $0.reflection.reflection.hkl == MillerIndex(2, 0, 0) })
        #expect(abs(r200.x) <= 1e-15 && abs(r200.y - -2 / a) <= 1e-14, "\(r200.x) \(r200.y)")
        let r020 = try #require(projected.first { $0.reflection.reflection.hkl == MillerIndex(0, 2, 0) })
        #expect(abs(r020.x - 2 / a) <= 1e-14 && abs(r020.y) <= 1e-15)
        let off = ZoneAxis.project(list, onto: [0, 0, 1], sgMax: 0.3)
        let r111 = try #require(off.first { $0.reflection.reflection.hkl == MillerIndex(1, 1, 1) })
        #expect(abs(r111.sg - 1 / a) <= 1e-14)
    }

    /// Curved sphere: s_g = g . n + lambda |g|^2 / 2 for every kept reflection, in input order, at 120 and 300 kV
    /// wavelengths; at the zone law s_g = lambda |g|^2 / 2 > 0, so the ZOLZ reflections leave the slab when sgMax is below it;
    /// g and -g are excited differently by lambda |g|^2 (s_g(g) + s_g(-g) = lambda |g|^2), where the flat sphere has them
    /// opposite (sum 0). Mutants: the sign of the curvature term, lambda without the half.
    @Test func testCurvedSphereAddsHalfLambdaGSquared() throws {
        let s = IndexPresets.betaDoublePrime
        let list = try kept(s)
        for lambda in [0.03349215272672374, 0.01968748900684880] {
            let n = try #require(ZoneAxis(0, 1, 0).cartesian(in: s.cell))
            let projected = ZoneAxis.project(list, onto: n, wavelengthAngstroms: lambda, sgMax: .infinity)
            #expect(projected.count == list.count)
            for (p, r) in zip(projected, list) {
                let g = r.reflection.g
                let want = simd_dot(g, n) + lambda * simd_length_squared(g) / 2
                #expect(abs(p.sg - want) <= 1e-15 * max(1, abs(want)), "\(r.reflection.hkl)")
                #expect(p.sg > simd_dot(g, n), "the curvature raises s_g")
            }
            var byHKL: [MillerIndex: Double] = [:]
            for p in projected { byHKL[p.reflection.reflection.hkl] = p.sg }
            for (h, sg) in byHKL {
                let mate = try #require(byHKL[h.negated])
                let g = s.cell.gVector(h)
                #expect(abs((sg + mate) - lambda * simd_length_squared(g)) <= 1e-15, "\(h)")
            }
            // the ZOLZ of [0 1 0]: s_g = lambda |g|^2 / 2; a slab below it drops them, a slab above keeps them
            let zolz = list.filter { $0.reflection.hkl.k == 0 }
            let g0 = zolz.map { simd_length_squared($0.reflection.g) }.min() ?? 0
            let slab = lambda * g0 / 2
            let narrow = ZoneAxis.project(zolz, onto: n, wavelengthAngstroms: lambda, sgMax: slab * 0.99)
            let wide = ZoneAxis.project(zolz, onto: n, wavelengthAngstroms: lambda, sgMax: slab * 1.01)
            #expect(narrow.isEmpty)
            #expect(!wide.isEmpty)
        }
        // nil wavelength is the flat sphere
        let n = try #require(ZoneAxis(1, 0, 0).cartesian(in: s.cell))
        let flat = ZoneAxis.project(list, onto: n, sgMax: 0.2)
        let none = ZoneAxis.project(list, onto: n, wavelengthAngstroms: nil, sgMax: 0.2)
        let zero = ZoneAxis.project(list, onto: n, wavelengthAngstroms: 0, sgMax: 0.2)
        #expect(flat == none && flat == zero)
    }

    /// The DEVIATION in sign, stated and tested: py4DSTEM f050d207 writes s_g = (2 g_z - lambda |g|^2)/(2 - 2 lambda g_z) with
    /// k0 = [0, 0, -1/lambda], and its lab z is the NEGATIVE of this beam direction, so g_z = -(g . n). To first order
    /// (the denominator dropped) its s_g is then -(g . n) - lambda |g|^2 / 2, exactly the negative of this kit's s_g. A
    /// mutation with the opposite curvature sign fails here and in testCurvedSphereAddsHalfLambdaGSquared.
    @Test func testCurvatureSignIsTheNegativeOfPy4DSTEMsInTheLabFrameOfTheBeam() throws {
        let s = IndexPresets.magnesium
        let n = try #require(ZoneAxis(0, 0, 1).cartesian(in: s.cell))
        let lambda = 0.02507934045054800
        for p in ZoneAxis.project(try kept(s), onto: n, wavelengthAngstroms: lambda, sgMax: .infinity) {
            let g = p.reflection.reflection.g
            let gz = -simd_dot(g, n)
            let py4dstem = (2 * gz - lambda * simd_length_squared(g)) / 2
            #expect(abs(p.sg + py4dstem) <= 1e-15, "\(p.reflection.reflection.hkl)")
        }
    }

    /// The slab is |s_g| <= sgMax, inclusive: with sgMax equal to the s_g of one reflection that reflection is kept, and with
    /// sgMax one ulp below it is not. Mutant: `<=` to `<`.
    @Test func testSlabIsInclusive() throws {
        let s = IndexPresets.gold
        let list = try kept(s)
        let lambda = 0.0197
        let all = ZoneAxis.project(list, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: .infinity)
        let probe = try #require(all.first { $0.reflection.reflection.hkl == MillerIndex(1, 1, 1) })
        let edge = abs(probe.sg)
        let at = ZoneAxis.project(list, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: edge)
        #expect(at.contains { $0.reflection.reflection.hkl == MillerIndex(1, 1, 1) }, "|s_g| = sgMax is kept")
        let below = ZoneAxis.project(list, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: edge.nextDown)
        #expect(!below.contains { $0.reflection.reflection.hkl == MillerIndex(1, 1, 1) }, "one ulp below it is not")
        // a negative s_g is cut by its absolute value
        let neg = try #require(all.first { $0.sg < 0 })
        #expect(ZoneAxis.project(list, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: abs(neg.sg))
            .contains { $0.reflection.reflection.hkl == neg.reflection.reflection.hkl })
        // a NaN excitation error is never kept
        #expect(ZoneAxis.project(list, onto: [0, 0, 1], wavelengthAngstroms: .nan, sgMax: .infinity).isEmpty)
    }
}
