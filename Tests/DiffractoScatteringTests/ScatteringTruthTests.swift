//
//  ScatteringTruthTests.swift
//  Planted truth for the structure factors, independent of every oracle: exact
//  extinctions, closed forms of F, the reference implementations' ring counts after the |F| cut, the weak cut, and f_e's own shape.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

enum Presets {
    static let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    static func structure(_ cell: UnitCell, _ sites: [AtomSite]) -> CrystalStructure {
        CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
    }
    static func fccElement(_ z: Int, a: Double) -> CrystalStructure {
        structure(.cubic(a: a), fcc.map { AtomSite(z: z, fractional: $0) })
    }
    static var gold: CrystalStructure { fccElement(79, a: 4.0782) }
    static var aluminium: CrystalStructure { fccElement(13, a: 4.0495) }
    static var silicon: CrystalStructure {
        structure(.cubic(a: 5.4309), (fcc + fcc.map { $0 + SIMD3(0.25, 0.25, 0.25) }).map { AtomSite(z: 14, fractional: $0) })
    }
    static func hcp(_ z: Int, a: Double, c: Double) -> CrystalStructure {
        structure(.hexagonal(a: a, c: c), [AtomSite(z: z, fractional: [0, 0, 0]),
                                            AtomSite(z: z, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5])])
    }
    static var magnesium: CrystalStructure { hcp(12, a: 3.2094, c: 5.2108) }
    static var idealHCP: CrystalStructure { hcp(12, a: 3.2094, c: 3.2094 * (8.0 / 3.0).squareRoot()) }
    /// the reference implementation's crystal-parity occupancy case: Au 0.6 / Cu 0.4 on the origin, Cu on the faces, a = 3.75 Å.
    static var cu3AuOcc: CrystalStructure {
        structure(.cubic(a: 3.75), [AtomSite(z: 79, fractional: [0, 0, 0], occupancy: 0.6),
                                    AtomSite(z: 29, fractional: [0, 0, 0], occupancy: 0.4)]
                  + fcc.dropFirst().map { AtomSite(z: 29, fractional: $0) })
    }
    /// The β″ CIF exactly as the generator writes it (the core fixture builder).
    static var betaDoublePrime: CrystalStructure {
        let sites: [(String, String, Double, Double, Double, Double)] = [
            ("Mg1", "Mg", 0.0, 0.0, 0.0, 0.5), ("Mg2", "Mg", 0.3459, 0.0, 0.089, 1.0), ("Mg3", "Mg", 0.430, 0.0, 0.652, 0.8),
            ("Si1", "Si", 0.0565, 0.0, 0.649, 1.1), ("Si2", "Si", 0.1885, 0.0, 0.224, 0.5), ("Si3", "Si", 0.2171, 0.0, 0.617, 2.5),
        ]
        var lines = ["data_beta_double_prime_Mg5Si6", "_cell_length_a 15.16", "_cell_length_b 4.05", "_cell_length_c 6.74",
                     "_cell_angle_alpha 90.0", "_cell_angle_beta 105.3", "_cell_angle_gamma 90.0",
                     "_symmetry_space_group_name_H-M 'C 2/m'", "_symmetry_Int_Tables_number 12",
                     "loop_", "_symmetry_equiv_pos_as_xyz"]
        lines += ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z",
                  "1/2+x,1/2-y,z"].map { "  '\($0)'" }
        lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
                  "_atom_site_fract_z"]
        lines += sites.map { "  \($0.0) \($0.1) \($0.2) \($0.3) \($0.4)" }
        // The CIF parse cannot fail on this fixed text; a failure is a CIF-parser regression the CIF tests name.
        return (try? CIFImport.read(lines.joined(separator: "\n")).structure)!
    }
}

private func fe(_ z: Int, _ g: Double) -> Double { ElectronScatteringFactors.f(z: z, gSquared: g * g)! }

private func relativelyClose(_ x: Double, _ y: Double, _ tol: Double) -> Bool {
    abs(x - y) <= tol * max(abs(x), abs(y), 1e-300)
}

struct ScatteringTruthTests {

    /// fcc: |F| ≤ 1e-12 for mixed parity, F = 4 f_e / V (real) otherwise. Mutants: no 1/V; occupancy ignored (no move
    /// here); the table row Z − 1 → Z.
    @Test func testFCCExtinctionsAndClosedForm() throws {
        for (name, s, z) in [("Au", Presets.gold, 79), ("Al", Presets.aluminium, 13)] {
            let list = try s.kinematicReflections(kMax: 1.2)
            for k in list.all {
                let h = k.reflection.hkl, f = k.structureFactor
                let parity = [h.h, h.k, h.l].map { (($0 % 2) + 2) % 2 }.reduce(0, +)
                if parity == 1 || parity == 2 {
                    #expect(f.magnitude <= 1e-12, "\(name) \(h) mixed parity |F| \(f.magnitude)")
                } else {
                    let want = 4 * fe(z, f.gLength) / s.cell.volume
                    #expect(relativelyClose(f.real, want, 1e-13), "\(name) \(h)")
                    #expect(abs(f.imaginary) <= 1e-13 * want, "\(name) \(h)")
                }
            }
        }
    }

    /// Diamond: also h + k + l = 4n + 2 (all even) absent; |F| = 4√2 f_e / V for h + k + l odd, 8 f_e / V for 4n.
    /// Mutant: the diamond basis dropped (the 4n + 2 reflections appear).
    @Test func testDiamondExtinctionsAndClosedForm() throws {
        let s = Presets.silicon
        for k in try s.kinematicReflections(kMax: 1.2).all {
            let h = k.reflection.hkl, f = k.structureFactor
            let p = [h.h, h.k, h.l].map { (($0 % 2) + 2) % 2 }.reduce(0, +)
            let sum = ((h.h + h.k + h.l) % 4 + 4) % 4
            let base = fe(14, f.gLength) / s.cell.volume
            if p == 1 || p == 2 || (p == 0 && sum == 2) {
                #expect(f.magnitude <= 1e-12, "Si \(h) |F| \(f.magnitude)")
            } else if p == 3 {
                #expect(relativelyClose(f.magnitude, 4 * 2.0.squareRoot() * base, 1e-13), "Si \(h)")
            } else {
                #expect(relativelyClose(f.magnitude, 8 * base, 1e-13), "Si \(h)")
            }
        }
    }

    /// hcp: |F| = 2 f_e |cos(π((h + 2k)/3 + l/2))| / V; 000l with l odd and h − k = 3n with l odd are absent.
    /// Mutant: the phase sign flipped (no move on |F|; Im F is checked by the pins).
    @Test func testHCPClosedForm() throws {
        let s = Presets.magnesium
        for k in try s.kinematicReflections(kMax: 1.2).all {
            let h = k.reflection.hkl, f = k.structureFactor
            let want = 2 * fe(12, f.gLength) * abs(cos(Double.pi * (Double(h.h + 2 * h.k) / 3 + Double(h.l) / 2)))
                / s.cell.volume
            #expect(abs(f.magnitude - want) <= 1e-13 * max(want, 1e-3), "Mg \(h)")
            if ((h.l % 2) + 2) % 2 == 1 && ((h.h - h.k) % 3 + 3) % 3 == 0 { #expect(f.magnitude <= 1e-12, "Mg \(h)") }
        }
    }

    /// The occupancy case: unmixed (0.6 f_Au + 0.4 f_Cu + 3 f_Cu)/V; mixed (0.6 f_Au + 0.4 f_Cu − f_Cu)/V, the
    /// superlattice reflections that exist only through the occupancy. Mutant: occupancy ignored.
    @Test func testOccupancyClosedForm() throws {
        let s = Presets.cu3AuOcc
        for k in try s.kinematicReflections(kMax: 1.2).all {
            let h = k.reflection.hkl, f = k.structureFactor
            let p = [h.h, h.k, h.l].map { (($0 % 2) + 2) % 2 }.reduce(0, +)
            let au = fe(79, f.gLength), cu = fe(29, f.gLength)
            let want = (p == 0 || p == 3 ? 0.6 * au + 0.4 * cu + 3 * cu : 0.6 * au + 0.4 * cu - cu) / s.cell.volume
            #expect(relativelyClose(f.real, want, 1e-13), "Cu3Au \(h)")
            #expect(abs(f.imaginary) <= 1e-13 * abs(want), "Cu3Au \(h)")
        }
    }

    /// The reference implementations' ring counts after the cut (k_max 1.2 Å⁻¹): Au 8, Al 8, Si 10, Mg 24, ideal hcp 23, β″ 417 rings and
    /// 1442 reflections. Mutants: families on the uncut list (Si 12, Mg 25); the cut at 1e-3; the cut ≥ instead of >.
    @Test func testRingCountsAfterTheCut() throws {
        let cases: [(String, CrystalStructure, Int)] = [
            ("Au", Presets.gold, 8), ("Al", Presets.aluminium, 8), ("Si", Presets.silicon, 10),
            ("Mg", Presets.magnesium, 24), ("ideal hcp", Presets.idealHCP, 23),
        ]
        for (name, s, rings) in cases {
            let list = try s.kinematicReflections(kMax: 1.2)
            #expect(list.families.count == rings, "\(name): \(list.families.map(\.family.name))")
        }
        let si = try Presets.silicon.kinematicReflections(kMax: 1.2).families.map(\.family.name)
        #expect(!si.contains("{2 2 2}") && !si.contains("{4 4 2}") && si.contains("{5 1 1} + {3 3 3}"), "\(si)")
        let mg = try Presets.magnesium.kinematicReflections(kMax: 1.2).families.map(\.family.name)
        #expect(!mg.contains("{3 0 -3 1}"), "\(mg)")

        let beta = try Presets.betaDoublePrime.kinematicReflections(kMax: 1.2)
        #expect(beta.kept.count == 1442)
        #expect(beta.families.count == 417)
        let dropped = Set(beta.all.filter { $0.reflection.isAllowed && $0.structureFactor.magnitude <= 1e-4 }
            .map { Presets.betaDoublePrime.laueGroup.canonical($0.reflection.hkl) })
        #expect(dropped == [MillerIndex(3, 3, 3), MillerIndex(14, 0, 1)], "\(dropped)")
    }

    /// β″ at 1.6 and 2.0 Å⁻¹: 953 and 1805 rings (the reference implementation's predictions).
    @Test func testBetaDoublePrimeRingsAtHigherKMax() throws {
        let s = Presets.betaDoublePrime
        #expect(try s.kinematicReflections(kMax: 1.6).families.count == 953)
        #expect(try s.kinematicReflections(kMax: 2.0).families.count == 1805)
    }

    /// Every operation-absent reflection has |F| ≤ 1e-12, and a single orbit is wholly kept or wholly dropped.
    @Test func testAbsencesVanishAndOrbitsAreCutWhole() throws {
        for s in [Presets.gold, Presets.silicon, Presets.magnesium, Presets.betaDoublePrime, Presets.cu3AuOcc] {
            let list = try s.kinematicReflections(kMax: 1.2)
            for k in list.all where !k.reflection.isAllowed { #expect(k.structureFactor.magnitude <= 1e-12) }
            let kept = Set(list.kept.map(\.reflection.hkl))
            for k in list.all {
                let orbit = s.laueGroup.orbit(of: k.reflection.hkl)
                let inside = orbit.filter { kept.contains($0) }.count
                #expect(inside == 0 || inside == orbit.count, "\(k.reflection.hkl): \(inside) of \(orbit.count) kept")
            }
            let sums = list.families.map { $0.intensity }
            #expect(sums.allSatisfy { $0 > 0 })
        }
    }

    /// The weak cut (the reference implementation's case): Cu on (0,0,0) and on (½,½,½) at occupancy q, q solved so (1 0 0) has |F| = 2.5e-4;
    /// the h + k + l odd reflections straddle 1e-4, each on the side (1 − q) f_e / V puts it. Mutant: the cut at 1e-3.
    @Test func testWeakCutStraddlesTheCut() throws {
        let a = 3.0
        let q = 1 - 2.5e-4 * a * a * a / fe(29, 1 / a)
        let s = Presets.structure(.cubic(a: a), [AtomSite(z: 29, fractional: [0, 0, 0]),
                                                 AtomSite(z: 29, fractional: [0.5, 0.5, 0.5], occupancy: q)])
        let list = try s.kinematicReflections(kMax: 1.2)
        let odd = list.all.filter { (($0.reflection.hkl.h + $0.reflection.hkl.k + $0.reflection.hkl.l) % 2 + 2) % 2 == 1 }
        let keptIn = odd.filter { $0.structureFactor.magnitude > 1e-4 && $0.structureFactor.magnitude < 1e-3 }.count
        let droppedIn = odd.filter { $0.structureFactor.magnitude < 1e-4 && $0.structureFactor.magnitude > 1e-5 }.count
        #expect(keptIn >= 3 && droppedIn >= 3, "kept \(keptIn), dropped \(droppedIn)")
        let kept = Set(list.kept.map(\.reflection.hkl))
        for k in odd {
            let want = (1 - q) * fe(29, k.structureFactor.gLength) / s.cell.volume
            #expect(relativelyClose(k.structureFactor.magnitude, want, 1e-12), "\(k.reflection.hkl)")
            #expect(kept.contains(k.reflection.hkl) == (want > 1e-4), "\(k.reflection.hkl) |F| \(want)")
        }
    }

    /// f_e is positive and strictly decreasing on [0, 1.2] Å⁻¹; f_e(Z, 0) = 2 Σ aᵢ; outside the table it is nil.
    @Test func testScatteringFactorShape() {
        for z in [1, 6, 12, 13, 14, 29, 74, 79, 103] {
            var last = Double.infinity
            for i in 0...240 {
                let g = 1.2 * Double(i) / 240
                let f = fe(z, g)
                #expect(f > 0 && f < last, "Z \(z) at g \(g)")
                last = f
            }
            let p = ElectronScatteringFactors.parameters[z - 1]
            #expect(relativelyClose(fe(z, 0), 2 * (p[0] + p[2] + p[4] + p[6] + p[8]), 1e-15), "Z \(z)")
        }
        #expect(ElectronScatteringFactors.f(z: 104, gSquared: 0) == nil)
        #expect(ElectronScatteringFactors.f(z: 0, gSquared: 0) == nil)
        #expect(ElectronScatteringFactors.parameters.count == 103)
        #expect(ElectronScatteringFactors.parameters.allSatisfy { $0.count == 10 })
    }

    /// Refusals: an element without a form factor, before any F; a bad cut; DiffractoCore's k_max refusal passed through.
    @Test func testRefusals() throws {
        let rf = Presets.structure(.cubic(a: 4), [AtomSite(z: 104, fractional: [0, 0, 0])])
        #expect(throws: KinematicError.unsupportedElements([104])) { try rf.kinematicReflections(kMax: 1.2) }
        // Refused before any F, even when no reflection lies within k_max (an empty list would hide the missing row).
        #expect(throws: KinematicError.unsupportedElements([104])) { try rf.kinematicReflections(kMax: 0.01) }
        #expect(throws: KinematicError.invalidCut(-1)) { try Presets.gold.kinematicReflections(kMax: 1.2, cut: -1) }
        #expect(throws: KinematicError.reflectionList(.invalidKMax(5))) { try Presets.gold.kinematicReflections(kMax: 5) }
        let one = try #require(StructureFactors.compute(MillerIndex(1, 1, 1), in: Presets.gold))
        #expect(relativelyClose(one.real, 4 * fe(79, one.gLength) / Presets.gold.cell.volume, 1e-13))
        #expect(StructureFactors.compute(MillerIndex(1, 1, 1), in: rf) == nil)
    }
}
