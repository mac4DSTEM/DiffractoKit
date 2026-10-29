//
//  RefuterRegressionTests.swift
//  An independent review's findings, each as a regression
//  test on the refuter's own CIF (Fixtures/refuter_*.cif). Truth: pymatgen-core 2026.9.23 / spglib 2.8.0 on the
//  UNSHIFTED structure (absences and families do not depend on the origin), computed independently of the refuter's
//  harness: Mg (c = 5.2112) P6₃/mmc 24 operations, 282 of 324 allowed at 1.2 Å⁻¹; one Fe atom in a 5 Å cube Pm-3m
//  48 operations, 924 of 924; the 7× Po supercell (c = 23.45) 112 operations, 122 of 782 allowed at 0.9.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct RefuterRegressionTests {

    private func read(_ name: String) throws -> CrystalStructure {
        try CIFImport.read(try Fixtures.text("refuter_\(name).cif"), fileBaseName: name).structure
    }

    /// Two structures that differ by an origin shift must agree on everything a pattern shows.
    private func assertSameDiffraction(_ a: CrystalStructure, _ b: CrystalStructure, kMax: Double, _ label: String,
                                       sourceLocation: SourceLocation = #_sourceLocation) throws {
        #expect(a.operations.count == b.operations.count, "\(label) operations", sourceLocation: sourceLocation)
        #expect(a.laueClass == b.laueClass, "\(label) Laue class", sourceLocation: sourceLocation)
        let la = try a.reflections(kMax: kMax), lb = try b.reflections(kMax: kMax)
        #expect(la.reflections.map(\.hkl) == lb.reflections.map(\.hkl), "\(label) list", sourceLocation: sourceLocation)
        let wrong = zip(la.reflections, lb.reflections).filter { $0.isAllowed != $1.isAllowed }.map(\.0.hkl)
        #expect(wrong.isEmpty, "\(label): \(wrong.count) allowed flags differ, e.g. \(wrong.prefix(4))", sourceLocation: sourceLocation)
        #expect(la.families().map(\.name) == lb.families().map(\.name), "\(label) families", sourceLocation: sourceLocation)
        #expect(la.families().map(\.multiplicity) == lb.families().map(\.multiplicity), "\(label) multiplicities", sourceLocation: sourceLocation)
    }

    /// F1(b): hcp Mg shifted by (0.1, 0.2, 0.3), 8 decimals. A translation of 0.20000001 failed the 1e-9 absence test
    /// for (2 −1 0): 262 allowed instead of 282, ring {1 1 −2 0} of 2 instead of 6. Mutant: the absence tolerance
    /// without the translation's uncertainty.
    @Test func testShiftedMagnesiumKeepsItsAbsences() throws {
        let shifted = try read("F2a_Mg_shifted")
        #expect(shifted.operations.count == 24)
        #expect(shifted.laueClass == .sixOverMmm)
        let list = try shifted.reflections(kMax: 1.2)
        #expect(list.reflections.count == 324)
        #expect(list.allowed.count == 282)
        let ring = try #require(list.families().first { $0.name == "{1 1 -2 0}" })
        #expect(ring.multiplicity == 6)
        let unshifted = CrystalStructure.detectingSymmetry(
            cell: shifted.cell, sites: [AtomSite(z: 12, fractional: [1.0 / 3, 2.0 / 3, 0.25]),
                                        AtomSite(z: 12, fractional: [2.0 / 3, 1.0 / 3, 0.75])])
        try assertSameDiffraction(shifted, unshifted, kMax: 1.2, "Mg shifted")
    }

    /// F1(b) with a file loop: the same Mg with 24 exact rational operations at another origin. The file's exact
    /// translations are used for the operations it names (the 8-decimal coordinates alone also give 282 now).
    @Test func testShiftedMagnesiumWithItsLoopUsesTheExactTranslations() throws {
        let s = try read("F2f_Mg_loop_0.37")
        #expect(s.fileOperations?.count == 24)
        #expect(s.operations.count == 24)
        #expect(s.fileOperationsAreSymmetries == true)
        #expect(try s.reflections(kMax: 1.2).allowed.count == 282)
    }

    /// F1(a): a lone atom at a general position is an origin shift of Pm-3m: 48 operations, every reflection allowed,
    /// ring {1 1 0} of 12. The fixed k/24 snap had moved translations by up to 5e-3: 43 operations, 80 wrong flags at 1.2.
    /// Mutant: the snap within the position tolerance instead of the translation's own uncertainty.
    @Test func testOneAtomAtAGeneralPositionIsPmThreeBarM() throws {
        let s = try read("P1_cubic_1fe")
        #expect(s.operations.count == 48)
        #expect(s.laueClass == .mThreeBarM)
        // At 1.19 Å⁻¹ (h² + k² + l² ≤ 35; 36 sits on the float edge of 1.2): every lattice point, all allowed.
        let list = try s.reflections(kMax: 1.19)
        var points = 0
        for h in -6...6 { for k in -6...6 { for l in -6...6 where h * h + k * k + l * l <= 35 && (h, k, l) != (0, 0, 0) {
            points += 1
        } } }
        #expect(list.reflections.count == points)
        #expect(list.allowed.count == points)
        let ring = try #require(list.families().first { $0.name == "{1 1 0}" })
        #expect(ring.multiplicity == 12)
    }

    /// F1(a): rock salt shifted by (0.9172, 0.6292, 0.5141), exact 4-decimal coordinates: Fm-3m's 192 operations and
    /// the unshifted structure's absences and rings (the snap had given 160 operations, no Laue class, 24 wrong flags).
    @Test func testShiftedRockSaltIsFmThreeBarM() throws {
        let shifted = try read("F2d_NaCl_shifted")
        #expect(shifted.operations.count == 192)
        #expect(shifted.laueClass == .mThreeBarM)
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let unshifted = CrystalStructure.detectingSymmetry(
            cell: shifted.cell,
            sites: fcc.map { AtomSite(z: 11, fractional: $0) }
                + fcc.map { AtomSite(z: 17, fractional: wrapFractional($0 + SIMD3(0.5, 0.5, 0.5))) })
        try assertSameDiffraction(shifted, unshifted, kMax: 1.2, "NaCl shifted")
    }

    /// F1(a): gold shifted by (0.2484, 0.3, 0.07): the inversion's translation is 0.4968, not 0.5.
    @Test func testShiftedGoldKeepsItsOffGridTranslations() throws {
        let s = try read("F2c_Au_shifted_exactdecimals")
        #expect(s.operations.count == 192)
        #expect(s.laueClass == .mThreeBarM)
        let inversions = s.operations.filter { $0.rotation == .inversion }
        #expect(inversions.count == 4)
        #expect(inversions.contains { abs($0.translationMod1.x - 0.4968) < 1e-9 }, "\(inversions.map(\.xyzString))")
    }

    /// F1(a): quartz shifted: P3₂21's 6 operations, Laue -3m (the snap had given 5 and no class).
    @Test func testShiftedQuartzIsThreeBarM() throws {
        let s = try read("F2e_quartz_shifted")
        #expect(s.operations.count == 6)
        #expect(s.laueClass == .threeBarM)
    }

    /// F1(b): a 7× supercell has translations of k/7, off the 1/24 grid, written to 8 decimals: 112 operations and
    /// 122 allowed of 782 at 0.9 Å⁻¹ (the 1e-9 test had flagged every l = ±7 sublattice reflection absent: 28).
    @Test func testSevenfoldSupercellKeepsItsSublatticeReflections() throws {
        let s = try read("F3_super7")
        #expect(s.operations.count == 112)
        #expect(s.laueClass == .fourOverMmm, "a 1 × 1 × 7 cell is tetragonal by its metric; spglib's 221 is the primitive cubic cell")
        let list = try s.reflections(kMax: 0.9)
        #expect(list.reflections.count == 782)
        #expect(list.allowed.count == 122)
        #expect(list.allowed.allSatisfy { $0.hkl.l % 7 == 0 })
    }

    /// F2: two species sharing one site (Fe2+ 0.5 and Fe3+ 0.5 at the origin) under F-centring keep both rows: 8 Fe
    /// sites of 0.5 (Fe 4.0), 4 O; the expansion had merged the rows (Fe 2.0). Mutant: images of any row merge.
    @Test func testCoSitedSpeciesAreNotMerged() throws {
        let s = try read("mixed_centred")
        #expect(s.sites.count == 12)
        let ironSites: [AtomSite] = s.sites.filter { $0.z == 26 }
        let ironOccupancy: Double = ironSites.reduce(0) { $0 + $1.occupancy }
        #expect(abs(ironOccupancy - 4.0) <= 1e-12)
        let oxygenSites: [AtomSite] = s.sites.filter { $0.z == 8 }
        #expect(oxygenSites.count == 4)
        let labels: Set<String> = Set(s.sites.compactMap(\.label))
        let expectedLabels: Set<String> = ["Fe2", "Fe3", "O1"]
        #expect(labels == expectedLabels)
    }

    /// F4: corundum in rhombohedral axes is -3c by its atoms (Laue -3m) and its rings are named in three indices; the
    /// Miller–Bravais form is for hexagonal axes only. Mutant: the four-index form by the Laue class alone.
    @Test func testRhombohedralAxesKeepThreeIndexNames() throws {
        let s = try read("corundum_rh_p1")
        #expect(s.sites.count == 10)
        #expect(s.laueClass == .threeBarM)
        #expect(s.operations.count == 12)
        let list = try s.reflections(kMax: 1.2)
        #expect(!list.usesMillerBravais)
        for family in list.families() {
            for orbit in family.name.components(separatedBy: " + ") {
                #expect(orbit.split(separator: " ").count == 3, "\(family.name)")
            }
        }
        #expect(try CIFImport.read(try Fixtures.text("calcite_R-3c_hex.cif")).structure.reflections(kMax: 1.2).usesMillerBravais)
    }

    /// F5: a finite but huge cell is refused, never trapped or left to hang; a large list is refused before it is built.
    @Test func testHugeCellsAreRefusedNotTrapped() {
        let cases: [(Double, Double)] = [(1e30, 1.2), (1e18, 1.2), (1e10, 1.2), (100, 4)]
        for (a, k) in cases {
            let s = CrystalStructure.detectingSymmetry(cell: .cubic(a: a), sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
            let error = #expect(throws: ReflectionListError.self, "a = \(a)") { try s.reflections(kMax: k) }
            guard case .cellTooLargeOrInvalid? = error else {
                Issue.record("a = \(a): \(String(describing: error))")
                continue
            }
        }
        let thin = CrystalStructure.detectingSymmetry(cell: .orthorhombic(a: 100, b: 3, c: 3),
                                                     sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
        #expect(throws: Never.self, "a long thin cell is a small list") { try thin.reflections(kMax: 4) }
    }

    /// F5: an operation whose rotation part is singular or a dilation is not a symmetry and is refused.
    @Test func testNonInvertibleOperationsAreRefused() {
        for raw in ["x,x,z", "2x,y,z", "x+1000y,y+1000x,z"] {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: raw) }
            #expect(error == .notInvertible(raw))
        }
        #expect(throws: Never.self) { try SymmetryOperation(xyz: "x-y,x,z") }
    }

    /// F1 fix, the grid: magnesium at the standard origin written to four decimals (⅓ as 0.3333) has every translation
    /// exactly on the k/24 grid (within its own uncertainty), so the absence test is exact. Mutant: the grid step dropped.
    @Test func testCoarseCoordinatesGiveExactGridTranslations() throws {
        let cif = """
        data_mg4
        _cell_length_a 3.2094
        _cell_length_b 3.2094
        _cell_length_c 5.2108
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 120
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Mg1 0.3333 0.6667 0.25
        Mg2 0.6667 0.3333 0.75
        """
        let s = try CIFImport.read(cif).structure
        #expect(s.operations.count == 24)
        for op in s.operations {
            let scaled = op.translationMod1 * 24
            let worst = max(abs(scaled.x - scaled.x.rounded()), abs(scaled.y - scaled.y.rounded()), abs(scaled.z - scaled.z.rounded()))
            #expect(abs(worst - 0) <= 1e-12, "\(op.xyzString)")
        }
        #expect(s.translationUncertainty == 1e-9)
    }

    /// F1 fix, the tolerance: the shifted magnesium written to SIX decimals, where the least-squares refinement does
    /// not cancel the rounding (two sites), keeps P6₃/mmc's 282 allowed reflections. Mutant: the absence tolerance
    /// without the translations' uncertainty.
    @Test func testShiftedMagnesiumAtSixDecimalsKeepsItsAbsences() throws {
        let cif = """
        data_mg6
        _cell_length_a 3.2094
        _cell_length_b 3.2094
        _cell_length_c 5.2112
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 120
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Mg1 0.433333 0.866667 0.55
        Mg2 0.766667 0.533333 0.05
        """
        let s = try CIFImport.read(cif).structure
        #expect(s.operations.count == 24)
        #expect(s.translationUncertainty > 1e-9, "the shifted translations are off the grid and carry the rounding")
        #expect(try s.reflections(kMax: 1.2).allowed.count == 282)
    }

    /// F1 fix, the tolerance where refinement cannot help: rock salt shifted by an arbitrary vector and written to FIVE
    /// decimals, so each of the 8 sites rounds differently and the refined translations keep about 1e-6 of error. Its
    /// absences and rings equal the unshifted structure's. Mutant: the absence tolerance without the translations'
    /// uncertainty.
    @Test func testShiftedRockSaltAtFiveDecimalsKeepsItsAbsences() throws {
        let shift = SIMD3<Double>(0.123456789, 0.234567891, 0.345678912)
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let na = fcc.map { wrapFractional($0 + shift) }, cl = fcc.map { wrapFractional($0 + SIMD3(0.5, 0.5, 0.5) + shift) }
        var rows = ""
        for (i, p) in na.enumerated() { rows += String(format: "Na%d Na %.5f %.5f %.5f\n", i + 1, p.x, p.y, p.z) }
        for (i, p) in cl.enumerated() { rows += String(format: "Cl%d Cl %.5f %.5f %.5f\n", i + 1, p.x, p.y, p.z) }
        let cif = """
        data_nacl5
        _cell_length_a 5.64
        _cell_length_b 5.64
        _cell_length_c 5.64
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        \(rows)
        """
        let shifted = try CIFImport.read(cif).structure
        #expect(shifted.operations.count == 192)
        let unshifted = CrystalStructure.detectingSymmetry(
            cell: shifted.cell, sites: fcc.map { AtomSite(z: 11, fractional: $0) }
                + fcc.map { AtomSite(z: 17, fractional: wrapFractional($0 + SIMD3(0.5, 0.5, 0.5))) })
        try assertSameDiffraction(shifted, unshifted, kMax: 1.2, "NaCl shifted, 5 decimals")
    }

    /// Every detected operation maps the site set onto itself within the position tolerance, on every fixture and the
    /// apps' split gold. Mutant: the final translation not re-verified after refinement and snapping.
    @Test func testEveryDetectedOperationIsASymmetryOfTheSites() throws {
        var structures: [(String, CrystalStructure)] = []
        let files = try FileManager.default.contentsOfDirectory(atPath: Fixtures.directory.path).filter { $0.hasSuffix(".cif") }
        for file in files { structures.append((file, try CIFImport.read(try Fixtures.text(file)).structure)) }
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let split = fcc.map { AtomSite(z: 79, fractional: $0) }
            + [AtomSite(z: 79, fractional: [0.004, 0, 0]), AtomSite(z: 79, fractional: [0, 0.004, 0])]
        structures.append(("split gold", CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: split)))
        for (name, s) in structures {
            for op in s.operations {
                for site in s.sites {
                    let image = wrapFractional(op.apply(site.fractional))
                    let hit = s.sites.contains { other in
                        other.z == site.z && fractionalAxisDelta(other.fractional.x, image.x) < 5e-3
                            && fractionalAxisDelta(other.fractional.y, image.y) < 5e-3
                            && fractionalAxisDelta(other.fractional.z, image.z) < 5e-3
                    }
                    #expect(hit, "\(name): \(op.xyzString) sends a site nowhere")
                    if !hit { break }
                }
            }
        }
    }

    /// F6: a cell 1e-4 from cubic (c/a − 1) with cubic fractional coordinates is tetragonal (4/mmm, 16 operations, as
    /// spglib), so no orbit is split across two rings; at the reference implementations' 1e-3 metric tolerance it read m-3m with every ring
    /// named twice.
    @Test func testNearlyCubicCellIsTetragonalAndNoRingIsNamedTwice() throws {
        let cell = UnitCell.tetragonal(a: 3.905, c: 3.905 * (1 + 1e-4))
        let sites = [AtomSite(z: 38, fractional: [0, 0, 0]), AtomSite(z: 22, fractional: [0.5, 0.5, 0.5]),
                     AtomSite(z: 8, fractional: [0.5, 0.5, 0]), AtomSite(z: 8, fractional: [0.5, 0, 0.5]),
                     AtomSite(z: 8, fractional: [0, 0.5, 0.5])]
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        #expect(s.laueClass == .fourOverMmm)
        #expect(s.operations.count == 16)
        let names = try s.reflections(kMax: 1.2).families().flatMap { $0.name.components(separatedBy: " + ") }
        #expect(names.count == Set(names).count, "an orbit named in two rings")
    }
}
