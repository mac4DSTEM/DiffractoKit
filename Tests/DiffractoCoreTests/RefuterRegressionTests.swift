//
//  RefuterRegressionTests.swift
//  The Gate B refuter's findings on b4cfe17 (registration section 8, amendment of 2026-10-15), each as a regression
//  test on the refuter's own CIF (Fixtures/refuter_*.cif). Truth: pymatgen-core 2026.9.23 / spglib 2.8.0 on the
//  UNSHIFTED structure (absences and families do not depend on the origin), computed independently of the refuter's
//  harness: Mg (c = 5.2112) P6₃/mmc 24 operations, 282 of 324 allowed at 1.2 Å⁻¹; one Fe atom in a 5 Å cube Pm-3m
//  48 operations, 924 of 924; the 7× Po supercell (c = 23.45) 112 operations, 122 of 782 allowed at 0.9.
//

import XCTest
import simd
import DiffractoCore

final class RefuterRegressionTests: XCTestCase {

    private func read(_ name: String) throws -> CrystalStructure {
        try CIFImport.read(try Fixtures.text("refuter_\(name).cif"), fileBaseName: name).structure
    }

    /// Two structures that differ by an origin shift must agree on everything a pattern shows.
    private func assertSameDiffraction(_ a: CrystalStructure, _ b: CrystalStructure, kMax: Double, _ label: String,
                                       file: StaticString = #filePath, line: UInt = #line) throws {
        XCTAssertEqual(a.operations.count, b.operations.count, "\(label) operations", file: file, line: line)
        XCTAssertEqual(a.laueClass, b.laueClass, "\(label) Laue class", file: file, line: line)
        let la = try a.reflections(kMax: kMax), lb = try b.reflections(kMax: kMax)
        XCTAssertEqual(la.reflections.map(\.hkl), lb.reflections.map(\.hkl), "\(label) list", file: file, line: line)
        let wrong = zip(la.reflections, lb.reflections).filter { $0.isAllowed != $1.isAllowed }.map(\.0.hkl)
        XCTAssertTrue(wrong.isEmpty, "\(label): \(wrong.count) allowed flags differ, e.g. \(wrong.prefix(4))", file: file, line: line)
        XCTAssertEqual(la.families().map(\.name), lb.families().map(\.name), "\(label) families", file: file, line: line)
        XCTAssertEqual(la.families().map(\.multiplicity), lb.families().map(\.multiplicity), "\(label) multiplicities",
                       file: file, line: line)
    }

    /// F1(b): hcp Mg shifted by (0.1, 0.2, 0.3), 8 decimals. A translation of 0.20000001 failed the 1e-9 absence test
    /// for (2 −1 0): 262 allowed instead of 282, ring {1 1 −2 0} of 2 instead of 6. Mutant: the absence tolerance
    /// without the translation's uncertainty.
    func testShiftedMagnesiumKeepsItsAbsences() throws {
        let shifted = try read("F2a_Mg_shifted")
        XCTAssertEqual(shifted.operations.count, 24)
        XCTAssertEqual(shifted.laueClass, .sixOverMmm)
        let list = try shifted.reflections(kMax: 1.2)
        XCTAssertEqual(list.reflections.count, 324)
        XCTAssertEqual(list.allowed.count, 282)
        XCTAssertEqual(try XCTUnwrap(list.families().first { $0.name == "{1 1 -2 0}" }).multiplicity, 6)
        let unshifted = CrystalStructure.detectingSymmetry(
            cell: shifted.cell, sites: [AtomSite(z: 12, fractional: [1.0 / 3, 2.0 / 3, 0.25]),
                                        AtomSite(z: 12, fractional: [2.0 / 3, 1.0 / 3, 0.75])])
        try assertSameDiffraction(shifted, unshifted, kMax: 1.2, "Mg shifted")
    }

    /// F1(b) with a file loop: the same Mg with 24 exact rational operations at another origin. The file's exact
    /// translations are used for the operations it names (the 8-decimal coordinates alone also give 282 now).
    func testShiftedMagnesiumWithItsLoopUsesTheExactTranslations() throws {
        let s = try read("F2f_Mg_loop_0.37")
        XCTAssertEqual(s.fileOperations?.count, 24)
        XCTAssertEqual(s.operations.count, 24)
        XCTAssertEqual(s.fileOperationsAreSymmetries, true)
        XCTAssertEqual(try s.reflections(kMax: 1.2).allowed.count, 282)
    }

    /// F1(a): a lone atom at a general position is an origin shift of Pm-3m: 48 operations, every reflection allowed,
    /// ring {1 1 0} of 12. The fixed k/24 snap had moved translations by up to 5e-3: 43 operations, 80 wrong flags at 1.2.
    /// Mutant: the snap within the position tolerance instead of the translation's own uncertainty.
    func testOneAtomAtAGeneralPositionIsPmThreeBarM() throws {
        let s = try read("P1_cubic_1fe")
        XCTAssertEqual(s.operations.count, 48)
        XCTAssertEqual(s.laueClass, .mThreeBarM)
        // At 1.19 Å⁻¹ (h² + k² + l² ≤ 35; 36 sits on the float edge of 1.2): every lattice point, all allowed.
        let list = try s.reflections(kMax: 1.19)
        var points = 0
        for h in -6...6 { for k in -6...6 { for l in -6...6 where h * h + k * k + l * l <= 35 && (h, k, l) != (0, 0, 0) {
            points += 1
        } } }
        XCTAssertEqual(list.reflections.count, points)
        XCTAssertEqual(list.allowed.count, points)
        XCTAssertEqual(try XCTUnwrap(list.families().first { $0.name == "{1 1 0}" }).multiplicity, 12)
    }

    /// F1(a): rock salt shifted by (0.9172, 0.6292, 0.5141), exact 4-decimal coordinates: Fm-3m's 192 operations and
    /// the unshifted structure's absences and rings (the snap had given 160 operations, no Laue class, 24 wrong flags).
    func testShiftedRockSaltIsFmThreeBarM() throws {
        let shifted = try read("F2d_NaCl_shifted")
        XCTAssertEqual(shifted.operations.count, 192)
        XCTAssertEqual(shifted.laueClass, .mThreeBarM)
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let unshifted = CrystalStructure.detectingSymmetry(
            cell: shifted.cell,
            sites: fcc.map { AtomSite(z: 11, fractional: $0) }
                + fcc.map { AtomSite(z: 17, fractional: wrapFractional($0 + SIMD3(0.5, 0.5, 0.5))) })
        try assertSameDiffraction(shifted, unshifted, kMax: 1.2, "NaCl shifted")
    }

    /// F1(a): gold shifted by (0.2484, 0.3, 0.07): the inversion's translation is 0.4968, not 0.5.
    func testShiftedGoldKeepsItsOffGridTranslations() throws {
        let s = try read("F2c_Au_shifted_exactdecimals")
        XCTAssertEqual(s.operations.count, 192)
        XCTAssertEqual(s.laueClass, .mThreeBarM)
        let inversions = s.operations.filter { $0.rotation == .inversion }
        XCTAssertEqual(inversions.count, 4)
        XCTAssertTrue(inversions.contains { abs($0.translationMod1.x - 0.4968) < 1e-9 }, "\(inversions.map(\.xyzString))")
    }

    /// F1(a): quartz shifted: P3₂21's 6 operations, Laue -3m (the snap had given 5 and no class).
    func testShiftedQuartzIsThreeBarM() throws {
        let s = try read("F2e_quartz_shifted")
        XCTAssertEqual(s.operations.count, 6)
        XCTAssertEqual(s.laueClass, .threeBarM)
    }

    /// F1(b): a 7× supercell has translations of k/7, off the 1/24 grid, written to 8 decimals: 112 operations and
    /// 122 allowed of 782 at 0.9 Å⁻¹ (the 1e-9 test had flagged every l = ±7 sublattice reflection absent: 28).
    func testSevenfoldSupercellKeepsItsSublatticeReflections() throws {
        let s = try read("F3_super7")
        XCTAssertEqual(s.operations.count, 112)
        XCTAssertEqual(s.laueClass, .fourOverMmm, "a 1 × 1 × 7 cell is tetragonal by its metric; spglib's 221 is the primitive cubic cell")
        let list = try s.reflections(kMax: 0.9)
        XCTAssertEqual(list.reflections.count, 782)
        XCTAssertEqual(list.allowed.count, 122)
        XCTAssertTrue(list.allowed.allSatisfy { $0.hkl.l % 7 == 0 })
    }

    /// F2: two species sharing one site (Fe2+ 0.5 and Fe3+ 0.5 at the origin) under F-centring keep both rows: 8 Fe
    /// sites of 0.5 (Fe 4.0), 4 O; the expansion had merged the rows (Fe 2.0). Mutant: images of any row merge.
    func testCoSitedSpeciesAreNotMerged() throws {
        let s = try read("mixed_centred")
        XCTAssertEqual(s.sites.count, 12)
        XCTAssertEqual(s.sites.filter { $0.z == 26 }.reduce(0) { $0 + $1.occupancy }, 4.0, accuracy: 1e-12)
        XCTAssertEqual(s.sites.filter { $0.z == 8 }.count, 4)
        XCTAssertEqual(Set(s.sites.compactMap(\.label)), ["Fe2", "Fe3", "O1"])
    }

    /// F4: corundum in rhombohedral axes is -3c by its atoms (Laue -3m) and its rings are named in three indices; the
    /// Miller–Bravais form is for hexagonal axes only. Mutant: the four-index form by the Laue class alone.
    func testRhombohedralAxesKeepThreeIndexNames() throws {
        let s = try read("corundum_rh_p1")
        XCTAssertEqual(s.sites.count, 10)
        XCTAssertEqual(s.laueClass, .threeBarM)
        XCTAssertEqual(s.operations.count, 12)
        let list = try s.reflections(kMax: 1.2)
        XCTAssertFalse(list.usesMillerBravais)
        for family in list.families() {
            for orbit in family.name.components(separatedBy: " + ") {
                XCTAssertEqual(orbit.split(separator: " ").count, 3, family.name)
            }
        }
        XCTAssertTrue(try CIFImport.read(try Fixtures.text("calcite_R-3c_hex.cif")).structure.reflections(kMax: 1.2).usesMillerBravais)
    }

    /// F5: a finite but huge cell is refused, never trapped or left to hang; a large list is refused before it is built.
    func testHugeCellsAreRefusedNotTrapped() {
        let cases: [(Double, Double)] = [(1e30, 1.2), (1e18, 1.2), (1e10, 1.2), (100, 4)]
        for (a, k) in cases {
            let s = CrystalStructure.detectingSymmetry(cell: .cubic(a: a), sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
            XCTAssertThrowsError(try s.reflections(kMax: k), "a = \(a)") { error in
                guard case .cellTooLargeOrInvalid = error as? ReflectionListError else {
                    return XCTFail("\(error)")
                }
            }
        }
        let thin = CrystalStructure.detectingSymmetry(cell: .orthorhombic(a: 100, b: 3, c: 3),
                                                     sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
        XCTAssertNoThrow(try thin.reflections(kMax: 4), "a long thin cell is a small list")
    }

    /// F5: an operation whose rotation part is singular or a dilation is not a symmetry and is refused.
    func testNonInvertibleOperationsAreRefused() {
        for raw in ["x,x,z", "2x,y,z", "x+1000y,y+1000x,z"] {
            XCTAssertThrowsError(try SymmetryOperation(xyz: raw)) {
                XCTAssertEqual($0 as? SymmetryOperationError, .notInvertible(raw))
            }
        }
        XCTAssertNoThrow(try SymmetryOperation(xyz: "x-y,x,z"))
    }

    /// F6: a cell 1e-4 from cubic (c/a − 1) with cubic fractional coordinates is tetragonal (4/mmm, 16 operations, as
    /// spglib), so no orbit is split across two rings; at the apps' 1e-3 metric tolerance it read m-3m with every ring
    /// named twice.
    func testNearlyCubicCellIsTetragonalAndNoRingIsNamedTwice() throws {
        let cell = UnitCell.tetragonal(a: 3.905, c: 3.905 * (1 + 1e-4))
        let sites = [AtomSite(z: 38, fractional: [0, 0, 0]), AtomSite(z: 22, fractional: [0.5, 0.5, 0.5]),
                     AtomSite(z: 8, fractional: [0.5, 0.5, 0]), AtomSite(z: 8, fractional: [0.5, 0, 0.5]),
                     AtomSite(z: 8, fractional: [0, 0.5, 0.5])]
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        XCTAssertEqual(s.laueClass, .fourOverMmm)
        XCTAssertEqual(s.operations.count, 16)
        let names = try s.reflections(kMax: 1.2).families().flatMap { $0.name.components(separatedBy: " + ") }
        XCTAssertEqual(names.count, Set(names).count, "an orbit named in two rings")
    }
}
