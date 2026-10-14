//
//  CIFImportTests.swift
//  Planted truth for the importer (registration 5.3, 5.4): cell contents the files' formulas state, Laue classes by
//  the structure, absences by the operations, the apps' family counts and names, and the refusals. β″ first.
//

import XCTest
import simd
import DiffractoCore

final class CIFImportTests: XCTestCase {

    // MARK: β″ (monoclinic first)

    func testBetaDoublePrimeImportsAs22AtomsC2OverM() throws {
        let doc = try CIFImport.read(Fixtures.betaDoublePrimeCIF, fileBaseName: "beta_double_prime_Mg5Si6")
        let s = doc.structure
        XCTAssertEqual(doc.dataBlockName, "beta_double_prime_Mg5Si6")
        XCTAssertEqual(doc.displayName, "beta_double_prime_Mg5Si6")
        XCTAssertEqual(doc.declaredSpaceGroupNumber, 12)
        XCTAssertEqual(doc.declaredSpaceGroupSymbol, "C 2/m")
        XCTAssertTrue(doc.hadSymmetryOperations)
        XCTAssertEqual(s.sites.count, 22)
        XCTAssertEqual(s.formula, "Mg10 Si12")
        XCTAssertEqual(s.operations.count, 8)
        XCTAssertEqual(s.fileOperations?.count, 8)
        XCTAssertEqual(s.fileOperationsAreSymmetries, true)
        XCTAssertEqual(s.symmetrySource, .detectedAfterFileExpansion)
        XCTAssertEqual(s.laueClass, .twoOverM)
        XCTAssertEqual(s.crystalSystem, .monoclinic)
        XCTAssertEqual(s.spaceGroupNumber, 12)
        XCTAssertEqual(s.cell.metricFamily(), .monoclinic(uniqueAxis: .b))
        XCTAssertTrue(s.isUsable)
        // C-centring: h + k odd is absent.
        let list = try s.reflections(kMax: 1.2)
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: s.operations))
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(0, 1, 0), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(2, 0, 0), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(1, 1, 0), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: s.operations))
        for r in list.reflections { XCTAssertEqual(r.isAllowed, (r.hkl.h + r.hkl.k) % 2 == 0, "\(r.hkl)") }
    }

    /// The registered predictions: 1448 allowed reflections and 419 families at k_max 1.2 Å⁻¹ (macTEM's S2c record:
    /// gemmi's 1448; the apps' 1442 is after the |F| cut, slice 2's), and the first ten families by name, d and
    /// multiplicity. Mutants: the family tolerance at 1e-4 or 0; the multiplicity as the orbit count.
    func testBetaDoublePrimePredictedListAndFamilies() throws {
        let s = try CIFImport.read(Fixtures.betaDoublePrimeCIF).structure
        let list = try s.reflections(kMax: 1.2)
        XCTAssertEqual(list.allowed.count, 1448)
        let families = list.families()
        XCTAssertEqual(families.count, 419)
        let first: [(String, Double, Int)] = [
            ("{2 0 0}", 7.3113, 2), ("{0 0 1}", 6.5011, 2), ("{2 0 -1}", 5.6555, 2), ("{2 0 1}", 4.3246, 2),
            ("{1 1 0}", 3.9031, 4), ("{4 0 0}", 3.6557, 2), ("{4 0 -1}", 3.6207, 2), ("{1 1 -1}", 3.4554, 4),
            ("{2 0 -2}", 3.3124, 2), ("{0 0 2}", 3.2506, 2),
        ]
        for (i, (name, d, multiplicity)) in first.enumerated() {
            XCTAssertEqual(families[i].name, name, "family \(i)")
            XCTAssertEqual(families[i].d, d, accuracy: 6e-5, "family \(i)")
            XCTAssertEqual(families[i].multiplicity, multiplicity, "family \(i)")
            XCTAssertEqual(families[i].orbits.count, 1, "family \(i): no accidental coincidence in β″")
        }
        XCTAssertTrue(families.allSatisfy { $0.orbits.count == 1 }, "β″ has no multi-orbit family at 1.2")
        XCTAssertTrue(families.allSatisfy { $0.multiplicity == 2 || $0.multiplicity == 4 })
        XCTAssertEqual(families.filter { $0.multiplicity == 2 }.count + families.filter { $0.multiplicity == 4 }.count, 419)
        // The lattice rings (absent ones included) are more.
        XCTAssertGreaterThan(list.families(allowedOnly: false).count, 419)
        XCTAssertEqual(list.reflections.count, list.allowed.count + list.reflections.filter { !$0.isAllowed }.count)
    }

    /// The same 22 atoms as a P1 list: the operations are detected and everything else follows.
    func testBetaDoublePrimeAsP1DetectsTheSameSymmetry() throws {
        let fromCIF = try CIFImport.read(Fixtures.betaDoublePrimeCIF).structure
        let p1 = CrystalStructure.detectingSymmetry(cell: fromCIF.cell, sites: fromCIF.sites)
        XCTAssertEqual(p1.symmetrySource, .detected)
        XCTAssertEqual(p1.operations.count, 8)
        XCTAssertEqual(p1.laueClass, .twoOverM)
        for op in fromCIF.operations {
            XCTAssertTrue(p1.operations.contains { $0.isEquivalent(to: op, translationTolerance: 1e-6) }, "\(op)")
        }
        let a = try fromCIF.reflections(kMax: 1.2).families(), b = try p1.reflections(kMax: 1.2).families()
        XCTAssertEqual(a.count, b.count)
        XCTAssertEqual(a.map(\.name), b.map(\.name))
        XCTAssertEqual(a.map(\.multiplicity), b.map(\.multiplicity))
    }

    // MARK: The apps' fixtures, general expectations

    private func read(_ file: String) throws -> CIFDocument {
        try CIFImport.read(try Fixtures.text(file), fileBaseName: file)
    }

    func testGoldP1DetectsFm3mAndHasEightFamilies() throws {
        let s = try read("gold_P1.cif").structure
        XCTAssertEqual(s.sites.count, 4)
        XCTAssertEqual(s.symmetrySource, .detected)
        XCTAssertEqual(s.operations.count, 192)
        XCTAssertEqual(s.laueClass, .mThreeBarM)
        let families = try s.reflections(kMax: 1.2).families()
        XCTAssertEqual(families.map(\.name), ["{1 1 1}", "{2 0 0}", "{2 2 0}", "{3 1 1}", "{2 2 2}", "{4 0 0}", "{3 3 1}", "{4 2 0}"])
        XCTAssertEqual(families[0].multiplicity, 8)
        XCTAssertEqual(families[1].multiplicity, 6)
        XCTAssertEqual(families[3].multiplicity, 24)
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: s.operations), "mixed parity is absent")
    }

    /// One atom and the four centring translations: a legal loop that declares less than the atoms have (Laue -1 by
    /// the loop, m-3m by the cell it builds). The structure's operations are detected from the expanded cell
    /// (amendment 2); the loop is kept and is a subset of them.
    func testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m() throws {
        let doc = try read("gold_centring_ops.cif")
        let s = doc.structure
        XCTAssertEqual(s.sites.count, 4)
        XCTAssertTrue(doc.hadSymmetryOperations)
        XCTAssertEqual(s.symmetrySource, .detectedAfterFileExpansion)
        XCTAssertEqual(s.fileOperations?.count, 4)
        XCTAssertEqual(s.fileOperationsAreSymmetries, true)
        XCTAssertEqual(s.operations.count, 192)
        XCTAssertEqual(s.laueClass, .mThreeBarM)
        XCTAssertEqual(PointSymmetry.laueGroup(of: s.fileOperations ?? []).laueClass, .oneBar, "the loop alone is -1")
    }

    /// Silicon: Fd-3m; {2 0 0}, {4 2 0}, {6 0 0} absent by the operations; {2 2 2} and {4 4 2} allowed by them (their
    /// structure factor is zero by the basis, slice 2's), so 12 allowed families at 1.2 where the apps list 10.
    func testSiliconAbsencesAndFamilies() throws {
        let s = try read("silicon_Fd-3m.cif").structure
        XCTAssertEqual(s.sites.count, 8)
        XCTAssertEqual(s.laueClass, .mThreeBarM)
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(2, 0, 0), by: s.operations))
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(4, 2, 0), by: s.operations))
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(6, 0, 0), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(2, 2, 2), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(4, 4, 2), by: s.operations))
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(1, 1, 1), by: s.operations))
        let names = try s.reflections(kMax: 1.2).families().map(\.name)
        XCTAssertEqual(names.count, 12, "\(names)")
        XCTAssertTrue(names.contains("{2 2 2}"))
        XCTAssertTrue(names.contains("{4 4 2}"))
        XCTAssertTrue(names.contains("{5 1 1} + {3 3 3}"), "one ring of two orbits, named in descending order: \(names)")
        XCTAssertFalse(names.contains("{2 0 0}"))
        let ring = try XCTUnwrap(try s.reflections(kMax: 1.2).families().first { $0.name == "{5 1 1} + {3 3 3}" })
        XCTAssertEqual(ring.multiplicity, 32)
        XCTAssertEqual(ring.orbits, [MillerIndex(5, 1, 1), MillerIndex(3, 3, 3)])
    }

    func testMagnesiumFamiliesMatchMacTEM() throws {
        let s = try read("magnesium_P63mmc.cif").structure
        XCTAssertEqual(s.sites.count, 2)
        XCTAssertEqual(s.laueClass, .sixOverMmm)
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: s.operations), "the 6₃ screw: 000l with l odd")
        // 25 symmetry-allowed rings: macTEM's 24 plus {3 0 -3 1}, allowed by P6₃/mmc's operations and zero by the basis
        // (amendment 5; slice 2's |F| cut removes it, as Si's {2 2 2}).
        let families = try s.reflections(kMax: 1.2).families()
        XCTAssertEqual(families.count, 25)
        XCTAssertTrue(families.contains { $0.name == "{3 0 -3 1}" })
        func family(_ name: String) throws -> ReflectionFamily { try XCTUnwrap(families.first { $0.name == name }, name) }
        XCTAssertEqual(try family("{1 0 -1 0}").multiplicity, 6)
        XCTAssertEqual(try family("{0 0 0 2}").multiplicity, 2)
        XCTAssertEqual(try family("{1 0 -1 1}").multiplicity, 12)
        XCTAssertEqual(try family("{1 1 -2 0}").multiplicity, 6)
        XCTAssertEqual(try family("{1 1 -2 0}").orbits.count, 1)
    }

    func testRockSaltWithItsFull192Operations() throws {
        let doc = try read("rocksalt_Fm-3m_192ops.cif")
        XCTAssertEqual(doc.structure.sites.count, 8)
        XCTAssertEqual(doc.structure.spaceGroupNumber, 225)
        XCTAssertEqual(doc.structure.operations.count, 192)
        XCTAssertEqual(doc.structure.laueClass, .mThreeBarM)
    }

    /// Trigonal structures in hexagonal axes import with their own Laue class, -3m, not 6/mmm (the apps refused them).
    func testTrigonalInHexagonalAxesIsThreeBarM() throws {
        let quartz = try read("alpha_quartz_P3221.cif").structure
        XCTAssertEqual(quartz.sites.count, 9)
        XCTAssertEqual(quartz.laueClass, .threeBarM)
        XCTAssertEqual(quartz.crystalSystem, .trigonal)
        XCTAssertEqual(quartz.cell.metricFamily(), .hexagonal, "the metric alone says hexagonal")
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: quartz.operations), "3₂: 000l with l = 3n only")
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(0, 0, 3), by: quartz.operations))

        let calcite = try read("calcite_R-3c_hex.cif").structure
        XCTAssertEqual(calcite.sites.count, 30)
        XCTAssertEqual(calcite.laueClass, .threeBarM)
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: calcite.operations), "R-centring: −h + k + l = 3n")
        XCTAssertFalse(ReflectionList.isAllowed(MillerIndex(0, 0, 3), by: calcite.operations), "the c-glide: 000l with l = 6n")
        XCTAssertTrue(ReflectionList.isAllowed(MillerIndex(0, 0, 6), by: calcite.operations))

        let cdi2 = try read("cdi2_P-3m1.cif").structure
        XCTAssertEqual(cdi2.sites.count, 3)
        XCTAssertEqual(cdi2.laueClass, .threeBarM)
    }

    /// Fluorapatite P6₃/m is 6/m, not 6/mmm: the general orbit has 12 members, hk0 6 (12 and 24 under 6/mmm).
    func testFluorapatiteIsSixOverM() throws {
        let s = try read("fluorapatite_P63m.cif").structure
        XCTAssertEqual(s.sites.count, 42)
        XCTAssertEqual(s.laueClass, .sixOverM)
        XCTAssertEqual(s.laueGroup.orbit(of: MillerIndex(2, 1, 3)).count, 12)
        XCTAssertEqual(s.laueGroup.orbit(of: MillerIndex(2, 1, 0)).count, 6)
        XCTAssertNotEqual(s.laueGroup.canonical(MillerIndex(2, 1, 0)), s.laueGroup.canonical(MillerIndex(1, 2, 0)))
    }

    /// Pyrite Pa-3 is m-3, not m-3m: {2 1 0} and {1 2 0} are two orbits of 12.
    func testPyriteIsMThreeBar() throws {
        let s = try read("pyrite_Pa-3.cif").structure
        XCTAssertEqual(s.sites.count, 12)
        XCTAssertEqual(s.laueClass, .mThreeBar)
        XCTAssertEqual(s.laueGroup.orbit(of: MillerIndex(2, 1, 0)).count, 12)
        XCTAssertNotEqual(s.laueGroup.canonical(MillerIndex(2, 1, 0)), s.laueGroup.canonical(MillerIndex(1, 2, 0)))
    }

    func testOtherFixtures() throws {
        XCTAssertEqual(try read("ws2_2H.cif").structure.sites.count, 6)
        XCTAssertEqual(try read("ws2_2H.cif").structure.laueClass, .sixOverMmm)
        let accidental = try read("tetragonal_a_equals_c.cif").structure
        XCTAssertEqual(accidental.cell.metricFamily(lengthTolerance: 1e-3, angleToleranceDeg: 0.05), .cubic, "the metric is cubic")
        XCTAssertNotEqual(accidental.crystalSystem, .cubic, "the atoms are not")
        let hcp = try read("ideal_hcp_P1.cif").structure
        XCTAssertEqual(hcp.laueClass, .sixOverMmm)
        let hcpFamilies = try hcp.reflections(kMax: 1.2).families()
        XCTAssertEqual(hcpFamilies.count, 24, "23 after the |F| cut (slice 2)")
        let ring = try XCTUnwrap(hcpFamilies.first { $0.orbits.count == 2 })
        XCTAssertEqual(ring.name, "{3 0 -3 2} + {0 0 0 6}")
        XCTAssertEqual(ring.multiplicity, 14)
        // The seeded quartz and pyrite texts (incomplete loops, kept as what they are): 12 and 8 atoms, the loop a
        // subset of the detected operations.
        let seededQuartz = try read("alpha_quartz_P3221_as_seeded.cif").structure
        XCTAssertEqual(seededQuartz.sites.count, 12)
        XCTAssertEqual(seededQuartz.fileOperationsAreSymmetries, true)
        let seededPyrite = try read("pyrite_Pa-3_as_seeded.cif").structure
        XCTAssertEqual(seededPyrite.sites.count, 8)
        XCTAssertEqual(seededPyrite.laueClass, .mThreeBar)
    }

    // MARK: Refusals and precision (5.4)

    static let goldP1 = """
    data_gold
    _cell_length_a 4.0782
    _cell_length_b 4.0782
    _cell_length_c 4.0782
    _cell_angle_alpha 90
    _cell_angle_beta 90
    _cell_angle_gamma 90
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    Au1 Au 0 0 0
    Au2 Au 0.5 0.5 0
    Au3 Au 0.5 0 0.5
    Au4 Au 0 0.5 0.5

    """

    private func refusal(_ cif: String, _ file: StaticString = #filePath, line: UInt = #line) -> CIFImportError? {
        do { _ = try CIFImport.read(cif, fileBaseName: "t"); XCTFail("admitted", file: file, line: line); return nil }
        catch { return error }
    }

    func testRefusalsNameTheProblem() {
        XCTAssertEqual(refusal(Self.goldP1.replacingOccurrences(of: "_cell_length_c 4.0782\n", with: "")), .missingTag("_cell_length_c"))
        XCTAssertEqual(refusal(Self.goldP1.replacingOccurrences(of: "_cell_angle_beta 90", with: "_cell_angle_beta ninety")),
                       .unparsableNumber(tag: "_cell_angle_beta", value: "ninety"))
        XCTAssertEqual(refusal(Self.goldP1.components(separatedBy: "loop_")[0]), .noAtomSites)
        XCTAssertEqual(refusal(Self.goldP1.replacingOccurrences(of: "Au1 Au 0 0 0", with: "Xx1 Xx 0 0 0")),
                       .unknownElementSymbol("Xx", atomSiteLabel: "Xx1"))
        let badOp = Self.goldP1 + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'banana'\n"
        XCTAssertEqual(refusal(badOp), .unsupportedSymmetryOperation("banana"))
        XCTAssertEqual(refusal("data_gold\n_symmetry_space_group_name_H-M 'F m -3 m'\n" + Self.goldP1.dropFirst("data_gold\n".count)),
                       .missingSymmetryOperations(spaceGroup: "F m -3 m"))
        XCTAssertEqual(refusal("data_gold\n_symmetry_Int_Tables_number 225\n" + Self.goldP1.dropFirst("data_gold\n".count)
                               + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'x+1,y,z'\n"),
                       .missingSymmetryOperations(spaceGroup: "No. 225"))
        XCTAssertEqual(refusal(Self.goldP1 + "Au5 Au 0.25 0.25\n"), .truncatedLoop(firstTag: "_atom_site_label", columns: 5, values: 24))
        XCTAssertEqual(refusal(Self.goldP1 + "\n" + Self.goldP1.replacingOccurrences(of: "data_gold", with: "data_again")),
                       .multipleStructures(blockNames: ["gold", "again"]))
        XCTAssertEqual(refusal(Self.goldP1 + "\ndata_symops\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n"),
                       .symmetryOperationsOutsideStructure(blockNames: ["symops"]))
        let garbled = Self.goldP1 + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'1000000000000000000000x,y,z'\n"
        XCTAssertEqual(refusal(garbled), .unsupportedSymmetryOperation("1000000000000000000000x,y,z"))
    }

    /// Fm-3m listing only its 48 point operations is refused naming the gap; in a primitive (rhombohedral) cell it is
    /// not; duplicated rows do not pad the list.
    func testIncompleteOperationListsAreRefused() throws {
        let pointOps = PointSymmetry.holohedry(of: .cubic).rotations.map { SymmetryOperation(rotation: $0).xyzString }
        func loop(_ ops: [String]) -> String { "loop_\n_symmetry_equiv_pos_as_xyz\n" + ops.map { "'\($0)'" }.joined(separator: "\n") + "\n" }
        let header = "data_MgO\n_symmetry_space_group_name_H-M 'F m -3 m'\n_symmetry_Int_Tables_number 225\n"
        let cell = "_cell_length_a 4.21\n_cell_length_b 4.21\n_cell_length_c 4.21\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        let atoms = "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nMg Mg 0 0 0\nO O 0.5 0.5 0.5\n"
        XCTAssertEqual(refusal(header + cell + loop(pointOps) + atoms),
                       .incompleteSymmetryOperations(spaceGroupNumber: 225, declaredName: "F m -3 m", listed: 48, required: 192))
        XCTAssertEqual(refusal(header + cell + loop(pointOps + pointOps) + atoms),
                       .incompleteSymmetryOperations(spaceGroupNumber: 225, declaredName: "F m -3 m", listed: 48, required: 192))
        let primitiveCell = "_cell_length_a 2.977\n_cell_length_b 2.977\n_cell_length_c 2.977\n_cell_angle_alpha 60\n_cell_angle_beta 60\n_cell_angle_gamma 60\n"
        let primitive = try CIFImport.read(header + primitiveCell + loop(pointOps) + atoms)
        XCTAssertEqual(primitive.structure.sites.count, 2)
        XCTAssertEqual(CIFImport.spaceGroupOrder(225), 192)
        XCTAssertEqual(CIFImport.spaceGroupOrder(12), 8)
        XCTAssertEqual(CIFImport.spaceGroupOrder(194), 24)
        XCTAssertEqual(CIFImport.spaceGroupOrder(166), 36)
        XCTAssertNil(CIFImport.spaceGroupOrder(231))
        XCTAssertEqual(CIFImport.requiredOperationCount(spaceGroup: 166, alphaDeg: 60, betaDeg: 60, gammaDeg: 60), 12)
        XCTAssertEqual(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 105, gammaDeg: 90), 8)
        XCTAssertEqual(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 71, betaDeg: 80, gammaDeg: 75), 4)
        XCTAssertEqual(CIFImport.requiredOperationCount(spaceGroup: 63, alphaDeg: 90, betaDeg: 90, gammaDeg: 104), 8)
    }

    func testDeclarationsThatDoNotRefuse() throws {
        XCTAssertEqual(try CIFImport.read("data_gold\n_symmetry_space_group_name_H-M 'P 1'\n_space_group_IT_number 1\n"
                                          + Self.goldP1.dropFirst("data_gold\n".count)).structure.sites.count, 4)
        XCTAssertEqual(try CIFImport.read("data_gold\n_symmetry_space_group_name_H-M ?\n"
                                          + Self.goldP1.dropFirst("data_gold\n".count)).structure.sites.count, 4)
        let journal = try CIFImport.read("data_global\n_journal_name_full 'Some Journal'\n\n" + Self.goldP1)
        XCTAssertEqual(journal.structure.sites.count, 4)
        XCTAssertEqual(journal.dataBlockName, "global")
    }

    func testUncertaintySuffixAndDisplayName() throws {
        let doc = try CIFImport.read(Self.goldP1.replacingOccurrences(of: "_cell_length_a 4.0782", with: "_cell_length_a 4.0782(3)")
            .replacingOccurrences(of: "data_gold", with: "data_VESTA_phase_1"), fileBaseName: "my_gold")
        XCTAssertEqual(doc.structure.cell.a, 4.0782)
        XCTAssertEqual(doc.displayName, "my_gold")
        XCTAssertEqual(CIFImport.displayName(dataBlockName: "gold", fileBaseName: "f"), "gold")
        XCTAssertEqual(CIFImport.displayName(dataBlockName: "", fileBaseName: "f"), "f")
        XCTAssertEqual(CIFImport.displayName(dataBlockName: "VESTA_phase_12", fileBaseName: "f"), "f")
    }

    /// Written precision: `0.2500` is exact, `0.333` is a half-step of 1e-3, `0.33` 1e-2, nothing beyond 1e-2.
    func testCoordinateHalfStep() {
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.2500"), 0)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.5"), 0)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.333"), 1e-3)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.6666"), 1e-4)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.33"), 1e-2)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.3"), 1e-2)
        XCTAssertEqual(CIFImport.coordinateHalfStep("0.1234(5)"), 1e-4)
        XCTAssertEqual(CIFImport.coordinateHalfStep("x"), 0)
    }

    /// Magnesium written to three decimals (⅓ as 0.333) still detects 6/mmm; two distinct sites 0.3 Å apart are kept
    /// and refused by the close-contact rule; a 0.5/0.5 split site is not.
    func testPrecisionAndCloseContact() throws {
        let mg = """
        data_magnesium
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
        Mg1 0 0 0
        Mg2 0.667 0.333 0.5
        """
        let doc = try CIFImport.read(mg)
        XCTAssertEqual(doc.coarsestHalfStep, 1e-3)
        XCTAssertEqual(doc.structure.laueClass, .sixOverMmm)
        XCTAssertEqual(doc.structure.operations.count, 24)

        let close = Self.goldP1 + "Au5 Au 0.05 0 0\n"
        XCTAssertEqual(refusal(close)?.description.contains("atomic_sites_too_close") ?? true, false)
        if case .invalidStructure(let issues)? = refusal(close) {
            XCTAssertTrue(issues.contains { $0.contains("closer than any real interatomic distance") }, "\(issues)")
        } else {
            XCTFail("expected invalidStructure")
        }
        let split = Self.goldP1 + "_atom_site_occupancy\n" // not a valid edit of the loop; the split case is below
        _ = split
        let splitCIF = """
        data_split
        _cell_length_a 5
        _cell_length_b 5
        _cell_length_c 5
        _cell_angle_alpha 90
        _cell_angle_beta 90
        _cell_angle_gamma 90
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy
        Si1 0.49 0.5 0.5 0.5
        Si2 0.51 0.5 0.5 0.5
        """
        XCTAssertEqual(try CIFImport.read(splitCIF).structure.sites.count, 2)
    }
}

