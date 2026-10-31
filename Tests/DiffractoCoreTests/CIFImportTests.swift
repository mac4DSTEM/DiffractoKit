//
//  CIFImportTests.swift
//  Planted truth for the importer: cell contents the files' formulas state, Laue classes by
//  the structure, absences by the operations, the reference implementations' family counts and names, and the refusals. β″ first.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct CIFImportTests {

    // MARK: β″ (monoclinic first)

    @Test func testBetaDoublePrimeImportsAs22AtomsC2OverM() throws {
        let doc = try CIFImport.read(Fixtures.betaDoublePrimeCIF, fileBaseName: "beta_double_prime_Mg5Si6")
        let s = doc.structure
        #expect(doc.dataBlockName == "beta_double_prime_Mg5Si6")
        #expect(doc.displayName == "beta_double_prime_Mg5Si6")
        #expect(doc.declaredSpaceGroupNumber == 12)
        #expect(doc.declaredSpaceGroupSymbol == "C 2/m")
        #expect(doc.hadSymmetryOperations)
        #expect(s.sites.count == 22)
        #expect(s.formula == "Mg10 Si12")
        #expect(s.operations.count == 8)
        #expect(s.fileOperations?.count == 8)
        #expect(s.fileOperationsAreSymmetries == true)
        #expect(s.symmetrySource == .detectedAfterFileExpansion)
        #expect(s.laueClass == .twoOverM)
        #expect(s.crystalSystem == .monoclinic)
        #expect(s.spaceGroupNumber == 12)
        #expect(s.cell.metricFamily() == .monoclinic(uniqueAxis: .b))
        #expect(s.isUsable)
        // C-centring: h + k odd is absent.
        let list = try s.reflections(kMax: 1.2)
        #expect(!ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: s.operations))
        #expect(!ReflectionList.isAllowed(MillerIndex(0, 1, 0), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(2, 0, 0), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(1, 1, 0), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: s.operations))
        for r in list.reflections { #expect(r.isAllowed == ((r.hkl.h + r.hkl.k) % 2 == 0), "\(r.hkl)") }
    }

    /// The predictions written before the code: 1448 allowed reflections and 419 families at k_max 1.2 Å⁻¹ (the reference implementation's
    /// record: gemmi's 1448; the reference implementations' 1442 is after DiffractoScattering's |F| cut), and the first ten families by name, d and
    /// multiplicity. Mutants: the family tolerance at 1e-4 or 0; the multiplicity as the orbit count.
    @Test func testBetaDoublePrimePredictedListAndFamilies() throws {
        let s = try CIFImport.read(Fixtures.betaDoublePrimeCIF).structure
        let list = try s.reflections(kMax: 1.2)
        #expect(list.allowed.count == 1448)
        let families = list.families()
        #expect(families.count == 419)
        let first: [(String, Double, Int)] = [
            ("{2 0 0}", 7.3113, 2), ("{0 0 1}", 6.5011, 2), ("{2 0 -1}", 5.6555, 2), ("{2 0 1}", 4.3246, 2),
            ("{1 1 0}", 3.9031, 4), ("{4 0 0}", 3.6557, 2), ("{4 0 -1}", 3.6207, 2), ("{1 1 -1}", 3.4554, 4),
            ("{2 0 -2}", 3.3124, 2), ("{0 0 2}", 3.2506, 2),
        ]
        for (i, (name, d, multiplicity)) in first.enumerated() {
            #expect(families[i].name == name, "family \(i)")
            #expect(abs(families[i].d - d) <= 6e-5, "family \(i)")
            #expect(families[i].multiplicity == multiplicity, "family \(i)")
            #expect(families[i].orbits.count == 1, "family \(i): no accidental coincidence in β″")
        }
        #expect(families.allSatisfy { $0.orbits.count == 1 }, "β″ has no multi-orbit family at 1.2")
        let multiplicityIsTwoOrFour: Bool = families.allSatisfy { (f: ReflectionFamily) -> Bool in f.multiplicity == 2 || f.multiplicity == 4 }
        #expect(multiplicityIsTwoOrFour)
        let twos: Int = families.filter { (f: ReflectionFamily) -> Bool in f.multiplicity == 2 }.count
        let fours: Int = families.filter { (f: ReflectionFamily) -> Bool in f.multiplicity == 4 }.count
        let twosAndFours: Int = twos + fours
        #expect(twosAndFours == 419)
        // The lattice rings (absent ones included) are more.
        #expect(list.families(allowedOnly: false).count > 419)
        #expect(list.reflections.count == list.allowed.count + list.reflections.filter { !$0.isAllowed }.count)
    }

    /// The same 22 atoms as a P1 list: the operations are detected and everything else follows.
    @Test func testBetaDoublePrimeAsP1DetectsTheSameSymmetry() throws {
        let fromCIF = try CIFImport.read(Fixtures.betaDoublePrimeCIF).structure
        let p1 = CrystalStructure.detectingSymmetry(cell: fromCIF.cell, sites: fromCIF.sites)
        #expect(p1.symmetrySource == .detected)
        #expect(p1.operations.count == 8)
        #expect(p1.laueClass == .twoOverM)
        for op in fromCIF.operations {
            #expect(p1.operations.contains { $0.isEquivalent(to: op, translationTolerance: 1e-6) }, "\(op)")
        }
        let a = try fromCIF.reflections(kMax: 1.2).families(), b = try p1.reflections(kMax: 1.2).families()
        #expect(a.count == b.count)
        #expect(a.map(\.name) == b.map(\.name))
        #expect(a.map(\.multiplicity) == b.map(\.multiplicity))
    }

    // MARK: The reference implementations' fixtures, general expectations

    private func read(_ file: String) throws -> CIFDocument {
        try CIFImport.read(try Fixtures.text(file), fileBaseName: file)
    }

    @Test func testGoldP1DetectsFm3mAndHasEightFamilies() throws {
        let s = try read("gold_P1.cif").structure
        #expect(s.sites.count == 4)
        #expect(s.symmetrySource == .detected)
        #expect(s.operations.count == 192)
        #expect(s.laueClass == .mThreeBarM)
        let families = try s.reflections(kMax: 1.2).families()
        #expect(families.map(\.name) == ["{1 1 1}", "{2 0 0}", "{2 2 0}", "{3 1 1}", "{2 2 2}", "{4 0 0}", "{3 3 1}", "{4 2 0}"])
        #expect(families[0].multiplicity == 8)
        #expect(families[1].multiplicity == 6)
        #expect(families[3].multiplicity == 24)
        #expect(!ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: s.operations), "mixed parity is absent")
    }

    /// One atom and the four centring translations: a legal loop that declares less than the atoms have (Laue -1 by
    /// the loop, m-3m by the cell it builds). The structure's operations are detected from the expanded cell; the loop is kept and is a subset of them.
    @Test func testGoldWithCentringOperationsExpandsToFourSitesAndDetectsFm3m() throws {
        let doc = try read("gold_centring_ops.cif")
        let s = doc.structure
        #expect(s.sites.count == 4)
        #expect(doc.hadSymmetryOperations)
        #expect(s.symmetrySource == .detectedAfterFileExpansion)
        #expect(s.fileOperations?.count == 4)
        #expect(s.fileOperationsAreSymmetries == true)
        #expect(s.operations.count == 192)
        #expect(s.laueClass == .mThreeBarM)
        #expect(PointSymmetry.laueGroup(of: s.fileOperations ?? []).laueClass == .oneBar, "the loop alone is -1")
    }

    /// Silicon: Fd-3m; {2 0 0}, {4 2 0}, {6 0 0} absent by the operations; {2 2 2} and {4 4 2} allowed by them (their
    /// structure factor is zero by the basis, removed by DiffractoScattering's cut), so 12 allowed families at 1.2 where the reference implementations list 10.
    @Test func testSiliconAbsencesAndFamilies() throws {
        let s = try read("silicon_Fd-3m.cif").structure
        #expect(s.sites.count == 8)
        #expect(s.laueClass == .mThreeBarM)
        #expect(!ReflectionList.isAllowed(MillerIndex(2, 0, 0), by: s.operations))
        #expect(!ReflectionList.isAllowed(MillerIndex(4, 2, 0), by: s.operations))
        #expect(!ReflectionList.isAllowed(MillerIndex(6, 0, 0), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(2, 2, 2), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(4, 4, 2), by: s.operations))
        #expect(ReflectionList.isAllowed(MillerIndex(1, 1, 1), by: s.operations))
        let names = try s.reflections(kMax: 1.2).families().map(\.name)
        #expect(names.count == 12, "\(names)")
        #expect(names.contains("{2 2 2}"))
        #expect(names.contains("{4 4 2}"))
        #expect(names.contains("{5 1 1} + {3 3 3}"), "one ring of two orbits, named in descending order: \(names)")
        #expect(!names.contains("{2 0 0}"))
        let ring = try #require(try s.reflections(kMax: 1.2).families().first { $0.name == "{5 1 1} + {3 3 3}" })
        #expect(ring.multiplicity == 32)
        #expect(ring.orbits == [MillerIndex(5, 1, 1), MillerIndex(3, 3, 3)])
    }

    @Test func testMagnesiumFamiliesMatchReference() throws {
        let s = try read("magnesium_P63mmc.cif").structure
        #expect(s.sites.count == 2)
        #expect(s.laueClass == .sixOverMmm)
        #expect(!ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: s.operations), "the 6₃ screw: 000l with l odd")
        // 25 symmetry-allowed rings: the reference implementation's 24 plus {3 0 -3 1}, allowed by P6₃/mmc's operations and zero by the basis
        // (DiffractoScattering's |F| cut removes it, as Si's {2 2 2}).
        let families = try s.reflections(kMax: 1.2).families()
        #expect(families.count == 25)
        #expect(families.contains { $0.name == "{3 0 -3 1}" })
        func family(_ name: String) throws -> ReflectionFamily { try #require(families.first { $0.name == name }, "\(name)") }
        #expect(try family("{1 0 -1 0}").multiplicity == 6)
        #expect(try family("{0 0 0 2}").multiplicity == 2)
        #expect(try family("{1 0 -1 1}").multiplicity == 12)
        #expect(try family("{1 1 -2 0}").multiplicity == 6)
        #expect(try family("{1 1 -2 0}").orbits.count == 1)
    }

    @Test func testRockSaltWithItsFull192Operations() throws {
        let doc = try read("rocksalt_Fm-3m_192ops.cif")
        #expect(doc.structure.sites.count == 8)
        #expect(doc.structure.spaceGroupNumber == 225)
        #expect(doc.structure.operations.count == 192)
        #expect(doc.structure.laueClass == .mThreeBarM)
    }

    /// Trigonal structures in hexagonal axes import with their own Laue class, -3m, not 6/mmm (the reference implementations refused them).
    @Test func testTrigonalInHexagonalAxesIsThreeBarM() throws {
        let quartz = try read("alpha_quartz_P3221.cif").structure
        #expect(quartz.sites.count == 9)
        #expect(quartz.laueClass == .threeBarM)
        #expect(quartz.crystalSystem == .trigonal)
        #expect(quartz.cell.metricFamily() == .hexagonal, "the metric alone says hexagonal")
        #expect(!ReflectionList.isAllowed(MillerIndex(0, 0, 1), by: quartz.operations), "3₂: 000l with l = 3n only")
        #expect(ReflectionList.isAllowed(MillerIndex(0, 0, 3), by: quartz.operations))

        let calcite = try read("calcite_R-3c_hex.cif").structure
        #expect(calcite.sites.count == 30)
        #expect(calcite.laueClass == .threeBarM)
        #expect(!ReflectionList.isAllowed(MillerIndex(1, 0, 0), by: calcite.operations), "R-centring: −h + k + l = 3n")
        #expect(!ReflectionList.isAllowed(MillerIndex(0, 0, 3), by: calcite.operations), "the c-glide: 000l with l = 6n")
        #expect(ReflectionList.isAllowed(MillerIndex(0, 0, 6), by: calcite.operations))

        let cdi2 = try read("cdi2_P-3m1.cif").structure
        #expect(cdi2.sites.count == 3)
        #expect(cdi2.laueClass == .threeBarM)
    }

    /// Fluorapatite P6₃/m is 6/m, not 6/mmm: the general orbit has 12 members, hk0 6 (12 and 24 under 6/mmm).
    @Test func testFluorapatiteIsSixOverM() throws {
        let s = try read("fluorapatite_P63m.cif").structure
        #expect(s.sites.count == 42)
        #expect(s.laueClass == .sixOverM)
        #expect(s.laueGroup.orbit(of: MillerIndex(2, 1, 3)).count == 12)
        #expect(s.laueGroup.orbit(of: MillerIndex(2, 1, 0)).count == 6)
        #expect(s.laueGroup.canonical(MillerIndex(2, 1, 0)) != s.laueGroup.canonical(MillerIndex(1, 2, 0)))
    }

    /// Pyrite Pa-3 is m-3, not m-3m: {2 1 0} and {1 2 0} are two orbits of 12.
    @Test func testPyriteIsMThreeBar() throws {
        let s = try read("pyrite_Pa-3.cif").structure
        #expect(s.sites.count == 12)
        #expect(s.laueClass == .mThreeBar)
        #expect(s.laueGroup.orbit(of: MillerIndex(2, 1, 0)).count == 12)
        #expect(s.laueGroup.canonical(MillerIndex(2, 1, 0)) != s.laueGroup.canonical(MillerIndex(1, 2, 0)))
    }

    @Test func testOtherFixtures() throws {
        #expect(try read("ws2_2H.cif").structure.sites.count == 6)
        #expect(try read("ws2_2H.cif").structure.laueClass == .sixOverMmm)
        let accidental = try read("tetragonal_a_equals_c.cif").structure
        #expect(accidental.cell.metricFamily(lengthTolerance: 1e-3, angleToleranceDeg: 0.05) == .cubic, "the metric is cubic")
        #expect(accidental.crystalSystem != .cubic, "the atoms are not")
        let hcp = try read("ideal_hcp_P1.cif").structure
        #expect(hcp.laueClass == .sixOverMmm)
        let hcpFamilies = try hcp.reflections(kMax: 1.2).families()
        #expect(hcpFamilies.count == 24, "23 after the |F| cut")
        let ring = try #require(hcpFamilies.first { $0.orbits.count == 2 })
        #expect(ring.name == "{3 0 -3 2} + {0 0 0 6}")
        #expect(ring.multiplicity == 14)
        // The seeded quartz and pyrite texts (incomplete loops, kept as what they are): 12 and 8 atoms, the loop a
        // subset of the detected operations.
        let seededQuartz = try read("alpha_quartz_P3221_as_seeded.cif").structure
        #expect(seededQuartz.sites.count == 12)
        #expect(seededQuartz.fileOperationsAreSymmetries == true)
        let seededPyrite = try read("pyrite_Pa-3_as_seeded.cif").structure
        #expect(seededPyrite.sites.count == 8)
        #expect(seededPyrite.laueClass == .mThreeBar)
    }

    // MARK: Refusals and precision

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
        do { _ = try CIFImport.read(cif, fileBaseName: "t"); Issue.record("admitted"); return nil }
        catch { return error }
    }

    @Test func testRefusalsNameTheProblem() {
        #expect(refusal(Self.goldP1.replacingOccurrences(of: "_cell_length_c 4.0782\n", with: "")) == .missingTag("_cell_length_c"))
        #expect(refusal(Self.goldP1.replacingOccurrences(of: "_cell_angle_beta 90", with: "_cell_angle_beta ninety")) == .unparsableNumber(tag: "_cell_angle_beta", value: "ninety"))
        #expect(refusal(Self.goldP1.components(separatedBy: "loop_")[0]) == .noAtomSites)
        #expect(refusal(Self.goldP1.replacingOccurrences(of: "Au1 Au 0 0 0", with: "Xx1 Xx 0 0 0")) == .unknownElementSymbol("Xx", atomSiteLabel: "Xx1"))
        let badOp = Self.goldP1 + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'banana'\n"
        #expect(refusal(badOp) == .unsupportedSymmetryOperation("banana"))
        #expect(refusal("data_gold\n_symmetry_space_group_name_H-M 'F m -3 m'\n" + Self.goldP1.dropFirst("data_gold\n".count)) == .missingSymmetryOperations(spaceGroup: "F m -3 m"))
        #expect(refusal("data_gold\n_symmetry_Int_Tables_number 225\n" + Self.goldP1.dropFirst("data_gold\n".count)
                               + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'x+1,y,z'\n") == .missingSymmetryOperations(spaceGroup: "No. 225"))
        #expect(refusal(Self.goldP1 + "Au5 Au 0.25 0.25\n") == .truncatedLoop(firstTag: "_atom_site_label", columns: 5, values: 24))
        #expect(refusal(Self.goldP1 + "\n" + Self.goldP1.replacingOccurrences(of: "data_gold", with: "data_again")) == .multipleStructures(blockNames: ["gold", "again"]))
        #expect(refusal(Self.goldP1 + "\ndata_symops\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n") == .symmetryOperationsOutsideStructure(blockNames: ["symops"]))
        let garbled = Self.goldP1 + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'1000000000000000000000x,y,z'\n"
        #expect(refusal(garbled) == .unsupportedSymmetryOperation("1000000000000000000000x,y,z"))
    }

    /// Fm-3m listing only its 48 point operations is refused naming the gap; duplicated rows do not pad the list. In a
    /// primitive (rhombohedral) cell 48 is enough, but these conventional-basis operations do not map that cell onto
    /// itself and are refused by name (`OperationListCheckTests` has the right primitive list).
    @Test func testIncompleteOperationListsAreRefused() throws {
        let pointOps = PointSymmetry.holohedry(of: .cubic).rotations.map { SymmetryOperation(rotation: $0).xyzString }
        func loop(_ ops: [String]) -> String { "loop_\n_symmetry_equiv_pos_as_xyz\n" + ops.map { "'\($0)'" }.joined(separator: "\n") + "\n" }
        let header = "data_MgO\n_symmetry_space_group_name_H-M 'F m -3 m'\n_symmetry_Int_Tables_number 225\n"
        let cell = "_cell_length_a 4.21\n_cell_length_b 4.21\n_cell_length_c 4.21\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        let atoms = "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nMg Mg 0 0 0\nO O 0.5 0.5 0.5\n"
        #expect(refusal(header + cell + loop(pointOps) + atoms) == .incompleteSymmetryOperations(spaceGroupNumber: 225, declaredName: "F m -3 m", listed: 48, required: 192))
        #expect(refusal(header + cell + loop(pointOps + pointOps) + atoms) == .incompleteSymmetryOperations(spaceGroupNumber: 225, declaredName: "F m -3 m", listed: 48, required: 192))
        let primitiveCell = "_cell_length_a 2.977\n_cell_length_b 2.977\n_cell_length_c 2.977\n_cell_angle_alpha 60\n_cell_angle_beta 60\n_cell_angle_gamma 60\n"
        guard case .inconsistentSymmetryOperations(.notACellSymmetry)? = refusal(header + primitiveCell + loop(pointOps) + atoms)
        else { Issue.record("conventional-basis operations admitted in a primitive cell"); return }
        #expect(CIFImport.spaceGroupOrder(225) == 192)
        #expect(CIFImport.spaceGroupOrder(12) == 8)
        #expect(CIFImport.spaceGroupOrder(194) == 24)
        #expect(CIFImport.spaceGroupOrder(166) == 36)
        #expect(CIFImport.spaceGroupOrder(231) == nil)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 166, alphaDeg: 60, betaDeg: 60, gammaDeg: 60) == 12)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 105, gammaDeg: 90) == 8)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 71, betaDeg: 80, gammaDeg: 75) == 4)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 63, alphaDeg: 90, betaDeg: 90, gammaDeg: 104) == 8)
    }

    @Test func testDeclarationsThatDoNotRefuse() throws {
        #expect(try CIFImport.read("data_gold\n_symmetry_space_group_name_H-M 'P 1'\n_space_group_IT_number 1\n"
                                          + Self.goldP1.dropFirst("data_gold\n".count)).structure.sites.count == 4)
        #expect(try CIFImport.read("data_gold\n_symmetry_space_group_name_H-M ?\n"
                                          + Self.goldP1.dropFirst("data_gold\n".count)).structure.sites.count == 4)
        let journal = try CIFImport.read("data_global\n_journal_name_full 'Some Journal'\n\n" + Self.goldP1)
        #expect(journal.structure.sites.count == 4)
        #expect(journal.dataBlockName == "global")
    }

    @Test func testUncertaintySuffixAndDisplayName() throws {
        let doc = try CIFImport.read(Self.goldP1.replacingOccurrences(of: "_cell_length_a 4.0782", with: "_cell_length_a 4.0782(3)")
            .replacingOccurrences(of: "data_gold", with: "data_VESTA_phase_1"), fileBaseName: "my_gold")
        #expect(doc.structure.cell.a == 4.0782)
        #expect(doc.displayName == "my_gold")
        #expect(CIFImport.displayName(dataBlockName: "gold", fileBaseName: "f") == "gold")
        #expect(CIFImport.displayName(dataBlockName: "", fileBaseName: "f") == "f")
        #expect(CIFImport.displayName(dataBlockName: "VESTA_phase_12", fileBaseName: "f") == "f")
    }

    /// Written precision: `0.2500` is exact, `0.333` is a half-step of 1e-3, `0.33` 1e-2, nothing beyond 1e-2.
    @Test func testCoordinateHalfStep() {
        #expect(CIFImport.coordinateHalfStep("0.2500") == 0)
        #expect(CIFImport.coordinateHalfStep("0.5") == 0)
        #expect(CIFImport.coordinateHalfStep("0.333") == 1e-3)
        #expect(CIFImport.coordinateHalfStep("0.6666") == 1e-4)
        #expect(CIFImport.coordinateHalfStep("0.33") == 1e-2)
        #expect(CIFImport.coordinateHalfStep("0.3") == 1e-2)
        #expect(CIFImport.coordinateHalfStep("0.1234(5)") == 1e-4)
        #expect(CIFImport.coordinateHalfStep("x") == 0)
    }

    /// Magnesium written to three decimals (⅓ as 0.333) still detects 6/mmm; two distinct sites 0.3 Å apart are kept
    /// and refused by the close-contact rule; a 0.5/0.5 split site is not.
    @Test func testPrecisionAndCloseContact() throws {
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
        #expect(doc.coarsestHalfStep == 1e-3)
        #expect(doc.structure.laueClass == .sixOverMmm)
        #expect(doc.structure.operations.count == 24)

        let close = Self.goldP1 + "Au5 Au 0.05 0 0\n"
        #expect(refusal(close)?.description.contains("atomic_sites_too_close") ?? true == false)
        if case .invalidStructure(let issues)? = refusal(close) {
            #expect(issues.contains { $0.contains("closer than any real interatomic distance") }, "\(issues)")
        } else {
            Issue.record("expected invalidStructure")
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
        #expect(try CIFImport.read(splitCIF).structure.sites.count == 2)
    }
}

// MARK: Completeness of the operation list (0.9.0: symbol-only declarations, primitive C-monoclinic with beta = 90)

struct CIFSymmetryCompletenessTests {
    private struct Pins: Decodable {
        struct Case: Decodable { let name: String; let text: String; let atoms: Int; let operations: Int }
        let symbols: [[JSONValue]]
        let cases: [Case]
        let conventionalOpsListed: Int
    }
    private enum JSONValue: Decodable {
        case int(Int), string(String)
        init(from decoder: Decoder) throws {
            let c = try decoder.singleValueContainer()
            if let i = try? c.decode(Int.self) { self = .int(i) } else { self = .string(try c.decode(String.self)) }
        }
    }
    private func pins() throws -> Pins {
        try JSONDecoder().decode(Pins.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("cif-symmetry-pins.json")))
    }

    /// Every spelling spglib 2.8.0 lists for the 230 groups (short, full, unspaced, every setting) reads to the order and the
    /// centring of the kit's own table: exact, no tolerance. RED UNDER: no symbol parser (the symbol-only file unchecked).
    @Test func everyHermannMauguinSpellingReadsToItsGroupsOrder() throws {
        var checked = 0
        for row in try pins().symbols {
            guard case .int(let number) = row[0], case .string(let symbol) = row[1] else { continue }
            let facts = try #require(CIFImport.hermannMauguinFacts(symbol), "unreadable symbol \(symbol)")
            #expect(facts.order * facts.centring == CIFImport.spaceGroupOrder(number), "\(symbol) (No. \(number))")
            checked += 1
        }
        #expect(checked == 1048)
        #expect(CIFImport.hermannMauguinFacts("P 1") != nil)
        #expect(CIFImport.hermannMauguinFacts("not a group") == nil)
        #expect(CIFImport.hermannMauguinFacts("") == nil)
        #expect(CIFImport.hermannMauguinFacts("P 7") == nil)
    }

    /// A group named by symbol alone is held to the same completeness as one named by number: rock salt listing its 48
    /// point operations is refused, in its conventional cell, naming the symbol. RED UNDER: the check ran only with an IT number.
    @Test func aGroupNamedBySymbolAloneIsChecked() throws {
        let pointOps = PointSymmetry.holohedry(of: .cubic).rotations.map { SymmetryOperation(rotation: $0).xyzString }
        let loop = "loop_\n_symmetry_equiv_pos_as_xyz\n" + pointOps.map { "'\($0)'" }.joined(separator: "\n") + "\n"
        let header = "data_MgO\n_symmetry_space_group_name_H-M 'F m -3 m'\n"
        let cell = "_cell_length_a 4.21\n_cell_length_b 4.21\n_cell_length_c 4.21\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        let atoms = "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nMg Mg 0 0 0\nO O 0.5 0.5 0.5\n"
        do {
            _ = try CIFImport.read(header + cell + loop + atoms)
            Issue.record("a 48-operation Fm-3m list was admitted")
        } catch {
            #expect(error == .incompleteSymmetryOperationsForSymbol(symbol: "F m -3 m", listed: 48, required: 192))
        }
        // A symbol the parser cannot read stays unchecked for its count (as before). In the primitive rhombohedral cell 48
        // passes the count (192 / 4); these conventional-basis operations are then refused as no symmetry of that cell.
        let primitiveCell = "_cell_length_a 2.977\n_cell_length_b 2.977\n_cell_length_c 2.977\n_cell_angle_alpha 60\n_cell_angle_beta 60\n_cell_angle_gamma 60\n"
        do {
            _ = try CIFImport.read(header + primitiveCell + loop + atoms)
            Issue.record("conventional-basis operations admitted in a primitive cell")
        } catch {
            guard case .inconsistentSymmetryOperations(.notACellSymmetry) = error else { Issue.record("\(error)"); return }
        }
        let unreadable = header.replacingOccurrences(of: "F m -3 m", with: "F banana")
        #expect(try CIFImport.read(unreadable + cell + loop + atoms).structure.sites.count > 0)
    }

    /// A primitive C-centred monoclinic cell with beta = 90 has two right angles like the conventional cell; the kit judged it
    /// conventional and refused its correct 4-operation list. The rotations settle it (not diagonal in a primitive basis).
    /// Pins: pymatgen 2026.9.24 primitive structure, spglib 2.8.0 operations in that cell. RED UNDER: angles alone.
    @Test func aPrimitiveCMonoclinicCellWithBeta90IsNotRefused() throws {
        for c in try pins().cases {
            let doc = try CIFImport.read(c.text, fileBaseName: "t")
            #expect(doc.structure.sites.count == c.atoms, "\(c.name)")
            #expect(c.operations == 4, "\(c.name)")
        }
        // The same ops cut to two are refused: the primitive list needs four.
        let cut = try pins().cases[0].text.split(separator: "\n", omittingEmptySubsequences: false)
        let ops = cut.enumerated().filter { $0.element.hasPrefix("'") }.map(\.offset)
        var kept = cut.map(String.init)
        for k in ops.dropFirst(2).reversed() { kept.remove(at: k) }
        do {
            _ = try CIFImport.read(kept.joined(separator: "\n"))
            Issue.record("a half list was admitted")
        } catch {
            if case .incompleteSymmetryOperations(let n, _, let listed, let required) = error {
                #expect(n == 12 && listed == 2 && required == 4)
            } else { Issue.record("wrong refusal \(error)") }
        }
        // The conventional cell with the centring dropped is still refused: diagonal rotations, two right angles.
        #expect(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 90, gammaDeg: 90, rotationsDiagonal: true) == 8)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 90, gammaDeg: 102.68, rotationsDiagonal: false) == 4)
        #expect(CIFImport.requiredOperationCount(spaceGroup: 12, alphaDeg: 90, betaDeg: 90, gammaDeg: 90) == 8)
    }
}
