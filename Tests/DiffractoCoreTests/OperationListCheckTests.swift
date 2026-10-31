//
//  OperationListCheckTests.swift
//  A listed operation list whose count is right but whose operations are wrong: not symmetries of the cell, not closed,
//  or of another crystal class than the IT number's. The groups are spglib 2.8.0's (`hall-symbol-pins.json`,
//  tools/pins/cif-symmetry/hall-pins.py). pymatgen-core 2026.10.16 `CifParser` admits every one of these lists and
//  expands it as given.
//

import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct OperationListCheckTests {
    /// The fixture's 12-character operation: rotation as - 0 +, translation in twelfths as hex.
    static func decode(_ code: String) -> SymmetryOperation {
        let c = Array(code)
        let entries = c[0..<9].map { $0 == "-" ? -1 : $0 == "+" ? 1 : 0 }
        let rotation = IntegerMatrix(SIMD3(entries[0], entries[1], entries[2]), SIMD3(entries[3], entries[4], entries[5]),
                                     SIMD3(entries[6], entries[7], entries[8]))
        let t = c[9..<12].map { Double(Int(String($0), radix: 16) ?? 0) / 12 }
        return SymmetryOperation(rotation: rotation, translation: SIMD3(t[0], t[1], t[2]))
    }

    static func group(_ symbol: String) throws -> [SymmetryOperation] {
        let pins = try HallSymbolTests.pins()
        return try #require(pins.groups.first { $0.symbol == symbol }).operations.map(decode)
    }

    static func cif(number: Int?, cell: [Double], ops: [String], sites: [(String, Double, Double, Double)]) -> String {
        var s = "data_t\n"
        if let number { s += "_symmetry_Int_Tables_number \(number)\n" }
        for (tag, v) in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"], cell) {
            s += "_cell_\(tag) \(v)\n"
        }
        s += "loop_\n_symmetry_equiv_pos_as_xyz\n" + ops.map { "'\($0)'" }.joined(separator: "\n") + "\n"
        s += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        for (i, site) in sites.enumerated() { s += "\(site.0)\(i) \(site.0) \(site.1) \(site.2) \(site.3)\n" }
        return s
    }

    static func problem(_ text: String) -> SymmetryListProblem? {
        do { _ = try CIFImport.read(text); return nil } catch {
            guard case .inconsistentSymmetryOperations(let p) = error else { return nil }
            return p
        }
    }

    static let general: [(String, Double, Double, Double)] = [("Si", 0.1, 0.2, 0.3)]

    /// The crystal class read from the rotations equals the class of the IT number for every one of spglib's 530
    /// settings (centred, rhombohedral and origin-shifted ones included), and every IT number has a class.
    @Test func theCrystalClassOfEverySpglibSettingIsItsNumbers() throws {
        let pins = try HallSymbolTests.pins()
        var failures: [String] = []
        for g in pins.groups {
            let listed = CIFImport.crystalClass(of: g.operations.map { Self.decode($0).rotation })
            let required = CIFImport.crystalClass(spaceGroup: g.number)
            if listed == nil || listed != required { failures.append("\(g.symbol) No. \(g.number): \(listed ?? "nil") vs \(required ?? "nil")") }
        }
        #expect(failures.isEmpty, "\(failures.prefix(10))")
        #expect((1...230).allSatisfy { CIFImport.crystalClass(spaceGroup: $0) != nil })
        #expect(CIFImport.crystalClass(spaceGroup: 231) == nil)
        #expect(Set((1...230).compactMap { CIFImport.crystalClass(spaceGroup: $0) }).count == 32)
    }

    /// Fm-3m in its primitive 60-degree cell: the 48 operations listed in the conventional basis pass the count (192 / 4)
    /// but do not map the cell onto itself, and are refused by name; spglib's operations brought into that basis are
    /// admitted with the 48 images of a general site. RED UNDER: no metric check (the wrong list was expanded).
    @Test func operationsOfAnotherBasisAreRefused() throws {
        let cell = [10.0, 10.0, 10.0, 60, 60, 60]    // large enough for 48 images of a general site
        let pointOps = try Self.group("-P 4 2 3").map(\.xyzString)
        guard case .notACellSymmetry? = Self.problem(Self.cif(number: 225, cell: cell, ops: pointOps, sites: Self.general))
        else { Issue.record("conventional-basis operations admitted in the primitive cell"); return }

        // x_prim = P^-1 x_conv with P = [[0,1,1],[1,0,1],[1,1,0]] / 2: R' = P^-1 R P, t' = P^-1 t.
        let p: [[Double]] = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
        let pInverse: [[Double]] = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        func mul(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
            (0..<3).map { i in (0..<3).map { j in (0..<3).reduce(0) { $0 + a[i][$1] * b[$1][j] } } }
        }
        var primitiveOps = Set<String>()
        for op in try Self.group("-F 4 2 3") {
            let r = (0..<3).map { i in (0..<3).map { j in Double(op.rotation[i, j]) } }
            let rp = mul(mul(pInverse, r), p).map { $0.map { Int($0.rounded()) } }
            let t = (0..<3).map { i in (0..<3).reduce(0.0) { $0 + pInverse[i][$1] * op.translation[$1] } }
            let prim = SymmetryOperation(rotation: IntegerMatrix(SIMD3(rp[0]), SIMD3(rp[1]), SIMD3(rp[2])),
                                         translation: wrapFractional(SIMD3(t[0], t[1], t[2])))
            primitiveOps.insert(prim.xyzString)
        }
        #expect(primitiveOps.count == 48)
        let doc = try CIFImport.read(Self.cif(number: 225, cell: cell, ops: primitiveOps.sorted(), sites: Self.general))
        #expect(doc.structure.sites.count == 48)
    }

    /// P4/mmm with one of its 16 operations swapped for an operation of another group: the count is right, the list is
    /// not closed, and the refusal names a missing product. RED UNDER: no closure check (16 atoms of a structure with no
    /// symmetry were built).
    @Test func aListThatIsNotClosedIsRefused() throws {
        let ops = try Self.group("-P 4 2").map(\.xyzString)
        #expect(ops.count == 16 && ops.contains("-y,x,z"))
        let swapped = ops.filter { $0 != "-y,x,z" } + ["y,x+1/2,z"]
        let cell = [4.0, 4.0, 5.0, 90, 90, 90]
        guard case .notClosed? = Self.problem(Self.cif(number: 123, cell: cell, ops: swapped, sites: Self.general))
        else { Issue.record("a list that is not a group was admitted"); return }
        #expect(try CIFImport.read(Self.cif(number: 123, cell: cell, ops: ops, sites: Self.general)).structure.sites.count == 16)
    }

    /// A hand-rounded list (P3_1 with z+0.33 and z+0.67, which compose to z+0.66) is still a group: translations within
    /// 0.004 of a twelfth (all of them) are read as twelfths. An a-b swap in an elongated cell is no symmetry of it.
    @Test func roundedTranslationsAreClosedAndElongatedCellsAreJudgedPerAxis() throws {
        let p31 = ["x,y,z", "-y,x-y,z+0.33", "-x+y,-x,z+0.67"]
        let hexagonal = [5.0, 5.0, 7.0, 90, 90, 120]
        #expect(try CIFImport.read(Self.cif(number: 144, cell: hexagonal, ops: p31, sites: Self.general)).structure.sites.count == 3)
        let swap = ["x,y,z", "y,x,-z"]
        guard case .notACellSymmetry? = Self.problem(Self.cif(number: nil, cell: [3.0, 4.0, 50.0, 90, 90, 90], ops: swap,
                                                             sites: Self.general))
        else { Issue.record("an a-b swap was admitted in a 3 x 4 x 50 cell"); return }
    }

    /// A file that declares P4/m (No. 83) and lists the 8 operations of Pmmm: closed, the right count, symmetries of the
    /// tetragonal cell, but of class mmm, not 4/m. RED UNDER: no class check (the Pmmm cell was built).
    @Test func aListOfAnotherCrystalClassIsRefused() throws {
        let cell = [4.0, 4.0, 5.0, 90, 90, 90]
        let pmmm = try Self.group("-P 2 2").map(\.xyzString)
        #expect(Self.problem(Self.cif(number: 83, cell: cell, ops: pmmm, sites: Self.general))
            == .wrongCrystalClass(spaceGroupNumber: 83, listed: "mmm", required: "4/m"))
        let p4m = try Self.group("-P 4").map(\.xyzString)
        #expect(try CIFImport.read(Self.cif(number: 83, cell: cell, ops: p4m, sites: Self.general)).structure.sites.count == 8)
    }

    /// A Hall symbol in hexagonal axes (`-R 3 2"`) on a cell written in rhombohedral axes: the group is closed, but its
    /// three-fold is no symmetry of that cell, so the file is refused rather than given R-centring images.
    @Test func aHallSettingThatDoesNotFitTheCellIsRefused() throws {
        var text = "data_t\n_space_group_name_Hall '-R 3 2\"'\n"
        for (tag, v) in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"],
                            ["5.1", "5.1", "5.1", "55", "55", "55"]) { text += "_cell_\(tag) \(v)\n" }
        text += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        text += "Al1 Al 0.35 0.35 0.35\n"
        guard case .notACellSymmetry? = Self.problem(text) else { Issue.record("a Hall setting of another cell was used"); return }
        let rhombohedral = text.replacingOccurrences(of: "-R 3 2\"", with: "-P 3* 2")
        #expect(try CIFImport.read(rhombohedral).structure.sites.count == 2)
    }

    @Test func theProblemsSayWhatIsWrong() {
        #expect(SymmetryListProblem.notACellSymmetry(operation: "-x,y,z").description.contains("-x,y,z"))
        #expect(SymmetryListProblem.notClosed(first: "a", second: "b", product: "c").description.contains("\"c\""))
        #expect(SymmetryListProblem.wrongCrystalClass(spaceGroupNumber: 83, listed: "mmm", required: "4/m").description
            .contains("4/m"))
        #expect(CIFImportError.inconsistentSymmetryOperations(.notACellSymmetry(operation: "-x,y,z")).description
            .contains("refused"))
    }
}
