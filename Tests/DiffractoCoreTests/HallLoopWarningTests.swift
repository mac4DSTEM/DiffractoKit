//
//  HallLoopWarningTests.swift
//  A CIF that lists a full operation loop and also states a Hall symbol of another group: the loop is used, as
//  pymatgen-core 2026.10.16 `CifParser` uses it (it reads no Hall symbol), and the document carries a warning that names
//  the symbol and the two orders. The groups are spglib 2.8.0's (`hall-symbol-pins.json`).
//

import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoCore

struct HallLoopWarningTests {
    static func cif(hall: String?, ops: [String]) -> String {
        var s = "data_t\n_symmetry_Int_Tables_number 123\n"
        if let hall { s += "_space_group_name_Hall '\(hall)'\n" }
        for (tag, v) in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"],
                            [10.0, 10.0, 12.0, 90, 90, 90]) {
            s += "_cell_\(tag) \(v)\n"
        }
        s += "loop_\n_symmetry_equiv_pos_as_xyz\n" + ops.map { "'\($0)'" }.joined(separator: "\n") + "\n"
        s += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        s += "Si1 Si 0.1 0.2 0.3\n"
        return s
    }

    /// A P4/mmm loop beside the Hall symbol of mmm (-P 2 2): the 16 images of a general site are read from the loop, and
    /// one warning names the symbol and the orders 16 and 8. RED UNDER: `CIFDocument` had no warnings.
    @Test func aLoopOfAnotherGroupThanTheHallSymbolIsUsedWithAWarning() throws {
        let loop = try OperationListCheckTests.group("-P 4 2").map(\.xyzString)
        #expect(loop.count == 16)
        let document = try CIFImport.read(Self.cif(hall: "-P 2 2", ops: loop))
        #expect(document.structure.sites.count == 16)
        #expect(document.warnings == [.operationLoopIsNotTheHallGroup(hallSymbol: "-P 2 2", loopOrder: 16, hallOrder: 8)])
        let text = try #require(document.warnings.first).description
        #expect(text.contains("-P 2 2") && text.contains("16") && text.contains("8") && text.contains("loop was used"))
    }

    /// The same loop beside its own Hall symbol (-P 4 2), and beside none: no warning, the same 16 sites.
    @Test func aLoopOfTheHallGroupOrWithoutAHallSymbolHasNoWarning() throws {
        let loop = try OperationListCheckTests.group("-P 4 2").map(\.xyzString)
        for hall in ["-P 4 2", nil] as [String?] {
            let document = try CIFImport.read(Self.cif(hall: hall, ops: loop))
            #expect(document.structure.sites.count == 16)
            #expect(document.warnings.isEmpty, "\(hall ?? "no Hall symbol")")
        }
    }

    /// The loop written in reverse order, with every translation moved by a lattice vector, is still the Hall group: the
    /// comparison is of sets, modulo lattice translations.
    @Test func theComparisonIgnoresOrderAndLatticeTranslations() throws {
        let ops = try OperationListCheckTests.group("-P 4 2")
        let shifted = ops.reversed().map { op in
            SymmetryOperation(rotation: op.rotation, translation: op.translation + SIMD3(1, 0, -1)).xyzString
        }
        let document = try CIFImport.read(Self.cif(hall: "-P 4 2", ops: shifted))
        #expect(document.warnings.isEmpty)
    }

    /// Every file read without a Hall symbol, or with only a Hall symbol, has no warning.
    @Test func aHallOnlyFileHasNoWarning() throws {
        var s = "data_t\n_space_group_name_Hall '-P 4 2'\n"
        for (tag, v) in zip(["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"],
                            [10.0, 10.0, 12.0, 90, 90, 90]) {
            s += "_cell_\(tag) \(v)\n"
        }
        s += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        s += "Si1 Si 0.1 0.2 0.3\n"
        let document = try CIFImport.read(s)
        #expect(document.structure.sites.count == 16)
        #expect(document.warnings.isEmpty)
    }
}
