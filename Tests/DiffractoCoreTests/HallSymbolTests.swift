//
//  HallSymbolTests.swift
//  The Hall-symbol reader against spglib 2.8.0's 530 settings, exact, and the CIF that names its group by Hall symbol
//  alone. Pins: tools/pins/cif-symmetry/hall-pins.py.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct HallSymbolTests {
    struct Pins: Decodable {
        struct Group: Decodable { let hall: Int; let number: Int; let symbol: String; let operations: [String] }
        struct PlantedCIF: Decodable {
            let name: String; let text: String; let hall: String
            let truthAtoms: Int; let pymatgenAtoms: Int; let listed: Int?; let groupOrder: Int
        }
        let spglib: String
        let groups: [Group]
        let cifs: [PlantedCIF]
    }

    static func pins() throws -> Pins {
        try JSONDecoder().decode(Pins.self, from: Data(contentsOf: Fixtures.directory.appendingPathComponent("hall-symbol-pins.json")))
    }

    /// The fixture's encoding: the rotation row by row as - 0 +, the translation in twelfths as three hex digits.
    static func encode(_ op: SymmetryOperation) -> String {
        var s = ""
        for row in op.rotation.rows { for v in [row.x, row.y, row.z] { s.append(v < 0 ? "-" : v > 0 ? "+" : "0") } }
        for t in [op.translation.x, op.translation.y, op.translation.z] {
            let twelfths = Int((t * 12).rounded())
            s.append(Character(String(((twelfths % 12) + 12) % 12, radix: 16)))
        }
        return s
    }

    /// Every one of the 530 Hall settings spglib lists reads to exactly spglib's operation set: same count, same
    /// rotations, same translations in twelfths, centring included. No tolerance (readers are tested exact).
    @Test func everySpglibHallSymbolReadsToItsOperations() throws {
        let pins = try Self.pins()
        #expect(pins.spglib == "2.8.0")
        #expect(pins.groups.count == 530)
        var failures: [String] = []
        for group in pins.groups {
            do {
                let ops = try HallSymbol.operations(group.symbol)
                let encoded = ops.map(Self.encode)
                if encoded.count != Set(encoded).count || Set(encoded) != Set(group.operations) {
                    failures.append("\(group.hall) \(group.symbol): \(encoded.count) vs \(group.operations.count)")
                }
                if ops.first != SymmetryOperation.identity { failures.append("\(group.symbol): identity not first") }
            } catch {
                failures.append("\(group.hall) \(group.symbol): \(error)")
            }
        }
        #expect(failures.isEmpty, "\(failures.count) settings differ: \(failures.prefix(20))")
    }

    /// A CIF that names its group by Hall symbol alone builds the cell spglib's operations give: with no loop it is
    /// expanded with the symbol's group (pymatgen builds the asymmetric unit as P1 here); with a partial loop it is
    /// refused naming the symbol and the group's order (pymatgen expands the partial loop silently).
    /// RED UNDER: the Hall symbol unread (no loop: refused as missing operations; partial loop: admitted, 2 atoms of 8).
    @Test func aCIFNamedByHallSymbolAloneIsReadAndChecked() throws {
        let pins = try Self.pins()
        #expect(pins.cifs.count == 6)
        for planted in pins.cifs {
            if let listed = planted.listed {
                #expect(throws: CIFImportError.incompleteSymmetryOperationsForSymbol(
                    symbol: planted.hall, listed: listed, required: planted.groupOrder), "\(planted.name)") {
                    try CIFImport.read(planted.text)
                }
                #expect(planted.pymatgenAtoms != planted.truthAtoms, "\(planted.name): pymatgen is the wrong side here")
            } else {
                let doc = try CIFImport.read(planted.text)
                #expect(doc.structure.sites.count == planted.truthAtoms, "\(planted.name)")
                #expect(doc.declaredHallSymbol == planted.hall, "\(planted.name)")
                #expect(doc.hadSymmetryOperations, "\(planted.name)")
                #expect(doc.structure.fileOperations?.count == planted.groupOrder, "\(planted.name)")
                #expect(planted.pymatgenAtoms < planted.truthAtoms, "\(planted.name): pymatgen expands as P1")
            }
        }
        // An unreadable Hall symbol with no loop is still refused; "P 1" with no loop is the cell as listed.
        let rockSalt = try #require(pins.cifs.first).text
        #expect(throws: CIFImportError.missingSymmetryOperations(spaceGroup: "-F 4 2 7")) {
            try CIFImport.read(rockSalt.replacingOccurrences(of: "-F 4 2 3", with: "-F 4 2 7"))
        }
        let p1 = try CIFImport.read(rockSalt.replacingOccurrences(of: "-F 4 2 3", with: "P 1"))
        #expect(p1.structure.sites.count == 2)
        #expect(!p1.hadSymmetryOperations)
        // A Hall symbol smaller than the IT number beside it is refused like a short loop (a mistyped "-P 1").
        let mistyped = rockSalt.replacingOccurrences(of: "-F 4 2 3", with: "-P 1") + "_space_group_IT_number 225\n"
        #expect(throws: CIFImportError.incompleteSymmetryOperations(
            spaceGroupNumber: 225, declaredName: nil, listed: 2, required: 192)) { try CIFImport.read(mistyped) }
        // A loop, even an identity-only one, wins over the Hall symbol: as before, the non-P1 symbol is refused.
        let identityLoop = rockSalt + "loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n"
        #expect(throws: CIFImportError.missingSymmetryOperations(spaceGroup: "-F 4 2 3")) {
            try CIFImport.read(identityLoop)
        }
        // A full loop beside the Hall symbol and the IT number reads as before (the loop is expanded).
        let full = try HallSymbol.operations("-F 4 2 3").map { "'\($0.xyzString)'" }.joined(separator: "\n")
        let doc = try CIFImport.read(rockSalt + "_space_group_IT_number 225\nloop_\n_symmetry_equiv_pos_as_xyz\n" + full + "\n")
        #expect(doc.structure.sites.count == 8)
    }

    @Test func spellingsAndRefusals() throws {
        // The lattice run into the first matrix, a lowercase lattice letter, extra spaces.
        let reference = Set(try HallSymbol.operations("-P 2ybc").map(Self.encode))
        #expect(Set(try HallSymbol.operations("-P2ybc").map(Self.encode)) == reference)
        #expect(Set(try HallSymbol.operations(" -p  2ybc ").map(Self.encode)) == reference)
        #expect(reference.count == 4)
        // An origin shift is read mod 12 twelfths.
        #expect(Set(try HallSymbol.operations("P 31 2 (0 0 16)").map(Self.encode))
            == Set(try HallSymbol.operations("P 31 2 (0 0 4)").map(Self.encode)))
        #expect(try HallSymbol.operations("P 1") == [.identity])
        for bad in ["", "Q 2", "P", "P 5", "P 2q", "P 4 4", "P 2 2 (0 0)", "P 2 2 (0 0 1", "P 3 2 1 1 1", "P 2yb3", "P 3 2 1", "ß 2", "P 2 2 (12 0 99999999999999999999)"] {
            #expect(throws: HallSymbolError.self, "\(bad)") { try HallSymbol.operations(bad) }
        }
        #expect(HallSymbolError.unreadable(symbol: "P 5", reason: "r").description.contains("P 5"))
    }
}
