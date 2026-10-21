//
//  TetragonalSettingTests.swift
//  Registration 0.0.4 section 6 (tier 1): a tetragonal structure is tetragonal whichever axis carries its 4-fold. The
//  pins are pymatgen 2026.9.24's SpacegroupAnalyzer (spglib 2.8.0) on rutile, anatase, beta-Sn, an orthorhombic control
//  with a = c by accident and a cubic one, each in three settings (tools/kit-0.0.4/tetragonal-pins.py, symprec 1e-4,
//  written before the fix). The truth beside them: the setting changes nothing, so the three settings give the same
//  Laue class, the same count of operations and the same rings (multiplicity, |g|).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct TetragonalSettingTests {

    struct Site: Decodable { let z: Int; let frac: [Double] }
    struct Entry: Decodable {
        let name: String
        let setting: String
        let cell: CellParams
        let sites: [Site]
        let spaceGroupNumber: Int
        let pointGroup: String
        let laueClass: String
        let pointGroupOperations: Int
        let laueGroupOrder: Int
        let spaceGroupOperations: Int

        var atoms: [AtomSite] {
            sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2])) }
        }
        var structure: CrystalStructure { CrystalStructure.detectingSymmetry(cell: cell.unitCell, sites: atoms) }
        var label: String { "\(name) (4-fold along \(setting))" }
    }
    struct File: Decodable {
        let versions: [String: String]
        let kMax: Double
        let entries: [Entry]
    }

    static let file: File = {
        do {
            let data = try Data(contentsOf: Fixtures.directory.appendingPathComponent("tetragonal-pins.json"))
            let decoder = JSONDecoder()
            decoder.keyDecodingStrategy = .convertFromSnakeCase
            return try decoder.decode(File.self, from: data)
        } catch { fatalError("tetragonal-pins.json: \(error)") }
    }()

    static let tetragonalNames = ["rutile_TiO2", "anatase_TiO2", "beta_Sn"]

    @Test func testPinsCarryTheirVersionsAndShape() {
        let v = Self.file.versions
        #expect(v["pymatgen"] == "2026.9.24")
        #expect(v["pymatgen_core"] == "2026.9.23")
        #expect(v["spglib"] == "2.8.0")
        #expect(v["numpy"] == "2.5.4")
        #expect(Self.file.entries.count == 15, "five structures in three settings")
        for name in Self.tetragonalNames {
            let ofThisName = Self.file.entries.filter { $0.name == name }
            let settings: Set<String> = Set(ofThisName.map(\.setting))
            let expectedSettings: Set<String> = ["a", "b", "c"]
            #expect(settings == expectedSettings)
        }
    }

    /// The kit's Laue class, its Laue group's order and its space-group operation count equal pymatgen's on every
    /// structure in every setting (tetragonal in all three, and the two controls untouched).
    @Test func testLaueClassAndOperationCountEqualPymatgenInEverySetting() {
        for entry in Self.file.entries {
            let s = entry.structure
            #expect(s.laueClass?.symbol == entry.laueClass, "\(entry.label): Laue class \(String(describing: s.laueClass))")
            #expect(s.laueGroup.order == entry.laueGroupOrder, "\(entry.label): Laue group order")
            #expect(s.operations.count == entry.spaceGroupOperations, "\(entry.label): space-group operations")
            #expect(LaueClass(pointGroupSymbol: entry.pointGroup)?.symbol == entry.laueClass, "\(entry.label): pin table")
        }
    }

    /// The setting changes nothing: the three settings of one structure give the same class and count of operations.
    @Test func testTheThreeSettingsAgree() {
        for name in Self.tetragonalNames + ["ortho_a_equals_c", "cubic_SrTiO3"] {
            let structures = Self.file.entries.filter { $0.name == name }.map { ($0, $0.structure) }
            let reference = structures[0].1
            for (entry, s) in structures {
                #expect(s.laueClass == reference.laueClass, "\(entry.label)")
                #expect(s.operations.count == reference.operations.count, "\(entry.label)")
                #expect(s.laueGroup.order == reference.laueGroup.order, "\(entry.label)")
            }
        }
    }

    /// The same rings in every setting: the same number of rings, and ring by ring the same |g|, the same multiplicity
    /// and the same number of orbits, with and without the systematic absences.
    @Test func testTheThreeSettingsGiveTheSameRingMultiplicities() throws {
        for name in Self.tetragonalNames + ["ortho_a_equals_c", "cubic_SrTiO3"] {
            let entries = Self.file.entries.filter { $0.name == name }
            let lists = try entries.map { try $0.structure.reflections(kMax: Self.file.kMax) }
            for allowedOnly in [true, false] {
                let rings = lists.map { $0.families(allowedOnly: allowedOnly) }
                for (index, entry) in entries.enumerated() where index > 0 {
                    #expect(rings[index].count == rings[0].count, "\(entry.label) allowedOnly \(allowedOnly): ring count")
                    for (got, want) in zip(rings[index], rings[0]) {
                        assertClose(got.g, want.g, relative: 1e-9, "\(entry.label) ring |g|")
                        #expect(got.multiplicity == want.multiplicity,
                                "\(entry.label) ring at |g| \(want.g): multiplicity \(got.multiplicity) vs \(want.multiplicity)")
                        #expect(got.orbits.count == want.orbits.count, "\(entry.label) ring at |g| \(want.g): orbits")
                    }
                }
            }
        }
    }

    // MARK: The metric family and its holohedry (the cause, not only the symptom)

    static let oddAxisCells: [(UnitCell.Axis, UnitCell)] = [
        (.c, .orthorhombic(a: 3.7, b: 3.7, c: 5.1)),
        (.b, .orthorhombic(a: 3.7, b: 5.1, c: 3.7)),
        (.a, .orthorhombic(a: 5.1, b: 3.7, c: 3.7)),
    ]

    /// Two equal lengths with 90 degrees and a third that differs read tetragonal about the odd axis, whichever it is;
    /// three different lengths stay orthorhombic and three equal ones cubic.
    @Test func testMetricFamilyNamesTheOddAxis() {
        for (axis, cell) in Self.oddAxisCells {
            #expect(cell.metricFamily() == .tetragonal(uniqueAxis: axis), "odd axis \(axis)")
        }
        #expect(UnitCell.orthorhombic(a: 3.7, b: 5.1, c: 6.2).metricFamily() == .orthorhombic)
        #expect(UnitCell.orthorhombic(a: 3.7, b: 3.7, c: 3.7).metricFamily() == .cubic)
    }

    /// The holohedry of each is 4/mmm (order 16), a symmetry of its own metric (every |g| preserved), with its 4-fold
    /// about the odd axis and about no other.
    @Test func testHolohedryIsBuiltAboutTheOddAxis() {
        let unit: [UnitCell.Axis: MillerIndex] = [.a: MillerIndex(1, 0, 0), .b: MillerIndex(0, 1, 0), .c: MillerIndex(0, 0, 1)]
        for (axis, cell) in Self.oddAxisCells {
            let group = PointSymmetry.holohedry(of: cell.metricFamily())
            #expect(group.order == 16, "odd axis \(axis)")
            #expect(group.laueClass == .fourOverMmm, "odd axis \(axis)")
            for r in group.rotations {
                for h in -2...2 { for k in -2...2 { for l in -2...2 {
                    let hkl = MillerIndex(h, k, l)
                    assertClose(cell.gSquaredFromMetric(hkl.applying(r)), cell.gSquaredFromMetric(hkl), relative: 1e-12,
                                "odd axis \(axis): \(r) moves |g| of \(hkl)")
                } } }
            }
            let fourFoldAxes = Set(UnitCell.Axis.allCases.filter { candidate in
                group.properRotations.contains { $0.order == 4 && unit[candidate]!.applying($0) == unit[candidate]! }
            })
            #expect(fourFoldAxes == [axis], "odd axis \(axis): 4-fold about \(fourFoldAxes)")
        }
    }
}
