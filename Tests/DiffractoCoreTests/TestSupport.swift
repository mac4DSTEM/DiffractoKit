//
//  TestSupport.swift
//  Role: the fixtures directory (read by #filePath, not as resources), the β″ CIF text rebuilt from the pinned
//        generator's values (tools/crystal-structures/make_beta_double_prime.py, verbatim from mac4DSTEM 7448401b; the
//        pins carry its sha256), and the pins file written by tools/slice-1/pins.py (registration section 4).
//

import Foundation
import XCTest
import DiffractoCore

enum Fixtures {
    static let directory: URL = URL(fileURLWithPath: #filePath)
        .deletingLastPathComponent().appendingPathComponent("Fixtures")

    static func text(_ name: String) throws -> String {
        try String(contentsOf: directory.appendingPathComponent(name), encoding: .utf8)
    }

    /// The β″ CIF exactly as the generator writes it (C2/m, Andersen et al. 1998, Table 3 set 3).
    static var betaDoublePrimeCIF: String {
        let sites: [(String, String, Double, Double, Double, Double)] = [
            ("Mg1", "Mg", 0.0, 0.0, 0.0, 0.5),
            ("Mg2", "Mg", 0.3459, 0.0, 0.089, 1.0),
            ("Mg3", "Mg", 0.430, 0.0, 0.652, 0.8),
            ("Si1", "Si", 0.0565, 0.0, 0.649, 1.1),
            ("Si2", "Si", 0.1885, 0.0, 0.224, 0.5),
            ("Si3", "Si", 0.2171, 0.0, 0.617, 2.5),
        ]
        let ops = ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z",
                   "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z", "1/2+x,1/2-y,z"]
        var lines = [
            "data_beta_double_prime_Mg5Si6",
            "_chemical_name_mineral            'beta-double-prime'",
            "_chemical_formula_sum             'Mg10 Si12'",
            "_chemical_formula_structural      'Mg5 Si6'",
            "_publ_author_name                 'Andersen, S. J.'",
            "_journal_name_full                'Acta Materialia'",
            "_journal_volume                   46",
            "_journal_year                     1998",
            "_journal_page_first               3283",
            "_cell_length_a                    15.16",
            "_cell_length_b                    4.05",
            "_cell_length_c                    6.74",
            "_cell_angle_alpha                 90.0",
            "_cell_angle_beta                  105.3",
            "_cell_angle_gamma                 90.0",
            "_symmetry_space_group_name_H-M    'C 2/m'",
            "_symmetry_Int_Tables_number       12",
            "loop_",
            "_symmetry_equiv_pos_as_xyz",
        ]
        lines += ops.map { "  '\($0)'" }
        lines += [
            "loop_",
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
            "_atom_site_B_iso_or_equiv",
        ]
        lines += sites.map { (lab, el, x, y, z, b) in
            String(format: "  %-5@ %-3@ %8.4f %8.4f %8.4f  1.0  %.1f", lab, el, x, y, z, b)
        }
        return lines.joined(separator: "\n") + "\n"
    }

    /// The β″ cell as the apps' `Crystal.betaDoublePrime` builds it in code: the six sites expanded by hand under C2/m
    /// (22 atoms), for the detection tests (no CIF, no operations).
    static var betaDoublePrimeSites: [AtomSite] {
        func fourI(_ z: Int, _ x: Double, _ zc: Double) -> [AtomSite] {
            [SIMD3(x, 0, zc), SIMD3(-x, 0, -zc), SIMD3(x + 0.5, 0.5, zc), SIMD3(-x + 0.5, 0.5, -zc)]
                .map { AtomSite(z: z, fractional: wrapFractional($0)) }
        }
        let mg = 12, si = 14
        var sites: [AtomSite] = [
            AtomSite(z: mg, fractional: [0, 0, 0]),
            AtomSite(z: mg, fractional: [0.5, 0.5, 0]),
        ]
        sites += fourI(mg, 0.3459, 0.089)
        sites += fourI(mg, 0.430, 0.652)
        sites += fourI(si, 0.0565, 0.649)
        sites += fourI(si, 0.1885, 0.224)
        sites += fourI(si, 0.2171, 0.617)
        return sites
    }

    static let betaDoublePrimeCell = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)
}

// MARK: - The pins (tools/slice-1/pins.py)

struct Pins: Decodable {
    let versions: [String: String]
    let generatedAt: String
    let cells: [CellPin]
    let structures: [StructurePin]

    static let fileName = "slice-1-pins.json"

    static func load() throws -> Pins {
        let data = try Data(contentsOf: Fixtures.directory.appendingPathComponent(fileName))
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        return try decoder.decode(Pins.self, from: data)
    }

    func structure(_ name: String) throws -> StructurePin {
        try XCTUnwrap(structures.first { $0.name == name }, "pins hold no structure \(name)")
    }
}

struct CellParams: Decodable {
    let a: Double, b: Double, c: Double
    let alpha: Double, beta: Double, gamma: Double

    var unitCell: UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma)
    }
}

struct CellPin: Decodable {
    let name: String
    let a: Double, b: Double, c: Double
    let alpha: Double, beta: Double, gamma: Double
    let matrix: [[Double]]
    let metric: [[Double]]
    let volume: Double
    let reciprocal: [[Double]]
    let d: [DPin]

    var unitCell: UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma)
    }
}

struct DPin: Decodable {
    let hkl: [Int]
    let g: Double
    let d: Double
}

struct SitePin: Decodable {
    let z: Int
    let frac: [Double]
    let occupancy: Double
}

struct OperationPin: Decodable {
    let rotation: [[Int]]
    let translation: [Double]
    let xyz: String

    var operation: SymmetryOperation {
        SymmetryOperation(rotation: IntegerMatrix(rows: rotation.map { SIMD3($0[0], $0[1], $0[2]) }),
                          translation: SIMD3(translation[0], translation[1], translation[2]))
    }
}

struct ReflectionPin: Decodable {
    let hkl: [Int]
    let g: Double
    let d: Double
    let allowed: Bool
    let orbit: [[Int]]
    let canonical: [Int]
}

struct ReflectionsPin: Decodable {
    let kMax: Double
    let count: Int
    let allowedCount: Int
    let items: [ReflectionPin]
}

struct FamilyPin: Decodable {
    let g: Double
    let d: Double
    let members: [[Int]]
    let multiplicity: Int
    let orbits: [[Int]]
    let name: String
}

struct StructurePin: Decodable {
    let name: String
    let cifFile: String?
    let cifSha256: String?
    let cifText: String?
    let cell: CellParams
    let matrix: [[Double]]
    let metric: [[Double]]
    let volume: Double
    let reciprocal: [[Double]]
    let sites: [SitePin]
    let declaredSpaceGroupNumber: Int?
    let spaceGroupNumber: Int
    let pointGroup: String
    let laueClass: String
    let operations: [OperationPin]
    let laueRotations: [[[Int]]]
    let reflections: ReflectionsPin
    let families: [FamilyPin]
}

extension MillerIndex {
    init(_ array: [Int]) { self.init(array[0], array[1], array[2]) }
}

extension IntegerMatrix {
    init(_ rows: [[Int]]) { self.init(rows: rows.map { SIMD3($0[0], $0[1], $0[2]) }) }
}

/// |got − want| ≤ tolerance · max(1, |want|) when `relative`, else absolute.
func assertClose(_ got: Double, _ want: Double, relative: Double, _ message: @autoclosure () -> String = "",
                 file: StaticString = #filePath, line: UInt = #line) {
    let bound = relative * max(abs(want), abs(got))
    XCTAssertLessThanOrEqual(abs(got - want), bound, "\(message()) got \(got), want \(want)", file: file, line: line)
}
