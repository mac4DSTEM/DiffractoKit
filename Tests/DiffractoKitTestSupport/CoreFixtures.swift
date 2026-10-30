//
//  TestSupport.swift
//  Role: the fixtures directory (read by #filePath, not as resources), the β″ CIF text rebuilt from the pinned
//        generator's values (tools/crystal-structures/make_beta_double_prime.py, verbatim from reference commit 7448401b; the
//        pins carry its sha256), and the pins file written by tools/pins/cells/pins.py.
//

import Foundation
import Testing
import DiffractoCore

package enum Fixtures {
    package static let directory: URL = FixtureDirs.fixtures(of: "DiffractoCoreTests")

    package static func text(_ name: String) throws -> String {
        try String(contentsOf: directory.appendingPathComponent(name), encoding: .utf8)
    }

    /// The β″ CIF exactly as the generator writes it (C2/m, Andersen et al. 1998, Table 3 set 3).
    package static var betaDoublePrimeCIF: String {
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
            "  " + lab.padding(toLength: 5, withPad: " ", startingAt: 0) + " "
                + el.padding(toLength: 3, withPad: " ", startingAt: 0)
                + String(format: " %8.4f %8.4f %8.4f  1.0  %.1f", x, y, z, b)
        }
        return lines.joined(separator: "\n") + "\n"
    }

    /// The β″ cell as the reference implementations' `Crystal.betaDoublePrime` builds it in code: the six sites expanded by hand under C2/m
    /// (22 atoms), for the detection tests (no CIF, no operations).
    package static var betaDoublePrimeSites: [AtomSite] {
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

    package static let betaDoublePrimeCell = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)
}

// MARK: - The pins (tools/pins/cells/pins.py)

package struct Pins: Decodable {
    package let versions: [String: String]
    package let generatedAt: String
    package let cells: [CellPin]
    package let structures: [StructurePin]

    package static let fileName = "pymatgen-cells-pins.json"

    package static func load() throws -> Pins {
        let data = try Data(contentsOf: Fixtures.directory.appendingPathComponent(fileName))
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        return try decoder.decode(Pins.self, from: data)
    }

    package func structure(_ name: String) throws -> StructurePin {
        try #require(structures.first { $0.name == name }, "pins hold no structure \(name)")
    }
}

package struct CellParams: Decodable {
    package let a: Double, b: Double, c: Double
    package let alpha: Double, beta: Double, gamma: Double

    package var unitCell: UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma)
    }
}

package struct CellPin: Decodable {
    package let name: String
    package let a: Double, b: Double, c: Double
    package let alpha: Double, beta: Double, gamma: Double
    package let matrix: [[Double]]
    package let metric: [[Double]]
    package let volume: Double
    package let reciprocal: [[Double]]
    package let d: [DPin]

    package var unitCell: UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: alpha, betaDeg: beta, gammaDeg: gamma)
    }
}

package struct DPin: Decodable {
    package let hkl: [Int]
    package let g: Double
    package let d: Double
}

package struct SitePin: Decodable {
    package let z: Int
    package let frac: [Double]
    package let occupancy: Double
}

package struct OperationPin: Decodable {
    package let rotation: [[Int]]
    package let translation: [Double]
    package let xyz: String

    package var operation: SymmetryOperation {
        SymmetryOperation(rotation: IntegerMatrix(rows: rotation.map { SIMD3($0[0], $0[1], $0[2]) }),
                          translation: SIMD3(translation[0], translation[1], translation[2]))
    }
}

package struct ReflectionPin: Decodable {
    package let hkl: [Int]
    package let g: Double
    package let d: Double
    package let allowed: Bool
    /// Carried by the first items of each list only (size).
    package let orbit: [[Int]]?
    package let canonical: [Int]
}

/// What the file's loop is against the detected set.
package struct FileOperationsPin: Decodable {
    package let count: Int
    package let equalsDetected: Bool
    package let subsetOfDetected: Bool
}

package struct ReflectionsPin: Decodable {
    package let kMax: Double
    package let count: Int
    package let allowedCount: Int
    package let items: [ReflectionPin]
}

package struct FamilyPin: Decodable {
    package let g: Double
    package let d: Double
    package let members: [[Int]]
    package let multiplicity: Int
    package let orbits: [[Int]]
    package let name: String
}

package struct StructurePin: Decodable {
    package let name: String
    package let cifFile: String?
    package let cifSha256: String?
    package let cifText: String?
    package let cell: CellParams
    package let matrix: [[Double]]
    package let metric: [[Double]]
    package let volume: Double
    package let reciprocal: [[Double]]
    /// As pymatgen expanded them (coordinates within 1e-4 of ⅓ and ⅔ snapped to the fraction).
    package let sites: [SitePin]
    /// The sites with the written coordinates, where pymatgen snapped any; what the port is compared with.
    package let sitesAsWritten: [SitePin]?
    package let siteSnapMax: Double?
    package let declaredSpaceGroupNumber: Int?
    package let spaceGroupNumber: Int
    package let pointGroup: String
    package let laueClass: String
    package let operations: [OperationPin]
    package let laueRotations: [[[Int]]]
    package let reflections: ReflectionsPin
    package let families: [FamilyPin]
    package let fileOperations: FileOperationsPin?

    package var expectedSites: [SitePin] { sitesAsWritten ?? sites }
}

package extension MillerIndex {
    init(_ array: [Int]) { self.init(array[0], array[1], array[2]) }
}

package extension IntegerMatrix {
    init(_ rows: [[Int]]) { self.init(rows: rows.map { SIMD3($0[0], $0[1], $0[2]) }) }
}

/// |got − want| ≤ tolerance · max(1, |want|) when `relative`, else absolute.
package func assertClose(_ got: Double, _ want: Double, relative: Double, _ message: @autoclosure () -> String = "",
                 sourceLocation: SourceLocation = #_sourceLocation) {
    let bound = relative * max(abs(want), abs(got))
    #expect(abs(got - want) <= bound, "\(message()) got \(got), want \(want)", sourceLocation: sourceLocation)
}
