//
//  IndexFixtures.swift
//  Role: structure builders for the indexing tests, shared by the unit and the parity tests.
//

import Foundation
import Testing
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

/// Structure builders shared by this target's indexing tests (a small copy of DiffractoScatteringTests' Presets: a test
/// target cannot import another).
package enum IndexPresets {
    package static let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    package static func structure(_ cell: UnitCell, _ sites: [AtomSite]) -> CrystalStructure {
        CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
    }
    package static func fccElement(_ z: Int, a: Double) -> CrystalStructure {
        structure(.cubic(a: a), fcc.map { AtomSite(z: z, fractional: $0) })
    }
    package static var gold: CrystalStructure { fccElement(79, a: 4.0782) }
    package static var aluminium: CrystalStructure { fccElement(13, a: 4.0495) }
    package static var silicon: CrystalStructure {
        structure(.cubic(a: 5.4309), (fcc + fcc.map { $0 + SIMD3(0.25, 0.25, 0.25) }).map { AtomSite(z: 14, fractional: $0) })
    }
    package static func hcp(_ z: Int, a: Double, c: Double) -> CrystalStructure {
        structure(.hexagonal(a: a, c: c), [AtomSite(z: z, fractional: [0, 0, 0]),
                                            AtomSite(z: z, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5])])
    }
    package static var magnesium: CrystalStructure { hcp(12, a: 3.2094, c: 5.2108) }
    package static var idealHCP: CrystalStructure { hcp(12, a: 3.2094, c: 3.2094 * (8.0 / 3.0).squareRoot()) }
    /// the reference implementation's Weak-cut cell: Cu on the origin and on (1/2, 1/2, 1/2) at occupancy 0.9981472011156598, a = 3 Å.
    package static var weakCut: CrystalStructure {
        structure(.cubic(a: 3), [AtomSite(z: 29, fractional: [0, 0, 0]),
                                 AtomSite(z: 29, fractional: [0.5, 0.5, 0.5], occupancy: 0.9981472011156598)])
    }
    /// The β″ CIF exactly as the generator writes it (the core fixture builder).
    package static var betaDoublePrime: CrystalStructure {
        let sites: [(String, String, Double, Double, Double, Double)] = [
            ("Mg1", "Mg", 0.0, 0.0, 0.0, 0.5), ("Mg2", "Mg", 0.3459, 0.0, 0.089, 1.0), ("Mg3", "Mg", 0.430, 0.0, 0.652, 0.8),
            ("Si1", "Si", 0.0565, 0.0, 0.649, 1.1), ("Si2", "Si", 0.1885, 0.0, 0.224, 0.5), ("Si3", "Si", 0.2171, 0.0, 0.617, 2.5),
        ]
        var lines = ["data_beta_double_prime_Mg5Si6", "_cell_length_a 15.16", "_cell_length_b 4.05", "_cell_length_c 6.74",
                     "_cell_angle_alpha 90.0", "_cell_angle_beta 105.3", "_cell_angle_gamma 90.0",
                     "_symmetry_space_group_name_H-M 'C 2/m'", "_symmetry_Int_Tables_number 12",
                     "loop_", "_symmetry_equiv_pos_as_xyz"]
        lines += ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z",
                  "1/2+x,1/2-y,z"].map { "  '\($0)'" }
        lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
                  "_atom_site_fract_z"]
        lines += sites.map { "  \($0.0) \($0.1) \($0.2) \($0.3) \($0.4)" }
        // The CIF parse cannot fail on this fixed text; a failure is a slice-1 regression the CIF tests name.
        return (try? CIFImport.read(lines.joined(separator: "\n")).structure)!
    }

    package static let alOrigin = [AtomSite(z: 13, fractional: [0, 0, 0])]
    package static let alCentred = [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.5, 0.5, 0])]
    /// the reference implementation's monoclinic cells: M-P (a 7.2, b 4.9, c 5.6 Å, β 104°, one Al), M-C (C-centred), M-sq (a = c = 5, b = 4,
    /// cos β = -1/4: accidental (h k l) ~ (l k h) coincidences).
    package static func mono(_ a: Double, _ b: Double, _ c: Double, _ betaDeg: Double, sites: [AtomSite] = alOrigin) -> CrystalStructure {
        structure(.monoclinic(a: a, b: b, c: c, betaDeg: betaDeg), sites)
    }
    package static var mP: CrystalStructure { mono(7.2, 4.9, 5.6, 104) }
    package static var mC: CrystalStructure { mono(7.2, 4.9, 5.6, 104, sites: alCentred) }
    /// a = c and cos β = −1/4 make this metric C-centred orthorhombic: full detection (`conventional: true`)
    /// reads mmm and relates {1 0 0} and {0 0 1} by symmetry. the reference implementation's accidental-coincidence case is the file setting's
    /// 2/m, so this preset is detected in the file's setting only (`conventional: false`).
    package static var mSq: CrystalStructure {
        CrystalStructure.detectingSymmetry(cell: .monoclinic(a: 5, b: 4, c: 5, betaDeg: acos(-0.25) * 180 / .pi),
                                           sites: alOrigin, conventional: false)
    }
    package static var mSqFull: CrystalStructure { mono(5, 4, 5, acos(-0.25) * 180 / .pi) }

    /// The rings of `s` after the |F| cut, DiffractoCore's families.
    package static func rings(_ s: CrystalStructure, kMax: Double = 1.2) throws -> [ReflectionFamily] {
        try s.kinematicReflections(kMax: kMax).families.map(\.family)
    }

    /// An orthorhombic cell whose rings at g = 1, 1.5, 2, 3, ... Å⁻¹ are exact doubles (g = h/2, k, 2 l).
    package static var dyadicCell: CrystalStructure {
        structure(.orthorhombic(a: 2, b: 1, c: 0.5), alOrigin)
    }
    /// The kit's rings at exactly these g (Å⁻¹), in the order given; a g with no ring is a failure of this fixture.
    package static func dyadic(_ gs: [Double]) throws -> [ReflectionFamily] {
        let all = try rings(dyadicCell, kMax: 4)
        return try gs.map { g in try #require(all.first { $0.g == g }, "the dyadic cell has no ring at exactly g = \(g)") }
    }
}
