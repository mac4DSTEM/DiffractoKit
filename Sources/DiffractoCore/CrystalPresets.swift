//
//  CrystalPresets.swift
//  Role: the preset crystals of the two applications as `CrystalStructure`s: the cubic and hexagonal single-element
//        families (fcc, bcc, simple cubic, diamond, hcp), the named elements (Au, Al, Ni, Cu, Fe, Si, Mg) at the lattice
//        constants the applications use, and two literature structures (beta'' Mg5Si6 and 2H-WS2). Ported from
//        mac4DSTEM `Crystal.swift` and macTEM `ReferenceStructures.swift` (their sites and constants agree; macTEM's
//        four are a subset). The symmetry is DETECTED from the sites (`CrystalStructure.detectingSymmetry`), as for a
//        P1 CIF. Checked against pymatgen and spglib in `CrystalPresetsTests` (volume, space group, d-spacings).
//

import Foundation
import simd

/// The preset crystals. A namespace: every member returns a `CrystalStructure` with its symmetry detected from the atoms.
public enum CrystalPresets {

    /// A cell and its atoms, symmetry detected from the atoms.
    public static func structure(a: Double, b: Double, c: Double,
                                 alphaDeg: Double = 90, betaDeg: Double = 90, gammaDeg: Double = 90,
                                 sites: [AtomSite], name: String? = nil) -> CrystalStructure {
        CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: a, b: b, c: c, alphaDeg: alphaDeg, betaDeg: betaDeg, gammaDeg: gammaDeg),
            sites: sites, name: name)
    }

    /// Face-centred cubic, one element (4 atoms), Fm-3m.
    public static func fcc(a: Double, z: Int) -> CrystalStructure {
        structure(a: a, b: a, c: a, sites: [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
            .map { AtomSite(z: z, fractional: $0) })
    }

    /// Body-centred cubic, one element (2 atoms), Im-3m.
    public static func bcc(a: Double, z: Int) -> CrystalStructure {
        structure(a: a, b: a, c: a, sites: [[0, 0, 0], [0.5, 0.5, 0.5]].map { AtomSite(z: z, fractional: $0) })
    }

    /// Simple cubic, one element (1 atom), Pm-3m.
    public static func simpleCubic(a: Double, z: Int) -> CrystalStructure {
        structure(a: a, b: a, c: a, sites: [AtomSite(z: z, fractional: [0, 0, 0])])
    }

    /// Diamond cubic, one element (8 atoms), Fd-3m: fcc plus the (1/4, 1/4, 1/4) shifted sublattice.
    public static func diamond(a: Double, z: Int) -> CrystalStructure {
        let base: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let shift = SIMD3<Double>(0.25, 0.25, 0.25)
        return structure(a: a, b: a, c: a, sites: (base + base.map { $0 + shift }).map { AtomSite(z: z, fractional: $0) })
    }

    /// Hexagonal close-packed, one element (2 atoms), P6_3/mmc, in the conventional 120 degree basis.
    public static func hcp(a: Double, c: Double, z: Int) -> CrystalStructure {
        structure(a: a, b: a, c: c, gammaDeg: 120, sites: [
            AtomSite(z: z, fractional: [0, 0, 0]),
            AtomSite(z: z, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5]),
        ])
    }

    /// The cubic families offered for user-defined crystals.
    public enum CubicFamily: String, Sendable, CaseIterable {
        case fcc = "FCC"
        case bcc = "BCC"
        case simpleCubic = "Simple cubic"
        case diamond = "Diamond"

        /// A stable lower-case token for provenance records ("fcc", "bcc", "simple_cubic", "diamond").
        public var provenanceID: String {
            switch self {
            case .fcc: "fcc"
            case .bcc: "bcc"
            case .simpleCubic: "simple_cubic"
            case .diamond: "diamond"
            }
        }
    }

    /// A single-element cubic crystal of the given family.
    public static func cubic(_ family: CubicFamily, a: Double, z: Int) -> CrystalStructure {
        switch family {
        case .fcc: return fcc(a: a, z: z)
        case .bcc: return bcc(a: a, z: z)
        case .simpleCubic: return simpleCubic(a: a, z: z)
        case .diamond: return diamond(a: a, z: z)
        }
    }

    // MARK: Named elements (lattice constants in Å)

    /// Aluminium, fcc, a = 4.0495 Å: pure Al at room temperature. A CIF carries its own a (a 4.04 Å CIF moves the {220}
    /// ring across a k_max of 0.70 Å⁻¹ on the dataset the choice was made for; mac4DSTEM ADR 040).
    public static var aluminum: CrystalStructure { fcc(a: 4.0495, z: 13) }
    /// Gold, fcc, a = 4.0782 Å.
    public static var gold: CrystalStructure { fcc(a: 4.0782, z: 79) }
    /// Nickel, fcc, a = 3.5240 Å.
    public static var nickel: CrystalStructure { fcc(a: 3.5240, z: 28) }
    /// Copper, fcc, a = 3.6149 Å.
    public static var copper: CrystalStructure { fcc(a: 3.6149, z: 29) }
    /// Iron (alpha), bcc, a = 2.8665 Å.
    public static var iron: CrystalStructure { bcc(a: 2.8665, z: 26) }
    /// Silicon, diamond cubic, a = 5.4309 Å.
    public static var silicon: CrystalStructure { diamond(a: 5.4309, z: 14) }
    /// Magnesium, hcp, a = 3.2094 Å, c = 5.2108 Å.
    public static var magnesium: CrystalStructure { hcp(a: 3.2094, c: 5.2108, z: 12) }

    // MARK: Literature structures

    /// beta'' (Mg5Si6), the hardening precipitate of Al-Mg-Si alloys: S. J. Andersen, H. W. Zandbergen, J. Jansen,
    /// C. Traeholt, U. Tundal, O. Reiso, Acta Materialia 46(9), 3283-3298 (1998), doi:10.1016/S1359-6454(97)00493-X.
    /// Cell from the abstract (a 15.16, b 4.05, c 6.74 Å, beta 105.3 degrees); coordinates from Table 3, set 3, the
    /// C2/m refinement. The six published sites are expanded here under C2/m (unique axis b) to Mg10Si12 = 22 atoms, the
    /// cell content the paper states. Deliberately not the DFT-relaxed Materials Project cell, whose volume is a few
    /// percent high (about 1 % on d-spacings). The needle direction is [010], coherent with Al <100>, b = 4.05 Å.
    public static var betaDoublePrime: CrystalStructure {
        func fourI(_ z: Int, _ x: Double, _ zc: Double) -> [AtomSite] {
            [SIMD3(x, 0, zc), SIMD3(-x, 0, -zc), SIMD3(x + 0.5, 0.5, zc), SIMD3(-x + 0.5, 0.5, -zc)]
                .map { AtomSite(z: z, fractional: $0) }
        }
        let mg = 12, si = 14
        var sites: [AtomSite] = [AtomSite(z: mg, fractional: [0, 0, 0]), AtomSite(z: mg, fractional: [0.5, 0.5, 0])]
        sites += fourI(mg, 0.3459, 0.089)
        sites += fourI(mg, 0.430, 0.652)
        sites += fourI(si, 0.0565, 0.649)
        sites += fourI(si, 0.1885, 0.224)
        sites += fourI(si, 0.2171, 0.617)
        return structure(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3, sites: sites)
    }

    /// 2H tungsten disulfide, P6_3/mmc (No. 194): W. J. Schutte, J. L. de Boer, F. Jellinek, J. Solid State Chem. 70,
    /// 207-209 (1987), doi:10.1016/0022-4596(87)90057-0, Table I: a = 3.1532 Å, c = 12.323 Å, z(S) = 0.6225; W at
    /// 2(c) +-(1/3, 2/3, 1/4), S at 4(f) +-(1/3, 2/3, z; 1/3, 2/3, 1/2 - z). The measured cell, not the DFT-PBE cell whose
    /// c is about 15 % long.
    public static var tungstenDisulfide: CrystalStructure {
        let zS = 0.6225
        return structure(a: 3.1532, b: 3.1532, c: 12.323, gammaDeg: 120, sites: [
            AtomSite(z: 74, fractional: [1.0 / 3.0, 2.0 / 3.0, 0.25]),
            AtomSite(z: 74, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.75]),
            AtomSite(z: 16, fractional: [1.0 / 3.0, 2.0 / 3.0, zS]),
            AtomSite(z: 16, fractional: [1.0 / 3.0, 2.0 / 3.0, 1.5 - zS]),
            AtomSite(z: 16, fractional: [2.0 / 3.0, 1.0 / 3.0, 1.0 - zS]),
            AtomSite(z: 16, fractional: [2.0 / 3.0, 1.0 / 3.0, zS - 0.5]),
        ])
    }
}
