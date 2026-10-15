//
//  LaueFollowsStructureTests.swift
//  Role: the Laue class follows the structure's symmetry (the operations a file gives, and the atoms), with the lattice
//        metric only bounding it: a primitive 60-degree cell of MgO with the 48 operations of m-3m written in that basis is
//        m-3m; hcp written in a gamma = 60 degree cell is 6/mmm, its true class.
//

import Testing
import DiffractoCore

@Suite struct LaueFollowsStructureTests {

    /// The 48 cubic rotations (with inversion) of the F-centred cube written in the primitive 60-degree basis, as xyz strings.
    private func primitiveCubicOperations() -> [SymmetryOperation] { primitiveCubicOperationTexts().map { try! SymmetryOperation(xyz: $0) } }

    private func primitiveCubicOperationTexts() -> [String] {
        let ainv: [[Int]] = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        let a: [[Double]] = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
        var ops: [String] = []
        for p in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
            for s in 0..<8 {
                var P = [[Double]](repeating: [0, 0, 0], count: 3)
                for i in 0..<3 { P[i][p[i]] = (s >> i) & 1 == 1 ? -1 : 1 }
                var pa = [[Double]](repeating: [0, 0, 0], count: 3)
                for i in 0..<3 { for j in 0..<3 { for k in 0..<3 { pa[i][j] += P[i][k] * a[k][j] } } }
                var m = [[Int]](repeating: [0, 0, 0], count: 3)
                for i in 0..<3 { for j in 0..<3 {
                    var v = 0.0
                    for k in 0..<3 { v += Double(ainv[i][k]) * pa[k][j] }
                    m[i][j] = Int(v.rounded())
                } }
                let text = m.map { row in
                    ["x", "y", "z"].enumerated().compactMap { (j, n) -> String? in
                        row[j] == 0 ? nil : (row[j] == 1 ? "+\(n)" : row[j] == -1 ? "-\(n)" : "\(row[j] > 0 ? "+" : "")\(row[j])\(n)")
                    }.joined()
                }.joined(separator: ",")
                ops.append(text)
            }
        }
        return ops
    }

    @Test func primitiveRhombohedralCellOfMgOKeepsItsCubicOperations() {
        let ac = 4.212 / 2.0.squareRoot()
        let cell = UnitCell(a: ac, b: ac, c: ac, alphaDeg: 60, betaDeg: 60, gammaDeg: 60)
        let sites = [AtomSite(z: 12, fractional: SIMD3(0, 0, 0)), AtomSite(z: 8, fractional: SIMD3(0.5, 0.5, 0.5))]
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, fileOperations: primitiveCubicOperations())
        #expect(s.laueClass == .mThreeBarM)
    }

    @Test func hcpInAGammaSixtyCellIsHexagonal() {
        let cell = UnitCell(a: 3.2094, b: 3.2094, c: 5.2108, alphaDeg: 90, betaDeg: 90, gammaDeg: 60)
        let sites = [AtomSite(z: 12, fractional: SIMD3(0, 0, 0)), AtomSite(z: 12, fractional: SIMD3(1.0 / 3, 1.0 / 3, 0.5))]
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        #expect(s.laueClass == .sixOverMmm)
    }

    private func cif(name: String, a: Double, b: Double, c: Double, alpha: Double = 90, beta: Double = 90, gamma: Double = 90,
                     ops: [String], atoms: [(String, Double, Double, Double)]) -> String {
        var s = "data_\(name)\n_cell_length_a \(a)\n_cell_length_b \(b)\n_cell_length_c \(c)\n"
        s += "_cell_angle_alpha \(alpha)\n_cell_angle_beta \(beta)\n_cell_angle_gamma \(gamma)\n"
        if !ops.isEmpty { s += "loop_\n_symmetry_equiv_pos_as_xyz\n" + ops.map { "'\($0)'" }.joined(separator: "\n") + "\n" }
        s += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
        for (i, at) in atoms.enumerated() { s += "\(at.0)\(i) \(at.0) \(at.1) \(at.2) \(at.3) 1\n" }
        return s
    }

    /// The same cell read from a CIF, the way an application reads it.
    @Test func primitiveRhombohedralMgOFromACIFKeepsItsCubicClass() throws {
        let ac = 4.212 / 2.0.squareRoot()
        let ops = primitiveCubicOperationTexts()
        let text = cif(name: "MgOprim", a: ac, b: ac, c: ac, alpha: 60, beta: 60, gamma: 60, ops: ops,
                       atoms: [("Mg", 0, 0, 0), ("O", 0.5, 0.5, 0.5)])
        let doc = try CIFImport.read(text, fileBaseName: "mgo_prim")
        #expect(doc.structure.laueClass == .mThreeBarM)
    }

    /// A file that lists the 48 signed permutations of x, y, z literally (written in the primitive basis) for the
    /// 60-degree cell: not all of them are symmetries of that lattice, but the atoms and the lattice (an fcc lattice)
    /// have the full cubic group, so the class is m-3m.
    @Test func literalSignedPermutationsInAPrimitiveCellStillGiveTheStructuresClass() throws {
        var ops: [String] = []
        for p in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
            for s in 0..<8 {
                let names = ["x", "y", "z"]
                ops.append((0..<3).map { axis in (s >> axis) & 1 == 1 ? "-\(names[p[axis]])" : names[p[axis]] }.joined(separator: ","))
            }
        }
        for a in [4.212, 4.212 / 2.0.squareRoot()] {
        let text = cif(name: "MgOliteral", a: a, b: a, c: a, alpha: 60, beta: 60, gamma: 60, ops: ops,
                       atoms: [("Mg", 0, 0, 0), ("O", 0.5, 0.5, 0.5)])
        let doc = try CIFImport.read(text, fileBaseName: "t")
        #expect(doc.structure.laueClass == .mThreeBarM, "a = \(a)")
        }
    }
}
