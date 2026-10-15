//
//  ConventionalChecks.swift
//  Role: the public API DiffractoKit slice 1b adds to DiffractoCore, used from outside the package with a real call and a
//        result a reader verifies by hand: the Niggli reduction and its change of basis (`niggliReduced`, `ReducedCell`),
//        the conventional cell and its Bravais lattice (`conventionalCell`, `ConventionalCell`, `BravaisLattice`), the
//        symmetry of a structure in any setting (`detectingSymmetry(conventional:tolerance:search:)`,
//        `CrystalStructure.setting`, `conventionalSetting`, `ConventionalSetting` and its four conversions) and the
//        detection options `SymmetryDetection.Tolerance` and `Search`. The cells are small and built in: fcc gold in its
//        primitive cell (192 operations in the conventional cell), a C-centred monoclinic cell in its primitive cell,
//        and rutile with its 4-fold axis along a. `runConventionalChecks` is called from main.swift.
//

import Foundation
import simd
import DiffractoCore

/// The largest absolute component of a vector.
private func maxAbs(_ v: SIMD3<Double>) -> Double { max(abs(v.x), abs(v.y), abs(v.z)) }

/// The primitive fcc cell of gold (a = 4.0782 Å): a_p = a / sqrt 2, all three angles 60 degrees, one atom.
private func goldPrimitive() -> CrystalStructure {
    CrystalStructure.detectingSymmetry(cell: .rhombohedral(a: 4.0782 / 2.0.squareRoot(), alphaDeg: 60),
                                       sites: [AtomSite(z: 79, fractional: [0, 0, 0])])
}

/// The primitive cell of a C-centred monoclinic lattice with conventional a, b, c and beta (unique axis b): the
/// vectors (a + b) / 2, (b - a) / 2 and c. Their lengths and angles are worked out here from the conventional cell.
private func monoclinicCPrimitive(a: Double, b: Double, c: Double, betaDeg: Double) -> UnitCell {
    let conventional = UnitCell.monoclinic(a: a, b: b, c: c, betaDeg: betaDeg)
    let v = conventional.realLattice
    let p1 = 0.5 * (v[0] + v[1]), p2 = 0.5 * (v[1] - v[0]), p3 = v[2]
    func angle(_ x: SIMD3<Double>, _ y: SIMD3<Double>) -> Double {
        acos((x * y).sum() / (length(x) * length(y))) * 180 / .pi
    }
    return UnitCell(a: length(p1), b: length(p2), c: length(p3),
                    alphaDeg: angle(p2, p3), betaDeg: angle(p1, p3), gammaDeg: angle(p1, p2))
}

/// Rutile (P4_2/mnm, a = 4.593, c = 2.959, O at u = 0.3053) with the 4-fold axis along the FILE's a: the file's
/// (x, y, z) is the conventional (z, x, y).
private func rutileAlongA() -> CrystalStructure {
    let u = 0.3053
    let conventional: [(Int, SIMD3<Double>)] = [
        (22, [0, 0, 0]), (22, [0.5, 0.5, 0.5]),
        (8, [u, u, 0]), (8, [-u, -u, 0]), (8, [0.5 + u, 0.5 - u, 0.5]), (8, [0.5 - u, 0.5 + u, 0.5]),
    ]
    return CrystalStructure.detectingSymmetry(
        cell: .orthorhombic(a: 2.959, b: 4.593, c: 4.593),
        sites: conventional.map { AtomSite(z: $0.0, fractional: SIMD3($0.1.z, $0.1.x, $0.1.y)) })
}

func runConventionalChecks(_ report: inout Report) {
    let aAu = 4.0782

    report.section("1b: Niggli reduction and its change of basis") { r in
        // fcc as a rhombohedral cell of a_p = a / sqrt 2 and 60 degrees is already Niggli-reduced.
        let primitive = UnitCell.rhombohedral(a: aAu / 2.0.squareRoot(), alphaDeg: 60)
        let reduced = primitive.niggliReduced()
        r.equal("default epsilon constant", UnitCell.niggliRelativeEpsilon, 1e-5)
        r.equal("the reduction converged", reduced.converged, true)
        r.equal("its steps are counted", reduced.steps >= 1, true)
        r.equal("P is unimodular with determinant +1", reduced.transformation.determinant, 1)
        r.near("epsilon is the relative epsilon times V^(2/3)", reduced.epsilon, 1e-5 * pow(primitive.volume, 2.0 / 3.0), 1e-12)
        r.near("the reduced volume is the lattice's", reduced.cell.volume, aAu * aAu * aAu / 4, 1e-9)
        r.near("the reduced lengths are a / sqrt 2", reduced.cell.a, aAu / 2.0.squareRoot(), 1e-9)
        r.near("(b)", reduced.cell.b, aAu / 2.0.squareRoot(), 1e-9)
        r.near("(c)", reduced.cell.c, aAu / 2.0.squareRoot(), 1e-9)
        r.near("the reduced angles are 60 degrees", reduced.cell.alphaDeg + reduced.cell.betaDeg + reduced.cell.gammaDeg, 180, 1e-7)

        // A deliberately skewed basis of the same lattice: a' = a, b' = a + b, c' = a + b + c (columns of P
        // in terms of the reduced basis). The reduction gets back the primitive lengths and a unimodular P.
        let v = primitive.realLattice
        let skew = [v[0], v[0] + v[1], v[0] + v[1] + v[2]]
        func angle(_ x: SIMD3<Double>, _ y: SIMD3<Double>) -> Double { acos((x * y).sum() / (length(x) * length(y))) * 180 / .pi }
        let skewed = UnitCell(a: length(skew[0]), b: length(skew[1]), c: length(skew[2]),
                              alphaDeg: angle(skew[1], skew[2]), betaDeg: angle(skew[0], skew[2]), gammaDeg: angle(skew[0], skew[1]))
        let back = skewed.niggliReduced()
        r.equal("a skewed basis of the same lattice converged", back.converged, true)
        r.near("has the same volume", back.cell.volume, skewed.volume, 1e-9)
        r.near("reduces to a = b = c = a / sqrt 2", back.cell.a + back.cell.b + back.cell.c, 3 * aAu / 2.0.squareRoot(), 1e-8)
        r.equal("with a unimodular change of basis", abs(back.transformation.determinant), 1)
        // x_input = P x_reduced: the reduced metric is P^T G P, so the reduced cell's c-axis length from P and the input.
        let cols = (0..<3).map { j in
            skewed.realLattice[0] * Double(back.transformation[0, j]) + skewed.realLattice[1] * Double(back.transformation[1, j])
                + skewed.realLattice[2] * Double(back.transformation[2, j])
        }
        r.near("P's columns are the reduced basis vectors in terms of the input", length(cols[0]), back.cell.a, 1e-9)
        r.near("(b axis)", length(cols[1]), back.cell.b, 1e-9)
        r.near("(c axis)", length(cols[2]), back.cell.c, 1e-9)
        let invalid = UnitCell(a: -1, b: 1, c: 1, alphaDeg: 90, betaDeg: 90, gammaDeg: 90).niggliReduced()
        r.equal("an invalid cell does not converge and keeps the identity", !invalid.converged && invalid.transformation == .identity
                && invalid.steps == 0, true)
        r.equal("a coarser epsilon is accepted", UnitCell.cubic(a: 4).niggliReduced(relativeEpsilon: 1e-3).converged, true)
        r.equal("ReducedCell is Equatable", reduced == primitive.niggliReduced() && reduced != skewed.niggliReduced(), true)
    }

    report.section("1b: the conventional cell of fcc gold") { r in
        let primitive = UnitCell.rhombohedral(a: aAu / 2.0.squareRoot(), alphaDeg: 60)
        let conventional = primitive.conventionalCell()
        r.equal("default tolerance constant", UnitCell.conventionalCellTolerance, 1e-9)
        r.equal("the lattice is cF", conventional.bravais, .cF)
        r.equal("its name", conventional.bravais.description, "cF")
        r.equal("its centring letter", conventional.bravais.centring, "F")
        r.equal("its holohedry has 48 operations", conventional.bravais.holohedryOrder, 48)
        r.equal("the cell holds four lattice points", conventional.bravais.centringMultiplicity, 4)
        r.equal("the transformation has determinant 4", conventional.transformation.determinant, 4)
        r.equal("and so does the multiplicity", conventional.centringMultiplicity, 4)
        r.near("the conventional cube edge is a", conventional.cell.a, aAu, 1e-9)
        r.near("(b)", conventional.cell.b, aAu, 1e-9)
        r.near("(c)", conventional.cell.c, aAu, 1e-9)
        r.equal("the angles are 90 degrees", conventional.cell.alphaDeg == 90 && conventional.cell.betaDeg == 90
                && conventional.cell.gammaDeg == 90, true)
        r.equal("the cell reads cubic", conventional.cell.metricFamily(), .cubic)
        r.equal("a cubic cell has no unique axis", conventional.uniqueAxis, nil)
        r.equal("the setting is stated", conventional.setting.isEmpty, false)
        r.near("the conventional volume is four primitive volumes", conventional.cell.volume, 4 * primitive.volume, 1e-8)
        r.equal("four centring translations, the origin first", conventional.centringTranslations.count == 4
                && conventional.centringTranslations[0] == SIMD3(0, 0, 0), true)
        let faceCentreNames: [String] = conventional.centringTranslations.map { (t: SIMD3<Double>) -> String in
            let coordinates: [Double] = [t.x, t.y, t.z]
            let doubled: [Int] = coordinates.map { (c: Double) -> Int in Int((c * 2).rounded()) }
            let digits: [String] = doubled.map { (i: Int) -> String in String(i) }
            return digits.joined(separator: "")
        }
        let faceCentres: Set<String> = Set(faceCentreNames)
        r.equal("they are the face centres of the cube", faceCentres, ["000", "011", "101", "110"])
        r.equal("the reduced cell it came from is the Niggli one", conventional.reduced == primitive.niggliReduced(), true)
        r.equal("an already conventional cube maps to itself", UnitCell.cubic(a: 4).conventionalCell().transformation == .identity
                && UnitCell.cubic(a: 4).conventionalCell().bravais == .cP, true)
        r.equal("a conventional bcc cube is cI", UnitCell.rhombohedral(a: 2.8665 * 3.0.squareRoot() / 2, alphaDeg: 109.47122063449069)
                .conventionalCell().bravais, .cI)
        r.equal("ConventionalCell is Equatable", conventional == primitive.conventionalCell(), true)
        r.equal("an invalid cell is returned as aP", UnitCell(a: 0, b: 1, c: 1, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
                .conventionalCell().bravais, .aP)
    }

    report.section("1b: BravaisLattice") { r in
        r.equal("fourteen lattices", BravaisLattice.allCases.count, 14)
        r.equal("raw values are the names", BravaisLattice.allCases.map(\.rawValue),
                ["aP", "mP", "mC", "oP", "oC", "oI", "oF", "tP", "tI", "hR", "hP", "cP", "cI", "cF"])
        r.equal("holohedry orders", BravaisLattice.allCases.map(\.holohedryOrder), [2, 4, 4, 8, 8, 8, 8, 16, 16, 12, 24, 48, 48, 48])
        r.equal("lattice points per cell", BravaisLattice.allCases.map(\.centringMultiplicity), [1, 1, 2, 1, 2, 2, 4, 1, 2, 3, 1, 1, 2, 4])
        r.equal("centring letters", String(BravaisLattice.allCases.map(\.centring)), "PPCPCIFPIRPPIF")
        r.equal("description prints the name", "\(BravaisLattice.hR)", "hR")
        r.equal("made from a raw value", BravaisLattice(rawValue: "oF"), .oF)
    }

    report.section("1b: a C-centred monoclinic primitive cell") { r in
        let a = 8.0, b = 5.0, c = 4.0, beta = 100.0
        let primitive = monoclinicCPrimitive(a: a, b: b, c: c, betaDeg: beta)
        let conventional = primitive.conventionalCell()
        r.equal("the lattice is mC", conventional.bravais, .mC)
        r.equal("the unique axis is b", conventional.uniqueAxis, .b)
        r.equal("two lattice points per cell", conventional.centringMultiplicity, 2)
        r.equal("the cell reads monoclinic along b", conventional.cell.metricFamily(), .monoclinic(uniqueAxis: .b))
        r.near("the conventional volume is twice the primitive one", conventional.cell.volume, 2 * primitive.volume, 1e-8)
        r.near("and a b c sin(beta) by the closed form", conventional.cell.volume, a * b * c * sin(beta * .pi / 180), 1e-8)
        r.near("b is the shortest lattice vector along the 2-fold axis", conventional.cell.b, b, 1e-9)
        r.near("the ac-plane lattice is a rectangle-like mesh with the same areas", conventional.cell.a * conventional.cell.c
               * sin(conventional.cell.betaDeg * .pi / 180), a * c * sin(beta * .pi / 180), 1e-8)
        r.equal("the angle beta is at least 90 degrees", conventional.cell.betaDeg >= 90, true)
        r.equal("the other two angles are 90 degrees", conventional.cell.alphaDeg == 90 && conventional.cell.gammaDeg == 90, true)
        r.equal("the second lattice point is the C-face centre", conventional.centringTranslations.count == 2
                && maxAbs(conventional.centringTranslations[1] - SIMD3(0.5, 0.5, 0)) < 1e-9, true)
        let lengths = [conventional.cell.a, conventional.cell.c].sorted()
        r.near("a and c are the mesh's two shortest vectors (4 and 8)", lengths[0] + lengths[1], a + c, 1e-9)
        // The reflections follow h_conv = h_input P; d does not change.
        let hkl = MillerIndex(1, 1, 1)
        let hConv = hkl.applying(conventional.transformation)
        r.near("d is the same in both cells", conventional.cell.d(hConv), primitive.d(hkl), 1e-9)

        // Structure: one Fe atom on every lattice point of the primitive cell.
        let s = CrystalStructure.detectingSymmetry(cell: primitive, sites: [AtomSite(z: 26, fractional: [0, 0, 0])])
        r.equal("the file setting is stated", s.setting, .file)
        r.equal("the conventional setting is carried", s.conventionalSetting?.bravais, .mC)
        r.equal("the structure is 2/m in the conventional cell", s.conventionalSetting?.laueClass, .twoOverM)
        r.equal("with 8 operations (4 point operations, 2 centrings)", s.conventionalSetting?.operations.count, 8)
        r.equal("and 4 in the file's primitive cell", s.operations.count, 4)
    }

    report.section("1b: fcc gold in its primitive cell, detected in both settings") { r in
        let s = goldPrimitive()
        r.equal("the structure keeps the file's cell", s.setting, .file)
        r.near("(cell)", s.cell.a, aAu / 2.0.squareRoot(), 1e-12)
        let c = try unwrap(s.conventionalSetting)
        r.equal("the conventional setting is cF", c.bravais, .cF)
        r.equal("192 operations in the conventional cell", c.operations.count, 192)
        r.equal("48 in the primitive file cell", s.operations.count, 48)
        r.equal("Laue class m-3m in both", s.laueClass == .mThreeBarM && c.laueClass == .mThreeBarM, true)
        r.equal("a Laue group of 48 rotations in both", s.laueGroup.order == 48 && c.laueGroup.order == 48, true)
        r.equal("the file's crystal system is cubic", s.crystalSystem, .cubic)
        r.near("the conventional cell is the cube", c.cell.a, aAu, 1e-9)
        r.equal("the shortcuts are the cell's own", c.cell == c.conventionalCell.cell && c.bravais == c.conventionalCell.bravais
                && c.transformation == c.conventionalCell.transformation, true)
        r.equal("four atoms in the conventional cell", c.sites.count, 4)
        r.equal("all of them gold", c.sites.allSatisfy { $0.z == 79 }, true)
        r.equal("the translation uncertainty is tiny", c.translationUncertainty < 1e-6, true)
        let conv = c.structure
        r.equal("the conventional structure says so", conv.setting, .conventional)
        r.equal("it holds the conventional cell, sites and operations", conv.cell == c.cell && conv.sites == c.sites
                && conv.operations == c.operations, true)
        r.equal("it has no conventional setting of its own", conv.conventionalSetting, nil)
        r.equal("with 192 operations and m-3m", conv.operations.count == 192 && conv.laueClass == .mThreeBarM, true)
        let list = try conv.reflections(kMax: 1.2)
        r.equal("its families are gold's eight", list.families().map(\.name),
                ["{1 1 1}", "{2 0 0}", "{2 2 0}", "{3 1 1}", "{2 2 2}", "{4 0 0}", "{3 3 1}", "{4 2 0}"])

        // Conversions between the settings.
        let identityBack = try unwrap(c.fileOperation(of: .identity))
        r.equal("the identity is the identity", identityBack.rotation.isIdentity && maxAbs(identityBack.translationMod1) < 1e-12, true)
        let centring = try unwrap(c.fileOperation(of: SymmetryOperation(rotation: .identity, translation: [0, 0.5, 0.5])))
        r.equal("a centring translation is a lattice translation of the file's cell", centring.rotation.isIdentity
                && maxAbs(centring.translationMod1) < 1e-12, true)
        let fourFold = try SymmetryOperation(xyz: "-y,x,z")
        let inFile = try unwrap(c.fileOperation(of: fourFold))
        r.equal("a 4-fold of the cube is an order-4 rotation of the file's cell", inFile.rotation.order, 4)
        r.equal("and is not the identity", inFile.rotation.isIdentity, false)
        let roundTrip = try unwrap(c.conventionalOperation(of: inFile))
        r.equal("converted back it is the cube's 4-fold", roundTrip.rotation == fourFold.rotation && maxAbs(roundTrip.translation) < 1e-12, true)
        r.equal("every detected conventional rotation has a file counterpart",
                c.operations.allSatisfy { c.fileOperation(of: $0) != nil }, true)

        let fileHKL = MillerIndex(1, 1, 1)
        let hConv = c.conventionalIndices(of: fileHKL)
        r.near("a reflection keeps its d in the other setting", c.cell.d(hConv), s.cell.d(fileHKL), 1e-10)
        r.equal("and returns to the file's indices", c.fileIndices(ofConventional: hConv), fileHKL)
        r.equal("a conventional (1 0 0) is a centring absence of the file's lattice", c.fileIndices(ofConventional: MillerIndex(1, 0, 0)), nil)
        r.equal("a conventional (1 1 1) is a lattice reflection", c.fileIndices(ofConventional: MillerIndex(1, 1, 1)) != nil, true)
        r.equal("a conventional (2 0 0) is too", c.fileIndices(ofConventional: MillerIndex(2, 0, 0)) != nil, true)
        r.equal("ConventionalSetting is Equatable", c == (goldPrimitive().conventionalSetting ?? c), true)
    }

    report.section("1b: rutile with the 4-fold axis along a") { r in
        let s = rutileAlongA()
        let c = try unwrap(s.conventionalSetting)
        r.equal("the file's cell is tetragonal along a", s.cell.metricFamily(), .tetragonal(uniqueAxis: .a))
        r.equal("its Laue class is 4/mmm", s.laueClass, .fourOverMmm)
        r.equal("with 16 operations in the file's setting", s.operations.count, 16)
        r.equal("and 16 in the conventional one (a primitive lattice, nothing added)", c.operations.count, 16)
        r.equal("the lattice is tP", c.bravais, .tP)
        r.equal("the conventional unique axis is c", c.conventionalCell.uniqueAxis, .c)
        r.equal("the conventional cell reads tetragonal along c", c.cell.metricFamily(), .tetragonal(uniqueAxis: .c))
        r.near("with a = 4.593", c.cell.a, 4.593, 1e-9)
        r.near("and c = 2.959", c.cell.c, 2.959, 1e-9)
        r.equal("the change of basis is a permutation of the axes", abs(c.transformation.determinant) == 1
                && c.transformation != .identity, true)
        let along = c.conventionalIndices(of: MillerIndex(1, 0, 0))
        r.equal("the file's (1 0 0) is the conventional (0 0 1) up to sign", along.h == 0 && along.k == 0 && abs(along.l) == 1, true)
        r.near("a plane keeps its spacing in the conventional cell", c.cell.d(c.conventionalIndices(of: MillerIndex(1, 1, 0))),
               s.cell.d(MillerIndex(1, 1, 0)), 1e-10)
        r.equal("the structure itself still says the file's setting", s.setting, .file)
        let conv = c.structure
        r.equal("the conventional structure says so", conv.setting, .conventional)
        let list = try conv.reflections(kMax: 0.5)
        r.near("the first allowed ring of rutile is {1 1 0}, d = a / sqrt 2", list.allowed.first?.d ?? 0, 4.593 / 2.0.squareRoot(), 1e-9)
        r.equal("with the operations of the file's setting the file's own reflections agree", try s.reflections(kMax: 0.5)
                .allowed.map { $0.d }.first.map { abs($0 - 4.593 / 2.0.squareRoot()) < 1e-9 }, true)
    }

    report.section("1b: detectingSymmetry(conventional:tolerance:search:) and the options") { r in
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let sites = fcc.map { AtomSite(z: 79, fractional: $0) }
        let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: aAu), sites: sites)
        r.equal("the default asks for the conventional cell", gold.conventionalSetting != nil, true)
        r.equal("the cube is its own conventional cell", gold.conventionalSetting?.transformation, .identity)
        r.equal("192 operations", gold.operations.count, 192)

        let file = CrystalStructure.detectingSymmetry(cell: .cubic(a: aAu), sites: sites, conventional: false)
        r.equal("conventional: false carries no conventional setting", file.conventionalSetting, nil)
        r.equal("and detects the same operations", file.operations, gold.operations)
        let brute = CrystalStructure.detectingSymmetry(cell: .cubic(a: aAu), sites: sites, search: .bruteForce)
        r.equal("the brute-force scan finds the cell list's operations bit for bit", brute.operations, gold.operations)
        let cartesian = CrystalStructure.detectingSymmetry(cell: .cubic(a: aAu), sites: sites, tolerance: .cartesian(1e-3))
        r.equal("a Cartesian tolerance of 1e-3 A finds them too", cartesian.operations.count, 192)
        let all = CrystalStructure.detectingSymmetry(cell: .cubic(a: aAu), sites: sites, tolerance: .fractional(5e-3),
                                                     conventional: false, search: .cellList)
        r.equal("every parameter spelled out", all.operations, gold.operations)

        // Tolerance: a body-centred cell of 30 A with the centre atom moved 0.0004 * 30 = 0.012 A. Fractionally that is
        // within 5e-3 (bcc, 96 operations); in spglib's Cartesian sense it is 12 times 1e-3 A, so the centring is lost.
        let big = UnitCell.cubic(a: 30)
        let displaced = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5004, 0.5, 0.5])]
        let loose = CrystalStructure.detectingSymmetry(cell: big, sites: displaced, tolerance: .fractional(5e-3), conventional: false)
        let tight = CrystalStructure.detectingSymmetry(cell: big, sites: displaced, tolerance: .cartesian(1e-3), conventional: false)
        r.equal("0.012 A in a 30 A cell is symmetric at the fractional 5e-3", loose.operations.count, 96)
        r.equal("and is not at the Cartesian 1e-3 A", tight.operations.count < 96, true)
        r.equal("the proposed Cartesian tolerance is 1e-3 A", SymmetryDetection.Tolerance.proposedCartesian, .cartesian(1e-3))
        r.equal("it differs from the fractional default", SymmetryDetection.Tolerance.proposedCartesian != .default, true)
        r.equal("the default is slice 1's fractional 5e-3", SymmetryDetection.Tolerance.default, .fractional(SymmetryDetection.positionTolerance))
        r.equal("(its value)", SymmetryDetection.Tolerance.default, .fractional(5e-3))
        let proposed = CrystalStructure.detectingSymmetry(cell: big, sites: displaced, tolerance: .proposedCartesian, conventional: false)
        r.equal("detecting with the proposed tolerance equals the .cartesian(1e-3) call", proposed.operations, tight.operations)

        // SymmetryDetection.detect takes the same two options and agrees with the structure's.
        let bruteDetect = SymmetryDetection.detect(sites: sites, cell: .cubic(a: aAu), tolerance: .default, search: .bruteForce)
        let listDetect = SymmetryDetection.detect(sites: sites, cell: .cubic(a: aAu), tolerance: .default, search: .cellList)
        r.equal("detect(tolerance:search:) is the same either way", bruteDetect.operations == listDetect.operations
                && bruteDetect.operations.count == 192, true)
        r.equal("the two searches are different values", SymmetryDetection.Search.cellList != SymmetryDetection.Search.bruteForce, true)
        r.equal("a Cartesian detect with a tolerance below the half cell works", SymmetryDetection.detect(
            sites: sites, cell: .cubic(a: aAu), tolerance: .cartesian(1e-3)).operations.count, 192)
        r.equal("Setting values differ", CrystalStructure.Setting.file != CrystalStructure.Setting.conventional, true)
    }
}
