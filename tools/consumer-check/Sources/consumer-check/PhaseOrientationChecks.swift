//
//  PhaseOrientationChecks.swift
//  Role: the 0.6.0 public API used from outside the package: the orix IPF key, disorientation, the py4DSTEM direction
//        family, the preset crystals, the Materials Project models and the phase model. `runPhaseOrientationChecks` is
//        called from main.swift. Values are ones a reader can verify by hand.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

func runPhaseOrientationChecks(_ report: inout Report) {
    report.section("orientation 0.6.0") { r in
        let cubic = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        let c001 = try cubic.ipfColor(direction: [0, 0, 1], key: .orix)
        r.check("the orix cubic key is red at [001]", c001.x > 0.99 && c001.y < 0.01 && c001.z < 0.01, "\(c001)")
        for laue in LaueClass.allCases {
            let s = CrystalOrientationSymmetry(laueClass: laue)
            let a = try s.ipfColor(direction: [0.3, 0.5, 0.8], key: .native)
            if ![.mThreeBarM, .sixOverMmm, .twoOverM].contains(laue) { r.equal("the native key of \(laue) is the orix key", a, try s.ipfColor(direction: [0.3, 0.5, 0.8], key: .orix)) }
        }
        let m = simd_double3x3(simd_quatd(angle: 0.7, axis: simd_normalize(SIMD3(1, 1, 0))))
        let d = try cubic.disorientation(from: matrix_identity_double3x3, to: m)
        r.near("disorientation of a rotation inside the fundamental angle range", d.angleRad, 0.7, 1e-12)
        let quarter = simd_double3x3(simd_quatd(angle: .pi / 2, axis: [0, 0, 1]))
        r.near("a fourfold rotation is no disorientation in m-3m", try cubic.disorientation(from: matrix_identity_double3x3, to: quarter).angleRad, 0, 1e-12)
        let family = try cubic.py4DSTEMDirectionFamily(matrix_identity_double3x3)
        r.equal("the family of the identity is a signed permutation of the axes", family.columns.0.x.magnitude + family.columns.0.y.magnitude + family.columns.0.z.magnitude, 1)
    }
    report.section("presets, Materials Project, phase model") { r in
        r.equal("gold has four atoms", CrystalPresets.gold.sites.count, 4)
        r.equal("gold is m-3m", CrystalPresets.gold.laueClass, .mThreeBarM)
        r.equal("beta'' has Mg10Si12", CrystalPresets.betaDoublePrime.sites.count, 22)
        r.equal("cubic families", CrystalPresets.CubicFamily.allCases.count, 4)
        r.equal("formula parser", MaterialsProjectImport.parseFormula("Mg5Si6"), ["Mg": 5, "Si": 6])
        let doc = MaterialsProjectDocument(materialID: "mp-134", formulaPretty: "Al", symmetry: MaterialsProjectSymmetry(symbol: "Fm-3m", number: 225))
        r.equal("verdict matches", MaterialsProjectImport.check(doc, against: PhaseExpectation(formula: "Al", spaceGroupNumber: 225)), .matches)
        r.equal("verdict unchecked", MaterialsProjectImport.check(doc, against: PhaseExpectation()), .unchecked)
        let model = try unwrap(PhaseModelLibrary.model(id: "mg_hcp"))
        r.equal("magnesium orientation class", model.orientationClass, .sixOverMmm)
        r.equal("magnesium is usable", model.isUsable, true)
        r.equal("model source", model.source, .builtIn)
        r.equal("fingerprint length", model.contentFingerprint.count, 16)
        let issue = PhaseModelIssue(code: "x", message: "y")
        r.equal("issue", issue.description, "y")
    }
    report.section("second -3m setting, family order, Metal backend") { r in
        r.equal("two -3m settings", ThreeBarMSetting.allCases.map(\.rawValue), ["-31m", "-3m1"])
        let between = CrystalOrientationSymmetry(laueClass: .threeBarM)
        let onAxes = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarMOne)
        r.equal("the default -3m setting is -31m", between.threeBarMSetting, .threeBarOneM)
        r.equal("-3m1 reports its setting", onAxes.threeBarMSetting, .threeBarMOne)
        r.equal("other classes have no setting", CrystalOrientationSymmetry(laueClass: .mThreeBarM, setting: .threeBarMOne).threeBarMSetting, nil)
        let twofoldOnX = simd_double3x3(diagonal: [1, -1, -1])
        func has(_ s: CrystalOrientationSymmetry, _ m: simd_double3x3) -> Bool {
            s.operators.contains { op in (0..<3).allSatisfy { c in (0..<3).allSatisfy { abs(op[c][$0] - m[c][$0]) < 1e-12 } } }
        }
        r.equal("-3m1 has the twofold on a (x)", has(onAxes, twofoldOnX), true)
        r.equal("-31m has no twofold on a", has(between, twofoldOnX), false)
        r.equal("both settings have six proper operators", between.operators.count == 6 && onAxes.operators.count == 6, true)
        // CdI2 type, P-3m1: Cd at the origin, I at (1/3, 2/3, z) and (2/3, 1/3, -z).
        let cdi2 = CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: 4.24, b: 4.24, c: 6.84, alphaDeg: 90, betaDeg: 90, gammaDeg: 120),
            sites: [AtomSite(z: 48, fractional: [0, 0, 0]), AtomSite(z: 53, fractional: [1.0 / 3, 2.0 / 3, 0.25]),
                    AtomSite(z: 53, fractional: [2.0 / 3, 1.0 / 3, -0.25])])
        r.equal("CdI2 (P-3m1) names the -3m1 setting", PhaseModel.orientationSymmetry(of: cdi2)?.threeBarMSetting, .threeBarMOne)
        r.equal("CdI2 has orientation class -3m", PhaseModel.orientationClass(of: cdi2), .threeBarM)
        let model = PhaseModel(id: "cdi2", displayName: "CdI2", structure: cdi2, source: .custom)
        r.equal("the model carries the setting", model.orientationSetting, .threeBarMOne)
        r.equal("the model's symmetry is the -3m1 one", model.orientationSymmetry?.operators == onAxes.operators, true)
        r.equal("the explicit initialiser takes the setting", PhaseModel(id: "x", displayName: "x", structure: cdi2, orientationClass: .threeBarM,
                                                                      orientationSetting: .threeBarMOne, source: .custom).validationIssues.map(\.code).contains("hexagonal_symmetry_mismatch"), false)
        let cubic = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        r.equal("two operator orders", Py4DSTEMOperatorOrder.allCases, [.pymatgen, .sortedFractional])
        for order in Py4DSTEMOperatorOrder.allCases {
            let f = try cubic.py4DSTEMDirectionFamily(matrix_identity_double3x3, order: order)
            let l1 = [f.columns.0, f.columns.1, f.columns.2].map { $0.x.magnitude + $0.y.magnitude + $0.z.magnitude }
            r.equal("each cubic family column of the identity is a signed axis for \(order)", l1, [1, 1, 1])
        }
        r.equal("the default order is pymatgen", try cubic.py4DSTEMDirectionFamily(matrix_identity_double3x3) == cubic.py4DSTEMDirectionFamily(matrix_identity_double3x3, order: .pymatgen), true)
        r.equal("two backends", ACOMBackend.allCases, [.cpu, .metal])
    }
}
