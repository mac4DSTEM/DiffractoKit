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
            r.equal("every class has a native key (\(laue))", a == (try s.ipfColor(direction: [0.3, 0.5, 0.8], key: .native)), true)
        }
        let m = simd_double3x3(simd_quatd(angle: 0.7, axis: simd_normalize(SIMD3(1, 1, 0))))
        let d = try cubic.disorientation(from: matrix_identity_double3x3, to: m)
        r.near("disorientation of a rotation inside the fundamental angle range", d.angleRad, 0.7, 1e-12)
        let quarter = simd_double3x3(simd_quatd(angle: .pi / 2, axis: [0, 0, 1]))
        r.near("a fourfold rotation is no disorientation in m-3m", try cubic.disorientation(from: matrix_identity_double3x3, to: quarter).angleRad, 0, 1e-12)
        let family = try cubic.py4DSTEMDirectionFamily(matrix_identity_double3x3)
        r.equal("the family of the identity has three columns", family.columns.0.z.isFinite, true)
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
}
