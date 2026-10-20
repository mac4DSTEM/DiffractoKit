import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// The phase model: validation, orientation-class detection, provenance and the Materials Project route. Detection is
/// checked against pymatgen / spglib point groups (`Fixtures` of DiffractoCoreTests/presets-pins.json, tools/presets/pins.py):
/// whenever `orientationClass(of:)` names a class it is spglib's Laue class (exact), and the groups that return nil are the
/// ones whose standard operators are not in the kit's setting. Provenance, fingerprint and refusals are exact logic.
private struct OrientationPins: Decodable {
    struct Site: Decodable { let z: Int; let abc: [Double] }
    struct Lattice: Decodable { let a, b, c, alpha, beta, gamma: Double }
    struct Case: Decodable { let lattice: Lattice; let laue: String; let spaceGroupNumber: Int; let sites: [Site] }
    let orientation: [String: Case]
}

struct PhaseModelTests {
    private func structure(_ c: OrientationPins.Case) -> CrystalStructure {
        let l = c.lattice
        return CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: l.a, b: l.b, c: l.c, alphaDeg: l.alpha, betaDeg: l.beta, gammaDeg: l.gamma),
            sites: c.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.abc[0], $0.abc[1], $0.abc[2])) })
    }

    @Test func detectedOrientationClassIsSpglibsLaueClassOrNil() throws {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .appendingPathComponent("DiffractoCoreTests/Fixtures/presets-pins.json")
        let pins = try JSONDecoder().decode(OrientationPins.self, from: Data(contentsOf: file))
        #expect(pins.orientation.count == 15)
        var nils: [String] = []
        for (symbol, c) in pins.orientation.sorted(by: { $0.key < $1.key }) {
            let s = structure(c)
            #expect(s.laueClass?.symbol == c.laue, "\(symbol): kit detects \(String(describing: s.laueClass)), spglib \(c.laue)")
            if let found = PhaseModel.orientationClass(of: s) {
                #expect(found.symbol == c.laue, "\(symbol)")
            } else { nils.append(symbol) }
        }
        // Since the second -3m frame (-3m1: twofolds on a, P-3m1 and R-3m in hexagonal axes) all fifteen pinned groups are in a
        // standard setting of the kit, so none returns nil (before, the two -3m1 groups did).
        #expect(nils.isEmpty, "unsupported: \(nils)")
        for symbol in ["P-3m1", "R-3m"] {
            let c = try #require(pins.orientation[symbol])
            #expect(PhaseModel.orientationSymmetry(of: structure(c))?.threeBarMSetting == .threeBarMOne, "\(symbol)")
        }
        #expect(PhaseModel.orientationSymmetry(of: structure(try #require(pins.orientation["P-31m"])))?.threeBarMSetting == .threeBarOneM)
        print("PHASEMODEL orientation class not in the kit's standard setting: \(nils)")
    }

    @Test func presetsAndLibraryModelsKnowTheirOrientationClass() {
        let expected: [String: LaueClass] = ["au_fcc": .mThreeBarM, "al_fcc": .mThreeBarM, "ni_fcc": .mThreeBarM, "cu_fcc": .mThreeBarM,
            "fe_bcc": .mThreeBarM, "si_diamond": .mThreeBarM, "mg_hcp": .sixOverMmm, "ws2_2h": .sixOverMmm]
        #expect(PhaseModelLibrary.models.count == expected.count)
        for m in PhaseModelLibrary.models {
            #expect(m.orientationClass == expected[m.id], "\(m.id)")
            #expect(m.isUsable && m.supportsOrientationMapping && m.orientationMappingIssue == nil, "\(m.id)")
            #expect(PhaseModelLibrary.model(id: m.id)?.id == m.id)
        }
        #expect(PhaseModelLibrary.model(id: "nope") == nil)
        // beta'' is monoclinic b unique.
        let beta = PhaseModel(id: "beta", displayName: "beta''", structure: CrystalPresets.betaDoublePrime, source: .custom)
        #expect(beta.orientationClass == .twoOverM && beta.orientationFamilyName == "monoclinic" && beta.isUsable)
        let custom = PhaseModelLibrary.customCubic(family: .bcc, latticeA: 2.87, atomicNumber: 26)
        #expect(custom.id == "custom_cubic_bcc_z26" && custom.displayName == "Fe BCC, a = 2.87 Å" && custom.source == .custom)
    }

    @Test func aFaceCentredPrimitiveCellIsNotReducedInTheWrongFrame() {
        // Fm-3m written in its rhombohedral primitive cell: Laue class m-3m, but the cell is not the cubic frame.
        let a = 4.05
        let cell = UnitCell(a: a / 2.0.squareRoot(), b: a / 2.0.squareRoot(), c: a / 2.0.squareRoot(), alphaDeg: 60, betaDeg: 60, gammaDeg: 60)
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 13, fractional: [0, 0, 0])], conventional: false)
        let model = PhaseModel(id: "p", displayName: "p", structure: s, source: .imported)
        #expect(model.orientationClass == nil && !model.supportsOrientationMapping && model.orientationMappingIssue != nil)
        #expect(model.orientationFamilyName == "identity")
        // Still a usable crystal for phase identification.
        #expect(model.isUsable)
    }

    @Test func validationNamesEachProblem() {
        // A declared class that the cell does not have.
        let hexCell = CrystalPresets.magnesium
        let wrong = PhaseModel(id: "x", displayName: "x", structure: hexCell, orientationClass: .mThreeBarM, source: .custom)
        #expect(wrong.validationIssues.map(\.code) == ["cubic_symmetry_mismatch"] && !wrong.isUsable)
        let hexDeclared = PhaseModel(id: "x", displayName: "x", structure: CrystalPresets.gold, orientationClass: .threeBarM, source: .custom)
        #expect(hexDeclared.validationIssues.map(\.code) == ["hexagonal_symmetry_mismatch"])
        // An element with no electron scattering factors (Z = 104 is in the kit's table of elements, not in the scattering table).
        let heavy = CrystalPresets.fcc(a: 5, z: 104)
        let unsupported = PhaseModel(id: "u", displayName: "u", structure: heavy, source: .custom)
        #expect(unsupported.validationIssues.map(\.code).contains("unsupported_element"))
        // A cell with two atoms on top of each other.
        let crowded = CrystalPresets.structure(a: 4, b: 4, c: 4, sites: [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.001, 0, 0])])
        let model = PhaseModel(id: "c", displayName: "c", structure: crowded, source: .custom)
        #expect(!model.isUsable && model.shortestCloseContact != nil)
        // The unreduced declaration is not a validation issue.
        #expect(PhaseModel(id: "n", displayName: "n", structure: CrystalPresets.gold, orientationClass: nil, source: .custom).isUsable)
    }

    @Test func provenanceAndFingerprint() {
        let a = PhaseModelLibrary.model(id: "au_fcc")!
        let p = a.provenance
        #expect(p["crystal_symmetry"] == "cubic" && p["laue_class"] == "m-3m" && p["atomic_site_count"] == "4")
        #expect(p["crystal_model_source"] == "built_in" && p["cell_a_angstrom"] == "4.0782")
        // The fingerprint depends on what the model is, not on its id or name.
        let renamed = PhaseModel(id: "other", displayName: "other", structure: a.structure, source: .custom)
        #expect(renamed.contentFingerprint == a.contentFingerprint && a.contentFingerprint.count == 16)
        #expect(PhaseModelLibrary.model(id: "al_fcc")!.contentFingerprint != a.contentFingerprint)
        #expect(renamed.revisionID != a.revisionID)
        // Pinned digest of the gold model (FNV-1a over the canonical provenance string, recomputed independently in Python):
        // a change here is a change to the fingerprints written into session files.
        #expect(a.contentFingerprint == "cad6fd3edb39d0e5")
    }

    @Test func materialsProjectDocumentsBecomeModels() throws {
        let dir = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .appendingPathComponent("DiffractoCoreTests/Fixtures/mp-import")
        let when = Date(timeIntervalSince1970: 1_790_000_000)
        for (file, family, id) in [("input-aluminium-v3.json", "cubic", "mp-134"), ("input-siliconJSON.json", "cubic", "mp-149"),
                                   ("input-magnesiumJSON.json", "hexagonal", "mp-153")] {
            let doc = try MaterialsProjectImport.firstDocument(
                in: MaterialsProjectImport.decode(Data(contentsOf: dir.appendingPathComponent(file))))
            let model = try PhaseModel(materialsProject: doc, fetchedAt: when)
            #expect(model.id == id && model.source == .materialsProject && model.orientationFamilyName == family)
            #expect(model.displayName == "\(doc.formulaPretty!) (\(id))")
            #expect(model.spaceGroupNumber == doc.symmetry?.number)
            #expect(model.extraProvenance["materials_project_id"] == id)
            #expect(model.extraProvenance["materials_project_fetched"] == "2026-09-21T14:13:20Z")
            #expect(model.extraProvenance["materials_project_cell"] == "as served")
            #expect(model.provenance["materials_project_symbol"] == doc.symmetry?.symbol)
        }
        // A document that fails validation is refused with the issue text.
        let bad = MaterialsProjectDocument(
            materialID: "mp-5", structure: MaterialsProjectStructure(
                lattice: MaterialsProjectLattice(a: 4, b: 4, c: 4, alpha: 90, beta: 90, gamma: 90),
                sites: [MaterialsProjectSite(species: [MaterialsProjectSpecies(element: "Al", occu: 1)], abc: [0, 0, 0]),
                        MaterialsProjectSite(species: [MaterialsProjectSpecies(element: "Al", occu: 1)], abc: [0.01, 0, 0])]))
        #expect(throws: MaterialsProjectImport.Failure.self) { _ = try PhaseModel(materialsProject: bad, fetchedAt: when) }
        // Deprecated entries still import and say so.
        var deprecated = try MaterialsProjectImport.firstDocument(in: MaterialsProjectImport.decode(
            Data(contentsOf: dir.appendingPathComponent("input-aluminium-v3.json"))))
        deprecated = MaterialsProjectDocument(materialID: deprecated.materialID, formulaPretty: deprecated.formulaPretty, deprecated: true,
                                              structure: deprecated.structure, symmetry: deprecated.symmetry)
        #expect(try PhaseModel(materialsProject: deprecated, fetchedAt: when).extraProvenance["materials_project_deprecated"] == "true")
    }
}
