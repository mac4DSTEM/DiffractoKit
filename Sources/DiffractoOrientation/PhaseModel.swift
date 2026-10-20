//
//  PhaseModel.swift
//  Role: a complete phase definition for orientation mapping: a `CrystalStructure` (lattice, atoms, structure factors)
//        together with the Laue class its orientation symmetry works in, and the checks that keep the two consistent.
//        Ported from mac4DSTEM `CrystalModel.swift` and the model route of `MaterialsProjectImport.crystalModel`. The
//        application's three symmetry families (cubic, hexagonal, monoclinic, identity) are generalised to the kit's
//        eleven Laue classes (`CrystalOrientationSymmetry`); the application's `ACOMCrystalSymmetry` is not here.
//

import Foundation
import DiffractoCore
import DiffractoScattering

/// Origin of a phase definition. This describes the model, not the dataset; no phase is inferred from a file name or
/// from the look of a diffraction pattern.
public enum PhaseModelSource: String, Sendable, CaseIterable {
    case builtIn = "built_in"
    case custom = "custom"
    case imported = "imported"
    case materialsProject = "materials_project"
}

/// One reason a phase model is not usable. `code` is stable; `message` is for people.
public struct PhaseModelIssue: Sendable, Equatable, CustomStringConvertible {
    public let code: String
    public let message: String
    public init(code: String, message: String) {
        self.code = code
        self.message = message
    }
    public var description: String { message }
}

/// A phase model: the structure supplies lattice, atomic basis and structure factors; `orientationClass` names the Laue
/// class in which orientations of this phase are reduced and coloured. Keeping them together avoids applying one
/// crystal's presentation to a different lattice.
public struct PhaseModel: Identifiable, Sendable {
    public let id: String
    public let displayName: String
    public let structure: CrystalStructure
    /// The Laue class of the orientation symmetry, in the kit's standard setting (`CrystalOrientationSymmetry`), or nil
    /// when the phase can be identified but cannot be orientation-mapped ("Unreduced").
    public let orientationClass: LaueClass?
    public let source: PhaseModelSource
    /// International Tables space-group number when known; nil means "not recorded", not "P1".
    public let spaceGroupNumber: Int?
    /// Source facts beyond the fixed provenance keys (a Materials Project id, the fetch date, the reported symbol...).
    public let extraProvenance: [String: String]

    public init(id: String, displayName: String, structure: CrystalStructure, orientationClass: LaueClass?,
                source: PhaseModelSource, spaceGroupNumber: Int? = nil, extraProvenance: [String: String] = [:]) {
        self.id = id
        self.displayName = displayName
        self.structure = structure
        self.orientationClass = orientationClass
        self.source = source
        self.spaceGroupNumber = spaceGroupNumber
        self.extraProvenance = extraProvenance
    }

    /// A model whose orientation class is detected from the structure (`orientationClass(of:)`).
    public init(id: String, displayName: String, structure: CrystalStructure, source: PhaseModelSource,
                spaceGroupNumber: Int? = nil, extraProvenance: [String: String] = [:]) {
        self.init(id: id, displayName: displayName, structure: structure,
                  orientationClass: Self.orientationClass(of: structure), source: source,
                  spaceGroupNumber: spaceGroupNumber, extraProvenance: extraProvenance)
    }

    /// The orientation symmetry of this phase, or nil when it is unreduced.
    public var orientationSymmetry: CrystalOrientationSymmetry? { orientationClass.map { CrystalOrientationSymmetry(laueClass: $0) } }

    /// The Laue class the orientation reduction can work in for `structure`: the class detected from the atoms, when the
    /// structure's own proper Laue operators, written in its cell basis, are exactly the class's standard operators
    /// (`CrystalOrientationSymmetry.operators` documents the settings: cubic and orthorhombic axes, monoclinic b unique,
    /// tetragonal c, hexagonal a = x with the pymatgen -3m twofold). A cell in another setting, for example a face-centred
    /// structure in its primitive rhombohedral cell or a -3m1 cell whose twofolds lie on a, returns nil rather than a
    /// fundamental zone drawn in the wrong frame. (The application it came from supported m-3m, 6/mmm and 2/m only.)
    public static func orientationClass(of structure: CrystalStructure) -> LaueClass? {
        guard let laue = structure.laueClass else { return nil }
        let proper = Set(structure.laueGroup.rotations.filter { $0.determinant == 1 })
        guard proper == Set(CrystalOrientationSymmetry(laueClass: laue).fractionalRotations) else { return nil }
        return laue
    }

    /// The application-era family name of the orientation class ("cubic", "hexagonal", "monoclinic", "identity" for
    /// unreduced), kept for provenance records written before the kit; the other classes use their symbol.
    public var orientationFamilyName: String {
        switch orientationClass {
        case .mThreeBarM?: "cubic"
        case .sixOverMmm?: "hexagonal"
        case .twoOverM?: "monoclinic"
        case nil: "identity"
        case let other?: other.symbol
        }
    }

    // MARK: Validation

    /// Why the model is not usable, empty when it is. The structure's own checks are the kit's (`CrystalStructure.validationIssues`:
    /// cell lengths, angles and volume, a non-empty basis, occupancy, finite coordinates, close contacts); on top of them
    /// a site whose element has no electron scattering factors is `unsupported_element` (Z outside 1 to 103), and a
    /// declared orientation class that the structure does not have in the standard setting is `<system>_symmetry_mismatch`.
    /// An over-populated cell is the one modelling error that produces no visible symptom (site multiplicity multiplies
    /// straight into the structure factors), which is why the close contact check exists.
    public var validationIssues: [PhaseModelIssue] {
        var issues: [PhaseModelIssue] = []
        for issue in structure.validationIssues where issue.code != "invalid_atomic_number" {
            issues.append(PhaseModelIssue(code: issue.code, message: issue.message))
        }
        for site in structure.sites where site.z < 1 || !ElectronScatteringFactors.isSupported(z: site.z) {
            issues.append(PhaseModelIssue(code: "unsupported_element",
                                          message: "No electron scattering factors are available for Z = \(site.z)."))
        }
        if let declared = orientationClass, Self.orientationClass(of: structure) != declared {
            let system: String
            switch declared.crystalSystem {
            case .trigonal, .hexagonal: system = "hexagonal"
            case let other: system = other.rawValue
            }
            issues.append(PhaseModelIssue(
                code: "\(system)_symmetry_mismatch",
                message: "Laue class \(declared.symbol) is declared, but the structure's Laue group in its cell is not that class's standard setting."))
        }
        return issues
    }

    public var isUsable: Bool { validationIssues.isEmpty }

    /// Shortest centre-to-centre distance (Å) between two sites that cannot both be occupied by a whole atom, when under
    /// `CrystalStructure.closeContactLimit`.
    public var shortestCloseContact: Double? { structure.shortestCloseContact }

    /// Whether this model may drive orientation mapping. Separate from `isUsable` on purpose: an unreduced model is a good
    /// crystal for phase identification, but mapping one anyway would publish an unreduced orientation and an IPF colour
    /// that is no key at all. It is not a validation issue, because a validation issue would refuse the structure at import.
    public var supportsOrientationMapping: Bool { orientationClass != nil }

    /// Why orientation mapping is unavailable, or nil when it is available.
    public var orientationMappingIssue: String? {
        supportsOrientationMapping ? nil
            : "This phase's symmetry could not be reduced in a standard setting (its Laue group is undetermined, or is not that class's "
              + "standard setting in this cell). The structure can still be used to identify the phase, but it cannot produce an "
              + "orientation map or an IPF colour."
    }

    // MARK: Identity and provenance

    /// In-memory identity, to reject a completion when editable values changed while a plan was being generated.
    public var revisionID: String {
        let cell = structure.cell
        let sites = structure.sites.map {
            "\($0.z):\($0.fractional.x),\($0.fractional.y),\($0.fractional.z):\($0.occupancy)"
        }.joined(separator: ";")
        return "\(id)|\(orientationFamilyName)|\(cell.a),\(cell.b),\(cell.c)|\(cell.alphaDeg),\(cell.betaDeg),\(cell.gammaDeg)|\(sites)"
    }

    /// Immutable description to attach to products made with this phase.
    public var provenance: [String: String] {
        let cell = structure.cell
        let atomicBasis = structure.sites.map {
            "Z\($0.z)@(\($0.fractional.x),\($0.fractional.y),\($0.fractional.z));occ=\($0.occupancy)"
        }.joined(separator: "|")
        var result: [String: String] = [
            "crystal_model_id": id,
            "crystal_model_name": displayName,
            "crystal_model_source": source.rawValue,
            "crystal_symmetry": orientationFamilyName,
            "cell_a_angstrom": String(cell.a),
            "cell_b_angstrom": String(cell.b),
            "cell_c_angstrom": String(cell.c),
            "cell_alpha_degrees": String(cell.alphaDeg),
            "cell_beta_degrees": String(cell.betaDeg),
            "cell_gamma_degrees": String(cell.gammaDeg),
            "cell_volume_angstrom_cubed": String(cell.volume),
            "atomic_site_count": String(structure.sites.count),
            "atomic_basis": atomicBasis,
        ]
        if let orientationClass { result["laue_class"] = orientationClass.symbol }
        if let spaceGroupNumber { result["space_group_number"] = String(spaceGroupNumber) }
        for (key, value) in extraProvenance { result[key] = value }
        return result
    }

    /// A stable digest of what the model IS (cell, symmetry, atomic basis), independent of its id and name: two different
    /// CIFs with the same file name share an id, and a recipe that resolved its material by id alone would replay against
    /// the wrong crystal. FNV-1a over the provenance strings, not `Hasher` (seeded per process), because the value is
    /// written into session files.
    public var contentFingerprint: String {
        let keys = ["crystal_symmetry", "cell_a_angstrom", "cell_b_angstrom", "cell_c_angstrom",
                    "cell_alpha_degrees", "cell_beta_degrees", "cell_gamma_degrees", "atomic_basis"]
        let provenance = self.provenance
        let canonical = keys.map { "\($0)=\(provenance[$0] ?? "")" }.joined(separator: "\n")
        var hash: UInt64 = 0xcbf29ce484222325
        for byte in canonical.utf8 {
            hash ^= UInt64(byte)
            hash = hash &* 0x100000001b3
        }
        return String(format: "%016llx", hash)
    }
}

// MARK: - Built-in models

/// The built-in phase models, resolvable by id: a recipe recorded against a built-in id must keep resolving.
public enum PhaseModelLibrary {
    public static let models: [PhaseModel] = [
        PhaseModel(id: "au_fcc", displayName: "Gold (FCC)", structure: CrystalPresets.gold, source: .builtIn),
        PhaseModel(id: "al_fcc", displayName: "Aluminium (FCC)", structure: CrystalPresets.aluminum, source: .builtIn),
        PhaseModel(id: "ni_fcc", displayName: "Nickel (FCC)", structure: CrystalPresets.nickel, source: .builtIn),
        PhaseModel(id: "cu_fcc", displayName: "Copper (FCC)", structure: CrystalPresets.copper, source: .builtIn),
        PhaseModel(id: "fe_bcc", displayName: "Iron (BCC)", structure: CrystalPresets.iron, source: .builtIn),
        PhaseModel(id: "si_diamond", displayName: "Silicon (diamond)", structure: CrystalPresets.silicon, source: .builtIn),
        PhaseModel(id: "mg_hcp", displayName: "Magnesium (HCP)", structure: CrystalPresets.magnesium, source: .builtIn),
        PhaseModel(id: "ws2_2h", displayName: "Tungsten disulfide (2H-WS₂)", structure: CrystalPresets.tungstenDisulfide, source: .builtIn),
    ]

    public static func model(id: String) -> PhaseModel? { models.first { $0.id == id } }

    /// A single-element cubic model for a user-defined crystal.
    public static func customCubic(family: CrystalPresets.CubicFamily, latticeA: Double, atomicNumber: Int) -> PhaseModel {
        let symbol = Element.symbol(forAtomicNumber: atomicNumber) ?? "Z\(atomicNumber)"
        return PhaseModel(
            id: "custom_cubic_\(family.provenanceID)_z\(atomicNumber)",
            displayName: "\(symbol) \(family.rawValue), a = \(String(format: "%.4g", latticeA)) Å",
            structure: CrystalPresets.cubic(family, a: latticeA, z: atomicNumber), source: .custom)
    }
}

// MARK: - Materials Project

extension PhaseModel {
    /// A validated model from one fetched Materials Project document (`MaterialsProjectImport.structure(from:)`): the
    /// orientation class is detected from the atoms, never taken from MP's stated symmetry, which is recorded in
    /// `extraProvenance` only. A deprecated entry still imports (`materials_project_deprecated` rides in the provenance).
    /// Any validation issue becomes `MaterialsProjectImport.Failure.malformed`.
    public init(materialsProject doc: MaterialsProjectDocument, fetchedAt: Date) throws {
        let imported = try MaterialsProjectImport.structure(from: doc)
        let label = doc.formulaPretty?.trimmingCharacters(in: .whitespaces)
        let displayName = (label?.isEmpty == false) ? "\(label!) (\(doc.materialID))" : doc.materialID
        let cellProvenance: String
        switch imported.cell {
        case .asServed: cellProvenance = "as served"
        case .conventional(let centring): cellProvenance = "conventional from \(centring)-centred primitive"
        }
        let formatter = ISO8601DateFormatter()
        formatter.timeZone = TimeZone(secondsFromGMT: 0)
        formatter.formatOptions = [.withInternetDateTime]
        var extra: [String: String] = [
            "materials_project_id": doc.materialID,
            "materials_project_fetched": formatter.string(from: fetchedAt),
            "materials_project_cell": cellProvenance,
        ]
        if doc.deprecated == true { extra["materials_project_deprecated"] = "true" }
        if let symbol = doc.symmetry?.symbol { extra["materials_project_symbol"] = symbol }
        if let system = doc.symmetry?.crystalSystem { extra["materials_project_crystal_system"] = system }
        self.init(id: doc.materialID, displayName: displayName, structure: imported.structure, source: .materialsProject,
                  spaceGroupNumber: doc.symmetry?.number, extraProvenance: extra)
        let issues = validationIssues
        guard issues.isEmpty else {
            throw MaterialsProjectImport.Failure.malformed(
                "\(doc.materialID) failed validation: \(issues.map(\.message).joined(separator: "; "))")
        }
    }
}
