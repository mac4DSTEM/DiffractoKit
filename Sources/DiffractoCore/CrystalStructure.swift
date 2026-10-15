//
//  CrystalStructure.swift
//  Role: a cell with its atoms, its symmetry operations (from a file, or detected from the atoms), its Laue group
//        and class, and the validation the reference implementations' `CrystalModel` applies (reference commit 7448401b `CrystalModel.swift`:
//        cell, basis, occupancies, finite coordinates, the 0.5 Å close contact). No app type: nothing here is a
//        "model", an id, a source or a provenance dictionary; those stay in the reference implementations.
//

import Foundation
import simd

/// One atom in the unit cell.
public struct AtomSite: Sendable, Equatable {
    public var z: Int
    public var fractional: SIMD3<Double>
    public var occupancy: Double
    /// The `_atom_site_label` when the site came from a CIF.
    public var label: String?

    public init(z: Int, fractional: SIMD3<Double>, occupancy: Double = 1, label: String? = nil) {
        self.z = z
        self.fractional = fractional
        self.occupancy = occupancy
        self.label = label
    }
}

public struct CrystalStructure: Sendable, Equatable {
    /// Where the operations came from.
    public enum SymmetrySource: Sendable, Equatable {
        /// Given to the initialiser as they are, not checked against the sites (a structure built in code).
        case operationsGiven
        /// Found from the sites by the holohedry test (registration 3.4b).
        case detected
        /// A file's loop expanded the asymmetric unit, then the operations were detected from the expanded cell
        /// (3.4 as amended 2026-10-14); the loop is `fileOperations`.
        case detectedAfterFileExpansion
    }

    /// Which cell the numbers of this value (`cell`, `sites`, `operations`, `laueGroup`, and the reflections and rings
    /// built from them) are in. A file's own setting is never relabelled silently (slice 1b).
    public enum Setting: Sendable, Equatable {
        /// The cell the structure was given in.
        case file
        /// The conventional (centred) cell of the lattice: `ConventionalSetting.structure`.
        case conventional
    }

    public let cell: UnitCell
    public let sites: [AtomSite]
    /// The cell `cell`, `sites` and `operations` are in: `.file` for everything a caller builds or imports.
    public let setting: Setting
    /// The same structure in the conventional cell of its lattice, with the transformation between the two cells
    /// (`detectingSymmetry(conventional: true)`, the default); nil when the operations were given, or when
    /// `conventional: false`. `operations`, `laueGroup` and `laueClass` above are the file setting's, conjugated back
    /// from the operations detected in the conventional cell.
    public let conventionalSetting: ConventionalSetting?
    /// The space-group operations modulo lattice translations (the identity included; for a centred cell, the
    /// centring translations as operations with the identity rotation).
    public let operations: [SymmetryOperation]
    public let symmetrySource: SymmetrySource
    /// The operation loop of the file the structure came from (deduplicated, translations mod 1); nil when there was
    /// none. It expanded the sites; it may declare less symmetry than the atoms have.
    public let fileOperations: [SymmetryOperation]?
    /// Whether every file operation is among the detected ones (nil without a file loop). False means the file's loop
    /// is not a symmetry of the cell it describes: the file is inconsistent, or the metric's holohedry did not hold
    /// one of its rotations.
    public let fileOperationsAreSymmetries: Bool?
    /// {±R} over the operations' rotation parts.
    public let laueGroup: PointSymmetry
    /// nil only when the operations do not form a Laue group (a file whose loop is not a group).
    public let laueClass: LaueClass?
    /// The IT number the file states (`_symmetry_int_tables_number` / `_space_group_it_number`); nil means "not
    /// recorded", not P1.
    public let spaceGroupNumber: Int?
    /// A name for the structure (a CIF's data block, a formula); nil when none was given.
    public let name: String?
    /// Candidate rotations of the holohedry that failed the fixed position tolerance but would pass at the file's own
    /// written precision (3.4b): reported, not taken. Empty when the operations were given.
    public let unresolvedRotations: [IntegerMatrix]
    /// The largest uncertainty of a detected translation that is neither exact (from the file, or on the k/24 grid):
    /// its least-squares residual plus twice the written half-step, at least 1e-9. Reported; the absence test does not
    /// need it (`ReflectionList.isAllowed`).
    public let translationUncertainty: Double
    /// Why the symmetry detection did not run: the position tolerance it was given is outside its domain
    /// (`SymmetryDetection.Tolerance`). The operations are then the identity alone, `validationIssues` reports it
    /// (`invalid_symmetry_tolerance`) and `isUsable` is false. nil for every structure that was detected or given.
    public let symmetryRefusal: String?

    /// A structure whose operations are given as they are (a cell built in code with a known group); nothing is
    /// checked against the sites. A CIF goes through `CIFImport`, which detects.
    public init(cell: UnitCell, sites: [AtomSite], operations: [SymmetryOperation],
                spaceGroupNumber: Int? = nil, name: String? = nil) {
        self.cell = cell
        self.sites = sites
        self.setting = .file
        self.conventionalSetting = nil
        self.operations = operations.isEmpty ? [.identity] : operations
        self.symmetrySource = .operationsGiven
        self.fileOperations = nil
        self.fileOperationsAreSymmetries = nil
        self.laueGroup = PointSymmetry.laueGroup(of: self.operations)
        self.laueClass = laueGroup.laueClass
        self.spaceGroupNumber = spaceGroupNumber
        self.name = name
        self.unresolvedRotations = []
        self.translationUncertainty = 1e-9
        self.symmetryRefusal = nil
    }

    /// A structure whose operations are detected from its sites (registration 3.4b): a P1 file, a cell built in
    /// code, or a file's expanded cell with its loop as `fileOperations` (whose rotations join the candidates).
    /// `coarsestHalfStep` is the written precision of the coordinates when known, used only to report an unresolved
    /// candidate.
    ///
    /// `metricLengthTolerance` (relative, default 1e-9) is the tolerance of the lattice-symmetry test that finds the
    /// conventional cell and its Bravais lattice; `metricAngleToleranceDeg` has no effect on that route (it applies
    /// to the metric family in the file's cell when `conventional` is false). A primitive bcc cell written with the angle
    /// 109.4712° is a rhombohedral lattice at 1e-9 and reads Im-3m only from `metricLengthTolerance` ≈ 1e-5 (spglib's
    /// default symprec reads it so). If the lattice cannot be reduced, `conventionalSetting.reductionConverged` is false
    /// and the operations are those of the file's cell (see `ConventionalCell.reductionConverged`). `tolerance`:
    /// `.cartesian(t)` finds operations only to about half of `t` (see `SymmetryDetection.Tolerance`). A `tolerance`
    /// outside its domain (NaN, zero, negative, infinite, below 1e-9 or as wide as half a cell) is refused: the structure
    /// has the identity alone, `symmetryRefusal` names the reason and `validationIssues` reports it (`isUsable` is false).
    /// The two metric tolerances read NaN, negative and infinite as 0.
    public static func detectingSymmetry(
        cell: UnitCell, sites: [AtomSite], fileOperations: [SymmetryOperation]? = nil,
        spaceGroupNumber: Int? = nil, name: String? = nil, coarsestHalfStep: Double = 0,
        metricLengthTolerance: Double = SymmetryDetection.defaultMetricLengthTolerance,
        metricAngleToleranceDeg: Double = SymmetryDetection.defaultMetricAngleToleranceDeg,
        tolerance: SymmetryDetection.Tolerance = .default,
        conventional: Bool = true,
        search: SymmetryDetection.Search = .cellList,
        policy: SymmetryDetection.Policy = .default
    ) -> CrystalStructure {
        if conventional {
            return detectingInConventionalCell(
                cell: cell, sites: sites, fileOperations: fileOperations, spaceGroupNumber: spaceGroupNumber, name: name,
                coarsestHalfStep: coarsestHalfStep, metricLengthTolerance: metricLengthTolerance,
                metricAngleToleranceDeg: metricAngleToleranceDeg, tolerance: tolerance, search: search,
                policy: policy)
        }
        let detection = SymmetryDetection.detect(
            sites: sites, cell: cell, coarsestHalfStep: coarsestHalfStep,
            extraCandidates: fileOperations?.map(\.rotation) ?? [], preferredOperations: fileOperations ?? [],
            metricLengthTolerance: metricLengthTolerance, metricAngleToleranceDeg: metricAngleToleranceDeg,
            tolerance: tolerance, search: search, policy: policy)
        return CrystalStructure(cell: cell, sites: sites, detection: detection, fileOperations: fileOperations,
                                spaceGroupNumber: spaceGroupNumber, name: name)
    }

    init(cell: UnitCell, sites: [AtomSite], detection: SymmetryDetection.Result,
         fileOperations: [SymmetryOperation]?, spaceGroupNumber: Int?, name: String?,
         setting: Setting = .file, conventionalSetting: ConventionalSetting? = nil) {
        self.cell = cell
        self.sites = sites
        self.setting = setting
        self.conventionalSetting = conventionalSetting
        self.operations = detection.operations
        self.symmetrySource = fileOperations == nil ? .detected : .detectedAfterFileExpansion
        self.fileOperations = fileOperations
        self.fileOperationsAreSymmetries = fileOperations.map { ops in
            ops.allSatisfy { op in detection.operations.contains { $0.isEquivalent(to: op, translationTolerance: 1e-6) } }
        }
        self.laueGroup = PointSymmetry.laueGroup(of: detection.operations)
        self.laueClass = laueGroup.laueClass
        self.spaceGroupNumber = spaceGroupNumber
        self.name = name
        self.unresolvedRotations = detection.unresolvedRotations
        self.translationUncertainty = detection.translationUncertainty
        self.symmetryRefusal = detection.refusal
    }

    public var crystalSystem: CrystalSystem? { laueClass?.crystalSystem }

    /// The cell content by element, "Mg10 Si12" (symbols in increasing Z, occupancies summed and printed as integers
    /// when they are).
    public var formula: String {
        var counts: [Int: Double] = [:]
        for site in sites { counts[site.z, default: 0] += site.occupancy }
        return counts.keys.sorted().map { z in
            let n = counts[z]!
            let count = abs(n - n.rounded()) < 1e-9 ? String(Int(n.rounded())) : String(format: "%.3g", n)
            return "\(Element.symbol(forAtomicNumber: z) ?? "Z\(z)")\(count)"
        }.joined(separator: " ")
    }

    // MARK: Validation (the reference implementations' CrystalModel.validationIssues, without the app-only checks)

    public struct ValidationIssue: Sendable, Equatable, CustomStringConvertible {
        public let code: String
        public let message: String
        public var description: String { message }
    }

    /// Below the shortest interatomic contact known in any crystal (H–H in solid hydrogen, 0.74 Å).
    public static let closeContactLimit = 0.5

    /// The thresholds of the structure checks (DiffractoKit 0.2.0): the suite's choices, with `ValidationPolicy.default`
    /// the set every release before used.
    public struct ValidationPolicy: Sendable, Equatable {
        /// Two sites closer than this (Å) are a close contact. Default 0.5. A limit that is not a positive number
        /// (NaN, zero, negative) finds no contact.
        public var closeContactLimit: Double
        /// Two sites whose occupancies sum to no more than 1 + this are a split site, never a contact. Default 0.02
        /// (sum 1.02). A tolerance that is NaN or negative is taken as 0.
        public var splitSiteOccupancyTolerance: Double

        public init(closeContactLimit: Double = CrystalStructure.closeContactLimit,
                    splitSiteOccupancyTolerance: Double = 0.02) {
            self.closeContactLimit = closeContactLimit
            self.splitSiteOccupancyTolerance = splitSiteOccupancyTolerance
        }

        /// The suite's values.
        public static let `default` = ValidationPolicy()
    }

    public var validationIssues: [ValidationIssue] { validationIssues(using: .default) }

    /// `validationIssues` under the thresholds of `policy`.
    public func validationIssues(using policy: ValidationPolicy) -> [ValidationIssue] {
        var issues: [ValidationIssue] = []
        func add(_ code: String, _ message: String) { issues.append(.init(code: code, message: message)) }
        for issue in cell.issues {
            switch issue {
            case .invalidLengths: add("invalid_cell_lengths", issue.description)
            case .invalidAngles: add("invalid_cell_angles", issue.description)
            case .degenerateVolume: add("invalid_cell_volume", issue.description)
            }
        }
        if let refusal = symmetryRefusal { add("invalid_symmetry_tolerance", refusal) }
        if sites.isEmpty {
            add("missing_atomic_basis", "At least one occupied atomic site is required.")
        }
        for site in sites {
            if site.z < 1 || site.z > Element.maximumAtomicNumber {
                add("invalid_atomic_number", "Atomic numbers must be 1 to \(Element.maximumAtomicNumber); Z = \(site.z) is not.")
            }
            if !site.occupancy.isFinite || site.occupancy <= 0 || site.occupancy > 1 {
                add("invalid_occupancy", "Atomic occupancies must be finite and in (0, 1].")
            }
            if !site.fractional.x.isFinite || !site.fractional.y.isFinite || !site.fractional.z.isFinite {
                add("invalid_fractional_coordinate", "Atomic fractional coordinates must be finite.")
            }
        }
        // An over-populated cell is the one modelling error that produces no visible symptom: site multiplicity
        // multiplies straight into the structure factors. The message states what was measured and stops there.
        if let closest = shortestCloseContact(using: policy) {
            add("atomic_sites_too_close",
                String(format: "Two fully-occupied atomic sites are %.4f Å apart — closer than any real interatomic "
                       + "distance, so the cell holds more atoms than the structure has.", closest))
        }
        return issues
    }

    public var isUsable: Bool { validationIssues.isEmpty }

    public func isUsable(under policy: ValidationPolicy) -> Bool { validationIssues(using: policy).isEmpty }

    /// Shortest centre-to-centre distance (Å) between two sites that cannot both be occupied by a whole atom, when
    /// under `closeContactLimit`; nil when there is no such pair. Pairs whose occupancies sum to ≤ 1.02 are a split
    /// site and are skipped (the reference implementations' rule). Minimum image over the 27 neighbouring cells behind a per-axis prefilter.
    public var shortestCloseContact: Double? { shortestCloseContact(using: .default) }

    /// `shortestCloseContact` under the limit and the split-site rule of `policy`.
    public func shortestCloseContact(using policy: ValidationPolicy) -> Double? {
        guard sites.count > 1 else { return nil }
        let lattice = cell.realLattice
        guard lattice.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }),
              cell.volume.isFinite, cell.volume > 0
        else { return nil }
        let faceAreas = [
            cross(lattice[1], lattice[2]), cross(lattice[2], lattice[0]), cross(lattice[0], lattice[1]),
        ].map { length($0) }
        let limit = policy.closeContactLimit.isFinite && policy.closeContactLimit > 0 ? policy.closeContactLimit : 0
        guard limit > 0 else { return nil }
        let splitBound = 1 + (policy.splitSiteOccupancyTolerance.isFinite && policy.splitSiteOccupancyTolerance >= 0
                              ? policy.splitSiteOccupancyTolerance : 0)
        let maxAxisDelta = faceAreas.map { area -> Double in
            let height = area > 0 ? cell.volume / area : 0
            return height > 0 ? min(0.5, limit / height) : 0.5
        }
        var best = Double.infinity
        for i in 0..<sites.count {
            for j in (i + 1)..<sites.count {
                guard sites[i].occupancy + sites[j].occupancy > splitBound else { continue }
                let p = sites[i].fractional, q = sites[j].fractional
                guard fractionalAxisDelta(p.x, q.x) <= maxAxisDelta[0],
                      fractionalAxisDelta(p.y, q.y) <= maxAxisDelta[1],
                      fractionalAxisDelta(p.z, q.z) <= maxAxisDelta[2]
                else { continue }
                best = min(best, CrystalStructure.minimumImageDistanceSquared(p, q, lattice))
            }
        }
        let closest = best.squareRoot()
        return closest < limit ? closest : nil
    }

    /// Minimum-image squared distance (Å²) of two fractional positions over the 27 neighbouring cells.
    static func minimumImageDistanceSquared(_ p: SIMD3<Double>, _ q: SIMD3<Double>,
                                            _ lattice: [SIMD3<Double>]) -> Double {
        let base = p - q
        var best = Double.infinity
        for dx in -1...1 {
            for dy in -1...1 {
                for dz in -1...1 {
                    let f = base + SIMD3<Double>(Double(dx), Double(dy), Double(dz))
                    let cartesian = f.x * lattice[0] + f.y * lattice[1] + f.z * lattice[2]
                    best = min(best, (cartesian * cartesian).sum())
                }
            }
        }
        return best
    }

    // MARK: Reflections

    /// Every reflection with 0 < |g| ≤ `kMax` (Å⁻¹), with its systematic-absence flag from this structure's
    /// operations, and its families by this structure's Laue group (registration 3.8, 3.9). A reflection with |g|
    /// exactly `kMax` is in or out by the last bit of |g|, which can differ by one ulp between two settings of the same
    /// lattice, so the lists of the file's cell and of its conventional cell may differ by such a reflection (Gate B on
    /// 4876b18).
    public func reflections(kMax: Double) throws(ReflectionListError) -> ReflectionList {
        try reflections(kMax: kMax, policy: .default)
    }

    /// The same under the refusal bounds and tolerances of `policy` (DiffractoKit 0.2.0).
    public func reflections(kMax: Double, policy: ReflectionList.Policy) throws(ReflectionListError) -> ReflectionList {
        try ReflectionList(cell: cell, operations: operations, laueGroup: laueGroup, laueClass: laueClass, kMax: kMax,
                           policy: policy)
    }
}
