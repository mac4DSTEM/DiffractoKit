//
//  ConventionalSetting.swift
//  Role: the symmetry of a structure in ANY cell setting. The
//        sites are carried into the conventional cell of the lattice (x_conventional = P⁻¹ x_file, wrapped, the lattice
//        centring added as images and merged within the position tolerance), the operations are detected there by
//        the conventional-cell rule (a conventional cell has its holohedry in the standard position), and the result is reported
//        in BOTH settings: `CrystalStructure` carries the file setting's operations (conjugated back, R' = P R P⁻¹,
//        t' = P t) and `ConventionalSetting` the conventional cell's own structure, operations, Laue class (so its
//        absences and ring names) and the transformation between the two.
//

import Foundation
import simd

/// A structure in the conventional cell of its lattice, and how that cell relates to the one the structure came in.
public struct ConventionalSetting: Sendable, Equatable {
    /// The conventional cell, its Bravais lattice and P (x_file = P x_conventional; reflections h_conv = h_file P).
    public let conventionalCell: ConventionalCell
    /// The sites in the conventional cell (the centring images included).
    public let sites: [AtomSite]
    /// The operations detected in the conventional cell, centring translations included.
    public let operations: [SymmetryOperation]
    public let laueGroup: PointSymmetry
    public let laueClass: LaueClass?
    public let translationUncertainty: Double

    public var cell: UnitCell { conventionalCell.cell }
    /// False when the Niggli reduction did not converge (`ConventionalCell.reductionConverged`): the operations here are
    /// then the detection in the file's own cell, and `bravais` (`aP`) says nothing about the lattice.
    public var reductionConverged: Bool { conventionalCell.reductionConverged }
    public var bravais: BravaisLattice { conventionalCell.bravais }
    public var transformation: IntegerMatrix { conventionalCell.transformation }

    /// The structure in the conventional cell, `setting == .conventional`: its reflections and rings are the
    /// conventional setting's.
    public var structure: CrystalStructure {
        CrystalStructure(conventionalCell: cell, sites: sites, operations: operations, laueGroup: laueGroup,
                         laueClass: laueClass, translationUncertainty: translationUncertainty)
    }

    /// h' = h P: the file's (h, k, l) in the conventional cell's indices.
    public func conventionalIndices(of fileHKL: MillerIndex) -> MillerIndex {
        fileHKL.applying(transformation)
    }

    /// The file's indices of a conventional (h, k, l): h = h' P⁻¹; nil when they are not integers (a conventional
    /// reflection that is not a reflection of the file's lattice: a centring absence).
    public func fileIndices(ofConventional hkl: MillerIndex) -> MillerIndex? {
        let det = transformation.determinant
        let m = hkl.applying(transformation.adjugate)
        guard m.h % det == 0, m.k % det == 0, m.l % det == 0 else { return nil }
        return MillerIndex(m.h / det, m.k / det, m.l / det)
    }

    /// An operation of the conventional cell in the file's cell: R' = P R P⁻¹, t' = P t (translation mod 1);
    /// nil when R' is not an integer matrix.
    public func fileOperation(of op: SymmetryOperation) -> SymmetryOperation? {
        guard let r = transformation.conjugating(op.rotation) else { return nil }
        return SymmetryOperation(rotation: r, translation: wrapFractional(transformation.apply(op.translation)))
    }

    /// An operation of the file's cell in the conventional cell: R = P⁻¹ R' P, t = P⁻¹ t'.
    public func conventionalOperation(of op: SymmetryOperation) -> SymmetryOperation? {
        guard let r = transformation.conjugatedBy(op.rotation) else { return nil }
        return SymmetryOperation(rotation: r, translation: transformation.inverseApplied(to: op.translation))
    }
}

extension CrystalStructure {
    init(conventionalCell cell: UnitCell, sites: [AtomSite], operations: [SymmetryOperation], laueGroup: PointSymmetry,
         laueClass: LaueClass?, translationUncertainty: Double) {
        self.cell = cell
        self.sites = sites
        self.setting = .conventional
        self.conventionalSetting = nil
        self.operations = operations
        self.symmetrySource = .detected
        self.fileOperations = nil
        self.fileOperationsAreSymmetries = nil
        self.laueGroup = laueGroup
        self.laueClass = laueClass
        self.spaceGroupNumber = nil
        self.name = nil
        self.unresolvedRotations = []
        self.translationUncertainty = translationUncertainty
        self.symmetryRefusal = nil
    }

    /// `detectingSymmetry(conventional: true)`.
    static func detectingInConventionalCell(
        cell: UnitCell, sites: [AtomSite], fileOperations: [SymmetryOperation]?, spaceGroupNumber: Int?, name: String?,
        coarsestHalfStep: Double, metricLengthTolerance: Double, metricAngleToleranceDeg: Double,
        tolerance: SymmetryDetection.Tolerance, search: SymmetryDetection.Search,
        policy: SymmetryDetection.Policy = .default
    ) -> CrystalStructure {
        // An unusable position tolerance is refused before anything is built from it (the file's cell is checked here, the
        // conventional cell's by `detect`): the identity alone, with the reason.
        if let reason = tolerance.refusal(for: cell) {
            return refusedStructure(cell: cell, sites: sites, fileOperations: fileOperations, spaceGroupNumber: spaceGroupNumber,
                                    name: name, reason: reason)
        }
        let conventional = cell.conventionalCell(tolerance: metricLengthTolerance)
        let p = conventional.transformation
        if p == .identity || sites.isEmpty {
            // The file's cell is the conventional cell: nothing moves, the plain detection is the detection.
            let detection = SymmetryDetection.detect(
                sites: sites, cell: cell, coarsestHalfStep: coarsestHalfStep,
                extraCandidates: fileOperations?.map(\.rotation) ?? [], preferredOperations: fileOperations ?? [],
                metricLengthTolerance: metricLengthTolerance, metricAngleToleranceDeg: metricAngleToleranceDeg,
                tolerance: tolerance, search: search, policy: policy)
            let setting = ConventionalSetting(
                conventionalCell: conventional, sites: sites, operations: detection.operations,
                laueGroup: PointSymmetry.laueGroup(of: detection.operations),
                laueClass: PointSymmetry.laueGroup(of: detection.operations).laueClass,
                translationUncertainty: detection.translationUncertainty)
            return CrystalStructure(cell: cell, sites: sites, detection: detection, fileOperations: fileOperations,
                                    spaceGroupNumber: spaceGroupNumber, name: name, conventionalSetting: setting)
        }

        // The sites in the conventional cell: x' = P⁻¹ x, each lattice centring added as an image, duplicates merged.
        let matcher = PositionMatcher(cell: conventional.cell, tolerance: tolerance)
        var convSites: [AtomSite] = []
        for site in sites {
            let x = wrapFractional(p.inverseApplied(to: site.fractional))
            convSites += mergingCoincident(images: conventional.centringTranslations.map { t in
                var image = site
                image.fractional = wrapFractional(x + t)
                return image
            }, matcher: matcher)
        }

        // The file's own operations, in the conventional cell (rotations that are not lattice symmetries are dropped).
        let probe = ConventionalSetting(conventionalCell: conventional, sites: [], operations: [],
                                        laueGroup: PointSymmetry(rotations: []), laueClass: nil, translationUncertainty: 0)
        let convFileOperations = (fileOperations ?? []).compactMap { probe.conventionalOperation(of: $0) }
        // The written half-step is a fractional error of the file's coordinates; in the conventional cell it grows by
        // the largest row sum of |P⁻¹|.
        // Typed locals: this call took 2.9 s to type-check in one expression (a CI runner times out at that).
        let det: Double = Double(p.determinant)
        var growth: Double = 1
        for row in p.adjugate.rows {
            let sum: Double = abs(Double(row.x)) + abs(Double(row.y)) + abs(Double(row.z))
            growth = max(growth, sum / det)
        }
        let halfStep: Double = coarsestHalfStep * growth
        let extraRotations: [IntegerMatrix] = convFileOperations.map(\.rotation)
        let detection = SymmetryDetection.detect(
            sites: convSites, cell: conventional.cell, coarsestHalfStep: halfStep,
            extraCandidates: extraRotations, preferredOperations: convFileOperations,
            metricLengthTolerance: metricLengthTolerance, metricAngleToleranceDeg: metricAngleToleranceDeg,
            tolerance: tolerance, search: search, policy: policy)
        if let reason = detection.refusal {
            return refusedStructure(cell: cell, sites: sites, fileOperations: fileOperations, spaceGroupNumber: spaceGroupNumber,
                                    name: name, reason: reason)
        }

        // Back to the file's cell: the centring translations become lattice translations and fall out as duplicates.
        var fileSetting: [SymmetryOperation] = []
        for op in detection.operations {
            guard let back = probe.fileOperation(of: op) else { continue }
            if !fileSetting.contains(where: { $0.isEquivalent(to: back, translationTolerance: 1e-6) }) {
                fileSetting.append(back)
            }
        }
        fileSetting.sort(by: operationOrder)
        let unresolved = detection.unresolvedRotations.compactMap { p.conjugating($0) }
        let fileDetection = SymmetryDetection.Result(
            candidates: PointSymmetry(rotations: detection.candidates.rotations.compactMap { p.conjugating($0) }),
            operations: fileSetting, unresolvedRotations: unresolved,
            translationUncertainty: detection.translationUncertainty)
        let convGroup = PointSymmetry.laueGroup(of: detection.operations)
        let setting = ConventionalSetting(
            conventionalCell: conventional, sites: convSites, operations: detection.operations, laueGroup: convGroup,
            laueClass: convGroup.laueClass, translationUncertainty: detection.translationUncertainty)
        return CrystalStructure(cell: cell, sites: sites, detection: fileDetection, fileOperations: fileOperations,
                                spaceGroupNumber: spaceGroupNumber, name: name, conventionalSetting: setting)
    }

    /// The structure of a detection that was refused (an unusable position tolerance): the identity alone, in the file's
    /// cell, with the reason (`symmetryRefusal`).
    static func refusedStructure(cell: UnitCell, sites: [AtomSite], fileOperations: [SymmetryOperation]?,
                                 spaceGroupNumber: Int?, name: String?, reason: String) -> CrystalStructure {
        let detection = SymmetryDetection.Result(
            candidates: PointSymmetry(rotations: [.identity]), operations: [.identity], unresolvedRotations: [],
            translationUncertainty: 1e-9, refusal: reason)
        return CrystalStructure(cell: cell, sites: sites, detection: detection, fileOperations: fileOperations,
                                spaceGroupNumber: spaceGroupNumber, name: name)
    }

    /// The identity first, then by rotation and translation: the order `SymmetryDetection` returns.
    static func operationOrder(_ lhs: SymmetryOperation, _ rhs: SymmetryOperation) -> Bool {
        if lhs.rotation.isIdentity != rhs.rotation.isIdentity { return lhs.rotation.isIdentity }
        let l = lhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }, r = rhs.rotation.rows.flatMap { [$0.x, $0.y, $0.z] }
        if l != r { return l.lexicographicallyPrecedes(r) }
        return (lhs.translation.x, lhs.translation.y, lhs.translation.z)
            < (rhs.translation.x, rhs.translation.y, rhs.translation.z)
    }

    /// The centring images of ONE site, without any that coincides with an earlier image of the same site (a centring
    /// translation listed twice). Rows are never merged with each other: two rows on one position (a mixed-occupancy
    /// site written as two rows of one element) are two sites in the conventional cell too.
    static func mergingCoincident(images: [AtomSite], matcher: PositionMatcher) -> [AtomSite] {
        var kept: [AtomSite] = []
        for image in images where !kept.contains(where: { matcher.matches($0.fractional, image.fractional) }) {
            kept.append(image)
        }
        return kept
    }
}
