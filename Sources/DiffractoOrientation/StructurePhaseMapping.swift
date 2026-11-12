//
//  StructurePhaseMapping.swift
//  Role: a scan's phases mapped from the structures alone, as one call (DiffractoKit#99). Nothing here is a new method: it is
//        the order of the steps of DiffractoKit#81 and the values they share, in one place.
//
//    1. the matrix orientation of every position, by template matching (`matrixOrientations`);
//    2. the grain rule, with the spot test (`regularizingMatrixOrientations(…explainedBy:)`);
//    3. the matrix's spots removed, a neighbouring grain's too (`removingMatrixSpots(…neighbouringGrains:)`);
//    4. the two passes on the whole pattern (`classifyInTwoPasses`): the guarded free search, the relationships learned
//       from its strong positions, every position classified against the variants.
//
//  The reach is one number for all of it: the templates and reflections of the matrix and of every candidate are built to
//  it, and it is the radius the chance guard counts stray vectors in.
//
//  Checked against planted truth (`StructurePhaseMappingTests`); the measurements on real data are in docs/parity.md under
//  `classifyInTwoPasses`, made with these steps in this order.
//

import Foundation
import simd
import DiffractoCore

/// The settings of ``PhaseVectorMatcher/mapFromStructures(scan:scanWidth:scanHeight:matrix:candidates:settings:originX:originY:invAngstromPerPixel:maximumWorkers:isCancelled:progress:)``.
/// The tolerances mean what they mean in `PhaseVectorSettings` and `PhaseOrientationSearchSettings`, with the same defaults.
public struct StructurePhaseMappingSettings: Sendable, Equatable {
    /// How far from the origin the detector records, Å⁻¹: the templates and reflections are built to it.
    public var reachInvAngstrom: Double
    /// Vectors longer than this are not used, Å⁻¹. 0 uses every vector.
    public var maximumVectorInvAngstrom: Double
    /// A vector within this distance of a predicted matrix spot is the matrix's, Å⁻¹ (the removal and the grain rule's
    /// spot test).
    public var matrixToleranceInvAngstrom: Double
    /// A vector within this distance of a predicted spot of a phase is matched, Å⁻¹.
    public var pairRadiusInvAngstrom: Double
    /// A winner whose score is above this is not indexed, Å⁻¹.
    public var residualCutoffInvAngstrom: Double
    /// Vectors inside it are the direct beam, Å⁻¹.
    public var directBeamRadiusInvAngstrom: Double
    /// The grain rule's and the removal's neighbourhood, in scan positions to each side.
    public var neighbourhoodRadius: Int
    /// The neighbour rule of the two passes (off by default, as there).
    public var usesNeighbourSupport: Bool
    /// The fewest strong positions a relationship is learned from.
    public var minimumPositions: Int

    public init(reachInvAngstrom: Double, maximumVectorInvAngstrom: Double = 0, matrixToleranceInvAngstrom: Double = 0.02,
                pairRadiusInvAngstrom: Double = 0.02, residualCutoffInvAngstrom: Double = 0.07,
                directBeamRadiusInvAngstrom: Double = 0.15, neighbourhoodRadius: Int = 2, usesNeighbourSupport: Bool = false,
                minimumPositions: Int = 20) {
        self.reachInvAngstrom = reachInvAngstrom
        self.maximumVectorInvAngstrom = maximumVectorInvAngstrom
        self.matrixToleranceInvAngstrom = matrixToleranceInvAngstrom
        self.pairRadiusInvAngstrom = pairRadiusInvAngstrom
        self.residualCutoffInvAngstrom = residualCutoffInvAngstrom
        self.directBeamRadiusInvAngstrom = directBeamRadiusInvAngstrom
        self.neighbourhoodRadius = neighbourhoodRadius
        self.usesNeighbourSupport = usesNeighbourSupport
        self.minimumPositions = minimumPositions
    }
}

/// A phase given by its structure, for ``PhaseVectorMatcher/mapFromStructures(scan:scanWidth:scanHeight:matrix:candidates:settings:originX:originY:invAngstromPerPixel:maximumWorkers:isCancelled:progress:)``.
public struct StructurePhase: Sendable {
    public let name: String
    public let crystal: CrystalStructure
    /// Half the thickness of the slab of reciprocal space an orientation brings into the pattern, Å⁻¹ (a candidate's;
    /// the matrix uses the removal's own slab).
    public var excitationSlabInvAngstrom: Double

    public init(name: String, crystal: CrystalStructure, excitationSlabInvAngstrom: Double = 0.05) {
        self.name = name
        self.crystal = crystal
        self.excitationSlabInvAngstrom = excitationSlabInvAngstrom
    }
}

/// What ``PhaseVectorMatcher/mapFromStructures(scan:scanWidth:scanHeight:matrix:candidates:settings:originX:originY:invAngstromPerPixel:maximumWorkers:isCancelled:progress:)`` returns.
public struct StructurePhaseMapping: Sendable {
    /// The verdict of every position. Phase 0 is the matrix, phase k + 1 the k-th candidate. No verdict comes from a
    /// reference-library entry, so every `entryIndex` and the map's `matrixEntryIndex` are −1; `score`, `matchedCount`
    /// and `survivingCount` are the final verdict's, `removedCount` the removal's.
    public let map: PhaseMap
    /// The matrix orientation field after the grain rule.
    public let matrixField: MatrixOrientationField
    /// How many positions the grain rule changed.
    public let regularizedPositionCount: Int
    /// The matrix orientation of every position with the removal's in-plane correction applied: the one the
    /// relationships and the variants are built on.
    public let matrixOrientations: [simd_double3x3?]
    /// The scan after the matrix removal.
    public let residual: MatrixResidualScan
    /// The two passes: final and free-search verdicts (their `phaseIndex` counts the candidates from 0), the learned
    /// relationships, the angles to them.
    public let classification: TwoPassClassification
    /// The symmetries the relationships were clustered by.
    public let matrixSymmetry: CrystalOrientationSymmetry
    public let candidateSymmetries: [CrystalOrientationSymmetry]
}

extension PhaseVectorMatcher {

    /// A scan's phases from the structures alone: the matrix orientation of every position, the grain rule, the matrix
    /// removal and the two passes (see the file header).
    ///
    /// - Parameters:
    ///   - scan: every position's peaks, row-major, in detector pixels.
    ///   - matrix: the matrix phase.
    ///   - candidates: the phases to look for, each structure once.
    ///   - progress: 0 to 1; the matrix step reports up to 0.6, the later steps at their ends.
    /// - Returns: nil when the scan's size does not fit, there is no candidate, or the run was cancelled.
    /// - Throws: ``ACOMError`` when a structure is refused or its Laue class cannot be determined.
    public static func mapFromStructures(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                         matrix: StructurePhase, candidates: [StructurePhase],
                                         settings: StructurePhaseMappingSettings,
                                         originX: Float, originY: Float, invAngstromPerPixel: Double,
                                         maximumWorkers: Int? = nil,
                                         isCancelled: (@Sendable () -> Bool)? = nil,
                                         progress: (@Sendable (Double) -> Void)? = nil) throws(ACOMError) -> StructurePhaseMapping? {
        let (positions, overflow) = scanWidth.multipliedReportingOverflow(by: scanHeight)
        guard scanWidth > 0, scanHeight > 0, !overflow, scan.count == positions, !candidates.isEmpty else { return nil }
        let reach = settings.reachInvAngstrom, cut = settings.maximumVectorInvAngstrom
        let beam = settings.directBeamRadiusInvAngstrom

        let matrixCandidate = try PhaseOrientationCandidate(name: matrix.name, crystal: matrix.crystal, kMax: reach)
        let matrixSymmetry = CrystalOrientationSymmetry(laueClass: matrixCandidate.matcher.library.laueClass)
        var phases: [TwoPassPhase] = []
        for c in candidates {
            let candidate = try PhaseOrientationCandidate(name: c.name, crystal: c.crystal, kMax: reach,
                                                          excitationSlabInvAngstrom: c.excitationSlabInvAngstrom)
            phases.append(TwoPassPhase(candidate: candidate,
                                       symmetry: CrystalOrientationSymmetry(laueClass: candidate.matcher.library.laueClass)))
        }
        if isCancelled?() == true { return nil }

        // 1. the matrix of every position.
        let matrixProgress: (@Sendable (Double) -> Void)?
        if let progress { matrixProgress = { fraction in progress(0.6 * fraction) } } else { matrixProgress = nil }
        guard let raw = matrixOrientations(scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, matcher: matrixCandidate.matcher,
                                           originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                           directBeamRadiusInvAngstrom: beam, maximumVectorInvAngstrom: cut,
                                           maximumWorkers: maximumWorkers, isCancelled: isCancelled,
                                           progress: matrixProgress) else { return nil }
        // 2. the grain rule with the spot test.
        guard let regular = regularizingMatrixOrientations(raw, scan: scan, matcher: matrixCandidate.matcher,
                                                           originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                                           directBeamRadiusInvAngstrom: beam, maximumVectorInvAngstrom: cut,
                                                           neighbourhoodRadius: settings.neighbourhoodRadius,
                                                           explainedBy: matrixCandidate.reflections,
                                                           toleranceInvAngstrom: settings.matrixToleranceInvAngstrom) else { return nil }
        if isCancelled?() == true { return nil }
        progress?(0.7)
        // 3. the removal.
        guard let residual = removingMatrixSpots(scan: scan, field: regular.field, reflections: matrixCandidate.reflections,
                                                 originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                                 toleranceInvAngstrom: settings.matrixToleranceInvAngstrom,
                                                 directBeamRadiusInvAngstrom: beam, maximumVectorInvAngstrom: cut,
                                                 neighbouringGrains: MatrixNeighbourhood(radius: settings.neighbourhoodRadius,
                                                                                         symmetry: matrixSymmetry)) else { return nil }
        if isCancelled?() == true { return nil }
        progress?(0.75)
        // 4. the two passes on the whole pattern.
        let corrected = regular.field.correctedOrientations(inPlaneCorrectionsRad: residual.inPlaneCorrectionsRad)
        let search = PhaseOrientationSearchSettings(reachInvAngstrom: reach, pairRadiusInvAngstrom: settings.pairRadiusInvAngstrom,
                                                    residualCutoffInvAngstrom: settings.residualCutoffInvAngstrom,
                                                    directBeamRadiusInvAngstrom: beam, maximumVectorInvAngstrom: cut)
        guard let two = classifyInTwoPasses(residual: residual, matrixOrientations: corrected, phases: phases,
                                            matrixSymmetry: matrixSymmetry, settings: search,
                                            originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                            matrixReflections: matrixCandidate.reflections, scan: scan, scanWidth: scanWidth,
                                            usesNeighbourSupport: settings.usesNeighbourSupport,
                                            minimumPositions: settings.minimumPositions,
                                            matrixToleranceInvAngstrom: settings.matrixToleranceInvAngstrom,
                                            maximumWorkers: maximumWorkers, isCancelled: isCancelled) else { return nil }

        var map = PhaseMap(width: scanWidth, height: scanHeight, matrixEntryIndex: -1,
                           phaseNames: [matrix.name] + candidates.map(\.name), matrixPhaseIndex: 0)
        for (i, v) in two.verdicts.enumerated() {
            var r = PhaseVectorResult()
            r.verdict = v.verdict
            switch v.verdict {
            case .indexed: r.phaseIndex = Int32(v.phaseIndex + 1)
            case .matrix: r.phaseIndex = 0
            case .notIndexed, .noData: r.phaseIndex = -1
            }
            r.score = Float(v.score)
            r.matchedCount = Int32(v.matchedCount)
            r.survivingCount = Int32(v.survivingCount)
            r.removedCount = Int32(residual.removedCounts[i])
            if v.runnerUpPhaseIndex >= 0 {
                r.runnerUpPhaseIndex = Int32(v.runnerUpPhaseIndex + 1)
                r.runnerUpScore = Float(v.runnerUpScore)
            }
            map.results[i] = r
        }
        progress?(1)
        return StructurePhaseMapping(map: map, matrixField: regular.field, regularizedPositionCount: regular.changedPositions.count,
                                     matrixOrientations: corrected, residual: residual, classification: two,
                                     matrixSymmetry: matrixSymmetry, candidateSymmetries: phases.map(\.symmetry))
    }
}
