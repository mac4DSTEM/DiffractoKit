//
//  TwoPassClassification.swift
//  Role: the two passes of "find phase orientations from the structure alone" (DiffractoKit#81) as one call, scored on the
//        whole pattern, with the free-search verdict kept beside the final one.
//
//  The pieces are in PhaseOrientationSearch.swift (pass one, the guarded free search) and OrientationRelationshipLearning.swift
//  (learning, refinement, classification against the variants). What this adds is their order and three decisions:
//
//    1. the free search on every position's residual (the pattern after the matrix removal);
//    2. per phase, the relationships learned from its strong positions and refined on them; two relationships of one phase
//       that lie within their own spread of each other are one relationship and are merged;
//    3. every position classified against the variants, scored on the WHOLE pattern when it is given (the phase's spots
//       together with the matrix's): the removal serves the free search, the learning and the in-plane fit, not the verdict.
//       A precipitate that shares spots with the matrix loses them to the removal;
//       The whole pattern is the position's own: what the removal took for a neighbouring grain is left out of it
//       (`patternWithoutNeighbouringGrainSpots`, DiffractoKit#103);
//    4. a phase is called where at least two vectors lie on no matrix spot (the "firm" rule). Optionally a position with one
//       such vector is called too when one of its eight neighbours is firm with the same phase at the same orientation
//       (precipitates are connected; applied once, so it cannot creep).
//
//  Measured on Thronsen et al. dataset A, stride 3: see docs/parity.md. The one-vector positions are the precipitate edges;
//  no rule on them reaches both numbers of the typed zone axes (positions right, matrix positions called a phase). Tried
//  2026-11-11 and not done, each without gain there: a second refinement of the relationships on everything the classification
//  assigned; refining the matrix orientation per grain on its own spots.
//
//  Checked against planted truth (`TwoPassClassificationTests`, `TwoPassOtherMaterialsTests`). Not wired into
//  `PhaseVectorMatcher.map`.
//

import Foundation
import simd
import DiffractoCore

extension MatrixOrientationField {
    /// The orientations with the removal's in-plane correction applied (``MatrixResidualScan/inPlaneCorrectionsRad``): where the
    /// removal found the matrix's spots, not where the template's azimuthal bin put them. The relationships and the variants
    /// are built on the matrix orientation, so they need this one.
    ///
    /// - Returns: empty when the list is not one correction per position.
    public func correctedOrientations(inPlaneCorrectionsRad: [Double]) -> [simd_double3x3?] {
        guard inPlaneCorrectionsRad.count == orientations.count else { return [] }
        return zip(orientations, inPlaneCorrectionsRad).map { m, delta in
            guard let m else { return nil }
            // The spots turned by +delta are g·(cos δ x − sin δ y) and g·(sin δ x + cos δ y).
            let c = cos(delta), s = sin(delta)
            return simd_double3x3(columns: (c * m.columns.0 - s * m.columns.1, s * m.columns.0 + c * m.columns.1, m.columns.2))
        }
    }
}

/// A phase for ``PhaseVectorMatcher/classifyInTwoPasses(residual:matrixOrientations:phases:matrixSymmetry:settings:originX:originY:invAngstromPerPixel:matrixReflections:scan:scanWidth:usesNeighbourSupport:minimumPositions:maximumWorkers:isCancelled:)``.
public struct TwoPassPhase: Sendable {
    public let candidate: PhaseOrientationCandidate
    /// The phase's symmetry, by which its relationships are clustered and its variants counted.
    public let symmetry: CrystalOrientationSymmetry

    public init(candidate: PhaseOrientationCandidate, symmetry: CrystalOrientationSymmetry) {
        self.candidate = candidate
        self.symmetry = symmetry
    }
}

/// The result of the two passes.
public struct TwoPassClassification: Sendable {
    /// The final verdict of every position; `phaseIndex` is an index into the phases given.
    public let verdicts: [PhaseOrientationVerdict]
    /// The guarded free search of every position (pass one), kept: its orientation is the one the phase has HERE, found
    /// without any relationship.
    public let freeSearch: [PhaseOrientationVerdict]
    /// The relationships learned for every phase (empty for a phase that taught none; such a phase keeps its free-search
    /// verdicts).
    public let relationships: [[LearnedOrientationRelationship]]
    /// For every position the free search indexed: the angle between the relation it shows (its matrix orientation and the
    /// free-search orientation) and the nearest relationship learned for that phase, up to both symmetries, radians. Small
    /// where the phase keeps its relationship; large where the crystal is there at another orientation. NaN where the free
    /// search indexed nothing, the position has no matrix orientation, or the phase taught no relationship. A sparse
    /// pattern (four spots) can be indexed at a wrong orientation by the free search, so a large angle alone is not proof:
    /// read it together with the final verdict's score.
    public let relationAnglesRad: [Double]
    /// True where the verdict is a phase only because a neighbouring position is firmly that phase (the neighbour rule).
    public let neighbourSupported: [Bool]
}

extension PhaseVectorMatcher {

    /// Relationships of one phase that are one relationship: two whose relations lie within the sum of their spreads of each
    /// other, up to both symmetries, are merged into their mean, weighted by the positions they were learned from (the larger
    /// first). The spread is each cluster's own (``LearnedOrientationRelationship/spreadRad``), so no angle is set here. Measured
    /// on Thronsen et al. dataset A: theta' is learned as two relationships 53 and 45 positions strong; merged, 17 more
    /// positions are right and no more matrix positions are called a phase.
    public static func mergingCloseRelationships(_ relationships: [LearnedOrientationRelationship],
                                                 matrixSymmetry: CrystalOrientationSymmetry,
                                                 precipitateSymmetry: CrystalOrientationSymmetry) -> [LearnedOrientationRelationship] {
        var merged: [LearnedOrientationRelationship] = []
        for r in relationships.sorted(by: { $0.positionCount > $1.positionCount }) {
            var absorbed = false
            for k in merged.indices {
                let near = nearestEquivalent(of: r.relation, to: merged[k].relation, matrixSymmetry: matrixSymmetry,
                                             precipitateSymmetry: precipitateSymmetry)
                guard near.angleRad <= merged[k].spreadRad + r.spreadRad else { continue }
                let total = merged[k].positionCount + r.positionCount
                let t = total > 0 ? Double(r.positionCount) / Double(total) : 0.5
                let mean = simd_slerp(simd_quatd(merged[k].relation), simd_quatd(near.relation), t)
                merged[k] = LearnedOrientationRelationship(relation: simd_double3x3(mean), positionCount: total,
                                                           spreadRad: max(merged[k].spreadRad, r.spreadRad, near.angleRad),
                                                           matrixSymmetry: matrixSymmetry, precipitateSymmetry: precipitateSymmetry)
                absorbed = true
                break
            }
            if !absorbed { merged.append(r) }
        }
        return merged
    }

    /// The pattern of one position as a phase may be scored on it: the vectors the removal left, and those it took that lie
    /// on a spot of the position's OWN matrix orientation (a precipitate shares spots with its matrix, and loses them to
    /// the removal). What the removal took and the own orientation does not predict was taken for a neighbouring grain
    /// (``removingMatrixSpots(scan:field:reflections:originX:originY:invAngstromPerPixel:toleranceInvAngstrom:excitationSlabInvAngstrom:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:neighbouringGrains:)``):
    /// it is another matrix crystal's spot and no phase's vector (DiffractoKit#103: at a grain boundary of the simulated
    /// three-grain demo a phase was called on two of the neighbouring grain's six spots).
    ///
    /// - Parameters:
    ///   - all: every peak of the position, in detector pixels.
    ///   - left: the peaks the removal left (a subset of `all`, the same values).
    ///   - matrixOrientation: the position's matrix orientation, with the removal's in-plane correction applied.
    /// - Returns: `all` when nothing was removed or there is no orientation.
    public static func patternWithoutNeighbouringGrainSpots(all: [SIMD2<Double>], left: [SIMD2<Double>],
                                                            matrixOrientation: simd_double3x3?, matrixReflections: [ACOMReflection],
                                                            originX: Float, originY: Float, invAngstromPerPixel: Double,
                                                            matrixToleranceInvAngstrom: Double = 0.02,
                                                            matrixExcitationSlabInvAngstrom: Double = matrixRemovalSlabInvAngstrom,
                                                            directBeamRadiusInvAngstrom: Double = 0.15) -> [SIMD2<Double>] {
        guard let matrixOrientation, left.count < all.count else { return all }
        let spots = matrixSpots(reflections: matrixReflections, orientation: matrixOrientation,
                                excitationSlabInvAngstrom: matrixExcitationSlabInvAngstrom)
        let kept = Set(left)
        return all.filter { peak in
            if kept.contains(peak) { return true }
            let vector = PeakCoordinates.detectorPixels(originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel)
                .vector(of: peak)
            // The direct beam is nobody's spot; it is dropped where the vectors are made.
            if simd_length(vector) < directBeamRadiusInvAngstrom { return true }
            return spots.contains { simd_length($0 - vector) <= matrixToleranceInvAngstrom }
        }
    }

    /// A scan classified from the structures alone, in two passes (see the file header for the steps).
    ///
    /// - Parameters:
    ///   - residual: the scan after the matrix removal.
    ///   - matrixOrientations: the matrix orientation of every position, with the removal's in-plane correction applied
    ///     (``MatrixOrientationField/correctedOrientations(inPlaneCorrectionsRad:)``).
    ///   - settings: the free search's settings. The classification against the variants uses the same pair radius, cut-off
    ///     and vector limits without the chance guard: the variants are few, and that is today's known-variants rule.
    ///   - matrixReflections: when given, the phase-specific-reflection guard of the classification against the variants is
    ///     applied (the matrix orientations are already corrected, so no further correction is passed).
    ///   - scan: every peak of every position, the matrix's included. When given together with `matrixReflections`, the
    ///     variants are scored on the whole pattern and the firm rule is applied (file header, steps 3 and 4); without it
    ///     they are scored on the residual.
    ///   - scanWidth: the scan's width, for the neighbour rule.
    ///   - usesNeighbourSupport: the neighbour rule of step 4 (needs `scan`, `matrixReflections` and `scanWidth`). Off by
    ///     default: it finds more of a precipitate's edge and calls more matrix positions a phase, and a false precipitate
    ///     is the worse error.
    ///   - minimumPositions: the fewest strong positions a relationship is learned from.
    ///   - matrixToleranceInvAngstrom: how near a vector must lie to a predicted matrix spot to be the matrix's, Å⁻¹: the
    ///     tolerance the removal was run with (DiffractoKit#106). It decides which removed vectors are the position's own
    ///     matrix spots and which vectors of the whole pattern are a phase's own.
    /// - Returns: nil when the lengths disagree or the run was cancelled.
    public static func classifyInTwoPasses(residual: MatrixResidualScan, matrixOrientations: [simd_double3x3?],
                                           phases: [TwoPassPhase], matrixSymmetry: CrystalOrientationSymmetry,
                                           settings: PhaseOrientationSearchSettings,
                                           originX: Float, originY: Float, invAngstromPerPixel: Double,
                                           matrixReflections: [ACOMReflection]? = nil,
                                           scan: [[SIMD2<Double>]]? = nil, scanWidth: Int = 0,
                                           usesNeighbourSupport: Bool = false,
                                           minimumPositions: Int = 20,
                                           matrixToleranceInvAngstrom: Double = 0.02,
                                           maximumWorkers: Int? = nil,
                                           isCancelled: (@Sendable () -> Bool)? = nil,
                                           progress: (@Sendable (Double) -> Void)? = nil) -> TwoPassClassification? {
        let total = residual.scan.count
        guard matrixOrientations.count == total, residual.removedCounts.count == total else { return nil }
        if let scan, scan.count != total { return nil }
        let candidates = phases.map(\.candidate)
        var searchProgress: (@Sendable (Double) -> Void)?
        if let progress { searchProgress = { @Sendable fraction in progress(0.8 * fraction) } }
        guard var free = searchOrientations(residual: residual, candidates: candidates, settings: settings,
                                            originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                            maximumWorkers: maximumWorkers, isCancelled: isCancelled,
                                            progress: searchProgress) else { return nil }
        // A displaced matrix pattern: the removal took none of the position's vectors, yet at least half of the matrix's
        // predicted spots (and at least two) have a vector within three pair radii. The residual there is the matrix pattern
        // itself, a little off (a shifted origin, a wrong calibration), and a phase indexed on it is the matrix's spots read as
        // a precipitate. Such a position teaches no relationship and is never indexed (DiffractoKit#94). Needs the matrix's
        // reflections; without them nothing is judged.
        let hasUnusableMatrix: [Bool] = (0..<total).map { i in
            guard let reflections = matrixReflections, let m = matrixOrientations[i], residual.removedCounts[i] == 0 else { return false }
            let spots = matrixSpots(reflections: reflections, orientation: m, excitationSlabInvAngstrom: matrixRemovalSlabInvAngstrom)
            guard spots.count >= 2 else { return false }
            let vectors = residual.scan[i].compactMap {
                experimentalVector(peak: $0, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                   directBeamRadiusInvAngstrom: 0.15, maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            }
            let near = 3 * settings.pairRadiusInvAngstrom
            let found = spots.filter { s in vectors.contains { simd_length($0 - s) <= near } }.count
            return found >= 2 && 2 * found >= spots.count
        }
        for i in 0..<total where hasUnusableMatrix[i] && free[i].verdict == .indexed { free[i].verdict = .notIndexed }
        var passTwo = settings
        passTwo.chanceLevel = .infinity

        // 2. learned from the strong positions and refined on them; relationships within their own spread merged, and a
        // merged one refined again.
        var relationships: [[LearnedOrientationRelationship]] = []
        for k in phases.indices {
            let rough = learnOrientationRelationships(matrixOrientations: matrixOrientations, verdicts: free, candidateIndex: k,
                                                      matrixSymmetry: matrixSymmetry, precipitateSymmetry: phases[k].symmetry,
                                                      minimumPositions: minimumPositions)
            let strong = free.indices.filter { free[$0].verdict == .indexed && free[$0].phaseIndex == k }
            let inputs = strong.compactMap { i in matrixOrientations[i].map { (peaks: residual.scan[i], matrixOrientation: $0) } }
            func refined(_ r: LearnedOrientationRelationship) -> LearnedOrientationRelationship {
                refiningOrientationRelationship(r, positions: inputs, candidate: candidates[k], settings: passTwo,
                                                matrixSymmetry: matrixSymmetry, precipitateSymmetry: phases[k].symmetry,
                                                originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel)
            }
            let first = rough.map(refined)
            let merged = mergingCloseRelationships(first, matrixSymmetry: matrixSymmetry, precipitateSymmetry: phases[k].symmetry)
            relationships.append(merged.count == first.count ? first : merged.map(refined))
        }
        if isCancelled?() == true { return nil }
        progress?(0.85)

        // 3 and 4. every position against the variants.
        let sets = zip(candidates, relationships).map { PhaseVariantSet(candidate: $0, relationships: $1) }
        let onAllSpots = scan != nil && matrixReflections != nil
        var verdicts = [PhaseOrientationVerdict](repeating: PhaseOrientationVerdict(), count: total)
        var single = [PhaseOrientationVerdict?](repeating: nil, count: total)   // indexed on one own vector: for the neighbour rule
        let step = max(1, total / 100)
        for i in 0..<total {
            if let progress, i % step == 0 { progress(0.85 + 0.15 * Double(i) / Double(total)) }
            // The whole pattern, without what the removal took for a neighbouring grain (DiffractoKit#103).
            var whole: [SIMD2<Double>]? = scan?[i]
            if let all = whole, let reflections = matrixReflections {
                whole = patternWithoutNeighbouringGrainSpots(all: all, left: residual.scan[i], matrixOrientation: matrixOrientations[i],
                                                             matrixReflections: reflections, originX: originX, originY: originY,
                                                             invAngstromPerPixel: invAngstromPerPixel,
                                                             matrixToleranceInvAngstrom: matrixToleranceInvAngstrom,
                                                             directBeamRadiusInvAngstrom: passTwo.directBeamRadiusInvAngstrom)
            }
            var v = classifyAgainstVariants(peaks: residual.scan[i], removedCount: residual.removedCounts[i],
                                            matrixOrientation: matrixOrientations[i], phases: sets, settings: passTwo,
                                            originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                            matrixReflections: matrixReflections, matrixToleranceInvAngstrom: matrixToleranceInvAngstrom,
                                            allPeaks: whole, countsSingleOwnVector: true)
            if onAllSpots, matrixOrientations[i] != nil, v.verdict == .indexed, v.survivingCount <= passTwo.directMatrixMaximumVectors {
                single[i] = v
                v.verdict = .matrix
                v.orientation = nil
            }
            // A phase with no learned relationship stays on the guarded free search.
            if v.verdict != .indexed, free[i].verdict == .indexed, relationships[free[i].phaseIndex].isEmpty { v = free[i] }
            if hasUnusableMatrix[i] {
                single[i] = nil
                if v.verdict == .indexed { v.verdict = .notIndexed; v.orientation = nil }
            }
            verdicts[i] = v
        }
        var supported = [Bool](repeating: false, count: total)
        if usesNeighbourSupport, onAllSpots, scanWidth > 0, total % scanWidth == 0 {
            let height = total / scanWidth
            let firm = verdicts
            for y in 0..<height { for x in 0..<scanWidth {
                let i = y * scanWidth + x
                guard let mine = single[i], firm[i].verdict != .indexed, let orientation = mine.orientation else { continue }
                search: for dy in -1...1 { for dx in -1...1 where dx != 0 || dy != 0 {
                    let xx = x + dx, yy = y + dy
                    guard xx >= 0, yy >= 0, xx < scanWidth, yy < height else { continue }
                    let n = firm[yy * scanWidth + xx]
                    guard n.verdict == .indexed, n.phaseIndex == mine.phaseIndex, let theirs = n.orientation else { continue }
                    // The same variant: nearer than the angle within which two variants are one. A plain rotation angle: a
                    // neighbour indexed through another relationship of the phase, or fitted more than that angle away in
                    // the plane, does not support (the rule then calls fewer positions, never more).
                    let r = orientation.transpose * theirs
                    let angle = acos(min(1, max(-1, (r.columns.0.x + r.columns.1.y + r.columns.2.z - 1) / 2)))
                    if angle <= LearnedOrientationRelationship.defaultVariantMergeAngleRad {
                        verdicts[i] = mine
                        supported[i] = true
                        break search
                    }
                } }
            } }
        }

        // The free search's own orientation against the learned relationships.
        let angles: [Double] = free.indices.map { i in
            let f = free[i]
            guard f.verdict == .indexed, let m = matrixOrientations[i], let p = f.orientation, !relationships[f.phaseIndex].isEmpty else { return .nan }
            let observed = OrientationVariantSet.observedRelation(matrixOrientation: m, precipitateOrientation: p)
            return relationships[f.phaseIndex].map {
                relationDistance($0.relation, observed, matrixSymmetry: matrixSymmetry, precipitateSymmetry: phases[f.phaseIndex].symmetry)
            }.min() ?? .nan
        }
        return TwoPassClassification(verdicts: verdicts, freeSearch: free, relationships: relationships, relationAnglesRad: angles,
                                     neighbourSupported: supported)
    }
}
