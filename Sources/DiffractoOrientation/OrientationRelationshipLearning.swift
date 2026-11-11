//
//  OrientationRelationshipLearning.swift
//  Role: the second pass of "find phase orientations from the structure alone" (DiffractoKit#81, step 3). The positions with
//        strong evidence (the guarded free search of PhaseOrientationSearch.swift) teach the orientation relationship of a
//        phase to its matrix; every variant of that relationship follows from the matrix's symmetry; the rest of the scan is
//        classified against those variants only.
//
//  Why two passes: a free search over the whole orientation space cannot decide a position that carries two to five matched
//  vectors (measured: 81.6 % with the chance guard and 94.5 % without, on Thronsen et al. dataset A, against 98.7 % for typed
//  zone axes). Typed zone axes are a handful of references; that is what makes little evidence enough. A precipitate has a
//  definite orientation relationship to its matrix, so the handful can be learned: this is the paper's recipe with the
//  references found from the data. A relationship learned from the variants that show rich patterns also gives the variants
//  that never do (θ′ face-on from θ′ edge-on), because they are the same relationship under the matrix's symmetry, the way
//  orix and MTEX generate variants.
//
//  No Python call does this as one; the variant enumeration is `OrientationRelationship.variants` (group theory, on orix's
//  `Misorientation.equivalent`). Checked against planted truth (`OrientationRelationshipLearningTests`).
//

import Foundation
import simd
import DiffractoCore

/// An orientation relationship found from a scan, with its variants.
public struct LearnedOrientationRelationship: Sendable {
    /// Precipitate Cartesian vectors to matrix Cartesian vectors (the convention of ``OrientationVariantSet/matrix``): the
    /// mean of the cluster's observed relations.
    public let relation: simd_double3x3
    /// How many strong positions share it.
    public let positionCount: Int
    /// The median distance of those positions' relations to `relation`, radians.
    public let spreadRad: Double
    /// One precipitate-to-matrix rotation per variant: ``OrientationRelationship/variants(of:matrix:precipitate:)``, with
    /// variants closer than `variantMergeAngleRad` (up to the precipitate's symmetry) counted once. A measured relation is
    /// never exactly on a symmetric one, so a relationship with twelve variants would otherwise come back as twenty-four in
    /// near-identical pairs.
    public let variants: [simd_double3x3]

    public init(relation: simd_double3x3, positionCount: Int, spreadRad: Double,
                matrixSymmetry: CrystalOrientationSymmetry, precipitateSymmetry: CrystalOrientationSymmetry,
                variantMergeAngleRad: Double = 3 * Double.pi / 180) {
        self.relation = relation
        self.positionCount = positionCount
        self.spreadRad = spreadRad
        var merged: [simd_double3x3] = []
        for v in OrientationRelationship.variants(of: relation, matrix: matrixSymmetry, precipitate: precipitateSymmetry) {
            let known = merged.contains { m in
                precipitateSymmetry.operators.contains { op in
                    let t = m.transpose * v * op
                    return acos(max(-1, min(1, (t[0][0] + t[1][1] + t[2][2] - 1) / 2))) < variantMergeAngleRad
                }
            }
            if !known { merged.append(v) }
        }
        self.variants = merged
    }
}

/// A phase and the relationships it was found to take: what the second pass looks for.
public struct PhaseVariantSet: Sendable {
    /// The phase: its reflections and its slab (its templates are not used by the second pass).
    public let candidate: PhaseOrientationCandidate
    /// The relationships learned for it; empty for a phase that taught none (it is then not looked for by the second pass).
    public let relationships: [LearnedOrientationRelationship]

    public init(candidate: PhaseOrientationCandidate, relationships: [LearnedOrientationRelationship]) {
        self.candidate = candidate
        self.relationships = relationships
    }
}

extension PhaseVectorMatcher {

    /// The angle between two precipitate-to-matrix relations up to both symmetries: the smallest rotation angle of
    /// `a⁻¹ g b s` over the matrix operators g and the precipitate operators s, radians. 0 for two statements of the same
    /// relationship or of two of its variants.
    public static func relationDistance(_ a: simd_double3x3, _ b: simd_double3x3,
                                        matrixSymmetry: CrystalOrientationSymmetry,
                                        precipitateSymmetry: CrystalOrientationSymmetry) -> Double {
        var best = Double.pi
        let at = a.transpose
        for g in matrixSymmetry.operators {
            let left = at * g * b
            for s in precipitateSymmetry.operators {
                let m = left * s
                let cosine = max(-1, min(1, (m[0][0] + m[1][1] + m[2][2] - 1) / 2))
                best = min(best, acos(cosine))
            }
        }
        return best
    }

    /// The statement of relation `b` that is nearest to `a` among its equivalents `g b s` (g a matrix operator, s a
    /// precipitate operator), and its angle to `a`, radians.
    static func nearestEquivalent(of b: simd_double3x3, to a: simd_double3x3,
                                  matrixSymmetry: CrystalOrientationSymmetry,
                                  precipitateSymmetry: CrystalOrientationSymmetry) -> (relation: simd_double3x3, angleRad: Double) {
        var best = b, bestAngle = Double.infinity
        let at = a.transpose
        for g in matrixSymmetry.operators {
            let gb = g * b
            for s in precipitateSymmetry.operators {
                let candidate = gb * s
                let m = at * candidate
                let angle = acos(max(-1, min(1, (m[0][0] + m[1][1] + m[2][2] - 1) / 2)))
                if angle < bestAngle { bestAngle = angle; best = candidate }
            }
        }
        return (best, bestAngle)
    }

    /// The orientation relationships a phase takes, from the positions the guarded search indexed as that phase.
    ///
    /// Every such position gives one observed relation, matrix orientation times precipitate orientation transposed
    /// (``OrientationVariantSet/observedRelation(matrixOrientation:precipitateOrientation:)``). The observations are taken in
    /// order of evidence (most matched vectors first, then the better score, then the position) and clustered greedily: an
    /// observation joins the first cluster whose relation is within `clusterAngleRad` of it up to both symmetries
    /// (``relationDistance(_:_:matrixSymmetry:precipitateSymmetry:)``), else it starts one. A cluster's relation is the mean
    /// of its members' (each in the statement nearest to the cluster's first member; the quaternion mean), which takes out
    /// the template grid's few degrees of scatter. A cluster with fewer than
    /// `minimumPositions` positions, or with less than `minimumShare` of the observations, teaches nothing.
    ///
    /// Only positions with at least `minimumMatchedVectors` matched vectors are observations: a pattern of a few spots fits
    /// several orientations, and a relationship read from one of them is not the phase's.
    ///
    /// A zero-order pattern does not fix an orientation completely (the half turn about the beam, and the mirror pass of the
    /// template match), so one physical relationship can come back as two or more clusters when the matrix zone has no
    /// two-fold axis along the beam. Their variants predict the same patterns; nothing is lost by keeping them.
    ///
    /// - Parameters:
    ///   - matrixOrientations: the matrix orientation of every position (``MatrixOrientationField/orientations``).
    ///   - verdicts: the first pass, one per position.
    ///   - candidateIndex: the phase, as an index into the first pass's candidates.
    /// - Returns: the relationships, most positions first; empty when the two arrays differ in length.
    public static func learnOrientationRelationships(matrixOrientations: [simd_double3x3?], verdicts: [PhaseOrientationVerdict],
                                                     candidateIndex: Int,
                                                     matrixSymmetry: CrystalOrientationSymmetry,
                                                     precipitateSymmetry: CrystalOrientationSymmetry,
                                                     clusterAngleRad: Double = 5 * Double.pi / 180,
                                                     minimumPositions: Int = 20, minimumMatchedVectors: Int = 8,
                                                     minimumShare: Double = 0.1) -> [LearnedOrientationRelationship] {
        guard matrixOrientations.count == verdicts.count else { return [] }
        var observations: [(position: Int, relation: simd_double3x3, matched: Int, score: Double)] = []
        for (i, v) in verdicts.enumerated() where v.verdict == .indexed && v.phaseIndex == candidateIndex
            && v.matchedCount >= minimumMatchedVectors {
            guard let m = matrixOrientations[i], let p = v.orientation else { continue }
            observations.append((i, OrientationVariantSet.observedRelation(matrixOrientation: m, precipitateOrientation: p), v.matchedCount, v.score))
        }
        observations.sort { a, b in
            if a.matched != b.matched { return a.matched > b.matched }
            if a.score != b.score { return a.score < b.score }
            return a.position < b.position
        }
        var clusters: [(first: simd_double3x3, members: [simd_double3x3])] = []
        for o in observations {
            var joined = false
            for k in clusters.indices {
                let nearest = nearestEquivalent(of: o.relation, to: clusters[k].first, matrixSymmetry: matrixSymmetry,
                                                precipitateSymmetry: precipitateSymmetry)
                if nearest.angleRad <= clusterAngleRad { clusters[k].members.append(nearest.relation); joined = true; break }
            }
            if !joined { clusters.append((o.relation, [o.relation])) }
        }
        let share = Int((minimumShare * Double(observations.count)).rounded(.up))
        return clusters.filter { $0.members.count >= max(1, minimumPositions, share) }
            .sorted { $0.members.count > $1.members.count }
            .map { c in
                // The quaternion mean of the members, each on the first member's side.
                let reference = simd_quatd(c.first)
                var sum = SIMD4<Double>(repeating: 0)
                for m in c.members {
                    let v = simd_quatd(m).vector
                    sum += simd_dot(v, reference.vector) < 0 ? -v : v
                }
                let mean = simd_length(sum) > 0 ? simd_double3x3(simd_quatd(vector: simd_normalize(sum))) : c.first
                let distances = c.members.map { m -> Double in
                    let t = mean.transpose * m
                    return acos(max(-1, min(1, (t[0][0] + t[1][1] + t[2][2] - 1) / 2)))
                }.sorted()
                return LearnedOrientationRelationship(relation: mean, positionCount: distances.count, spreadRad: distances[distances.count / 2],
                                                      matrixSymmetry: matrixSymmetry, precipitateSymmetry: precipitateSymmetry)
            }
    }

    /// The known-variants score of one position's vectors against the best variant of a relation (the score of
    /// ``classifyAgainstVariants(peaks:removedCount:matrixOrientation:phases:settings:originX:originY:invAngstromPerPixel:)``),
    /// capped at the cut-off; the cut-off for a position no variant has spots for.
    static func bestVariantScore(vectors: [SIMD2<Double>], matrixOrientation: simd_double3x3, variants: [simd_double3x3],
                                 candidate: PhaseOrientationCandidate, settings: PhaseOrientationSearchSettings) -> Double {
        var best = settings.residualCutoffInvAngstrom
        for variant in variants {
            let predicted = matrixSpots(reflections: candidate.reflections, orientation: variant.transpose * matrixOrientation,
                                        excitationSlabInvAngstrom: candidate.excitationSlabInvAngstrom)
            guard !predicted.isEmpty else { continue }
            var sum = 0.0
            var hits = Set<Int>()
            for u in vectors {
                var nearest = 0, distance = Double.infinity
                for (k, s) in predicted.enumerated() { let d = simd_length(u - s); if d < distance { distance = d; nearest = k } }
                sum += distance
                hits.insert(nearest)
            }
            best = min(best, sum / Double(hits.count))
        }
        return best
    }

    /// A learned relation turned, by a few degrees at most, to where its variants fit the strong positions' spots best.
    ///
    /// The relation read from the template match carries the template grid's error in the beam direction (a few degrees),
    /// and a mean over positions does not remove an error they share. The spots do: the variants of one relation are seen
    /// in different projections at different positions, so the sum over the positions of the best variant's score (without
    /// the in-plane fit, which would hide the error) has its minimum at the true relation in all three angles. The
    /// minimum is searched by coordinate descent over the three components of a small rotation, in steps of 2°, 1°, 0.5°
    /// and 0.25°, within `limitRad` of the input. A step is taken only when it lowers the sum by at least `minimumGain`
    /// of its current value: where the positions show the variants in one kind of projection only, some angle is not fixed
    /// by the spots, and the relation must not drift along it on rounding noise. Measured 2026-11-11: with 2 % the
    /// synthetic demo's relation stays at 0.04° from the planted one and Thronsen et al. dataset A classifies 98.3 to
    /// 98.4 %; with 0.5 % the demo drifts to 0.5° and Thronsen falls to 98.29 %; with no threshold the demo drifts to 2.4°.
    ///
    /// - Parameters:
    ///   - positions: residual peaks and matrix orientation of the positions that taught the relation (at most
    ///     `maximumPositions` are used, evenly spread over the array).
    /// - Returns: the relation with the same position count and spread, its variants regenerated; the input when there is
    ///   nothing to refine against.
    public static func refiningOrientationRelationship(_ relationship: LearnedOrientationRelationship,
                                                       positions: [(peaks: [SIMD2<Double>], matrixOrientation: simd_double3x3)],
                                                       candidate: PhaseOrientationCandidate,
                                                       settings: PhaseOrientationSearchSettings,
                                                       matrixSymmetry: CrystalOrientationSymmetry,
                                                       precipitateSymmetry: CrystalOrientationSymmetry,
                                                       originX: Float, originY: Float, invAngstromPerPixel: Double,
                                                       maximumPositions: Int = 300,
                                                       limitRad: Double = 8 * Double.pi / 180,
                                                       minimumGain: Double = 0.02) -> LearnedOrientationRelationship {
        let stride = max(1, (positions.count + max(1, maximumPositions) - 1) / max(1, maximumPositions))
        var sample: [(vectors: [SIMD2<Double>], matrix: simd_double3x3)] = []
        for k in Swift.stride(from: 0, to: positions.count, by: stride) {
            let v = experimentalVectors(peaks: positions[k].peaks, originX: originX, originY: originY,
                                        invAngstromPerPixel: invAngstromPerPixel,
                                        directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                                        maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            if v.count > settings.directMatrixMaximumVectors { sample.append((v, positions[k].matrixOrientation)) }
        }
        guard !sample.isEmpty else { return relationship }
        func rotation(_ w: SIMD3<Double>) -> simd_double3x3 {
            let angle = simd_length(w)
            return angle > 0 ? simd_double3x3(simd_quatd(angle: angle, axis: w / angle)) : matrix_identity_double3x3
        }
        func cost(_ w: SIMD3<Double>) -> Double {
            let relation = rotation(w) * relationship.relation
            let variants = OrientationRelationship.variants(of: relation, matrix: matrixSymmetry, precipitate: precipitateSymmetry)
            return sample.reduce(0) { $0 + bestVariantScore(vectors: $1.vectors, matrixOrientation: $1.matrix, variants: variants,
                                                            candidate: candidate, settings: settings) }
        }
        var w = SIMD3<Double>(repeating: 0)
        var best = cost(w)
        for stepDeg in [2.0, 1.0, 0.5, 0.25] {
            let step = stepDeg * .pi / 180
            var improved = true
            while improved {
                improved = false
                for axis in 0..<3 {
                    for sign in [1.0, -1.0] {
                        var trial = w
                        trial[axis] += sign * step
                        guard simd_length(trial) <= limitRad else { continue }
                        let c = cost(trial)
                        if c < best * (1 - minimumGain) { best = c; w = trial; improved = true }
                    }
                }
            }
        }
        return LearnedOrientationRelationship(relation: rotation(w) * relationship.relation, positionCount: relationship.positionCount,
                                              spreadRad: relationship.spreadRad, matrixSymmetry: matrixSymmetry,
                                              precipitateSymmetry: precipitateSymmetry)
    }

    /// One position decided against the variants of the learned relationships.
    ///
    /// For every phase, relationship and variant: the precipitate's orientation here is the variant transposed times the
    /// position's matrix orientation; its spots are predicted, turned by the in-plane fit, and scored with the known-variants
    /// score of ``searchOrientations(peaks:removedCount:candidates:settings:originX:originY:invAngstromPerPixel:)``. The lowest
    /// score wins; it is indexed within the cut-off and the chance guard, and the guard counts the variants tried here
    /// (`orientationsTried`), not a template library: that is what lets a few matched vectors be evidence.
    ///
    /// - Parameters:
    ///   - matrixOrientation: this position's matrix orientation; nil (no matrix found here) leaves anything searchable
    ///     not indexed.
    ///   - phases: the result's `phaseIndex` is an index into this array.
    ///   - matrixReflections: the matrix crystal's reflections. When given, an indexed winner must have at least
    ///     `minimumSpecificReflections` matched spots that are its own: predicted spots farther than
    ///     `matrixToleranceInvAngstrom` from every spot the matrix has at this position. A winner without one is the matrix
    ///     (the evidence guard of `PhaseVectorSettings.knownVariantsMinimumSpecificReflections`: a precipitate that shares
    ///     its spot positions with the matrix is not shown by a vector the matrix removal just missed). Such a position's
    ///     verdict is `.matrix` whatever `removedCount` is; the refused winner's phase, score and counts stay for diagnosis
    ///     and its orientation is cleared.
    ///   - matrixInPlaneCorrectionRad: the in-plane correction the removal fitted here
    ///     (``MatrixResidualScan/inPlaneCorrectionsRad``), so the matrix's spots are where the removal had them.
    ///   - matrixExcitationSlabInvAngstrom: the slab the removal used.
    public static func classifyAgainstVariants(peaks: [SIMD2<Double>], removedCount: Int, matrixOrientation: simd_double3x3?,
                                               phases: [PhaseVariantSet], settings: PhaseOrientationSearchSettings,
                                               originX: Float, originY: Float, invAngstromPerPixel: Double,
                                               matrixReflections: [ACOMReflection]? = nil,
                                               matrixInPlaneCorrectionRad: Double = 0,
                                               matrixExcitationSlabInvAngstrom: Double = matrixRemovalSlabInvAngstrom,
                                               matrixToleranceInvAngstrom: Double = 0.02,
                                               minimumSpecificReflections: Int = 1) -> PhaseOrientationVerdict {
        var out = PhaseOrientationVerdict()
        let vectors = experimentalVectors(peaks: peaks, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                          directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                                          maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
        out.survivingCount = vectors.count
        out.orientationsTried = phases.reduce(0) { $0 + $1.relationships.reduce(0) { $0 + $1.variants.count } }
        if vectors.count <= settings.directMatrixMaximumVectors {
            out.verdict = removedCount > 0 ? .matrix : (vectors.isEmpty ? .noData : .notIndexed)
            return out
        }
        guard let matrixOrientation else { out.verdict = .notIndexed; return out }

        struct Scored { var phase: Int; var orientation: simd_double3x3; var score: Double; var matched: Int; var spots: Int; var specific: Int }
        let matrixSpotsHere: [SIMD2<Double>]? = matrixReflections.map {
            let c = cos(matrixInPlaneCorrectionRad), s = sin(matrixInPlaneCorrectionRad)
            return matrixSpots(reflections: $0, orientation: matrixOrientation, excitationSlabInvAngstrom: matrixExcitationSlabInvAngstrom)
                .map { SIMD2(c * $0.x - s * $0.y, s * $0.x + c * $0.y) }
        }
        var bestPerPhase: [Scored] = []
        for (phase, set) in phases.enumerated() {
            var best: Scored?
            for relationship in set.relationships {
                for variant in relationship.variants {
                    let m = variant.transpose * matrixOrientation
                    let predicted = matrixSpots(reflections: set.candidate.reflections, orientation: m,
                                                excitationSlabInvAngstrom: set.candidate.excitationSlabInvAngstrom)
                    guard !predicted.isEmpty else { continue }
                    let fitted = fittingInPlane(spots: predicted, vectors: vectors, toleranceInvAngstrom: settings.pairRadiusInvAngstrom)
                    var sum = 0.0, matched = 0
                    var hits = Set<Int>()
                    for u in vectors {
                        var nearest = 0, distance = Double.infinity
                        for (k, s) in fitted.spots.enumerated() { let d = simd_length(u - s); if d < distance { distance = d; nearest = k } }
                        sum += distance
                        hits.insert(nearest)
                        if distance <= settings.pairRadiusInvAngstrom { matched += 1 }
                    }
                    let score = sum / Double(hits.count)
                    if best == nil || score < best!.score {
                        let c = cos(fitted.rotationRad), s = sin(fitted.rotationRad)
                        let turned = simd_double3x3(columns: (c * m.columns.0 - s * m.columns.1, s * m.columns.0 + c * m.columns.1, m.columns.2))
                        // The distinct predicted spots hit within the pair radius that are not on a matrix spot.
                        var specific = 0
                        if let matrixSpotsHere {
                            var own = Set<Int>()
                            for u in vectors {
                                var nearest = 0, distance = Double.infinity
                                for (k, sp) in fitted.spots.enumerated() { let d = simd_length(u - sp); if d < distance { distance = d; nearest = k } }
                                guard distance <= settings.pairRadiusInvAngstrom else { continue }
                                if !matrixSpotsHere.contains(where: { simd_length($0 - fitted.spots[nearest]) <= matrixToleranceInvAngstrom }) { own.insert(nearest) }
                            }
                            specific = own.count
                        }
                        best = Scored(phase: phase, orientation: turned, score: score, matched: matched, spots: fitted.spots.count, specific: specific)
                    }
                }
            }
            if let best { bestPerPhase.append(best) }
        }
        guard var winner = bestPerPhase.first else { out.verdict = .notIndexed; return out }
        for s in bestPerPhase.dropFirst() where s.score < winner.score { winner = s }
        var runnerUp: Scored?
        for s in bestPerPhase where s.phase != winner.phase { if runnerUp == nil || s.score < runnerUp!.score { runnerUp = s } }
        out.phaseIndex = winner.phase
        out.orientation = winner.orientation
        out.score = winner.score
        out.matchedCount = winner.matched
        if let runnerUp { out.runnerUpPhaseIndex = runnerUp.phase; out.runnerUpScore = runnerUp.score }
        out.requiredMatches = requiredMatches(vectors: vectors.count, spots: winner.spots, orientationsTried: out.orientationsTried,
                                              settings: settings)
        out.verdict = winner.score <= settings.residualCutoffInvAngstrom && winner.matched >= out.requiredMatches ? .indexed : .notIndexed
        if out.verdict == .indexed, matrixSpotsHere != nil, minimumSpecificReflections > 0, winner.specific < minimumSpecificReflections {
            out.verdict = .matrix
            out.orientation = nil
        }
        return out
    }
}
