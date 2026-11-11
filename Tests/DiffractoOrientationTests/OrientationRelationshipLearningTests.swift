//
//  OrientationRelationshipLearningTests.swift
//  Step 3 of "find phase orientations from the structure alone" (DiffractoKit#81): the positions with strong evidence teach the
//  orientation relationship of a phase to its matrix; every variant follows from the matrix's symmetry; the rest of the scan is
//  classified against those variants, where two or three matched vectors are evidence again. Planted truth: beta'' in aluminium
//  with the relationship of Andersen et al. 1998, [010]beta'' || [001]Al and [001]beta'' || [310]Al.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct OrientationRelationshipLearningTests {
    static let q = 0.012
    static let origin: Float = 63.5
    static let reach = 0.75

    struct Setup {
        let candidate: PhaseOrientationCandidate
        let matrixSymmetry: CrystalOrientationSymmetry
        let precipitateSymmetry: CrystalOrientationSymmetry
        let relation: simd_double3x3              // planted, precipitate to matrix
        let variants: [simd_double3x3]
        let matrixOrientation: simd_double3x3     // Al [001] at 12 degrees
    }

    static func setup() throws -> Setup {
        let beta = CrystalPresets.betaDoublePrime, al = CrystalPresets.aluminum
        let relation = try OrientationRelationship.rotation(
            OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
            OrientationRelationship(candidate: .direction(SIMD3(0, 0, 1)), matrix: .direction(SIMD3(3, 1, 0))),
            precipitate: beta, matrix: al)
        let ms = CrystalOrientationSymmetry(laueClass: try #require(al.laueClass))
        let ps = CrystalOrientationSymmetry(laueClass: try #require(beta.laueClass))
        let variants = OrientationRelationship.variants(of: relation, matrix: ms, precipitate: ps)
        let m = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1)), inPlaneAngle: Float(12 * Double.pi / 180))
        return Setup(candidate: try PhaseOrientationCandidate(name: "beta''", crystal: beta, kMax: reach),
                     matrixSymmetry: ms, precipitateSymmetry: ps, relation: relation, variants: variants, matrixOrientation: m)
    }

    /// The residual peaks of a precipitate at one variant under the matrix orientation: its spots in a thin slab, the strongest
    /// `count` of them by nearness to the centre (nil: all).
    static func peaks(_ s: Setup, variant: simd_double3x3, count: Int? = nil, turnedDeg: Double = 0) -> [SIMD2<Double>] {
        let orientation = variant.transpose * s.matrixOrientation
        var spots = PhaseVectorMatcher.matrixSpots(reflections: s.candidate.reflections, orientation: orientation, excitationSlabInvAngstrom: 0.02)
            .filter { simd_length($0) > 0.15 && simd_length($0) < reach }
            .sorted { simd_length($0) < simd_length($1) }
        if let count { spots = Array(spots.prefix(count)) }
        let c = cos(turnedDeg * .pi / 180), sn = sin(turnedDeg * .pi / 180)
        return spots.map { SIMD2(c * $0.x - sn * $0.y, sn * $0.x + c * $0.y) }
            .map { SIMD2(Double(origin) + $0.x / q, Double(origin) + $0.y / q) }
    }

    /// The variants a [001] matrix grain shows with rich patterns: those whose pattern holds at least 12 spots.
    static func richVariants(_ s: Setup) -> [simd_double3x3] {
        s.variants.filter { peaks(s, variant: $0).count >= 12 }
    }

    /// The relationship is learned from the strong positions, and it is the planted one: one learned relationship, within 3
    /// degrees of Andersen's up to both symmetries, with all twelve variants. Fewer strong positions than the minimum teach
    /// nothing. Mutation: the precipitate's orientation taken for the relationship without the matrix's (the distance is large).
    @Test func theRelationshipIsLearnedFromTheStrongPositions() throws {
        let s = try Self.setup()
        let rich = Self.richVariants(s)
        #expect(rich.count >= 2, "rich variants \(rich.count) of \(s.variants.count)")
        var scan: [[SIMD2<Double>]] = []
        for k in 0..<24 { scan.append(Self.peaks(s, variant: rich[k % rich.count])) }
        let residual = MatrixResidualScan(scan: scan, removedCounts: [Int](repeating: 8, count: scan.count),
                                          inPlaneCorrectionsRad: [Double](repeating: 0, count: scan.count))
        let settings = PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach)
        let verdicts = try #require(PhaseVectorMatcher.searchOrientations(
            residual: residual, candidates: [s.candidate], settings: settings,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(verdicts.allSatisfy { $0.verdict == .indexed })
        let matrix = [simd_double3x3?](repeating: s.matrixOrientation, count: scan.count)
        let learned = PhaseVectorMatcher.learnOrientationRelationships(
            matrixOrientations: matrix, verdicts: verdicts, candidateIndex: 0,
            matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry, minimumPositions: 10)
        #expect(learned.count >= 1)
        let best = try #require(learned.first)
        let distance = PhaseVectorMatcher.relationDistance(best.relation, s.relation, matrixSymmetry: s.matrixSymmetry,
                                                           precipitateSymmetry: s.precipitateSymmetry) * 180 / .pi
        #expect(distance < 3, "the learned relationship is \(distance) deg from the planted one")
        #expect(best.variants.count == 12)
        #expect(learned.reduce(0) { $0 + $1.positionCount } <= scan.count)
        let none = PhaseVectorMatcher.learnOrientationRelationships(
            matrixOrientations: matrix, verdicts: verdicts, candidateIndex: 0,
            matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry, minimumPositions: 100)
        #expect(none.isEmpty)
    }

    /// A relation 4 degrees off the planted one (the scatter the template grid leaves) is refined against the strong positions'
    /// spots to within 0.5 degrees: the variants are seen in different projections, which fixes all three angles.
    /// Mutation: the refinement returning its input (the first expectation).
    @Test func theRelationshipIsRefinedAgainstTheSpots() throws {
        let s = try Self.setup()
        let rich = Self.richVariants(s)
        let positions = (0..<24).map { k in (peaks: Self.peaks(s, variant: rich[k % rich.count]), matrixOrientation: s.matrixOrientation) }
        // Off by 4 degrees about an oblique axis.
        let axis = simd_normalize(SIMD3<Double>(1, 2, 3))
        let off = simd_double3x3(simd_quatd(angle: 4 * .pi / 180, axis: axis)) * s.relation
        let rough = LearnedOrientationRelationship(relation: off, positionCount: 24, spreadRad: 0,
                                                   matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry)
        let before = PhaseVectorMatcher.relationDistance(rough.relation, s.relation, matrixSymmetry: s.matrixSymmetry,
                                                         precipitateSymmetry: s.precipitateSymmetry) * 180 / .pi
        #expect(before > 3.5)
        let refined = PhaseVectorMatcher.refiningOrientationRelationship(
            rough, positions: positions, candidate: s.candidate, settings: PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach),
            matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        let after = PhaseVectorMatcher.relationDistance(refined.relation, s.relation, matrixSymmetry: s.matrixSymmetry,
                                                        precipitateSymmetry: s.precipitateSymmetry) * 180 / .pi
        #expect(after < 0.5, "refined to \(after) deg from the planted relation (from \(before))")
        #expect(refined.positionCount == rough.positionCount && refined.variants.count == 12)
        // Nothing to refine against: the input comes back.
        let same = PhaseVectorMatcher.refiningOrientationRelationship(
            rough, positions: [], candidate: s.candidate, settings: PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach),
            matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        #expect(PhaseVectorMatcher.relationDistance(same.relation, rough.relation, matrixSymmetry: s.matrixSymmetry,
                                                    precipitateSymmetry: s.precipitateSymmetry) < 1e-9)
    }

    /// Against the learned variants a weak position is decided: three spots of a variant are too few for the free search
    /// (not indexed there) and are beta'' against the variants. Three spots of the same crystal at an orientation that is no
    /// variant (turned 45 degrees) stay not indexed. Another variant of the same relationship is found too.
    /// Mutations: the variants built without the matrix orientation of the position; the chance guard counting the whole
    /// template search instead of the variants tried (the weak position is refused).
    @Test func weakPositionsAreDecidedAgainstTheVariants() throws {
        let s = try Self.setup()
        let rich = Self.richVariants(s)
        // Learned from the planted relationship itself here (the learning has its own test).
        let relationship = LearnedOrientationRelationship(relation: s.relation, positionCount: 24, spreadRad: 0,
                                                          matrixSymmetry: s.matrixSymmetry, precipitateSymmetry: s.precipitateSymmetry)
        #expect(relationship.variants.count == 12)
        let settings = PhaseOrientationSearchSettings(reachInvAngstrom: Self.reach)
        let phases = [PhaseVariantSet(candidate: s.candidate, relationships: [relationship])]
        func same(_ a: simd_double3x3, _ b: simd_double3x3) -> Bool {
            (0..<3).allSatisfy { c in (0..<3).allSatisfy { r in abs(a[c][r] - b[c][r]) < 1e-9 } }
        }
        // Another variant than the first: every variant of the relationship is looked for. (The relationship is planted
        // here; that a variant never seen strongly is found is measured on Thronsen dataset A, theta' face-on.)
        let unseen = try #require(s.variants.last { v in !same(rich[0], v) && Self.peaks(s, variant: v).count >= 3 })
        for variant in [rich[0], unseen] {
            let weak = Self.peaks(s, variant: variant, count: 3)
            #expect(weak.count == 3)
            let free = PhaseVectorMatcher.searchOrientations(peaks: weak, removedCount: 8, candidates: [s.candidate], settings: settings,
                                                             originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
            #expect(free.verdict == .notIndexed, "the free search must refuse three spots: \(free)")
            let v = PhaseVectorMatcher.classifyAgainstVariants(peaks: weak, removedCount: 8, matrixOrientation: s.matrixOrientation,
                                                               phases: phases, settings: settings,
                                                               originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
            #expect(v.verdict == .indexed && v.phaseIndex == 0, "\(v)")
            #expect(v.matchedCount == 3 && v.score < 0.005)
            #expect(v.orientationsTried == 12)
        }
        let stray = Self.peaks(s, variant: rich[0], count: 3, turnedDeg: 45)
        let refused = PhaseVectorMatcher.classifyAgainstVariants(peaks: stray, removedCount: 8, matrixOrientation: s.matrixOrientation,
                                                                 phases: phases, settings: settings,
                                                                 originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        #expect(refused.verdict == .notIndexed, "\(refused)")
        // The phase-specific-reflection guard. A "matrix" whose spots here are exactly these three vectors: the winner has no
        // spot of its own and the position is the matrix. With the matrix turned 3 degrees and the removal's correction of
        // -3 degrees given, the same; with the correction left out the matrix's spots are 0.03 Å⁻¹ away at 0.6 Å⁻¹ and the
        // guard would miss. Mutation: the correction ignored (the second expectation).
        do {
            let vectors = PhaseVectorMatcher.experimentalVectors(
                peaks: Self.peaks(s, variant: rich[0]), originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                directBeamRadiusInvAngstrom: 0.15).sorted { simd_length($0) > simd_length($1) }
            let far = Array(vectors.prefix(3))
            #expect(far.allSatisfy { simd_length($0) > 0.5 })
            let farPeaks = far.map { SIMD2(Double(Self.origin) + $0.x / Self.q, Double(Self.origin) + $0.y / Self.q) }
            func fakeMatrix(turnedDeg: Double) -> [ACOMReflection] {
                let c = cos(turnedDeg * .pi / 180), sn = sin(turnedDeg * .pi / 180)
                return far.flatMap { u -> [ACOMReflection] in
                    let v = SIMD2(c * u.x - sn * u.y, sn * u.x + c * u.y)
                    let g = v.x * s.matrixOrientation.columns.0 + v.y * s.matrixOrientation.columns.1
                    return [ACOMReflection(g: g, intensity: 1), ACOMReflection(g: -g, intensity: 1)]
                }
            }
            func decide(_ reflections: [ACOMReflection]?, correctionDeg: Double = 0) -> PhaseOrientationVerdict {
                var open = settings
                open.chanceLevel = .infinity
                return PhaseVectorMatcher.classifyAgainstVariants(
                    peaks: farPeaks, removedCount: 8, matrixOrientation: s.matrixOrientation, phases: phases, settings: open,
                    originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                    matrixReflections: reflections, matrixInPlaneCorrectionRad: correctionDeg * .pi / 180)
            }
            #expect(decide(nil).verdict == .indexed)
            let onTheMatrix = decide(fakeMatrix(turnedDeg: 0))
            #expect(onTheMatrix.verdict == .matrix && onTheMatrix.orientation == nil && onTheMatrix.phaseIndex == 0)
            #expect(decide(fakeMatrix(turnedDeg: -3), correctionDeg: 3).verdict == .matrix)
            #expect(decide(fakeMatrix(turnedDeg: -3)).verdict == .indexed, "without the correction the guard looks in the wrong place")
        }
        // Few vectors and no matrix orientation.
        let empty = PhaseVectorMatcher.classifyAgainstVariants(peaks: [], removedCount: 8, matrixOrientation: s.matrixOrientation,
                                                               phases: phases, settings: settings,
                                                               originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        #expect(empty.verdict == .matrix)
        let noMatrix = PhaseVectorMatcher.classifyAgainstVariants(peaks: Self.peaks(s, variant: rich[0], count: 3), removedCount: 0,
                                                                  matrixOrientation: nil, phases: phases, settings: settings,
                                                                  originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q)
        #expect(noMatrix.verdict == .notIndexed)
    }
}
