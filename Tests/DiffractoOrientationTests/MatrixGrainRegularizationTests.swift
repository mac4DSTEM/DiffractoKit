//
//  MatrixGrainRegularizationTests.swift
//  The matrix orientation held to the grain (DiffractoKit#81): a position pulled to another zone takes its neighbourhood's zone
//  when its own pattern fits that zone, and only then. Planted truth: two aluminium grains meeting at a boundary, and a grain
//  of four positions inside another. Neither may be smeared.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct MatrixGrainRegularizationTests {
    typealias Base = MatrixOrientationFieldTests
    typealias Planted = MatrixOrientationFieldTests.Planted

    static func setup() throws -> (reflections: [ACOMReflection], matcher: ACOMMatcher) {
        let reflections = try Base.reflections()
        let library = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Base.reach)
        let matcher = try #require(ACOMMatcher(library: library))
        return (reflections, matcher)
    }

    static func family(_ field: MatrixOrientationField, _ i: Int) -> [Double]? {
        field.orientations[i].map { Base.family($0.columns.2) }
    }

    static func isZone(_ f: [Double]?, _ zone: SIMD3<Double>) -> Bool {
        guard let f else { return false }
        let t = Base.family(zone)
        return acos(min(1, f[0] * t[0] + f[1] * t[1] + f[2] * t[2])) * 180 / .pi < 8
    }

    static func field(_ scan: [[SIMD2<Double>]], width: Int, height: Int, _ matcher: ACOMMatcher) throws -> MatrixOrientationField {
        try #require(PhaseVectorMatcher.matrixOrientations(
            scan: scan, scanWidth: width, scanHeight: height, matcher: matcher,
            originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q))
    }

    static func regularized(_ field: MatrixOrientationField, _ scan: [[SIMD2<Double>]], _ matcher: ACOMMatcher) throws -> MatrixOrientationRegularization {
        try #require(PhaseVectorMatcher.regularizingMatrixOrientations(
            field, scan: scan, matcher: matcher, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q))
    }

    /// Two grains, [001] in columns 0...3 and [011] in columns 4...7, six rows. The boundary is not moved: no position changes.
    /// Then one position inside the [001] grain is given another zone, as a pattern full of precipitate spots can do: it, and
    /// only it, takes the grain's zone back, because its own pattern fits that zone.
    /// Mutations: the vote alone decides (the boundary columns 3 and 4 have a majority of their own grain, so this needs the
    /// small-grain test below); the score check reversed (the pulled position is not repaired).
    @Test func aBoundaryStaysAndAPulledPositionIsRepaired() throws {
        let (reflections, matcher) = try Self.setup()
        let width = 8, height = 6
        let a = Base.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections).peaks
        let b = Base.pattern(Planted(zone: SIMD3(0, 1, 1), inPlaneDeg: 40), reflections).peaks
        let scan = (0..<width * height).map { $0 % width < 4 ? a : b }
        let truth = try Self.field(scan, width: width, height: height, matcher)
        for i in 0..<scan.count {
            #expect(Self.isZone(Self.family(truth, i), i % width < 4 ? SIMD3(0, 0, 1) : SIMD3(0, 1, 1)), "position \(i)")
        }
        let untouched = try Self.regularized(truth, scan, matcher)
        #expect(untouched.changedPositions.isEmpty, "the boundary moved: \(untouched.changedPositions)")
        #expect(untouched.field == truth)

        // Position (row 2, column 1) pulled to [111]: the best [111]-like template for its [001] pattern, with that template's score.
        let pulledIndex = 2 * width + 1
        let zone111 = try #require(matcher.library.zoneAxes.max { simd_dot($0, simd_normalize(SIMD3(1, 1, 1))) < simd_dot($1, simd_normalize(SIMD3(1, 1, 1))) })
        let peaks = a.dropFirst().map { ACOMPeak(x: Float($0.x), y: Float($0.y), intensity: 1) }
        let wrong = matcher.match(peaks: peaks, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
                                  near: zone111, withinRadians: 5 * .pi / 180)
        #expect(wrong.templateIndex >= 0)
        var results = truth.results, orientations = truth.orientations
        results[pulledIndex] = wrong
        orientations[pulledIndex] = matcher.library.orientationMatrix(templateIndex: wrong.templateIndex, inPlaneAngle: wrong.inPlaneAngle, mirrored: wrong.mirrored)
        let pulled = MatrixOrientationField(width: width, height: height, results: results, orientations: orientations)
        #expect(!Self.isZone(Self.family(pulled, pulledIndex), SIMD3(0, 0, 1)))
        let repaired = try Self.regularized(pulled, scan, matcher)
        #expect(repaired.changedPositions == [pulledIndex])
        #expect(Self.isZone(Self.family(repaired.field, pulledIndex), SIMD3(0, 0, 1)))
        // And its spots go again at the repaired orientation.
        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: scan, field: repaired.field, reflections: reflections,
            originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q))
        #expect(residual.removedCounts[pulledIndex] == 8)
    }

    /// A grain of four positions ([111], 2 x 2) inside a [001] grain of 7 x 7: every one of its positions is outvoted in its
    /// 5 x 5 window (at most 3 of 24 neighbours share its zone), and every one keeps its zone, because its pattern does not
    /// fit the neighbours' zone. Mutation: the vote alone decides (the four positions are overwritten).
    @Test func aGrainOfAFewPositionsIsNotSmeared() throws {
        let (reflections, matcher) = try Self.setup()
        let side = 7
        let a = Base.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections).peaks
        let c = Base.pattern(Planted(zone: SIMD3(1, 1, 1), inPlaneDeg: 70), reflections).peaks
        func small(_ i: Int) -> Bool { (3...4).contains(i / side) && (3...4).contains(i % side) }
        let scan = (0..<side * side).map { small($0) ? c : a }
        let truth = try Self.field(scan, width: side, height: side, matcher)
        for i in 0..<scan.count {
            #expect(Self.isZone(Self.family(truth, i), small(i) ? SIMD3(1, 1, 1) : SIMD3(0, 0, 1)), "position \(i)")
        }
        let out = try Self.regularized(truth, scan, matcher)
        #expect(out.changedPositions.isEmpty, "smeared: \(out.changedPositions)")
        #expect(out.field == truth)
    }

    /// Malformed input is refused, and a position without an orientation is left alone.
    @Test func malformedInputIsRefused() throws {
        let (reflections, matcher) = try Self.setup()
        let a = Base.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections).peaks
        let centre = SIMD2<Double>(Double(Base.origin), Double(Base.origin))
        let scan = [a, a, [centre], a, a, a]
        let f = try Self.field(scan, width: 3, height: 2, matcher)
        let out = try Self.regularized(f, scan, matcher)
        #expect(out.field.orientations[2] == nil && out.changedPositions.isEmpty)
        #expect(PhaseVectorMatcher.regularizingMatrixOrientations(
            f, scan: Array(scan.dropLast()), matcher: matcher, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q) == nil)
        #expect(PhaseVectorMatcher.regularizingMatrixOrientations(
            f, scan: scan, matcher: matcher, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
            neighbourhoodRadius: 0) == nil)
    }
}
