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

    /// A precipitate-rich position: its own (wrong) zone scores clearly BETTER in the template match than the grain's zone, so the
    /// score test keeps it. Judged by the spots a zone explains, the grain's zone takes all eight matrix spots and the wrong zone
    /// fewer, so with the matrix's reflections given the position is repaired. The boundary and the small grain still do not
    /// move with the reflections given. Mutation: the explained-spot test left out (the pulled position keeps its zone).
    @Test func aPositionWhoseWrongZoneScoresBetterIsRepairedByTheSpotsTheGrainsZoneExplains() throws {
        let (reflections, matcher) = try Self.setup()
        let side = 5
        let a = Base.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections).peaks
        let scan = [[SIMD2<Double>]](repeating: a, count: side * side)
        let truth = try Self.field(scan, width: side, height: side, matcher)
        let centre = 2 * side + 2
        let zone111 = try #require(matcher.library.zoneAxes.max { simd_dot($0, simd_normalize(SIMD3(1, 1, 1))) < simd_dot($1, simd_normalize(SIMD3(1, 1, 1))) })
        let peaks = a.dropFirst().map { ACOMPeak(x: Float($0.x), y: Float($0.y), intensity: 1) }
        var wrong = matcher.match(peaks: peaks, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
                                  near: zone111, withinRadians: 5 * .pi / 180)
        wrong.score = 2          // clearly better than anything the grain's zone can score here
        var results = truth.results, orientations = truth.orientations
        results[centre] = wrong
        orientations[centre] = matcher.library.orientationMatrix(templateIndex: wrong.templateIndex, inPlaneAngle: wrong.inPlaneAngle, mirrored: wrong.mirrored)
        let pulled = MatrixOrientationField(width: side, height: side, results: results, orientations: orientations)

        let byScore = try Self.regularized(pulled, scan, matcher)
        #expect(byScore.changedPositions.isEmpty, "the score test alone must keep a zone that scores clearly better")
        let bySpots = try #require(PhaseVectorMatcher.regularizingMatrixOrientations(
            pulled, scan: scan, matcher: matcher, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
            explainedBy: reflections))
        #expect(bySpots.changedPositions == [centre])
        #expect(Self.isZone(Self.family(bySpots.field, centre), SIMD3(0, 0, 1)))

        // Two grains at a boundary and a 2 x 2 grain, with the reflections given: still nothing moves.
        let b = Base.pattern(Planted(zone: SIMD3(0, 1, 1), inPlaneDeg: 40), reflections).peaks
        let c = Base.pattern(Planted(zone: SIMD3(1, 1, 1), inPlaneDeg: 70), reflections).peaks
        let boundary = (0..<48).map { $0 % 8 < 4 ? a : b }
        let boundaryField = try Self.field(boundary, width: 8, height: 6, matcher)
        let small = (0..<49).map { (3...4).contains($0 / 7) && (3...4).contains($0 % 7) ? c : a }
        let smallField = try Self.field(small, width: 7, height: 7, matcher)
        for (f, sc) in [(boundaryField, boundary), (smallField, small)] {
            let out = try #require(PhaseVectorMatcher.regularizingMatrixOrientations(
                f, scan: sc, matcher: matcher, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
                explainedBy: reflections))
            #expect(out.changedPositions.isEmpty, "smeared with the reflections given: \(out.changedPositions)")
        }
    }

    /// Two matrix grains overlapping along the beam (an inclined boundary): the position holds both grains' spots. Its own
    /// orientation explains one grain; with the neighbouring grains taken into the removal the other grain's spots go too, and
    /// the position ends as matrix. Spots that are no grain's stay. Without the neighbours the second grain's spots are left,
    /// which is what a candidate phase was then matched to.
    /// Mutation: a neighbour's orientation removing whatever lies near any of its spots without the share test (the two spots
    /// planted on the neighbouring grain's positions at a single-grain position would go).
    @Test func aPositionWhereTwoGrainsOverlapEndsAsMatrix() throws {
        let (reflections, matcher) = try Self.setup()
        let a = Base.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections)
        let b = Base.pattern(Planted(zone: SIMD3(0, 1, 1), inPlaneDeg: 40), reflections)
        let both = a.peaks + b.peaks.dropFirst()
        let withExtras = both + Base.extras.map { SIMD2(Double(Base.origin) + $0.x / Base.q, Double(Base.origin) + $0.y / Base.q) }
        // Seven columns, three rows: grain A, the overlap column, grain B; the middle row's overlap position carries extras.
        let width = 7, height = 3
        let scan: [[SIMD2<Double>]] = (0..<width * height).map { i in
            let column = i % width
            if column < 3 { return a.peaks }
            if column > 3 { return b.peaks }
            return i / width == 1 ? withExtras : both
        }
        let field = try Self.regularized(try Self.field(scan, width: width, height: height, matcher), scan, matcher).field
        let symmetry = CrystalOrientationSymmetry(laueClass: matcher.library.laueClass)
        let alone = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: scan, field: field, reflections: reflections, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q))
        let together = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: scan, field: field, reflections: reflections, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
            neighbouringGrains: MatrixNeighbourhood(radius: 2, symmetry: symmetry)))
        func left(_ r: MatrixResidualScan, _ i: Int) -> Int {
            PhaseVectorMatcher.experimentalVectors(peaks: r.scan[i], originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
                                                   directBeamRadiusInvAngstrom: 0.15).count
        }
        for row in 0..<height {
            let i = row * width + 3
            #expect(left(alone, i) >= 6, "row \(row): without the neighbours the other grain's spots are left (\(left(alone, i)))")
            #expect(left(together, i) == (row == 1 ? Base.extras.count : 0), "row \(row): \(left(together, i)) vectors left")
            #expect(together.removedCounts[i] == a.matrixSpots + b.matrixSpots)
        }
        // A grain that is merely next door takes nothing: a grain-A position two columns from the overlap, with two spots
        // planted exactly on two of grain B's spot positions. Grain B is in its window, but only two of B's eight spots are
        // measured there (a quarter, under the half that is asked for), so they stay.
        let bVectors = PhaseVectorMatcher.experimentalVectors(peaks: b.peaks, originX: Base.origin, originY: Base.origin,
                                                              invAngstromPerPixel: Base.q, directBeamRadiusInvAngstrom: 0.15)
        var nextDoor = scan
        let position = 1 * width + 2
        nextDoor[position] = a.peaks + bVectors.prefix(2).map { SIMD2(Double(Base.origin) + $0.x / Base.q, Double(Base.origin) + $0.y / Base.q) }
        let nextField = try Self.regularized(try Self.field(nextDoor, width: width, height: height, matcher), nextDoor, matcher).field
        let kept = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: nextDoor, field: nextField, reflections: reflections, originX: Base.origin, originY: Base.origin, invAngstromPerPixel: Base.q,
            neighbouringGrains: MatrixNeighbourhood(radius: 2, symmetry: symmetry)))
        #expect(left(kept, position) == 2, "two spots on a neighbouring grain's positions were taken at a single-grain position")
        #expect(kept.removedCounts[position] == a.matrixSpots)
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
