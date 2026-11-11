//
//  MatrixOrientationFieldTests.swift
//  Step 1 of "find phase orientations from the CIF alone" (DiffractoKit#81): the matrix orientation per scan position from the
//  templates, and the removal of the spots that orientation explains. Planted truth only (no Python reference does this step as
//  one call): an aluminium scan of three grains, [001], [011] and [111], as the synthetic Al-Mg-Si demo has them, with extra
//  spots that are not the matrix's at some positions.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct MatrixOrientationFieldTests {
    static let q = 0.012                      // Å⁻¹ per pixel, the demo's
    static let origin: Float = 63.5
    static let reach = 0.75                   // Å⁻¹ to the detector's edge

    struct Planted {
        let zone: SIMD3<Double>
        let inPlaneDeg: Double
        var extra: [SIMD2<Double>] = []       // not the matrix's, Å⁻¹
    }

    static func reflections() throws -> [ACOMReflection] {
        try CrystalPresets.aluminum.kinematicReflections(kMax: reach, retainingAll: false).kept
            .map { ACOMReflection(g: $0.reflection.g, intensity: $0.structureFactor.intensity) }
    }

    /// The zero-order pattern of a zone at an in-plane angle, as pixel positions with the direct beam: the reflections with
    /// g·n = 0 in the detector basis of the zone, turned by the angle. Written out here, not taken from the library under test.
    static func pattern(_ p: Planted, _ reflections: [ACOMReflection]) -> (peaks: [SIMD2<Double>], matrixSpots: Int) {
        let n = simd_normalize(p.zone)
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: n)
        let c = cos(p.inPlaneDeg * .pi / 180), s = sin(p.inPlaneDeg * .pi / 180)
        var peaks = [SIMD2<Double>(Double(origin), Double(origin))]
        var count = 0
        for r in reflections where abs(simd_dot(r.g, n)) < 1e-6 {
            let x = simd_dot(r.g, basis.columns.0), y = simd_dot(r.g, basis.columns.1)
            peaks.append(SIMD2(Double(origin) + (c * x - s * y) / q, Double(origin) + (s * x + c * y) / q))
            count += 1
        }
        for e in p.extra { peaks.append(SIMD2(Double(origin) + e.x / q, Double(origin) + e.y / q)) }
        return (peaks, count)
    }

    static func family(_ v: SIMD3<Double>) -> [Double] {
        let n = simd_normalize(v)
        return [abs(n.x), abs(n.y), abs(n.z)].sorted()
    }

    static let extras: [SIMD2<Double>] = [SIMD2(0.31, 0.05), SIMD2(-0.31, -0.05), SIMD2(0.12, 0.27), SIMD2(-0.12, -0.27)]
    static let scan: [Planted] = [
        Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12, extra: extras),
        Planted(zone: SIMD3(0, 1, 1), inPlaneDeg: 40), Planted(zone: SIMD3(0, 1, 1), inPlaneDeg: 40, extra: extras),
        Planted(zone: SIMD3(1, 1, 1), inPlaneDeg: 70), Planted(zone: SIMD3(1, 1, 1), inPlaneDeg: 70, extra: extras),
    ]

    /// The matrix is found where it is, grain by grain, and its spots go: nothing is left of a pure matrix pattern in any of
    /// the three grains, and exactly the spots that are not the matrix's are left where some were planted.
    /// Mutations: one orientation for the whole scan (the [011] and [111] positions keep their spots); the in-plane angle
    /// taken with the wrong sign or left at the template's bin without the fit (spots survive at 40 and 70 degrees); the
    /// removal at any reflection of the crystal regardless of orientation (the planted extras stay, but so would a spot on
    /// another zone's position: the last expectation).
    @Test func theMatrixIsFoundPerPositionAndItsSpotsAreRemoved() throws {
        let reflections = try Self.reflections()
        let library = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)
        let matcher = try #require(ACOMMatcher(library: library))
        let planted = Self.scan.map { Self.pattern($0, reflections) }
        let peaks = planted.map(\.peaks)

        let field = try #require(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: peaks.count, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(field.width == peaks.count && field.height == 1)
        for (i, p) in Self.scan.enumerated() {
            let m = try #require(field.orientations[i], "position \(i) has no orientation")
            // The beam direction in crystal coordinates is the planted zone, up to the cubic symmetry, within the template spacing.
            let found = Self.family(m.columns.2), truth = Self.family(p.zone)
            let angle = acos(min(1, found[0] * truth[0] + found[1] * truth[1] + found[2] * truth[2])) * 180 / .pi
            // A pure matrix pattern lands on the nearest template (within 4 degrees at 600 directions). Spots that are not the
            // matrix's can move the winner to a neighbouring template (5.7 degrees measured for the [011] position); the
            // removal has to hold that, which the counts below check.
            #expect(angle < (p.extra.isEmpty ? 4 : 8), "position \(i): beam \(m.columns.2) is \(angle) deg from \(p.zone)")
        }

        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: peaks, field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        for (i, p) in Self.scan.enumerated() {
            let left = PhaseVectorMatcher.experimentalVectors(
                peaks: residual.scan[i], originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
                directBeamRadiusInvAngstrom: 0.15)
            #expect(residual.removedCounts[i] == planted[i].matrixSpots, "position \(i): \(residual.removedCounts[i]) of \(planted[i].matrixSpots) matrix spots removed")
            #expect(left.count == p.extra.count, "position \(i): \(left.count) vectors left, \(p.extra.count) planted beside the matrix")
            for e in p.extra {
                #expect(left.contains { simd_length($0 - e) < 1e-3 }, "position \(i): the planted \(e) was removed")
            }
        }
        // The direct beam is never a matrix spot: it stays in the residual scan for the caller's own direct-beam rule.
        #expect(residual.scan.allSatisfy { $0.contains { abs($0.x - Double(Self.origin)) < 1e-9 && abs($0.y - Double(Self.origin)) < 1e-9 } })
    }

    /// The case this step is for (the synthetic Al-Mg-Si demo's second grain): with the matrix listed at [001] only, an Al [011]
    /// position keeps all its spots and the phase map does not call it matrix; a candidate is asked to explain a matrix pattern.
    /// With the matrix found and removed position by position first, the same map calls all three grains matrix and counts
    /// the spots the matrix took. Mutation: `completing` left out (the pure positions read "no data", not matrix).
    @Test func aSecondMatrixGrainIsMatrixOnceItsOwnOrientationExplainsIt() throws {
        let reflections = try Self.reflections()
        let pure = [Self.scan[0], Self.scan[2], Self.scan[4]]
        let planted = pure.map { Self.pattern($0, reflections) }
        let peaks = planted.map(\.peaks)
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = Self.reach
        reference.minimumIntensityFraction = 0
        var matching = PhaseVectorSettings()
        matching.classificationRule = .knownVariants
        let phases = [
            PhaseDefinition(id: "al", displayName: "Al", crystal: CrystalPresets.aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "b001", displayName: "beta'' [001]", crystal: CrystalPresets.betaDoublePrime, role: .candidate, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "b010", displayName: "beta'' [010]", crystal: CrystalPresets.betaDoublePrime, role: .candidate, zoneAxes: [SIMD3(0, 1, 0)]),
        ]
        let library = try PhaseReferenceLibrary.build(phases: phases, settings: reference)
        let before = try #require(PhaseVectorMatcher.map(
            scan: peaks, scanWidth: 3, scanHeight: 1, library: library, settings: matching,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(before.results[0].verdict == .matrix)
        #expect(before.results[1].verdict != .matrix, "the [011] grain is not matrix to a list that holds Al at [001] only")
        #expect(before.results[1].removedCount == 0)

        let matcher = try #require(ACOMMatcher(library: try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)))
        let field = try #require(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: 3, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: peaks, field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        let onResidual = try #require(PhaseVectorMatcher.map(
            scan: residual.scan, scanWidth: 3, scanHeight: 1, library: library, settings: matching,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, matrixEntryIndex: before.matrixEntryIndex))
        let after = try #require(residual.completing(onResidual))
        for i in 0..<3 {
            #expect(after.results[i].verdict == .matrix, "position \(i): \(after.results[i].verdict)")
            #expect(after.results[i].phaseIndex == Int32(library.matrixPhaseIndex))
            #expect(Int(after.results[i].removedCount) == planted[i].matrixSpots)
        }
        #expect(residual.completing(PhaseMap(width: 2, height: 1, matrixEntryIndex: -1, phaseNames: [], matrixPhaseIndex: 0)) == nil)
    }

    /// Removal follows the orientation found, not the crystal: a spot that sits on a reflection of ANOTHER zone of the matrix
    /// is not the matrix's here and stays. On the [001] grain the {111} reflections are not in the zone; a spot planted at the
    /// {111} length survives.
    @Test func aSpotOnAReflectionOfAnotherZoneStays() throws {
        let reflections = try Self.reflections()
        let library = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)
        let matcher = try #require(ACOMMatcher(library: library))
        let g111 = 3.0.squareRoot() / 4.0495
        let stray = SIMD2(g111 * cos(0.3), g111 * sin(0.3)), opposite = SIMD2(-stray.x, -stray.y)
        let peaks = [Self.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12, extra: [stray, opposite]), reflections).peaks]
        let field = try #require(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: 1, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: peaks, field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        let left = PhaseVectorMatcher.experimentalVectors(
            peaks: residual.scan[0], originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q,
            directBeamRadiusInvAngstrom: 0.15)
        #expect(left.count == 2)
        #expect(residual.removedCounts[0] == 8, "the eight [001] spots within reach go, the two strays do not")
    }

    /// The in-plane fit, on a field built by hand so the angle error is known. The matrix at the planted angle plus 3 degrees:
    /// the fit takes the 3 degrees back and every spot goes (0.494 Å⁻¹ x 3 degrees is 0.026 Å⁻¹, more than the tolerance, so
    /// without the fit the spots stay). At plus 8 degrees no spot is within pairing reach: nothing is fitted and nothing removed.
    /// Mutations: the fit left out (red at 3 degrees); the fit's sign reversed (red at 3 degrees).
    @Test func theInPlaneAngleIsFittedBeforeTheRemoval() throws {
        let reflections = try Self.reflections()
        let planted = Self.pattern(Planted(zone: SIMD3(0, 0, 1), inPlaneDeg: 12), reflections)
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1))
        func field(offsetDeg: Double) -> MatrixOrientationField {
            let m = ZoneAxisFrame.matrix(detectorBasis: basis, inPlaneAngle: Float((12 + offsetDeg) * .pi / 180))
            return MatrixOrientationField(width: 1, height: 1, results: [.empty], orientations: [m])
        }
        let near = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: [planted.peaks], field: field(offsetDeg: 3), reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(abs(near.inPlaneCorrectionsRad[0] * 180 / .pi + 3) < 0.01, "correction \(near.inPlaneCorrectionsRad[0] * 180 / .pi) deg")
        #expect(near.removedCounts[0] == planted.matrixSpots)
        let far = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: [planted.peaks], field: field(offsetDeg: 8), reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(far.inPlaneCorrectionsRad[0] == 0 && far.removedCounts[0] == 0)
        #expect(far.scan[0] == planted.peaks)
    }

    /// What the matcher would not use is neither matched nor removed: a non-finite peak and a peak past the maximum vector stay
    /// and do not break the match, and the direct beam with ONE spot is not enough for an orientation.
    @Test func unusablePeaksAreNeitherMatchedNorRemoved() throws {
        let reflections = try Self.reflections()
        let matcher = try #require(ACOMMatcher(library: try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)))
        let pure = Self.pattern(Self.scan[0], reflections)
        let centre = SIMD2<Double>(Double(Self.origin), Double(Self.origin))
        let g200 = 2 / 4.0495
        let peaks: [[SIMD2<Double>]] = [
            pure.peaks + [SIMD2(.nan, 10), SIMD2(1e300, 5)],
            [centre, SIMD2(Double(Self.origin) + g200 / Self.q, Double(Self.origin))],
        ]
        let field = try #require(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: 2, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, maximumVectorInvAngstrom: 0.6))
        #expect(field.orientations[0] != nil)
        #expect(field.orientations[1] == nil, "the direct beam and one spot are one usable peak")
        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: peaks, field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q, maximumVectorInvAngstrom: 0.6))
        // Within 0.6 Å⁻¹ the [001] pattern holds the four {200} spots; the four {220} at 0.698 are past the maximum and stay.
        #expect(residual.removedCounts[0] == 4)
        #expect(residual.scan[0].count == pure.peaks.count + 2 - 4)
        #expect(residual.scan[0].contains { $0.x.isNaN } && residual.scan[0].contains { $0.x == 1e300 })
        #expect(residual.removedCounts[1] == 0 && residual.scan[1] == peaks[1])
    }

    /// `completing` changes a verdict only where the matrix took spots: a position with no peaks stays "no data", and an indexed
    /// position keeps its verdict and gains the count. Mutation: every "no data" turned into matrix (the first expectation).
    @Test func completingTouchesOnlyWhatTheMatrixExplained() throws {
        var map = PhaseMap(width: 3, height: 1, matrixEntryIndex: 0, phaseNames: ["Al", "x"], matrixPhaseIndex: 0)
        map.results[1].verdict = .indexed
        map.results[1].phaseIndex = 1
        map.results[1].removedCount = 2
        let residual = MatrixResidualScan(scan: [[], [], []], removedCounts: [0, 5, 7], inPlaneCorrectionsRad: [0, 0, 0])
        let out = try #require(residual.completing(map))
        #expect(out.results[0].verdict == .noData && out.results[0].removedCount == 0)
        #expect(out.results[1].verdict == .indexed && out.results[1].phaseIndex == 1 && out.results[1].removedCount == 7)
        #expect(out.results[2].verdict == .matrix && out.results[2].phaseIndex == 0 && out.results[2].removedCount == 7)
    }

    /// A position without two peaks has no orientation and loses nothing; a malformed scan is refused.
    @Test func positionsWithoutPeaksAndMalformedScans() throws {
        let reflections = try Self.reflections()
        let library = try ACOMTemplateLibrary.generate(crystal: CrystalPresets.aluminum, kMax: Self.reach)
        let matcher = try #require(ACOMMatcher(library: library))
        let centre = SIMD2<Double>(Double(Self.origin), Double(Self.origin))
        let peaks = [[centre], Self.pattern(Self.scan[0], reflections).peaks]
        let field = try #require(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: 2, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(field.orientations[0] == nil && field.orientations[1] != nil)
        let residual = try #require(PhaseVectorMatcher.removingMatrixSpots(
            scan: peaks, field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q))
        #expect(residual.scan[0] == [centre] && residual.removedCounts[0] == 0)
        #expect(PhaseVectorMatcher.matrixOrientations(
            scan: peaks, scanWidth: 3, scanHeight: 1, matcher: matcher,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q) == nil)
        #expect(PhaseVectorMatcher.removingMatrixSpots(
            scan: [peaks[0]], field: field, reflections: reflections,
            originX: Self.origin, originY: Self.origin, invAngstromPerPixel: Self.q) == nil)
    }
}
