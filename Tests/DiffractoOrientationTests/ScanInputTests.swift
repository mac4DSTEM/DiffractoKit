//
//  ScanInputTests.swift
//  Role: the scan as the kit takes it (registration slice 4, section 2.2 to 2.4): the two input forms, the refusal of a
//        malformed scan, the worker count that moves no number, progress and cancellation, and the documented defaults of
//        the two settings structs.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

/// Records the progress values a run reports (from several threads).
private final class ProgressLog: @unchecked Sendable {
    private let lock = NSLock()
    private var values: [Double] = []
    func add(_ v: Double) { lock.lock(); values.append(v); lock.unlock() }
    var all: [Double] { lock.lock(); defer { lock.unlock() }; return values }
}

private func sameResult(_ a: PhaseVectorResult, _ b: PhaseVectorResult) -> Bool {
    a.verdict == b.verdict && a.phaseIndex == b.phaseIndex && a.entryIndex == b.entryIndex
        && a.score.bitPattern == b.score.bitPattern && a.matchedCount == b.matchedCount
        && a.survivingCount == b.survivingCount && a.removedCount == b.removedCount
        && a.runnerUpPhaseIndex == b.runnerUpPhaseIndex && a.runnerUpScore.bitPattern == b.runnerUpScore.bitPattern
}

@Suite struct ScanInputTests {

    private func library() throws -> PhaseReferenceLibrary {
        var s = PhaseReferenceSettings()
        s.kMaxInvAngstrom = 1.0
        s.inPlaneStepDeg = 30
        return try PhaseReferenceLibrary.build(phases: [
            PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime, role: .candidate,
                            zoneAxes: [SIMD3(0, 1, 0)]),
        ], settings: s)
    }

    /// A 4 × 4 scan of calibrated vectors: matrix positions, candidate positions, noise, and empty positions.
    private func plantedScan(_ library: PhaseReferenceLibrary) -> [[SIMD2<Double>]] {
        let matrix = library.entries[library.matrixEntryIndices[0]].vectors.map(\.q)
        let candidate = library.entries[library.candidateEntryIndices[5]].vectors.map(\.q)
        var state: UInt64 = 0x9E37_79B9_7F4A_7C15
        func uniform() -> Double {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return Double(state >> 11) / Double(UInt64(1) << 53)
        }
        return (0..<16).map { i in
            switch i % 4 {
            case 0: return matrix
            case 1: return matrix + candidate
            case 2: return (0..<9).map { _ in SIMD2(0.3 * uniform() - 0.15 + 0.4, 0.5 * uniform() - 0.25) }
            default: return []
            }
        }
    }

    /// Mutation 16 of the registration: the experimental vector computed in Double before the subtraction. The pixel form
    /// subtracts the origin in Float, converts, then scales. Positions are chosen so that the Float difference is inexact
    /// (a large coordinate and a small origin), where the two orders differ.
    @Test func thePixelFormSubtractsInFloatAndTheCalibratedFormIsExact() {
        let x: Float = 3000.3, y: Float = -2000.7, ox: Float = 0.07, oy: Float = 0.03
        let scale = 0.0457
        let peak = SIMD2<Double>(Double(x), Double(y))
        let got = PhaseVectorMatcher.experimentalVector(
            peak: peak, originX: ox, originY: oy, invAngstromPerPixel: scale, directBeamRadiusInvAngstrom: 0)
        let expected = SIMD2(Double(x - ox) * scale, Double(y - oy) * scale)
        #expect(got == expected)
        let inDouble = SIMD2((Double(x) - Double(ox)) * scale, (Double(y) - Double(oy)) * scale)
        #expect(got != inDouble, "the fixture does not separate the Float order from the Double order")
        // The same through the enum.
        let viaEnum = PhaseVectorMatcher.experimentalVector(
            peak: peak, coordinates: .detectorPixels(originX: ox, originY: oy, invAngstromPerPixel: scale),
            directBeamRadiusInvAngstrom: 0)
        #expect(viaEnum == expected)

        // Calibrated vectors are used as given: no Float rounding, no scale.
        let q = SIMD2<Double>(0.123456789012345, -0.2)
        #expect(PhaseVectorMatcher.experimentalVector(
            peak: q, coordinates: .calibrated, directBeamRadiusInvAngstrom: 0.05) == q)
        let throughPixels = PhaseVectorMatcher.experimentalVector(
            peak: q, originX: 0, originY: 0, invAngstromPerPixel: 1, directBeamRadiusInvAngstrom: 0.05)
        #expect(throughPixels != q, "a double-precision vector survived a Float round trip")
        #expect(abs(throughPixels!.x - q.x) < 1e-7)
    }

    /// The direct beam is dropped at the radius (strictly above kept), the outer reach is exclusive, and a vector that is
    /// not finite is dropped, for both forms.
    @Test func bothFormsApplyTheSameUsabilityRule() {
        for coordinates in [PeakCoordinates.calibrated,
                            .detectorPixels(originX: 0, originY: 0, invAngstromPerPixel: 1)] {
            func vector(_ x: Double, reach: Double = 0) -> SIMD2<Double>? {
                PhaseVectorMatcher.experimentalVector(
                    peak: SIMD2(x, 0), coordinates: coordinates, directBeamRadiusInvAngstrom: 0.5,
                    maximumVectorInvAngstrom: reach)
            }
            #expect(vector(0.5) == nil, "a vector AT the direct-beam radius was kept")
            #expect(vector(0.75) != nil)
            #expect(vector(2, reach: 2) == nil, "a vector AT the outer reach was kept")
            #expect(vector(1.75, reach: 2) != nil)
            #expect(vector(.nan) == nil)
            #expect(vector(.infinity) == nil)
        }
    }

    /// The worker count moves no number (registration gate G1 runs 1 and 8): each position's result depends only on its
    /// own vectors. Compared bit for bit, NaN scores included.
    @Test func theWorkerCountMovesNoNumber() throws {
        let library = try library()
        let scan = plantedScan(library)
        let settings = PhaseVectorSettings()
        let reference = try #require(PhaseVectorMatcher.map(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: settings,
            coordinates: .calibrated, maximumWorkers: 1))
        #expect(reference.count(of: .matrix) + reference.count(of: .indexed) + reference.count(of: .notIndexed)
                + reference.count(of: .noData) == 16)
        #expect(reference.count(of: .noData) == 4, "the empty positions are 'no data'")
        #expect(reference.count(of: .indexed) > 0, "the planted candidate positions were not indexed")
        for workers in [2, 3, 8, nil] as [Int?] {
            let other = try #require(PhaseVectorMatcher.map(
                scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: settings,
                coordinates: .calibrated, maximumWorkers: workers))
            #expect(other.matrixEntryIndex == reference.matrixEntryIndex)
            #expect(other.results.count == reference.results.count)
            for (a, b) in zip(other.results, reference.results) {
                #expect(sameResult(a, b), "workers \(String(describing: workers)) moved a result")
            }
        }
    }

    /// Progress ends at 1 and never exceeds it; a cancel that arrives mid-run returns nil, never a partial map.
    @Test func progressEndsAtOneAndACancelledRunIsNil() throws {
        let library = try library()
        let scan = plantedScan(library)
        let log = ProgressLog()
        let map = PhaseVectorMatcher.map(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: PhaseVectorSettings(),
            coordinates: .calibrated, maximumWorkers: 3, progress: { log.add($0) })
        #expect(map != nil)
        #expect(log.all.max() == 1)
        #expect(log.all.allSatisfy { $0 > 0 && $0 <= 1 })

        final class Counter: @unchecked Sendable {
            let lock = NSLock(); var n = 0
            func next() -> Int { lock.lock(); defer { lock.unlock() }; n += 1; return n }
        }
        let polls = Counter()
        #expect(PhaseVectorMatcher.map(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: PhaseVectorSettings(),
            coordinates: .calibrated, matrixEntryIndex: 0, maximumWorkers: 1,
            isCancelled: { polls.next() > 6 }) == nil, "a run cancelled mid-way returned a map")
    }

    /// A malformed scan, an out-of-range matrix entry and an overflowing size are refused, never a trap.
    @Test func aMalformedScanOrEntryIsRefused() throws {
        let library = try library()
        let scan = plantedScan(library)
        let settings = PhaseVectorSettings()
        func map(_ s: [[SIMD2<Double>]], _ w: Int, _ h: Int, entry: Int? = nil) -> PhaseMap? {
            PhaseVectorMatcher.map(scan: s, scanWidth: w, scanHeight: h, library: library, settings: settings,
                                   coordinates: .calibrated, matrixEntryIndex: entry)
        }
        #expect(map(scan, 4, 4) != nil)
        #expect(map(scan, 4, 5) == nil, "16 positions are not 4 × 5")
        #expect(map(scan, 0, 16) == nil)
        #expect(map(scan, -4, -4) == nil)
        #expect(map(scan, Int.max, 2) == nil, "a width × height that overflows is not a scan")
        #expect(map([], 0, 0) == nil)
        #expect(map(scan, 4, 4, entry: library.entries.count) == nil)
        #expect(map(scan, 4, 4, entry: -1) == nil)

        #expect(PhaseVectorMatcher.fitMatrixOrientation(
            scan: scan, scanWidth: 3, scanHeight: 4, library: library, settings: settings,
            coordinates: .calibrated) == nil)
        #expect(PhaseVectorMatcher.fitMatrixOrientation(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: settings,
            coordinates: .calibrated) != nil)
        #expect(try PhaseVectorMatcher.fitZoneAxis(
            scan: scan, scanWidth: 3, scanHeight: 4, crystal: .aluminum, referenceSettings: PhaseReferenceSettings(),
            settings: settings, coordinates: .calibrated).isEmpty)
    }

    /// The pixel overloads are the enum overloads with the origin and scale in the enum: the two give the same map.
    @Test func thePixelOverloadsAreTheEnumOverloads() throws {
        let library = try library()
        let scale = 0.01
        let scan: [[SIMD2<Double>]] = plantedScan(library).map { $0.map { SIMD2(Double(Float($0.x / scale)) + 64,
                                                                                 Double(Float($0.y / scale)) + 64) } }
        let settings = PhaseVectorSettings()
        let a = try #require(PhaseVectorMatcher.map(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: settings,
            originX: 64, originY: 64, invAngstromPerPixel: scale))
        let b = try #require(PhaseVectorMatcher.map(
            scan: scan, scanWidth: 4, scanHeight: 4, library: library, settings: settings,
            coordinates: .detectorPixels(originX: 64, originY: 64, invAngstromPerPixel: scale)))
        for (x, y) in zip(a.results, b.results) { #expect(sameResult(x, y)) }
        #expect(a.matrixEntryIndex == b.matrixEntryIndex)
    }

    /// Registration section 3: the 24 settings and their documented defaults, and `default` is the empty initialiser.
    @Test func theDocumentedDefaultsAreTheTunedValues() {
        let v = PhaseVectorSettings.default
        #expect(v == PhaseVectorSettings())
        #expect(v.directBeamRadiusInvAngstrom == 0.15)
        #expect(v.maximumVectorInvAngstrom == 0)
        #expect(v.matrixToleranceInvAngstrom == 0.02)
        #expect(v.pairRadiusInvAngstrom == 0.02)
        #expect(v.minimumVectors == 2)
        #expect(v.matrixFallbackExplainedFraction == 0)
        #expect(v.minimumMatchedVectors == 3)
        #expect(v.friedelPairMinimumMatchedVectors == 2)
        #expect(v.chanceMatchMultiple == 5)
        #expect(v.notIndexedAboveInvAngstrom == 0.015)
        #expect(v.minimumPhaseContrastInvAngstrom == 0)
        #expect(v.orientationRelationshipToleranceDeg == 10)
        #expect(v.completenessAwareCrossPhaseRanking == false)
        #expect(v.classificationRule == .search)
        #expect(v.directMatrixMaximumVectors == 1)
        #expect(v.residualCutoffInvAngstrom == 0.07)
        #expect(v.knownVariantsMinimumSpecificReflections == 1)
        #expect(v.minimumMatrixEvidence == 0)

        let r = PhaseReferenceSettings.default
        #expect(r == PhaseReferenceSettings())
        #expect(r.kMaxInvAngstrom == 1.6)
        #expect(r.excitationSlabInvAngstrom == 0.05)
        #expect(r.excitationWidthInvAngstrom == 0.03)
        #expect(r.minimumIntensityFraction == 0.05)
        #expect(r.maximumVectorsPerEntry == 48)
        #expect(r.inPlaneStepDeg == 2)
        #expect(r.maximumEntries == 4000)
    }
}
