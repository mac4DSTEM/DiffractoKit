//
//  PolicyFixesTests.swift
//  DiffractoKit 0.2.0, the Gate B findings on the ring-fit policies (tier 2): F4 (no trap for any bin count), F1 (the three
//  bin counts and what `isValid` requires of them), F2 (the calibration policy's evidence replaces the ellipse policy's;
//  loosening coverage through `fitRing` needs the three bounds), F3 (`RingIndexing.reach` with an invalid argument).
//  Every test names under "RED UNDER" the production change that turns it red.
//

import Testing
import Foundation
import DiffractoGeometry

struct PolicyFixesTests {
    typealias T = RingFitAmendmentTests
    typealias Policy = EllipseCalibration.Policy

    // MARK: F4

    /// No Int traps `Policy(angularBinCount:)`: the derived bounds are computed without overflow, a count outside 6...3600 is
    /// not valid, and for every count that did not overflow before the derived values are the old expressions' (5/6 and 1/3,
    /// truncating towards zero).
    /// RED UNDER: the derivation `angularBinCount * 5 / 6` restored (the first extreme value traps).
    @Test func noBinCountTrapsTheInitialiser() {
        let extremes: [Int] = [Int.max, Int.min, Int.max - 1, Int.min + 1, Int.max / 5 + 1, Int.min / 5 - 1, Int.max / 5,
                               Int.min / 5, Int.max / 2, Int.min / 2, 1 << 62, -(1 << 62), 1 << 40, -(1 << 40), 0, -1, 1]
        for n in extremes {
            let p = Policy(angularBinCount: n)
            #expect(p.angularBinCount == n)
            #expect(!p.isValid, "\(n)")
        }
        for n in -2000...5000 {
            let p = Policy(angularBinCount: n)
            #expect(p.degeneracyBoundBins == n * 5 / 6, "\(n)")
            #expect(p.sparseFloorBins == n / 3, "\(n)")
            #expect(p.isValid == (n >= 6 && n <= 3600 && n / 3 >= 1), "\(n)")
        }
        for n in [Int.max / 5, Int.min / 5, 1 << 40, -(1 << 40), 123_456_789_012, -123_456_789_012] {
            #expect(Policy(angularBinCount: n).degeneracyBoundBins == n * 5 / 6)
        }
        #expect(Policy.default.degeneracyBoundBins == 30 && Policy.default.sparseFloorBins == 12)
    }

    // MARK: F1

    /// The bounds are bin counts: setting `angularBinCount` alone on `.default` keeps 30 and 12 (the initialiser derives
    /// them, a mutation does not); `isValid` requires 1 <= sparse floor <= degeneracy bound <= bin count, so a violated
    /// order is refused and a count mutated alone stays valid with the bounds it had.
    /// RED UNDER: the order check of `isValid` dropped (any of the four cases below becomes valid).
    @Test func theThreeBinCountsAreCheckedAgainstEachOther() {
        var mutated: Policy = .default
        mutated.angularBinCount = 72
        #expect(mutated.degeneracyBoundBins == 30 && mutated.sparseFloorBins == 12)
        #expect(mutated.isValid)
        #expect(mutated != Policy(angularBinCount: 72))
        #expect(Policy(angularBinCount: 72).degeneracyBoundBins == 60 && Policy(angularBinCount: 72).sparseFloorBins == 24)
        var small: Policy = .default
        small.angularBinCount = 24                       // the bound 30 now exceeds the bin count
        #expect(!small.isValid)
        #expect(!Policy(angularBinCount: 36, degeneracyBoundBins: 37).isValid)
        #expect(!Policy(angularBinCount: 36, degeneracyBoundBins: 10, sparseFloorBins: 11).isValid)
        #expect(!Policy(angularBinCount: 36, degeneracyBoundBins: 0, sparseFloorBins: 0).isValid)
        #expect(!Policy(angularBinCount: 36, degeneracyBoundBins: -3, sparseFloorBins: -5).isValid)
        #expect(!Policy(angularBinCount: 36, sparseFloorBins: -1).isValid)
        #expect(Policy(angularBinCount: 36, degeneracyBoundBins: 36, sparseFloorBins: 36).isValid)
        #expect(Policy(angularBinCount: 36, degeneracyBoundBins: 1, sparseFloorBins: 1).isValid)
    }

    // MARK: F2

    static let scene = T.Scene()
    static func gapped() -> DiffractionPattern {
        T.render(Self.scene) { dr, dc in T.inGap(dr, dc, start: 200, width: 120) ? 0 : T.ring(dr, dc, radius: 42, epsilon: 0.01, thetaDeg: 15) }
    }
    static func fitRing(_ pattern: DiffractionPattern, _ policy: RingCalibration.Policy, _ ellipse: Policy) throws -> RingRadiusFit {
        try RingCalibration.fitRing(pattern: pattern, centerQX: 64.47, centerQY: 60.89, innerRadius: 28, outerRadius: 56,
                                    policy: policy, ellipsePolicy: ellipse)
    }

    /// A ring with a 120 degree gap (about 25 of 36 bins): `fitRing` refuses it at the defaults; loosening only the
    /// calibration policy's `minimumOccupiedAngularBins` (and its evidence bound) changes nothing, because the ellipse
    /// fit's own `degeneracyBoundBins` still refuses before any fit; with the ellipse policy's bound loosened too the ring
    /// is fitted and no coverage refusal remains. The calibration policy's evidence replaces the ellipse policy's.
    /// RED UNDER: the ellipse policy's `degeneracyBoundBins` ignored by `fitRing`, or the calibration policy's bound ignored.
    @Test func looseningCoverageThroughFitRingNeedsBothPolicies() throws {
        let pattern = Self.gapped()
        let atDefault = try Self.fitRing(pattern, .default, .default)
        #expect(atDefault.fit == nil && atDefault.refusals.map(\.name) == ["insufficientAngularCoverage"])

        var calibration: RingCalibration.Policy = .default
        calibration.minimumOccupiedAngularBins = 20
        calibration.evidence.minimumOccupiedBins = 20
        let calibrationOnly = try Self.fitRing(pattern, calibration, .default)
        #expect(calibrationOnly.fit == nil && calibrationOnly.refusals.map(\.name) == ["insufficientAngularCoverage"])

        var ellipse: Policy = .default
        ellipse.degeneracyBoundBins = 20
        let ellipseOnly = try Self.fitRing(pattern, .default, ellipse)
        #expect(ellipseOnly.fit != nil)
        #expect(ellipseOnly.refusals.map(\.name).contains("insufficientAngularCoverage"))      // the calibration bound still refuses

        let both = try Self.fitRing(pattern, calibration, ellipse)
        #expect(both.fit != nil)
        #expect(!both.refusals.map(\.name).contains("insufficientAngularCoverage"), "\(both.refusals.map(\.name))")
    }

    /// The ellipse policy's `evidence` field is replaced by the calibration policy's: an ellipse policy whose evidence is
    /// invalid or strict gives the same ring as one with the default evidence. Every other invalid field of the ellipse
    /// policy is still refused (`FitError.invalidInput`).
    /// RED UNDER: the replacement `fitPolicy.evidence = policy.evidence` removed (the strict evidence then reaches the
    /// fit); the validity guard of the ellipse policy removed.
    @Test func theCalibrationEvidenceReplacesTheEllipsePolicysAndOtherFieldsAreChecked() throws {
        let pattern = T.render(Self.scene) { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.02, thetaDeg: 30) }
        let plain = try Self.fitRing(pattern, .default, .default)
        #expect(plain.isUsable)
        for evidence in [RingEvidence.Policy(minimumSectorPixels: 2), RingEvidence.Policy(minimumSectorPixels: 100_000),
                         RingEvidence.Policy(minimumDrop: .nan), RingEvidence.Policy(maximumShift: -1)] {
            var ellipse: Policy = .default
            ellipse.evidence = evidence
            #expect(try Self.fitRing(pattern, .default, ellipse) == plain)
        }
        var strictCalibration: RingCalibration.Policy = .default
        strictCalibration.evidence.minimumSectorPixels = 100_000
        let demanding = try Self.fitRing(pattern, strictCalibration, .default)
        #expect(demanding.report?.ringEvidence?.occupiedBins == 0)
        for bad in [Policy(maximumConicResidual: .nan), Policy(angularBinCount: 5), Policy(maximumAxisRatio: 0.5),
                    Policy(angularBinCount: Int.max)] {
            var refused = false
            do { _ = try Self.fitRing(pattern, .default, bad) }
            catch let caught as EllipseCalibration.FitError { guard case .invalidInput = caught else { throw caught }; refused = true }
            #expect(refused, "\(bad)")
        }
    }

    // MARK: F3 (RingIndexing.reach)

    /// `reach` answers NaN for the coverage or margin that `index` refuses (a coverage that is not a positive number, a
    /// margin that is negative or not finite), and the closed form otherwise.
    /// RED UNDER: the argument guard of `reach` removed (a negative margin gives a number or NaN by accident, an infinite
    /// margin gives infinity).
    @Test func reachAnswersNaNWhereIndexRefuses() {
        #expect(RingIndexing.reach(coverage: 3) == (9 + 2 * log(19.0)).squareRoot())
        #expect(RingIndexing.reach(coverage: 3, minimumMargin: 0) == 3)
        for (k, m) in [(3.0, -1.0), (3.0, Double.nan), (3.0, Double.infinity), (3.0, -Double.infinity), (0.0, 1.0), (-3.0, 1.0),
                       (Double.nan, 1.0), (Double.infinity, 1.0)] {
            #expect(RingIndexing.reach(coverage: k, minimumMargin: m).isNaN, "\(k) \(m)")
        }
    }
}
