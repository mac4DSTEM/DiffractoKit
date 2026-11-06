import Testing
@testable import DiffractoGeometry

@Suite("LinearScale.estimate")
struct LinearScaleEstimateTests {
    // Hand-computed: two shells at 20 and 23.094 px (ratio 1.1547, fcc {200}/{111}), reference 0.4247 Å⁻¹.
    @Test func clusterMedianAndShellRatio() throws {
        let p1 = [0.5, 19.8, 20.2, 19.9, 23.0, 23.2, 1.0]          // beam dropped; cluster mean (19.8+19.9+20.2)/3 = 19.9667
        let p2 = [20.0, 20.0, 23.1]                                  // r1 = 20, r2 = 23.1
        let p3 = [20.1, 23.0]                                        // r1 = 20.1, r2 = 23
        let e = try #require(LinearScale.estimate(radiiPerPosition: [p1, p2, p3],
                                                  referenceRadiusInvAngstrom: 0.4247,
                                                  secondShellRadiusInvAngstrom: 0.4247 * 1.1547))
        let r1s = [(19.8 + 19.9 + 20.2) / 3, 20.0, 20.1]
        #expect(e.observedRadiusPixels == 20.0)                       // median of (19.9667, 20, 20.1)
        #expect(abs(e.invAngstromPerPixel - 0.4247 / 20.0) < 1e-15)
        #expect(e.sampleCount == 3)
        #expect(abs(e.medianAbsoluteDeviationPixels - abs(r1s[0] - 20.0)) < 1e-12)  // |dev| = (0.0333, 0, 0.1): median 0.0333
        #expect(e.sameShellPeaksPerPosition == 2)                     // (3, 2, 1) sorted: median 2
        guard case .measured(let ratio, let expected, let n) = e.shellCheck else { Issue.record("not measured"); return }
        #expect(n == 3)
        #expect(abs(expected - 1.1547) < 1e-12)
        let r2 = [(23.0 + 23.2) / 2, 23.1, 23.0].sorted()[1]          // median 23.1
        #expect(abs(ratio - r2 / 20.0) < 1e-12)
        #expect(abs(try #require(e.shellCheck.mismatch) - abs(r2 / 20.0 / 1.1547 - 1)) < 1e-12)
    }

    @Test func oneShellIsReportedUnchecked() throws {
        let e = try #require(LinearScale.estimate(radiiPerPosition: [[10, 10.1], [10.05]],
                                                  referenceRadiusInvAngstrom: 0.5, secondShellRadiusInvAngstrom: 0.7))
        #expect(e.shellCheck == .notSelfChecked("only one shell is detectable"))
        let f = try #require(LinearScale.estimate(radiiPerPosition: [[10]], referenceRadiusInvAngstrom: 0.5,
                                                  secondShellRadiusInvAngstrom: nil))
        #expect(f.shellCheck == .notSelfChecked("this crystal has only one allowed shell within range"))
        #expect(f.shellCheck.mismatch == nil)
    }

    @Test func badReferenceAndNonIncreasingSecondShell() throws {
        #expect(LinearScale.estimate(radiiPerPosition: [[10]], referenceRadiusInvAngstrom: .nan, secondShellRadiusInvAngstrom: nil) == nil)
        #expect(LinearScale.estimate(radiiPerPosition: [[10]], referenceRadiusInvAngstrom: .infinity, secondShellRadiusInvAngstrom: nil) == nil)
        let e = try #require(LinearScale.estimate(radiiPerPosition: [[10]], referenceRadiusInvAngstrom: 0.5, secondShellRadiusInvAngstrom: 0.5))
        #expect(e.shellCheck == .notSelfChecked("this crystal has only one allowed shell within range"))
    }

    @Test func refusals() {
        #expect(LinearScale.estimate(radiiPerPosition: [[10]], referenceRadiusInvAngstrom: 0, secondShellRadiusInvAngstrom: nil) == nil)
        #expect(LinearScale.estimate(radiiPerPosition: [[1, .nan, 2]], referenceRadiusInvAngstrom: 0.5, secondShellRadiusInvAngstrom: nil) == nil)
        #expect(LinearScale.estimate(radiiPerPosition: [], referenceRadiusInvAngstrom: 0.5, secondShellRadiusInvAngstrom: nil) == nil)
    }
}
