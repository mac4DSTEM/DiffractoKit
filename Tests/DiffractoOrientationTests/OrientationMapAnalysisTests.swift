import DiffractoKitTestSupport
import Testing
import simd
@testable import DiffractoOrientation

@Suite("OrientationMapAnalysis")
struct OrientationMapAnalysisTests {
    private func r(_ t: Int, score: Float = 1, second: Float = 0, euler: EulerAngles = .zero) -> OrientationResult {
        OrientationResult(templateIndex: t, euler: euler, score: score, secondScore: second, phaseID: 0)
    }

    @Test func eulerTextAndBounds() {
        let e = EulerAngles(phi1: Float.pi / 2, Phi: Float.pi / 4, phi2: Float.pi)
        let res = [r(0, euler: e), .empty]
        #expect(OrientationMapAnalysis.eulerText(results: res, width: 2, height: 1, x: 0, y: 0) == "90.0°, 45.0°, 180.0°")
        #expect(OrientationMapAnalysis.eulerText(results: res, width: 2, height: 1, x: 1, y: 0) == nil)
        #expect(OrientationMapAnalysis.eulerText(results: res, width: 2, height: 1, x: 2, y: 0) == nil)
        #expect(OrientationMapAnalysis.eulerText(results: res, width: 3, height: 1, x: 0, y: 0) == nil)
    }

    @Test func thresholdAndFraction() {
        // reliability = 1 - second/score: 0, .25, .5, .75 on four matched positions; one unmatched is ignored
        let res = [r(0, score: 4, second: 4), r(1, score: 4, second: 3), r(2, score: 4, second: 2), r(3, score: 4, second: 1), .empty]
        #expect(OrientationMapAnalysis.reliabilityThreshold(results: res, percentile: 0) == 0)
        #expect(OrientationMapAnalysis.reliabilityThreshold(results: res, percentile: 0.5) == 0.5)   // index round(1.5) = 2
        #expect(OrientationMapAnalysis.reliabilityThreshold(results: res, percentile: 2) == 0.75)
        #expect(OrientationMapAnalysis.fractionOfMatched(results: res, reliabilityAtLeast: 0.5) == 0.5)
        #expect(OrientationMapAnalysis.reliabilityThreshold(results: [.empty], percentile: 0.1) == nil)
        #expect(OrientationMapAnalysis.fractionOfMatched(results: [], reliabilityAtLeast: 0) == nil)
    }

    @Test func planesAndWrap() {
        let res = [r(0, euler: EulerAngles(phi1: -Float.pi / 2, Phi: 2 * Float.pi, phi2: 3 * Float.pi))]
        #expect(abs(OrientationMapAnalysis.phi1Plane(res)[0] - 0.75) < 1e-6)
        #expect(abs(OrientationMapAnalysis.phi2Plane(res)[0] - 0.5) < 1e-6)
        #expect(OrientationMapAnalysis.PhiPlane(res) == [1])
    }

    @Test func ipfBytesMaskAndUnmatched() {
        let res = [r(0, score: 4, second: 3), r(0, score: 4, second: 0), .empty]
        let rgba = OrientationMapAnalysis.ipfZRGBA(results: res, maskingReliabilityBelow: 0.5) { _ in [1, 0.5, 0] }
        #expect(Array(rgba[0..<4]) == [97, 97, 97, 255])        // reliability .25 < .5: grey .38 -> round(96.9)
        #expect(Array(rgba[4..<8]) == [255, 128, 0, 255])       // coloured
        #expect(Array(rgba[8..<12]) == [0, 0, 0, 255])          // unmatched: black
        let open = OrientationMapAnalysis.ipfZRGBA(results: res, maskingReliabilityBelow: nil) { _ in [1, 0, 0] }
        #expect(Array(open[0..<4]) == [255, 0, 0, 255])
    }
}
