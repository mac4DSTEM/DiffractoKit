//
//  RingOverlayTests.swift
//  `RingOverlay` (ported from macTEM's RingOverlay.swift; the app source is the oracle, no Python reference). Pins
//  (Fixtures/ring-overlay-pins.json) come from tools/mactem-pieces/ring-overlay-pins.py, a pure-Python transcription of the
//  app source, NOT from this Swift. Closed forms (circle, axis vertex) are hand-derived and checked beside them.
//

import Foundation
import Testing
import DiffractoGeometry

private struct OverlayPins: Decodable {
    struct Ring: Decodable { let x, y, a, b, theta: Double }
    struct Dist: Decodable { let px, py, d: Double }
    struct Case: Decodable { let ring: Ring; let polyline8: [[Double]]; let distances: [Dist] }
    let cases: [Case]
}

struct RingOverlayTests {

    private static func pins() throws -> OverlayPins {
        let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/ring-overlay-pins.json")
        return try JSONDecoder().decode(OverlayPins.self, from: Data(contentsOf: url))
    }

    @Test func polylineAndDistancesMatchThePythonTranscription() throws {
        for c in try Self.pins().cases {
            let e = RingEllipse(x: c.ring.x, y: c.ring.y, a: c.ring.a, b: c.ring.b, theta: c.ring.theta)
            let line = RingOverlay.polyline(e, count: 8)
            #expect(line.count == 8)
            for (p, q) in zip(line, c.polyline8) {
                #expect(abs(p.x - q[0]) <= 1e-12 * max(1, abs(q[0])) && abs(p.y - q[1]) <= 1e-12 * max(1, abs(q[1])))
            }
            for d in c.distances {
                let got = RingOverlay.distance(fromX: d.px, y: d.py, to: e)
                #expect(abs(got - d.d) <= 1e-9 * max(1, d.d), "ring \(e) point (\(d.px), \(d.py)): \(got) vs \(d.d)")
            }
        }
    }

    @Test func circleClosedForms() {
        let c = RingEllipse(x: 100, y: 80, a: 50, b: 50, theta: 0.3)
        #expect(c.radius == 50)
        // A point at distance s from the centre is |s - r| from the circle.
        for (px, py, s) in [(100.0, 80.0, 0.0), (130.0, 120.0, 50.0), (103.0, 84.0, 5.0), (190.0, 80.0, 90.0)] {
            #expect(abs(RingOverlay.distance(fromX: px, y: py, to: c) - abs(s - 50)) < 1e-9)
        }
        // Every polyline point is at the radius, and the first is the end of the first semi-axis.
        for p in RingOverlay.polyline(c, count: 37) { #expect(abs(hypot(p.x - 100, p.y - 80) - 50) < 1e-12) }
        let first = RingOverlay.polyline(c)[0]
        #expect(RingOverlay.polyline(c).count == RingOverlay.polylinePoints && RingOverlay.polylinePoints == 180)
        #expect(abs(first.x - (100 + 50 * sin(0.3))) < 1e-12 && abs(first.y - (80 + 50 * cos(0.3))) < 1e-12)
    }

    @Test func majorAxisVertexAndImageAxes() {
        // theta = 0: the first semi-axis lies along the ROW axis (y). Vertices: (x, y +- a) and (x +- b, y).
        let e = RingEllipse(x: 64.5, y: 70.25, a: 40, b: 30, theta: 0)
        let line = RingOverlay.polyline(e, count: 4)            // phi = 0, 90, 180, 270 degrees
        #expect(abs(line[0].x - 64.5) < 1e-12 && abs(line[0].y - 110.25) < 1e-12)
        #expect(abs(line[1].x - 94.5) < 1e-12 && abs(line[1].y - 70.25) < 1e-12)
        #expect(abs(line[2].y - 30.25) < 1e-12 && abs(line[3].x - 34.5) < 1e-12)
        // A point on the major axis outside the ellipse is nearest the vertex.
        #expect(abs(RingOverlay.distance(fromX: 64.5, y: 115.25, to: e) - 5) < 1e-9)
        #expect(abs(RingOverlay.distance(fromX: 64.5, y: 70.25, to: e) - 30) < 1e-9)   // the centre: the minor semi-axis
    }

    @Test func fitConversionSwapsRowAndColumn() {
        let fit = EllipseCalibrationFit(centerQX: 70.25, centerQY: 64.5, a: 40, b: 30, theta: 0.7, normalizedResidual: 0,
                                        conicResidual: 0, sampleCount: 0, occupiedAngularBins: 0, model: .conic, profile: nil,
                                        profileFallbackReason: nil)
        let e = RingOverlay.ellipse(of: fit)
        #expect(e == RingEllipse(x: 64.5, y: 70.25, a: 40, b: 30, theta: 0.7))
    }

    @Test func predictedRingScalesBothAxes() throws {
        let e = RingEllipse(x: 10, y: 20, a: 40, b: 30, theta: 0.7)          // mean radius 35
        let p = try #require(RingOverlay.predicted(from: e, radius: 70))
        #expect(p == RingEllipse(x: 10, y: 20, a: 80, b: 60, theta: 0.7))
        for bad in [0.0, -1, .nan, .infinity] { #expect(RingOverlay.predicted(from: e, radius: bad) == nil) }
        #expect(RingOverlay.predicted(from: RingEllipse(x: 0, y: 0, a: .nan, b: 1, theta: 0), radius: 1) == nil)
        #expect(RingOverlay.predicted(from: RingEllipse(x: 0, y: 0, a: 0, b: 0, theta: 0), radius: 1) == nil)
    }

    @Test func guardsGiveNaNAndEmpty() {
        let e = RingEllipse(x: 0, y: 0, a: 5, b: 3, theta: 0)
        #expect(RingOverlay.distance(fromX: .nan, y: 0, to: e).isNaN)
        #expect(RingOverlay.distance(fromX: 0, y: .infinity, to: e).isNaN)
        #expect(RingOverlay.distance(fromX: 0, y: 0, to: RingEllipse(x: 0, y: 0, a: 0, b: 3, theta: 0)).isNaN)
        #expect(RingOverlay.distance(fromX: 0, y: 0, to: RingEllipse(x: 0, y: 0, a: 5, b: -3, theta: 0)).isNaN)
        #expect(RingOverlay.distance(fromX: 0, y: 0, to: RingEllipse(x: 0, y: 0, a: 5, b: 3, theta: .nan)).isNaN)
        #expect(RingOverlay.polyline(e, count: 0).isEmpty)
        #expect(RingOverlay.polyline(e, count: -3).isEmpty)
    }
}
