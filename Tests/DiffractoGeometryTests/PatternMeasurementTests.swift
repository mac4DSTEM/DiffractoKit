import Foundation
import Testing
@testable import DiffractoGeometry

/// Kit asks 1, 2, 3 and 9 of macTEM's diffraction room (DiffractoKit 0.10.0), against planted truth. No Python reference exists for these (see the file
/// headers of PatternMeasurement.swift and LinearScaleFit.swift); `validation` is "none".
struct PatternMeasurementTests {
    struct Lcg { var state: UInt64
        mutating func uniform() -> Double { state = state &* 6364136223846793005 &+ 1442695040888963407; return Double(state >> 11) / Double(1 << 53) } }

    static func pattern(n: Int = 512, noise: Double = 0.01, seed: UInt64 = 7, _ intensity: (Double, Double) -> Double) -> DiffractionPattern {
        var g = Lcg(state: seed)
        var px = [Float](repeating: 0, count: n * n)
        for y in 0 ..< n { for x in 0 ..< n {
            px[y * n + x] = Float(intensity(Double(y), Double(x)) * (1 + noise * (g.uniform() - 0.5) * 2))
        } }
        return DiffractionPattern(width: n, height: n, pixels: px)
    }

    // MARK: centre finder (ask 1)

    @Test func theBeamCentreOfAGaussianBeamWithRingsAndNoise() throws {
        let cy = 255.43, cx = 257.81
        let p = Self.pattern { y, x in
            let r = ((y - cy) * (y - cy) + (x - cx) * (x - cx)).squareRoot()
            return 60 + 20000 * exp(-r * r / 32) + 3000 * exp(-(r - 34) * (r - 34) / 2.88) + 2200 * exp(-(r - 39.2) * (r - 39.2) / 2.88)
        }
        let fit = try #require(PatternCentre.find(pattern: p))
        #expect(abs(fit.row - cy) < 0.02 && abs(fit.column - cx) < 0.02, Comment(rawValue: "\(fit)"))
        #expect(abs(fit.beamRadiusPixels - 4 * (2 * log(2.0)).squareRoot()) < 0.4, Comment(rawValue: "\(fit.beamRadiusPixels)"))      // the half-height disc of sigma 4: 4.71 px
        #expect(fit.validation == "none")
    }

    @Test func theBeamCentreOfAFlatTopDiscAtSubpixelPositions() throws {
        for (cy, cx) in [(200.25, 300.75), (256.5, 256.5), (100.9, 411.1)] {
            let p = Self.pattern(noise: 0.02) { y, x in
                let r = ((y - cy) * (y - cy) + (x - cx) * (x - cx)).squareRoot()
                return 100 + 8000 / (1 + exp((r - 18) / 0.8))
            }
            let fit = try #require(PatternCentre.find(pattern: p))
            #expect(abs(fit.row - cy) < 0.05 && abs(fit.column - cx) < 0.05, Comment(rawValue: "\(cy),\(cx): \(fit)"))
            #expect(abs(fit.beamRadiusPixels - 18) < 0.5)
        }
    }

    @Test func noBeamNoCentre() {
        // pure noise, a beam at the edge, and a pattern whose brightest region is a ring (the beam stop case)
        #expect(PatternCentre.find(pattern: Self.pattern(noise: 0.05) { _, _ in 100 }) == nil)
        let edge = Self.pattern { y, x in 100 + 9000 * exp(-(y * y + (x - 200) * (x - 200)) / 40) }
        #expect(PatternCentre.find(pattern: edge) == nil)
        let tiny = DiffractionPattern(width: 4, height: 4, pixels: [Float](repeating: 1, count: 16))
        #expect(PatternCentre.find(pattern: tiny) == nil)
        #expect(PatternCentre.find(pattern: DiffractionPattern(width: 16, height: 16, pixels: [Float](repeating: .nan, count: 256))) == nil)
    }

    // MARK: Δq from several rings (ask 2)

    @Test func deltaQFromSeveralRingsIsExactWithoutNoiseAndCoversTheTruthWithIt() throws {
        let a = 4.0782, dq = 0.0125
        let gs = [3.0, 4.0, 8.0, 11.0].map { $0.squareRoot() / a }
        let exact = try #require(LinearScale.fit(ringRadiiPixels: gs.map { $0 / dq }, indexedInvAngstrom: gs))
        #expect(abs(exact.invAngstromPerPixel / dq - 1) < 1e-12 && (exact.sigma ?? 1) < 1e-12 && exact.ringCount == 4 && exact.validation == "none")
        // 400 noisy draws (0.05 px per radius): the truth is inside one sigma in about 68 % of them
        var g = Lcg(state: 99)
        var inside = 0, worst = 0.0
        for _ in 0 ..< 400 {
            let radii = gs.map { $0 / dq + 0.05 * (g.uniform() - 0.5) * 2 * 3.0.squareRoot() }
            let f = try #require(LinearScale.fit(ringRadiiPixels: radii, indexedInvAngstrom: gs))
            if abs(f.invAngstromPerPixel - dq) <= (f.sigma ?? 0) { inside += 1 }
            worst = max(worst, abs(f.invAngstromPerPixel / dq - 1))
        }
        #expect(inside > 400 * 55 / 100 && inside < 400 * 82 / 100, Comment(rawValue: "inside \(inside) of 400"))
        #expect(worst < 0.002)
    }

    @Test func oneRingHasNoSigmaAndBadInputIsRefused() throws {
        let one = try #require(LinearScale.fit(ringRadiiPixels: [50], indexedInvAngstrom: [0.5]))
        #expect(one.invAngstromPerPixel == 0.01 && one.sigma == nil)
        let weighted = try #require(LinearScale.fit(ringRadiiPixels: [50, 100], indexedInvAngstrom: [0.5, 1.1], weights: [1, 1e-9]))
        #expect(abs(weighted.invAngstromPerPixel - 0.01) < 1e-6)
        #expect(LinearScale.fit(ringRadiiPixels: [], indexedInvAngstrom: []) == nil)
        #expect(LinearScale.fit(ringRadiiPixels: [1, 2], indexedInvAngstrom: [1]) == nil)
        #expect(LinearScale.fit(ringRadiiPixels: [1, -2], indexedInvAngstrom: [1, 2]) == nil)
        #expect(LinearScale.fit(ringRadiiPixels: [1, 2], indexedInvAngstrom: [1, .nan]) == nil)
        #expect(LinearScale.fit(ringRadiiPixels: [1, 2], indexedInvAngstrom: [1, 2], weights: [1, 0]) == nil)
    }

    // MARK: ellipse-corrected radial profile (ask 3)

    @Test func theProfileOfAnEllipticalRingIsSharpOnceCorrected() throws {
        // a ring on an ellipse (semi-axes 40 and 36, major axis 0.5 rad from the row axis) about a known centre
        let a = 40.0, b = 36.0, theta = 0.5, cy = 250.0, cx = 260.0
        let p = Self.pattern(noise: 0) { y, x in
            let dy = y - cy, dx = x - cx
            // py4DSTEM's frame: theta is the angle of the major axis from the first (row) axis toward the second (column) axis
            let u = dy * cos(theta) + dx * sin(theta), v = -dy * sin(theta) + dx * cos(theta)
            let e = ((u / a) * (u / a) + (v / b) * (v / b)).squareRoot()     // 1 on the ellipse
            return 100 + 5000 * exp(-(e - 1) * (e - 1) * a * a / (2 * 1.5 * 1.5))
        }
        func width(_ prof: PatternRadialProfile) -> Double {
            let top = prof.mean.compactMap { $0.isNaN ? nil : $0 }.max()!
            return Double(prof.mean.filter { !$0.isNaN && $0 > 100 + 0.5 * (top - 100) }.count) * prof.binWidthPixels
        }
        let plain = try #require(PatternRadialProfile.about(pattern: p, centerRow: cy, centerColumn: cx, maximumRadius: 60))
        let corrected = try #require(PatternRadialProfile.about(pattern: p, centerRow: cy, centerColumn: cx, ellipse: (a, b, theta), maximumRadius: 60))
        #expect(corrected.mean.count == 60 && corrected.count.reduce(0, +) > 0 && corrected.validation == "none")
        #expect(width(corrected) <= width(plain) - 2, Comment(rawValue: "plain \(width(plain)) corrected \(width(corrected))"))
        // the corrected peak sits at one radius: the half-height width of a sigma 1.5 px ring is about 3.5 px, up to two bins of discretisation
        #expect(width(corrected) <= 5)
    }

    @Test func theProfileOfACircularRingIsTheSameWithAndWithoutACircularEllipse() throws {
        let p = Self.pattern(noise: 0) { y, x in
            let r = ((y - 256) * (y - 256) + (x - 256) * (x - 256)).squareRoot()
            return 50 + 1000 * exp(-(r - 70) * (r - 70) / 8)
        }
        let plain = try #require(PatternRadialProfile.about(pattern: p, centerRow: 256, centerColumn: 256, maximumRadius: 100))
        let circle = try #require(PatternRadialProfile.about(pattern: p, centerRow: 256, centerColumn: 256, ellipse: (30, 30, 0.7), maximumRadius: 100))
        for k in 0 ..< 100 { #expect(abs(plain.mean[k] - circle.mean[k]) < 1e-3 * max(1, abs(plain.mean[k]))) }
        let peak = plain.mean.indices.max { plain.mean[$0] < plain.mean[$1] }!
        #expect(abs(plain.radiusPixels[peak] - 70) <= 1)
        #expect(PatternRadialProfile.about(pattern: p, centerRow: 600, centerColumn: 0) == nil)
        #expect(PatternRadialProfile.about(pattern: p, centerRow: 256, centerColumn: 256, binWidth: 0) == nil)
    }

    // MARK: spot refiner (ask 9)

    /// A window narrower than the spot cuts its tail and pulls the centroid toward the click; a wide window does not (the doc comment's claim).
    @Test func aWindowThatCutsTheTailBiasesTowardTheClick() throws {
        let cy = 100.4, cx = 150.4
        let p = Self.pattern(noise: 0) { y, x in 80 + 4000 * exp(-((y - cy) * (y - cy) + (x - cx) * (x - cx)) / (2 * 3.0 * 3.0)) }
        let narrow = try #require(PatternSpot.refine(pattern: p, row: 103, column: 150, half: 2))
        let wide = try #require(PatternSpot.refine(pattern: p, row: 103, column: 150, half: 12))
        #expect(abs(wide.row - cy) < 0.01)
        #expect(narrow.row > cy + 0.01 && narrow.row < 103, Comment(rawValue: "narrow window: \(narrow.row)"))     // pulled toward the click at 103
    }

    @Test func aSpotIsRefinedToASubpixelCentre() throws {
        for (cy, cx) in [(100.37, 150.62), (200.5, 300.5), (64.2, 400.8)] {
            let p = Self.pattern(noise: 0.01) { y, x in
                80 + 4000 * exp(-((y - cy) * (y - cy) + (x - cx) * (x - cx)) / (2 * 2.0 * 2.0))
            }
            // a click up to three pixels off
            for (oy, ox) in [(0.0, 0.0), (2.0, -1.0), (-3.0, 3.0)] {
                let r = try #require(PatternSpot.refine(pattern: p, row: cy.rounded() + oy, column: cx.rounded() + ox))
                #expect(abs(r.row - cy) < 0.05 && abs(r.column - cx) < 0.05, Comment(rawValue: "\(cy),\(cx) click offset \(oy),\(ox): \(r)"))
            }
        }
        // refusals: a window off the pattern, no signal, no finite pixels, a bad window
        let flat = DiffractionPattern(width: 32, height: 32, pixels: [Float](repeating: 5, count: 1024))
        #expect(PatternSpot.refine(pattern: flat, row: 16, column: 16) == nil)
        #expect(PatternSpot.refine(pattern: flat, row: -50, column: 16) == nil)
        #expect(PatternSpot.refine(pattern: flat, row: .nan, column: 16) == nil)
        #expect(PatternSpot.refine(pattern: flat, row: 16, column: 16, half: 0) == nil)
        #expect(PatternSpot.refine(pattern: flat, row: 1e300, column: 16) == nil && PatternSpot.refine(pattern: flat, row: 16, column: -1e300) == nil)
        #expect(PatternRadialProfile.about(pattern: flat, centerRow: 16, centerColumn: 16, ellipse: (.nan, 1, 0)) == nil)
        #expect(PatternRadialProfile.about(pattern: flat, centerRow: 16, centerColumn: 16, maximumRadius: 1e12)?.mean.count ?? 0 > 0)
    }
}
