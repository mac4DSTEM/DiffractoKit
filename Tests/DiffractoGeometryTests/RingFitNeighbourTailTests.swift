import Foundation
import Testing
@testable import DiffractoGeometry

/// The {111} ring of a planted Au pattern read 0.1 % low through `RingCalibration.fitRing` (macTEM's open item, 2026-10-22: R 33.977 planted, 33.942 fitted).
/// DIAGNOSIS (0.10.0): the cause is the tail of the NEIGHBOURING ring inside the annulus, not the fit's ring model, the noise, the central beam, the
/// background or the ring radius. A single planted ring is recovered to 3e-8 px for sigma 0.6, 1.2, 2.4 and R 34, 68, 136; with 1 % noise to 4e-5 px.
/// With the {200} ring 5.3 px further out (amplitude 2200 against 3000) the bias is -0.041 px for the annulus 29 to 37, -0.003 px for 29 to 36 and
/// -3e-6 px for 29 to 35; it grows linearly with the neighbour's amplitude (-0.010, -0.021, -0.041, -0.068 px at 550, 1100, 2200, 4400). What would
/// have refuted the cause: a bias that stays with the neighbour removed, or that does not shrink with the annulus edge. The profile model (one
/// asymmetric Gaussian on a smooth background) has no term for a second ring, so it absorbs the tail. Not changed: a number that agrees with the
/// model's definition; an annulus that ends before the neighbour's tail (the rule below) removes the bias.
struct RingFitNeighbourTailTests {
    static let r111 = 3.0.squareRoot() / (4.0782 * 0.0125)
    static let r200 = 4.0.squareRoot() / (4.0782 * 0.0125)

    static func pattern(rings: [(Double, Double)], sigma: Double = 1.2, noise: Double = 0, n: Int = 512) -> DiffractionPattern {
        var state: UInt64 = 20261022
        func u() -> Double { state = state &* 6364136223846793005 &+ 1442695040888963407; return Double(state >> 11) / Double(1 << 53) }
        let c = Double(n / 2)
        var px = [Float](repeating: 0, count: n * n)
        for y in 0 ..< n { for x in 0 ..< n {
            let dx = Double(x) - c, dy = Double(y) - c
            let r = (dx * dx + dy * dy).squareRoot()
            var v = 60.0
            for (rr, a) in rings { v += a * exp(-(r - rr) * (r - rr) / (2 * sigma * sigma)) }
            if noise > 0 { v *= 1 + noise * (u() - 0.5) * 2 }
            px[y * n + x] = Float(v)
        } }
        return DiffractionPattern(width: n, height: n, pixels: px)
    }

    static func radius(_ p: DiffractionPattern, _ inner: Double, _ outer: Double) throws -> Double {
        let c = Double(p.width / 2)
        let fit = try #require(try RingCalibration.fitRing(pattern: p, centerQX: c, centerQY: c, innerRadius: inner, outerRadius: outer).fit)
        return (fit.a + fit.b) / 2
    }

    @Test func aSingleRingIsRecoveredWhateverItsWidthAndRadius() throws {
        for (R, sigma, inner, outer) in [(Self.r111, 0.6, 27.0, 41.0), (Self.r111, 1.2, 25.0, 43.0), (Self.r111, 2.4, 27.0, 41.0), (68.0, 1.2, 60.0, 76.0), (136.0, 1.2, 128.0, 144.0)] {
            let p = Self.pattern(rings: [(R, 3000)], sigma: sigma)
            #expect(abs(try Self.radius(p, inner, outer) - R) < 1e-6, Comment(rawValue: "R \(R) sigma \(sigma)"))
        }
        let noisy = Self.pattern(rings: [(Self.r111, 3000)], noise: 0.01)
        #expect(abs(try Self.radius(noisy, 29, 37) - Self.r111) < 5e-4)
    }

    @Test func aNeighbouringRingsTailInTheAnnulusPullsTheRadiusLow() throws {
        let p = Self.pattern(rings: [(Self.r111, 3000), (Self.r200, 2200)])
        let bias = try [35.0, 36.0, 37.0, 38.0].map { try Self.radius(p, 29, $0) - Self.r111 }
        #expect(abs(bias[0]) < 1e-4)                                 // the annulus ends before the tail: exact
        #expect(bias[1] < 0 && bias[2] < bias[1] && bias[3] < bias[2])     // monotone in how much tail is inside
        #expect(abs(bias[2] + 0.041) < 0.002, Comment(rawValue: "demo annulus bias \(bias[2])"))      // the macTEM case: -0.035 px with 1 % noise
        let amplitudes = try [0.0, 550.0, 1100.0, 2200.0].map { amp -> Double in
            try Self.radius(Self.pattern(rings: [(Self.r111, 3000), (Self.r200, amp)]), 29, 37) - Self.r111
        }
        #expect(abs(amplitudes[0]) < 1e-6)
        #expect(abs(amplitudes[2] / amplitudes[1] - 2) < 0.1 && abs(amplitudes[3] / amplitudes[2] - 2) < 0.1)    // linear in the neighbour's amplitude
        // the other ring: the {200} annulus 36 to 43 has no {111} tail inside, so it is exact
        #expect(abs(try Self.radius(p, 36, 43) - Self.r200) < 0.02)
    }
}
