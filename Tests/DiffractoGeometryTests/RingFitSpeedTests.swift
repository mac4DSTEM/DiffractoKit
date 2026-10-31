import Foundation
import Testing
@testable import DiffractoGeometry

/// `fitRing` at R 600 px on a 2048² pattern with a band of ±90 px took 24 s (#12; py4DSTEM f050d207 `fit_ellipse_amorphous_ring`
/// on one thread: 13.3 s). DIAGNOSIS: the profile optimiser ran most of its fits to the 80-pass cap, for two reasons.
/// (1) The profile fit refines two starts: py4DSTEM's, both ring sigmas at a quarter of the band (45 px for a 2 px ring), and one at the
/// ring's half-maximum width. From the first, both sigmas fell to the 0.25 px floor (the ring drained into the background, R 578 px)
/// within 8 passes and then crept along the bound, each accepted step lowering the cost by 1e-9 to 1e-5 relative, above the 1e-10 stop,
/// to the cap; the robust refit from there did the same; `refine` then refused it (`ringWidthAtFloor`). That was 166 of 174 passes
/// (151 of 163 with 1 % noise). The second start converged in 6 + 2.
/// (2) A refused step never stopped the optimiser: at the optimum every step is refused, the damping grows tenfold a pass, and a
/// converged fit ran on to the cap (14 accepted and 66 refused passes on a clean R 200 px ring). `fit1D` got this stop on 2026-10-12.
/// What would have refuted the causes: a pass count that stays high with these runs removed, or a time that does not fall with it.
/// FIX: a refinement whose two ring sigmas sit at the floor for 5 accepted steps in a row stops (the robust refit still runs from there,
/// and may lift a sigma off the floor, as before), and a refused step below the accepted-step stop (1e-7) ends a run that has accepted a
/// step. Not done: refusing a double-floor fit before its robust refit (that refit lifted a sigma at R 296 px for a ring at 300 and ran
/// 80 passes before the ring evidence refused it), because it could refuse a fit that the refit rescues when only one start runs. The returned radius is unchanged bit for bit on the measured rings (R 600 ±90 clean and with 1 % noise, R 600 ±20,
/// R 300 ±20, R 100 ±20 clean and noisy, R 200 ±60); R 600 ±90 now takes 4.7 s clean and 5.4 s with 1 % noise. The remaining cost is about 145 ms a pass at
/// 680 000 pixels: 22 model evaluations for the finite-difference Jacobian and the normal matrix.
struct RingFitSpeedTests {
    final class PassCounter: @unchecked Sendable {
        private let lock = NSLock()
        private var count = 0
        func tick() { lock.lock(); count += 1; lock.unlock() }
        var value: Int { lock.lock(); defer { lock.unlock() }; return count }
    }

    static func pattern(n: Int, radius: Double, sigma: Double) -> DiffractionPattern {
        let c = Double(n / 2)
        var px = [Float](repeating: 0, count: n * n)
        for y in 0 ..< n { for x in 0 ..< n {
            let dx = Double(x) - c, dy = Double(y) - c
            let r = (dx * dx + dy * dy).squareRoot()
            px[y * n + x] = Float(60 + 3000 * exp(-(r - radius) * (r - radius) / (2 * sigma * sigma)))
        } }
        return DiffractionPattern(width: n, height: n, pixels: px)
    }

    /// Progress reports through `fitRing`: one per optimiser pass of the conic and of every profile refinement, plus a few.
    static func fitCountingPasses(n: Int, radius: Double, half: Double, sigma: Double) throws -> (RingRadiusFit, Int) {
        let p = Self.pattern(n: n, radius: radius, sigma: sigma)
        let passes = PassCounter()
        let c = Double(n / 2)
        let fit = try RingCalibration.fitRing(pattern: p, centerQX: c, centerQY: c, innerRadius: radius - half,
                                              outerRadius: radius + half, progress: { _ in passes.tick() })
        return (fit, passes.value)
    }

    /// The py4DSTEM start ends with both ring sigmas at the floor (refused); the half-maximum start gives the ring.
    @Test func aStartWhoseRingSigmasSitAtTheFloorStopsEarly() throws {
        let (fit, passes) = try Self.fitCountingPasses(n: 660, radius: 300, half: 20, sigma: 2)
        let ring = try #require(fit.fit)
        #expect(abs((ring.a + ring.b) / 2 - 300) < 1e-6)
        #expect(fit.isUsable)
        #expect(passes < 120, Comment(rawValue: "progress reports: \(passes)"))    // 174 before the stop, 104 after
    }

    /// Both starts converge; the runs ended only at the 80-pass cap after the optimum was reached.
    @Test func aConvergedRefinementStopsOnARefusedStep() throws {
        let (fit, passes) = try Self.fitCountingPasses(n: 540, radius: 200, half: 60, sigma: 1.5)
        let ring = try #require(fit.fit)
        #expect(abs((ring.a + ring.b) / 2 - 200) < 1e-6)
        #expect(fit.isUsable)
        #expect(passes < 60, Comment(rawValue: "progress reports: \(passes)"))     // 90 before the stop, 43 after
    }
}
