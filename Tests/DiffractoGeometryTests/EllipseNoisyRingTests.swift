//
//  EllipseNoisyRingTests.swift
//  DiffractoKit slice 3b (docs/registration-slice-3.md section 3): the profile ring fit on the three noisy plants of reference commit
//  7d6fdbd1 `EllipseNoisyRingTests` (same names, same plants, same bands), with the plant renderer (SplitMix64 + Box-Muller,
//  Float pixels) ported into this file. The renderer is also the one tools/slice-3b/apps-parity/main.swift carries a copy of;
//  EllipseAppsParityTests proves the two render the same pixels (FNV-1a hash of every case) before it compares any fit.
//
//  The profile ring fit (`EllipseCalibration.fitAmorphousRing`, behind `fitBestAvailable`) is refined from py4DSTEM's start, both
//  ring sigmas at max(1, ringWidth / 4), 10-17 px for a +-20 px annulus around a 3 px ring. On a noisy ring the 11-parameter LM
//  from there can fall into a basin that carries no ring information, and the residual cap does not refuse it: noise inflates the
//  normalised residual of any model (reference `archive/v5/ring-fit-two-start-gateD-2026-10-14.md`, and
//  `archive/v5/ring-sigma-start-gateD-2026-10-13.md` for the first step). A second start at the ring's half-maximum width escapes
//  those basins, but alone it ends in a worse one on a 1/r^2 tail with noise. The fit therefore refines from both starts and keeps
//  the lower-cost fit (ADR 069).
//
//  The first two rings are noisy 192 x 192 plants on which py4DSTEM's start was accepted wrong, found by the reference implementations' scan of 2568
//  plants (128-256 px square, R 40-100, ellipticity 0-3 %, a 60 degree gap or none, noise 0.03 and 0.10, a few seeds each; no data
//  file) with this same renderer: it accepted a profile fit off by more than 1 % in (a + b) / 2 on 221 of them, all at noise 0.10;
//  87 of the 221 are the dead ring (a ring sigma at its 0.25 px floor, I1 at 0), 130 the drifting kind (ring sigmas of 1-7 px, the
//  centre 10-40 px off). The third ring has a 1/r^2 background (0.8 (R / r)^2 around the pattern centre, the central beam's tail)
//  under noise 0.10, and guards the other side: py4DSTEM's start is right on it and a half-maximum start alone is not (10 of 626
//  scanned tail plants ended more than 1 % off from the half-maximum start alone). These are the reference implementations' scans; this file
//  re-runs only the three plants.
//
//  Neither fit may be refused or fall back to the conic: the profile model finds the ring on all three, and the tests assert the
//  ring's mean radius, its b / a (what a caller applies to every Bragg vector) and its centre, against the plant.
//
//  RED UNDER (measured by the break-once run in tools/slice-3b/mutants.py): the second start dropped turns the dead ring and the
//  drifting ring red; the first start dropped turns the tail ring red.
//

import Testing
import Foundation
import DiffractoGeometry

/// I = exp(-((rho - 1) R)^2 / 2 s^2), rho the conic's elliptical radius in py4DSTEM's (A, B, C) convention with qx the row,
/// s = 3 px, R = (a + b) / 2, a / b = 1 + eps, a along theta; optional 60 degree gap starting at (47 seed) mod 360 degrees;
/// optional background tail * (R / max(r, 1))^2 with r the distance from the plant's centre (added after the gap); additive
/// Gaussian noise (SplitMix64 seeded seed + 1, Box-Muller pairs over even and odd linear index), Float pixels.
struct EllipsePlant {
    let size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, noise: Double, seed: Int
    /// Amplitude of the 1/r^2 background, 0 for none.
    var tail = 0.0
    var b: Double { 2 * radius / (2 + eps) }
    var a: Double { b * (1 + eps) }
    var centerQX: Double { Double(size) / 2 + 0.37 + 2.5 * sin(1.3 * Double(seed) + 0.2) }
    var centerQY: Double { Double(size) / 2 - 0.41 + 2.5 * cos(1.9 * Double(seed) + 0.7) }

    func render() -> DiffractionPattern {
        let theta = thetaDegrees * .pi / 180
        let sin2 = sin(theta) * sin(theta), cos2 = cos(theta) * cos(theta)
        let A = sin2 / (b * b) + cos2 / (a * a)
        let C = cos2 / (b * b) + sin2 / (a * a)
        let B = 2 * (b * b - a * a) * sin(theta) * cos(theta) / (a * a * b * b)
        let mean = (a + b) / 2
        let gapStart = Double((47 * seed) % 360)
        var state = UInt64(seed) &+ 1
        func uniform() -> Double {
            state &+= 0x9E37_79B9_7F4A_7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return (Double((z ^ (z >> 31)) >> 11) + 0.5) / 9_007_199_254_740_992.0
        }
        var pixels = [Float](repeating: 0, count: size * size)
        var spare = 0.0
        for qx in 0..<size {
            for qy in 0..<size {
                let dx = Double(qx) - centerQX, dy = Double(qy) - centerQY
                let rho = max(0, A * dx * dx + B * dx * dy + C * dy * dy).squareRoot()
                let distance = (rho - 1) * mean
                var value = exp(-(distance * distance) / (2 * 3.0 * 3.0))
                if gap > 0 {
                    var angle = atan2(dy, dx) * 180 / .pi
                    angle = (angle - gapStart).truncatingRemainder(dividingBy: 360)
                    if angle < 0 { angle += 360 }
                    if angle < gap { value = 0 }
                }
                if tail > 0 {
                    let r = max(hypot(dx, dy), 1)
                    value += tail * (radius / r) * (radius / r)
                }
                if noise > 0 {
                    if (qx * size + qy) % 2 == 0 {
                        let u1 = uniform(), u2 = uniform()
                        let r = (-2 * log(u1)).squareRoot()
                        value += noise * r * cos(2 * .pi * u2)
                        spare = noise * r * sin(2 * .pi * u2)
                    } else {
                        value += spare
                    }
                }
                pixels[qx * size + qy] = Float(value)
            }
        }
        return DiffractionPattern(width: size, height: size, pixels: pixels)
    }
}

struct EllipseNoisyRingTests {

    /// A click's error from the plant (1.3, -0.9 px), the annulus (b - 20, a + 20).
    private func assertRingIsFound(_ plant: EllipsePlant, _ name: String) {
        let fit: EllipseCalibrationFit
        do {
            fit = try EllipseCalibration.fitBestAvailable(
                pattern: plant.render(),
                centerQX: plant.centerQX + 1.3, centerQY: plant.centerQY - 0.9,
                innerRadius: plant.b - 20, outerRadius: plant.a + 20)
        } catch {
            Issue.record("\(name) was refused: \(error)")
            return
        }
        #expect(fit.model == .asymmetricGaussian,
                "\(name): fell back to \(fit.model) (\(fit.profileFallbackReason ?? "no reason"))")
        #expect(abs((fit.a + fit.b) / 2 / plant.radius - 1) <= 0.01,
                "\(name): mean radius \((fit.a + fit.b) / 2) against \(plant.radius)")
        #expect(abs(fit.b / fit.a - plant.b / plant.a) <= 0.01,
                "\(name): b / a \(fit.b / fit.a) against \(plant.b / plant.a)")
        #expect(hypot(fit.centerQX - plant.centerQX, fit.centerQY - plant.centerQY) <= 1.0,
                "\(name): centre (\(fit.centerQX), \(fit.centerQY)) against (\(plant.centerQX), \(plant.centerQY))")
    }

    /// The dead ring: R 70, eps 1 %, theta 0, a 60 degree gap, noise 0.10, seed 3.
    /// From py4DSTEM's start alone it ends at I1 0 with a ring sigma at 0.25 px and a mean radius 16.7 % low.
    @Test func testNoisyGappedRingIsNotDrainedIntoTheBackground() {
        assertRingIsFound(EllipsePlant(size: 192, radius: 70, eps: 0.01, thetaDegrees: 0, gap: 60, noise: 0.10, seed: 3),
                          "dead ring n192-R70-e10-t0.0-g60-s100-k3")
    }

    /// The drifting ring: R 70, eps 1 %, theta 112.5, no gap, noise 0.10, seed 2.
    /// From py4DSTEM's start alone it keeps ring sigmas of 2.9 and 3.9 px but ends 39 px off centre and 41 % low.
    @Test func testNoisyRingDoesNotDriftOffCentre() {
        assertRingIsFound(EllipsePlant(size: 192, radius: 70, eps: 0.01, thetaDegrees: 112.5, gap: 0, noise: 0.10, seed: 2),
                          "drifting ring n192-R70-e10-t112.5-g0-s100-k2")
    }

    /// The tail ring: R 65, eps 1 %, theta 112.5, no gap, noise 0.10, 1/r^2 tail 0.8, seed 5.
    /// From py4DSTEM's start it ends 0.002 % from the true radius; from the half-maximum start alone
    /// ends 2.3 % high with an inner ring sigma of 132 px.
    @Test func testNoisyRingOnAnInverseSquareTailKeepsItsRadius() {
        assertRingIsFound(EllipsePlant(size: 192, radius: 65, eps: 0.01, thetaDegrees: 112.5, gap: 0, noise: 0.10, seed: 5, tail: 0.8),
                          "tail ring n192-R65-e10-t112.5-g0-s100-k5-T80")
    }

    /// A ring-less fit on a real pattern: the mean diffraction pattern of the Ni65Cu35 ROI6 cube
    /// (`COPL_Ni65Cu35_C_ROI6_240911_..._bin_4_20240912.h5`, 200 x 200 positions of 32 x 32 px; float32 mean as the reference implementations compute
    /// it), annulus [11.125, 14.125] px about the centroid of its brightest 1 % (15.001983, 15.162129). Without the floor guard
    /// the fit accepts a ring-less profile fit there (I1 0, a ring sigma at the 0.25 px floor, b/a 0.947, where the conic gives
    /// 0.989 and the wider annulus 0.970); with it, the fit falls back to the conic and says why.
    ///
    /// DISABLED, not skipped silently: the pattern is derived from the owner's data, which this kit never copies, and
    /// `../data/README.md` holds no manifest line for it yet. When the line exists (the file under `$MACEM_DATA`, default
    /// `../data`, at the path below), remove the trait and set the path. The floor guard itself is held by
    /// EllipseCalibrationTests (`testARingFittedWithinPointThreePxOfTheFloorIsRefusedAndFallsBackToTheConic`).
    @Test(.disabled("owner's data, manifest line pending"))
    func testRinglessFitOnARealMeanPatternFallsBackToTheConic() throws {
        let root = ProcessInfo.processInfo.environment["MACEM_DATA"]
            ?? URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
                .deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("data").path
        // The path under data/ is fixed by the manifest line, once there is one.
        let url = URL(fileURLWithPath: root).appendingPathComponent("4dstem/ni65cu35-roi6-mean-32x32.f32")
        let data = try Data(contentsOf: url)
        let pixels = data.withUnsafeBytes { Array($0.bindMemory(to: Float.self)) }
        #expect(pixels.count == 32 * 32)
        let fit = try EllipseCalibration.fitBestAvailable(
            pattern: DiffractionPattern(width: 32, height: 32, pixels: pixels),
            centerQX: 15.001983, centerQY: 15.162129, innerRadius: 11.125, outerRadius: 14.125)
        #expect(fit.model == .conic, "a ring-less profile fit was accepted")
        #expect(fit.profile == nil)
        #expect(fit.profileFallbackReason?.contains("0.25 px floor") ?? false,
                "fallback reason was \(fit.profileFallbackReason ?? "nil")")
        #expect(abs(fit.b / fit.a - 0.989) <= 0.001)
    }
}
