//
//  RingScaleTests.swift
//  `RingScale` (macTEM's RingScale.swift, whose beam-energy form is the kit's closed form). The app source is the oracle; no
//  Python reference. Checks: (1) the app's forms, written out here from its header (q = 2 sin(1/2 atan(R/L))/lambda,
//  R = L tan(2 asin(lambda q/2)), q = dq R), agree with the kit on a grid, to 1e-12 relative (the kit measured 8.2e-16 for q and
//  3.3e-15 for R, tools/kit-0.0.4/q-pins.py); (2) the linear helper and the file comparison, closed forms.
//

import Foundation
import Testing
import DiffractoGeometry

struct RingScaleTests {

    private static func rel(_ a: Double, _ b: Double) -> Double { abs(a - b) / abs(b) }

    @Test func appFormsEqualTheKitsOnAGrid() throws {
        var cases = 0
        for volts in [60e3, 80e3, 120e3, 200e3, 300e3] {
            for d in [1.2, 2.0391, 2.3546, 3.1355] {
                for radius in [20.0, 150.0, 400.0, 1100.0, 2000.0] {
                    let cal = try RingCalibration.calibrate(radius: radius, d: d, highTensionVolts: volts)
                    let L = try #require(cal.lengthPixels), lambda = try #require(cal.wavelengthAngstroms)
                    for r in [3.0, 77.7, radius, 1.5 * radius, 3000.0] {
                        let appQ = 2 * sin(0.5 * atan(r / L)) / lambda               // the app's header form
                        let q = RingScale.q(radius: r, calibration: cal)
                        #expect(Self.rel(q, appQ) < 1e-12, "q: V \(volts) d \(d) R \(radius) r \(r)")
                        #expect(q == RingCalibration.q(radius: r, cameraLengthPixels: L, wavelengthAngstroms: lambda))
                        let appR = L * tan(2 * asin(lambda * appQ / 2))
                        let back = RingScale.radius(q: q, calibration: cal)
                        #expect(Self.rel(back, appR) < 1e-12 && Self.rel(back, r) < 1e-12, "R: V \(volts) r \(r)")
                        cases += 1
                    }
                    // The calibrating ring is 1/d, and the linear form is the same multiplication Δq R.
                    #expect(Self.rel(RingScale.q(radius: radius, calibration: cal), 1 / d) < 1e-12)
                    #expect(RingScale.linear(radius: 7, deltaQ: cal.deltaQ) == cal.deltaQ * 7)
                }
            }
        }
        #expect(cases == 5 * 4 * 5 * 5)
    }

    @Test func withoutABeamEnergyTheScaleIsLinear() throws {
        let cal = try RingCalibration.calibrate(radius: 200, d: 2.0, highTensionVolts: nil)
        #expect(cal.lengthPixels == nil)
        #expect(cal.deltaQ == 0.5 / 200)                           // (1/d)/R = 0.0025 exactly
        #expect(RingScale.q(radius: 100, calibration: cal) == 0.25)
        #expect(RingScale.radius(q: 0.25, calibration: cal) == 100)
        #expect(RingScale.q(radius: 200, calibration: cal) == 0.5)
    }

    @Test func linearIsLowAtLargeRadiusOnAFlatDetector() throws {
        // Au {111} calibrates (V 200 kV, R 300 px); {200} then sits where the exact form puts it. The linear Δq·R of that
        // radius differs from the true 1/d by a small, nonzero relative amount (the app's header: about 1e-5 at 200 kV).
        let cal = try RingCalibration.calibrate(radius: 300, d: 4.0782 / 3.0.squareRoot(), highTensionVolts: 200e3)
        let r200 = RingScale.radius(q: 2 / 4.0782, calibration: cal)
        let exact = RingScale.q(radius: r200, calibration: cal)
        let linear = RingScale.linear(radius: r200, deltaQ: cal.deltaQ)
        #expect(Self.rel(exact, 2 / 4.0782) < 1e-12)
        let diff = abs(linear - exact) / exact
        #expect(diff > 1e-7 && diff < 1e-3, "linear differs from exact by \(diff)")
    }

    @Test func guardsGiveNaN() throws {
        let cal = try RingCalibration.calibrate(radius: 200, d: 2.0, highTensionVolts: 200e3)
        for bad in [0.0, -1, .nan, .infinity] {
            #expect(RingScale.q(radius: bad, calibration: cal).isNaN)
            #expect(RingScale.radius(q: bad, calibration: cal).isNaN)
            #expect(RingScale.linear(radius: bad, deltaQ: 0.01).isNaN)
            #expect(RingScale.linear(radius: 5, deltaQ: bad).isNaN)
        }
        #expect(RingScale.linear(radius: .greatestFiniteMagnitude, deltaQ: 2).isNaN)          // the product overflows
        // Past the Bragg limit (lambda q / 2 >= 1/sqrt 2) the exact inverse is NaN, never a negative radius.
        let lambda = try #require(cal.wavelengthAngstroms)
        #expect(RingScale.radius(q: 0.75 * 2 / lambda * 1.0, calibration: cal).isNaN)
    }

    @Test func fitVersusFileIsTheRelativeDifference() {
        #expect(abs(RingScale.fitVersusFile(fitDeltaQ: 0.0105, fileDeltaQ: 0.01) - 0.05) < 1e-15)
        #expect(RingScale.fitVersusFile(fitDeltaQ: 0.01, fileDeltaQ: 0.01) == 0)
        #expect(RingScale.fitVersusFile(fitDeltaQ: 0.005, fileDeltaQ: 0.01) == -0.5)
        for bad in [0.0, -1, .nan, .infinity] {
            #expect(RingScale.fitVersusFile(fitDeltaQ: bad, fileDeltaQ: 0.01).isNaN)
            #expect(RingScale.fitVersusFile(fitDeltaQ: 0.01, fileDeltaQ: bad).isNaN)
        }
    }
}
