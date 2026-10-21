//
//  EllipseStopRuleTests.swift
//  DiffractoKit slice 3b: `EllipseCalibration.fit1D`'s Levenberg-Marquardt loop (the reference implementation 2026-10-12) once tested convergence only
//  after an ACCEPTED step. At the optimum the cost change of a step is roundoff, every step was refused, damping grew until the
//  candidate equalled the point, and the 100-pass cap threw `didNotConverge` on about one clean ring in seven
//  (reference `archive/v5/ellipse-stop-rule-gateD-2026-10-12.md`). Ported from reference commit 7d6fdbd1
//  `EllipseStopRuleTests.swift` (one test, same name, same four rings, same bands).
//
//  These four rings are clean 96 x 96 plants the old rule refused (the reference implementations scanned 675 such plants, 111 of them refused, with
//  this same renderer; the refusal depends on roundoff, so the set is a property of this renderer and this arithmetic, and the
//  check below is that they FIT). Ground truth is analytic: the plant's centre and ellipse, and the objective's own minimum on a
//  ring of width s, R (1 + 1.75 (s/R)^2) for the mean radius (an r dr measure, the conic's r^2, the squared intensity; first
//  order, so it sits 0.007 px under the fit at R = 30, which is the tolerance's origin).
//
//  RED UNDER: the stop test applied only after an accepted step (the `else` branch's `relativeStep < 1e-10` removed): one or more
//  of the four rings throws `didNotConverge`.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import DiffractoGeometry

struct EllipseStopRuleTests {
    private let size = 96
    private let radius = 30.0
    private let sigma = 3.0

    /// I = exp(-((rho - 1) R)^2 / 2 s^2), rho the conic's elliptical radius,
    /// py4DSTEM's (A, B, C) convention with qx the row.
    private func ring(centerQX: Double, centerQY: Double, a: Double, b: Double, theta: Double) -> DiffractionPattern {
        let sin2 = sin(theta) * sin(theta), cos2 = cos(theta) * cos(theta)
        let A = sin2 / (b * b) + cos2 / (a * a)
        let C = cos2 / (b * b) + sin2 / (a * a)
        let B = 2 * (b * b - a * a) * sin(theta) * cos(theta) / (a * a * b * b)
        let mean = (a + b) / 2
        var pixels = [Float](repeating: 0, count: size * size)
        for qx in 0..<size {
            for qy in 0..<size {
                let dx = Double(qx) - centerQX, dy = Double(qy) - centerQY
                let rho = max(0, A * dx * dx + B * dx * dy + C * dy * dy).squareRoot()
                let distance = (rho - 1) * mean
                pixels[qx * size + qy] = Float(exp(-(distance * distance) / (2 * sigma * sigma)))
            }
        }
        return DiffractionPattern(width: size, height: size, pixels: pixels)
    }

    @Test func testCleanRingsTheOldStopRuleRefusedAreFittedAtTheirOptimum() {
        // (index, ellipticity, theta in degrees). The index only places the
        // centre off the pixel grid, deterministically.
        for (index, eps, thetaDegrees) in [(1, 0.0, 0.0), (4, 0.01, 20.0), (37, 0.03, 112.5), (43, 0.03, 160.0)] {
            let name = "ring \(index), eps \(eps), theta \(thetaDegrees)"
            let centerQX = Double(size) / 2 + 0.37 + 2.5 * sin(1.3 * Double(index) + 0.2)
            let centerQY = Double(size) / 2 - 0.41 + 2.5 * cos(1.9 * Double(index) + 0.7)
            let b = 2 * radius / (2 + eps), a = b * (1 + eps)
            let theta = thetaDegrees * .pi / 180
            let pattern = ring(centerQX: centerQX, centerQY: centerQY, a: a, b: b, theta: theta)
            let fit: EllipseCalibrationFit
            do {
                // A click's error from the plant, the annulus (b - 20, a + 20).
                fit = try EllipseCalibration.fit1D(
                    pattern: pattern, centerQX: centerQX + 1.3, centerQY: centerQY - 0.9,
                    innerRadius: b - 20, outerRadius: a + 20)
            } catch {
                Issue.record("\(name) was refused: \(error)")
                continue
            }
            #expect(abs(fit.centerQX - centerQX) <= 1e-6, "\(name)")
            #expect(abs(fit.centerQY - centerQY) <= 1e-6, "\(name)")
            #expect(abs((fit.a + fit.b) / 2 - radius * (1 + 1.75 * pow(sigma / radius, 2))) <= 0.02, "\(name)")
            #expect(abs(fit.a / fit.b - (1 + eps)) <= 1e-8, "\(name)")
            if eps > 0 {
                let error = abs(fit.theta - theta).truncatingRemainder(dividingBy: .pi)
                #expect(min(error, .pi - error) < 1e-6, "\(name)")
            }
        }
    }
}
