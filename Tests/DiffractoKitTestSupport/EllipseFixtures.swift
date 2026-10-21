//
//  EllipseFixtures.swift
//  Role: the planted elliptical ring both the unit tests and the parity tests render.
//

import Foundation
import DiffractoGeometry

/// I = exp(-((rho - 1) R)^2 / 2 s^2), rho the conic's elliptical radius in py4DSTEM's (A, B, C) convention with qx the row,
/// s = 3 px, R = (a + b) / 2, a / b = 1 + eps, a along theta; optional 60 degree gap starting at (47 seed) mod 360 degrees;
/// optional background tail * (R / max(r, 1))^2 with r the distance from the plant's centre (added after the gap); additive
/// Gaussian noise (SplitMix64 seeded seed + 1, Box-Muller pairs over even and odd linear index), Float pixels.
package struct EllipsePlant {
    package let size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, noise: Double, seed: Int
    /// Amplitude of the 1/r^2 background, 0 for none.
    package var tail = 0.0

    package init(size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, noise: Double, seed: Int,
                 tail: Double = 0.0) {
        self.size = size; self.radius = radius; self.eps = eps; self.thetaDegrees = thetaDegrees
        self.gap = gap; self.noise = noise; self.seed = seed; self.tail = tail
    }
    package var b: Double { 2 * radius / (2 + eps) }
    package var a: Double { b * (1 + eps) }
    package var centerQX: Double { Double(size) / 2 + 0.37 + 2.5 * sin(1.3 * Double(seed) + 0.2) }
    package var centerQY: Double { Double(size) / 2 - 0.41 + 2.5 * cos(1.9 * Double(seed) + 0.7) }

    package func render() -> DiffractionPattern {
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
