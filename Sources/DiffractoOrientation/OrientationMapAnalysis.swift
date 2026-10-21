//
//  OrientationMapAnalysis.swift
//  Role: the numbers a scientist reads off a scan's `[OrientationResult]` (row-major, width × height): the reliability gate and the
//        matched fraction, the Euler read-out text, the per-position scalar planes (reliability, score, in-plane angle, Euler
//        angles, symmetry disorientation) and the IPF-Z colour bytes with the low-reliability mask. Moved from mac4DSTEM's
//        `OrientationMap` (its container and presentation types stay in the app); same arithmetic, on plain arrays.
//        A position is MATCHED when `templateIndex >= 0`; an unmatched one has Euler (0, 0, 0) and is never given a number.
//

import Foundation
import simd

public enum OrientationMapAnalysis {

    /// The Euler angles at (x, y) as "φ1°, Φ°, φ2°" (one decimal), or nil outside the map or at an unmatched position.
    public static func eulerText(results: [OrientationResult], width: Int, height: Int, x: Int, y: Int) -> String? {
        guard width > 0, height > 0, results.count == width * height,
              x >= 0, x < width, y >= 0, y < height else { return nil }
        let result = results[y * width + x]
        guard result.templateIndex >= 0 else { return nil }
        let d = result.euler.degrees
        return String(format: "%.1f°, %.1f°, %.1f°", d.0, d.1, d.2)
    }

    /// The reliability below which `fraction` (0 to 1, clamped) of the MATCHED positions fall: the sorted matched values at
    /// index round((n − 1) · fraction). Nil when nothing matched.
    public static func reliabilityThreshold(results: [OrientationResult], percentile fraction: Float) -> Float? {
        let matched = results.filter { $0.templateIndex >= 0 }.map(\.reliability).sorted()
        guard !matched.isEmpty else { return nil }
        let index = Int((Float(matched.count - 1) * max(0, min(1, fraction))).rounded())
        return matched[index]
    }

    /// The fraction of matched positions with reliability at or above `threshold`; nil when nothing matched.
    public static func fractionOfMatched(results: [OrientationResult], reliabilityAtLeast threshold: Float) -> Float? {
        let matched = results.filter { $0.templateIndex >= 0 }
        guard !matched.isEmpty else { return nil }
        return Float(matched.filter { $0.reliability >= threshold }.count) / Float(matched.count)
    }

    public static func reliabilityPlane(_ results: [OrientationResult]) -> [Float] { results.map(\.reliability) }
    public static func scorePlane(_ results: [OrientationResult]) -> [Float] { results.map(\.score) }
    public static func symmetryDisorientationPlane(_ results: [OrientationResult]) -> [Float] {
        results.map(\.symmetryDisorientationRad)
    }

    /// In-plane rotation angle over 2π (radians / 2π), for a cyclic colour map.
    public static func inPlaneAnglePlane(_ results: [OrientationResult]) -> [Float] {
        let twoPi = Float(2 * Double.pi)
        return results.map { $0.inPlaneAngle / twoPi }
    }

    /// Which positions the conjugated (mirror) pass won: 1 for a mirrored match, 0 otherwise. The in-plane angle of a mirrored
    /// match carries py4DSTEM's extra half turn (`inPlaneAnglePlane` publishes it raw), so a map of angles is read together
    /// with this plane.
    public static func mirroredPlane(_ results: [OrientationResult]) -> [UInt8] { results.map { $0.mirrored ? 1 : 0 } }

    /// The in-plane angle over 2π with a mirrored match's added half turn taken off, i.e. the angle of the pattern itself
    /// (the same plane as `inPlaneAnglePlane` for a forward match), wrapped to [0, 1).
    public static func patternInPlaneAnglePlane(_ results: [OrientationResult]) -> [Float] {
        let twoPi = Float(2 * Double.pi)
        return results.map { r in
            var a = r.inPlaneAngle / twoPi - (r.mirrored ? 0.5 : 0)
            if a < 0 { a += 1 }
            return a
        }
    }

    /// φ1 wrapped to [0, 1) (radians / 2π).
    public static func phi1Plane(_ results: [OrientationResult]) -> [Float] { cyclic(results) { $0.euler.phi1 } }
    /// φ2 wrapped to [0, 1) (radians / 2π).
    public static func phi2Plane(_ results: [OrientationResult]) -> [Float] { cyclic(results) { $0.euler.phi2 } }
    /// The middle Euler angle Φ over π, clamped to [0, 1].
    public static func PhiPlane(_ results: [OrientationResult]) -> [Float] {
        results.map { max(0, min(1, $0.euler.Phi / Float.pi)) }
    }

    /// IPF-Z colours as RGBA bytes (4 per position, alpha 255). The pole is the sample-Z direction in the crystal frame, the
    /// third column of `euler.py4DSTEMOrientationMatrix`, coloured by `color`. A matched position with reliability below
    /// `threshold` is neutral grey 0.38 (a low-confidence match must not read as texture); an unmatched one is black.
    public static func ipfZRGBA(results: [OrientationResult], maskingReliabilityBelow threshold: Float?,
                                color: (SIMD3<Double>) -> SIMD3<Float>) -> [UInt8] {
        var rgba = [UInt8]()
        rgba.reserveCapacity(results.count * 4)
        for result in results {
            let matched = result.templateIndex >= 0
            let masked = threshold.map { matched && result.reliability < $0 } ?? false
            let rgb: SIMD3<Float> = masked
                ? SIMD3<Float>(repeating: 0.38)
                : matched ? color(result.euler.py4DSTEMOrientationMatrix.columns.2) : SIMD3<Float>(repeating: 0)
            rgba.append(UInt8((max(0, min(1, rgb.x)) * 255).rounded()))
            rgba.append(UInt8((max(0, min(1, rgb.y)) * 255).rounded()))
            rgba.append(UInt8((max(0, min(1, rgb.z)) * 255).rounded()))
            rgba.append(255)
        }
        return rgba
    }

    private static func cyclic(_ results: [OrientationResult], _ angle: (OrientationResult) -> Float) -> [Float] {
        let twoPi = 2 * Float.pi
        return results.map {
            var value = angle($0).truncatingRemainder(dividingBy: twoPi)
            if value < 0 { value += twoPi }
            return value / twoPi
        }
    }
}
