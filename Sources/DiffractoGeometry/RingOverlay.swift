//
//  RingOverlay.swift
//  Ported from macTEM `TEMCore/Diffraction/RingOverlay.swift` (the app's source is the oracle; pins in P22). The ring geometry on the image: a fitted ring as an ellipse in IMAGE
//        coordinates (x = column, y = row), the polyline the room draws, the predicted ring of another family scaled from the
//        calibrating ring, and the distance of a point to an ellipse (the nearest-ring hit, and the harness's own check).
//
//  The kit's `EllipseCalibrationFit` holds its centre in py4DSTEM detector order, `centerQX` the ROW and `centerQY` the COLUMN;
//  `theta` is in [0, pi), measured from the row axis towards the column axis, the major axis `a` along (cos theta, sin theta)
//  in (row, column). This file is the one place that turns that into (x, y): an ellipse point is
//      x = column = cx + a cos(phi) sin(theta) + b sin(phi) cos(theta)
//      y = row    = cy + a cos(phi) cos(theta) - b sin(phi) sin(theta).
//  No state; no kit symbol beyond the fit's fields.
//

import Foundation

/// A ring ellipse in image coordinates (pixels, the image's top-left pixel centre at (0, 0)).
public struct RingEllipse: Sendable, Equatable {
    /// The centre's column.
    public var x: Double
    /// The centre's row.
    public var y: Double
    /// The semi-axes (pixels) and the angle of the first from the row axis towards the column axis (radians).
    public var a: Double
    public var b: Double
    public var theta: Double

    public init(x: Double, y: Double, a: Double, b: Double, theta: Double) {
        self.x = x; self.y = y; self.a = a; self.b = b; self.theta = theta
    }

    /// The mean radius (a + b)/2, the quantity the calibration is made from.
    public var radius: Double { (a + b) / 2 }
}

/// Overlay geometry of a fitted ring: image-coordinate ellipse, drawing polyline, predicted rings, point-to-ellipse distance.
public enum RingOverlay {

    /// The points the room draws per ring.
    public static let polylinePoints = 180

    /// The fit as an image-coordinate ellipse: x is the kit's `centerQY` (column), y its `centerQX` (row).
    public static func ellipse(of fit: EllipseCalibrationFit) -> RingEllipse {
        RingEllipse(x: fit.centerQY, y: fit.centerQX, a: fit.a, b: fit.b, theta: fit.theta)
    }

    /// `count` points on the ellipse (x = column, y = row), the first at phi = 0 (the end of the first semi-axis); the polygon is
    /// not closed (the last point is not the first).
    public static func polyline(_ ring: RingEllipse, count: Int = polylinePoints) -> [SIMD2<Double>] {
        guard count > 0 else { return [] }
        let st = sin(ring.theta), ct = cos(ring.theta)
        return (0..<count).map { i in
            let phi = 2 * Double.pi * Double(i) / Double(count)
            let u = ring.a * cos(phi), v = ring.b * sin(phi)
            return SIMD2(ring.x + u * st + v * ct, ring.y + u * ct - v * st)
        }
    }

    /// The ellipse of a ring of another family, predicted from the calibrating ring: the same centre, b/a and angle, both
    /// semi-axes scaled by `radius / ring.radius`. nil unless every number is finite and positive.
    public static func predicted(from ring: RingEllipse, radius: Double) -> RingEllipse? {
        guard radius.isFinite, radius > 0, ring.radius.isFinite, ring.radius > 0 else { return nil }
        let scale = radius / ring.radius
        return RingEllipse(x: ring.x, y: ring.y, a: ring.a * scale, b: ring.b * scale, theta: ring.theta)
    }

    /// The orthogonal distance (pixels) of the point (x = column, y = row) to the ellipse: Newton on the nearest-point
    /// condition from the best of 3600 samples, 60 iterations at most. NaN for a non-finite number or a non-positive axis.
    public static func distance(fromX px: Double, y py: Double, to ring: RingEllipse) -> Double {
        guard px.isFinite, py.isFinite, ring.x.isFinite, ring.y.isFinite, ring.theta.isFinite,
              ring.a.isFinite, ring.a > 0, ring.b.isFinite, ring.b > 0 else { return .nan }
        let a = ring.a, b = ring.b
        let dr = py - ring.y, dc = px - ring.x
        let ct = cos(ring.theta), st = sin(ring.theta)
        let u = dr * ct + dc * st, v = -dr * st + dc * ct
        let samples = 3600
        var t = 0.0, best = Double.infinity
        for i in 0..<samples {
            let phi = 2 * Double.pi * Double(i) / Double(samples)
            let ex = a * cos(phi) - u, ey = b * sin(phi) - v
            let d2 = ex * ex + ey * ey
            if d2 < best { best = d2; t = phi }
        }
        for _ in 0..<60 {
            let x = a * cos(t), y = b * sin(t)
            let dx = -a * sin(t), dy = b * cos(t)
            let g1 = (x - u) * dx + (y - v) * dy
            var g2 = dx * dx + dy * dy + (x - u) * (-x) + (y - v) * (-y)
            if g2 <= 0 { g2 = dx * dx + dy * dy }
            let step = g1 / g2
            t -= step
            if abs(step) < 1e-15 { break }
        }
        return hypot(a * cos(t) - u, b * sin(t) - v)
    }
}
