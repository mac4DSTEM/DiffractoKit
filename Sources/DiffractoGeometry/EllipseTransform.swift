//
//  EllipseTransform.swift
//  Role: the detector-distortion correction a fitted ellipse (a, b, theta) defines, applied to a detector offset, and its
//        exact inverse.
//
//  SEED: reference commit 7d6fdbd1 `Core/Data/Calibration.swift`, `hasEllipse`, `acceptedEllipse`, `ellipseTransform`,
//  `ellipseCorrectedOffset` and `ellipseUncorrectedOffset` only (about 50 lines; none of that type's provenance machinery).
//  The formulas are py4DSTEM `BraggVectors.cal._transform(..., ellipse=True)`'s matrix in the native (qx = row, qy = column)
//  frame:
//    e = b / a,  s = sin(theta - pi/2),  c = cos(theta - pi/2)
//    t00 = e s^2 + c^2,  t01 = s c (1 - e),  t11 = s^2 + e c^2
//    qx' = t00 dy + t01 dx,  qy' = t01 dy + t11 dx,  returned as (x: qy', y: qx')
//  The axis swap is the seed's: py4DSTEM's [qx, qy] is an offset's [dy, dx] (dx along the columns, dy along the rows).
//  `a >= b` is NOT required: py4DSTEM's `_transform` uses e = b/a with no ordering, so an a < b ellipse is exactly corrected.
//  The one predicate (`accepted`): theta finite, both semi-axes finite and positive; anything else corrects nothing
//  (identity), because e = b/a <= 0 would otherwise be applied to every vector.
//

import Foundation

public struct EllipseTransform: Sendable, Equatable {
    public let a: Double
    public let b: Double
    public let theta: Double
    let t00: Double
    let t01: Double
    let t11: Double

    /// The ellipse triple when it is one, nil otherwise: finite, a > 0, b > 0 (the seed's `acceptedEllipse`).
    public static func accepted(a: Double?, b: Double?, theta: Double?) -> (a: Double, b: Double, theta: Double)? {
        guard let a, let b, let theta,
              a.isFinite, b.isFinite, theta.isFinite, a > 0, b > 0 else { return nil }
        return (a, b, theta)
    }

    /// nil unless (a, b, theta) is an ellipse by `accepted`.
    public init?(a: Double?, b: Double?, theta: Double?) {
        guard let ellipse = Self.accepted(a: a, b: b, theta: theta) else { return nil }
        self.a = ellipse.a
        self.b = ellipse.b
        self.theta = ellipse.theta
        let e = ellipse.b / ellipse.a
        let sine = sin(ellipse.theta - .pi / 2)
        let cosine = cos(ellipse.theta - .pi / 2)
        t00 = e * sine * sine + cosine * cosine
        t01 = sine * cosine * (1 - e)
        t11 = sine * sine + e * cosine * cosine
    }

    /// py4DSTEM `_transform` applied to a detector offset (dx along the columns, dy along the rows), as the seed's
    /// `ellipseCorrectedOffset`.
    public func correctedOffset(dx: Float, dy: Float) -> (x: Float, y: Float) {
        let qx = t00 * Double(dy) + t01 * Double(dx)
        let qy = t01 * Double(dy) + t11 * Double(dx)
        return (Float(qy), Float(qx))
    }

    /// The exact inverse of `correctedOffset`: a calibrated offset back to the raw (distorted) detector frame. The seed's
    /// `ellipseUncorrectedOffset`; the offset is returned unchanged when the matrix is singular.
    public func uncorrectedOffset(dx: Float, dy: Float) -> (x: Float, y: Float) {
        let det = t00 * t11 - t01 * t01
        guard abs(det) > .leastNormalMagnitude else { return (dx, dy) }
        // The forward map returned (x: qy, y: qx); undo the same axis bookkeeping.
        let qx = Double(dy)
        let qy = Double(dx)
        let rawDY = ( t11 * qx - t01 * qy) / det
        let rawDX = (-t01 * qx + t00 * qy) / det
        return (Float(rawDX), Float(rawDY))
    }

    /// `correctedOffset` where (a, b, theta) may not be an ellipse: the offset unchanged then (the seed's behaviour with
    /// no accepted ellipse).
    public static func correctedOffset(dx: Float, dy: Float, a: Double?, b: Double?, theta: Double?) -> (x: Float, y: Float) {
        guard let transform = EllipseTransform(a: a, b: b, theta: theta) else { return (dx, dy) }
        return transform.correctedOffset(dx: dx, dy: dy)
    }
}
