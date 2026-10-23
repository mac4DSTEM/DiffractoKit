//
//  LinearScaleFit.swift
//  Role: Δq (Å⁻¹ per pixel) from SEVERAL indexed rings or spots by least squares (macTEM kit ask 2; py4DSTEM `get_dq_from_indexed_peaks` is the
//        same idea in the cubic case: q_k = Δq r_k with no intercept). `LinearScale.estimate` uses one shell; this uses every ring that is indexed,
//        and says how well they agree.
//
//  Model: g_k = Δq r_k + e_k, r_k the observed radius (pixels, after any ellipse correction), g_k the indexed |g| (Å⁻¹), through the origin.
//  Weighted least squares with weights w_k (default 1): Δq = Σ w r g / Σ w r². σ(Δq) is the standard error from the residuals,
//  sqrt(Σ w e² / (n - 1) / Σ w r²), nil with one ring. No Python reference computes this with an uncertainty; planted truth only, `validation: "none"`.
//

import Foundation

public struct LinearScaleFit: Sendable, Equatable {
    public let invAngstromPerPixel: Double
    /// Standard error of Δq; nil with fewer than two rings.
    public let sigma: Double?
    /// g_k - Δq r_k, Å⁻¹, in the order given.
    public let residualsInvAngstrom: [Double]
    /// Root mean square of the relative residuals e_k / g_k.
    public let rmsRelativeResidual: Double
    public let ringCount: Int
    public var validation: String { "none" }
}

extension LinearScale {
    /// Δq from `ringRadiiPixels` against `indexedInvAngstrom` (same length, at least one ring, every radius and |g| finite and positive).
    /// `weights` (positive, same length) default to 1. Nil when the input does not meet that.
    public static func fit(ringRadiiPixels: [Double], indexedInvAngstrom: [Double], weights: [Double]? = nil) -> LinearScaleFit? {
        let n = ringRadiiPixels.count
        guard n >= 1, indexedInvAngstrom.count == n, weights.map({ $0.count == n && $0.allSatisfy { $0.isFinite && $0 > 0 } }) ?? true,
              ringRadiiPixels.allSatisfy({ $0.isFinite && $0 > 0 }), indexedInvAngstrom.allSatisfy({ $0.isFinite && $0 > 0 }) else { return nil }
        var rg = 0.0, rr = 0.0
        for k in 0 ..< n {
            let w = weights?[k] ?? 1
            rg += w * ringRadiiPixels[k] * indexedInvAngstrom[k]
            rr += w * ringRadiiPixels[k] * ringRadiiPixels[k]
        }
        let dq = rg / rr
        let residuals = (0 ..< n).map { indexedInvAngstrom[$0] - dq * ringRadiiPixels[$0] }
        var ss = 0.0, rel = 0.0
        for k in 0 ..< n {
            ss += (weights?[k] ?? 1) * residuals[k] * residuals[k]
            rel += pow(residuals[k] / indexedInvAngstrom[k], 2)
        }
        return LinearScaleFit(invAngstromPerPixel: dq, sigma: n > 1 ? (ss / Double(n - 1) / rr).squareRoot() : nil, residualsInvAngstrom: residuals,
                              rmsRelativeResidual: (rel / Double(n)).squareRoot(), ringCount: n)
    }
}
