//
//  RingCalibration+Q.swift
//  Role: the exact scattering vector q of a pixel radius on a flat detector, and its inverse. DiffractoKit 0.1.0,
//        added to the closed forms of `RingCalibration.swift` without touching them.
//
//  THE FORMULAS (exact; nothing here is a small-angle form):
//    θ  = ½ atan(R / L)               the Bragg angle of a ring at radius R on a flat detector perpendicular to the beam
//                                     (R·p = L·tan 2θ, so tan 2θ = R/L_px with R and L both in pixels)
//    q  = 2 sin θ / λ = 2 sin(½ atan(R/L)) / λ          Å⁻¹, no 2π; for the calibrating ring q = 1/d (Bragg's law)
//    R  = L tan(2 asin(qλ/2))         the inverse
//  The small-angle form q = R/(λL) is high by 3 (R/L)² / 8 (about 4e-5 at the 200 kV Au {111} ring, as `cameraLengthPixels`'s
//  header says of λL = Rd); these forms have no such error.
//
//  Total, as every closed form here: NaN outside the domain, never a plausible number and never an infinity.
//

import Foundation

extension RingCalibration {
    /// q = 2 sin(½ atan(R/L)) / λ in Å⁻¹ (no 2π), for a ring at radius R (pixels) on a flat detector perpendicular to the beam
    /// at the camera length L (pixels, `cameraLengthPixels`) and wavelength λ (Å). For the calibrating ring's own (R, L_px, λ)
    /// it is 1/d. The angle is `atan2(R, L)`, so R/L is never formed and cannot overflow or lose its quotient to underflow.
    /// NaN unless R, L and λ are finite and positive, and NaN where the result is not a finite positive number: an overflow (a
    /// subnormal λ) or an underflow to 0 (a positive R so small against L that atan(R/L) is 0) is never given as infinity or 0.
    /// R may be as large as you like against L: 2θ stays below 180° and q below 2/λ, so there is no other limit.
    /// VALIDITY: against 60-digit arithmetic (tools/kit-0.0.4/q-pins.py, 100 cases, V 60 to 300 kV, L 500 px to 1e5 px, R 1 px
    /// to 2000 px) the worst relative error is 8.2e-16: atan2 and sin each round once and neither amplifies (the form has no
    /// subtraction). The 1/d identity on the 320 ring-calibration pin cases holds to 5.5e-16. The trig-free form
    /// 2 sin(½ atan(R/L)) = (R/h) / √((1 + L/h)/2), h = hypot(R, L), was measured too and is no better (both 2.2e-16 in
    /// NumPy on the grid), so the atan2/sin form stays.
    public static func q(radius: Double, cameraLengthPixels: Double, wavelengthAngstroms: Double) -> Double {
        guard radius.isFinite, radius > 0, cameraLengthPixels.isFinite, cameraLengthPixels > 0,
              wavelengthAngstroms.isFinite, wavelengthAngstroms > 0 else { return .nan }
        let q = 2 * sin(atan2(radius, cameraLengthPixels) / 2) / wavelengthAngstroms
        return q.isFinite && q > 0 ? q : .nan
    }

    /// R = L tan(2 asin(qλ/2)) in pixels, the inverse of `q(radius:cameraLengthPixels:wavelengthAngstroms:)`, for q in Å⁻¹
    /// (no 2π), L in pixels and λ in Å. NaN unless q, L and λ are finite and positive and qλ/2 is below
    /// `maximumLambdaOver2d` (1/√2: 2θ below 90°, as `twoTheta`: tan 2θ is negative past it and a flat detector perpendicular
    /// to the beam shows no ring there), and NaN where the result is not a finite positive number (an overflow, or an
    /// underflow to 0), never an infinity.
    /// VALIDITY: the form is exact, but qλ/2 is a double and tan 2θ diverges at 90°. The rounding of qλ/2 (1.1e-16 relative)
    /// reaches R multiplied by 1/((1 − s²)(1 − 2s²)), s = qλ/2 and 1 − 2s² = cos 2θ, so against 60-digit arithmetic
    /// (tools/kit-0.0.4/q-pins.py, 120 cases) the worst relative error is 3.3e-15 on the 100 grid cases (V 60 to 300 kV, L
    /// 500 px to 1e5 px, R 1 px to 2000 px, 2θ up to 76°) and 4.5e-14 at 0.99 of the limit (the 20 cases at s√2 = 0.99, the
    /// factor about 100, 2θ = 88.9°), about 1e-15/cos 2θ in general. That passes the 1e-12 band for 2θ below about 89.9°;
    /// closer to the limit the number is right only to that error, and the `maximumLambdaOver2d` guard is the only refusal (one
    /// double below it still gives a finite positive R, not compared to anything).
    public static func radius(q: Double, cameraLengthPixels: Double, wavelengthAngstroms: Double) -> Double {
        guard q.isFinite, q > 0, cameraLengthPixels.isFinite, cameraLengthPixels > 0,
              wavelengthAngstroms.isFinite, wavelengthAngstroms > 0 else { return .nan }
        let sinTheta = q * (wavelengthAngstroms / 2)           // λ/2 is exact in binary: one rounding in all
        guard sinTheta < maximumLambdaOver2d else { return .nan }
        let radius = cameraLengthPixels * tan(2 * asin(sinTheta))
        // Not `> 0`: the sign is the Bragg guard's to keep (tan 2θ is negative past 90°), so each guard answers for itself
        // and a test that removes one is seen. Finite and nonzero: an overflow and an underflow to 0 are NaN.
        return radius.isFinite && radius != 0 ? radius : .nan
    }
}
