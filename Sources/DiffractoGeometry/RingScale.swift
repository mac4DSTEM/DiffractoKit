//
//  RingScale.swift
//  Ported from macTEM TEMCore/Diffraction/RingScale.swift (the app's source is the oracle; P24). |q| of a ring radius and the
//  ring radius of a |q| for a calibration, the linear Δq·R form, and the fit's scale against a file's.
//
//  With a beam energy the calibration holds a camera length and a wavelength, and the kit's closed forms apply:
//      q = 2 sin(1/2 atan(R / L_px)) / lambda        R = L_px tan(2 asin(lambda q / 2))      (Å⁻¹, no 2π)
//  (`RingCalibration.q(radius:cameraLengthPixels:wavelengthAngstroms:)` and `radius(q:...)`, which this calls: the app's
//  forms agree with them to 1e-12 relative; `q` is bit for bit). Without one the calibration has only the linear scale Δq (Å⁻¹ per pixel, 1/d over the
//  ring's radius) and q = Δq R is all that can be said. For a flat detector the two differ by 1.4e-5 relative at Au {200}
//  from {111} at 200 kV (8.7e-6 at 300 kV), growing with R. No state.
//

import Foundation

/// The reciprocal scale of a ring under a `RingCalibration`.
public enum RingScale {

    /// |q| (Å⁻¹) of a ring of radius `radius` px under `calibration`: the exact `RingCalibration.q` when the calibration has a
    /// camera length and a wavelength (a beam energy), else Δq · R. NaN unless R is finite and positive (and the calibration
    /// usable): never 0 or infinity.
    public static func q(radius: Double, calibration: RingCalibration) -> Double {
        guard radius.isFinite, radius > 0 else { return .nan }
        if let length = calibration.lengthPixels, let wavelength = calibration.wavelengthAngstroms {
            return RingCalibration.q(radius: radius, cameraLengthPixels: length, wavelengthAngstroms: wavelength)
        }
        return linear(radius: radius, deltaQ: calibration.deltaQ)
    }

    /// The ring radius (px) of |q| (Å⁻¹) under `calibration`: the inverse of `q(radius:calibration:)`. NaN unless q is finite
    /// and positive; past the detector's Bragg limit (2θ at 90°) the exact form is NaN.
    public static func radius(q: Double, calibration: RingCalibration) -> Double {
        guard q.isFinite, q > 0 else { return .nan }
        if let length = calibration.lengthPixels, let wavelength = calibration.wavelengthAngstroms {
            return RingCalibration.radius(q: q, cameraLengthPixels: length, wavelengthAngstroms: wavelength)
        }
        guard calibration.deltaQ.isFinite, calibration.deltaQ > 0 else { return .nan }
        let r = q / calibration.deltaQ
        return r.isFinite && r > 0 ? r : .nan
    }

    /// The linear form Δq · R (Å⁻¹), one multiplication; NaN unless both are finite and positive and the product is.
    public static func linear(radius: Double, deltaQ: Double) -> Double {
        guard radius.isFinite, radius > 0, deltaQ.isFinite, deltaQ > 0 else { return .nan }
        let q = deltaQ * radius
        return q.isFinite && q > 0 ? q : .nan
    }

    /// Δq_fit / Δq_file − 1: how far the fitted scale is from the file's (0.05 means the fit's is 5 % larger). NaN unless both
    /// are finite and positive.
    public static func fitVersusFile(fitDeltaQ: Double, fileDeltaQ: Double) -> Double {
        guard fitDeltaQ.isFinite, fitDeltaQ > 0, fileDeltaQ.isFinite, fileDeltaQ > 0 else { return .nan }
        return fitDeltaQ / fileDeltaQ - 1
    }
}
