//
//  RingCalibration+Fit.swift
//  Role: the fit-driven calibration: a ring of known d fitted in a pattern, the fit's numbers checked against every
//        registered refusal, and Δq and the camera length computed from the fitted radius only when nothing refuses them.
//        Slice 3b of DiffractoKit (docs/registration-slice-3.md section 2).
//
//  SEED: macTEM 9482831 `Sources/TEMCore/Analysis/RingCalibration+Fit.swift`, logic unchanged. The one difference is the
//  fit it runs: macTEM's own `EllipseCalibration` copy there is the pre-ADR-069 file; this runs the kit's
//  `EllipseCalibration` (mac4DSTEM 7d6fdbd1, with ADR 069's two starts and the 0.30 px ring-width floor), so the numbers
//  macTEM registered on its copy move when it imports the kit (stated in the consumers note).
//
//  THE PATH: `EllipseCalibration.fitBestAvailable` (the conic `fit1D` seeds py4DSTEM's 11-parameter asymmetric-Gaussian ring
//  and falls back to the conic when the profile does not fit), with the coverage bound held (`acceptSparseCoverage:
//  false`); its result becomes a `FitReport`; `refusals` decides; `calibrate` runs on R = (a + b)/2 of the PROFILE fit. A
//  conic fallback is refused by `refusals` (`modelNotAsymmetricGaussian`): its radius is biased high by about 1.75 (s/R)²,
//  so it is never offered.
//
//  WHAT IS A REFUSAL AND WHAT IS AN ERROR. A ring the fit finds but that must not be calibrated from is a refusal, listed
//  by name in the result. Too little of the ring (fewer than 30 of 36 azimuthal bins) is one of those: the fit throws
//  `insufficientAngularCoverage` before fitting, and here it becomes the refusal of the same name, so a caller reads one
//  list. Every other fit error (no signal, no convergence, an impossible ellipse, the conic's own residual cap) means no
//  ring was fitted at all, and reaches the caller as thrown. An impossible d, V or pixel pitch is the caller's mistake
//  and throws `Refusal.invalidInput` before any fit runs.
//
//  ROW AND COLUMN: the fit's "qx" is the pattern's ROW and "qy" its COLUMN (py4DSTEM's order); a `DiffractionPattern` is
//  built as `DiffractionPattern(width: columns, height: rows, pixels:)`.
//

import Foundation

/// What a fit-driven calibration of one ring gives: the fit (nil when the ring was refused before fitting), the numbers
/// the decision was made on, every refusal that applies, and the calibration — nil when any refusal withholds Δq.
public struct RingFitCalibration: Sendable, Equatable {
    public let fit: EllipseCalibrationFit?
    public let report: RingCalibration.FitReport?
    public let refusals: [RingCalibration.Refusal]
    public let calibration: RingCalibration?

    // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
    public init(fit: EllipseCalibrationFit?, report: RingCalibration.FitReport?, refusals: [RingCalibration.Refusal],
                calibration: RingCalibration?) {
        self.fit = fit
        self.report = report
        self.refusals = refusals
        self.calibration = calibration
    }
}

extension RingCalibration {

    /// Fits the ring in the annulus (innerRadius, outerRadius) about the start centre (row `centerQX`, column `centerQY`)
    /// and calibrates from it with d (Å), V (volts; nil when not stated, never defaulted) and the detector pixel pitch
    /// (mm; nil when not stated). Throws `Refusal.invalidInput` for an impossible d, V or pitch, and the fit's own error
    /// when no ring was fitted; a ring that was fitted but must not be calibrated from comes back with its refusals.
    /// Synchronous and CPU-bound (a few hundred ms for a ring of R ≈ 100 px): a caller on the main actor runs it from a
    /// `@concurrent` function.
    public static func fit(pattern: DiffractionPattern, centerQX: Double, centerQY: Double,
                            innerRadius: Double, outerRadius: Double,
                            d: Double, highTensionVolts: Double?, pixelPitchMM: Double?) throws -> RingFitCalibration {
        guard d.isFinite, d > 0 else { throw Refusal.invalidInput("the ring spacing must be a positive number of Å") }
        if let volts = highTensionVolts, !(volts.isFinite && volts > 0) {
            throw Refusal.invalidInput("the high tension must be a positive number of volts")
        }
        if let pitch = pixelPitchMM, !(pitch.isFinite && pitch > 0) {
            throw Refusal.invalidInput("the pixel pitch must be a positive number of millimetres")
        }

        let fit: EllipseCalibrationFit
        do {
            fit = try EllipseCalibration.fitBestAvailable(
                pattern: pattern, centerQX: centerQX, centerQY: centerQY,
                innerRadius: innerRadius, outerRadius: outerRadius, acceptSparseCoverage: false)
        } catch EllipseCalibration.FitError.insufficientAngularCoverage(let bins) {
            var refused: [Refusal] = [.insufficientAngularCoverage(occupiedBins: bins)]
            if highTensionVolts == nil { refused.append(.noHighTension) }
            return RingFitCalibration(fit: nil, report: nil, refusals: refused, calibration: nil)
        }

        let report = FitReport(
            model: fit.model, semiAxisA: fit.a, semiAxisB: fit.b,
            occupiedAngularBins: fit.occupiedAngularBins, sparseCoverage: fit.sparseCoverage,
            normalizedResidual: fit.normalizedResidual,
            innerSigma: fit.profile?.innerSigma, outerSigma: fit.profile?.outerSigma,
            annulusWidth: outerRadius - innerRadius)
        let refused = refusals(report, highTensionVolts: highTensionVolts)
        guard !refused.contains(where: \.withholdsDeltaQ) else {
            return RingFitCalibration(fit: fit, report: report, refusals: refused, calibration: nil)
        }
        let calibration = try calibrate(radius: report.radius, d: d, highTensionVolts: highTensionVolts,
                                        pixelPitchMM: pixelPitchMM)
        return RingFitCalibration(fit: fit, report: report, refusals: refused, calibration: calibration)
    }
}
