//
//  RingCalibration+Fit.swift
//  Role: the fit-driven calibration: a ring of known d fitted in a pattern, the fit's numbers checked against every
//        refusal, and Δq and the camera length computed from the fitted radius only when nothing refuses them.
//
//  SEED: reference commit 9482831 `Sources/TEMCore/Analysis/RingCalibration+Fit.swift`, logic unchanged. The one difference is the
//  fit it runs: the reference implementation's own `EllipseCalibration` copy there is the pre-ADR-069 file; this runs the kit's
//  `EllipseCalibration` (reference commit 7d6fdbd1, with ADR 069's two starts and the 0.30 px ring-width floor), so the numbers
//  the reference implementation measured on its copy move when it imports the kit.
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
//  THE SPLIT (DiffractoKit 0.1.0): `fitRing` fits the ring and decides the fit's
//  own refusals with no spacing, so a caller can measure a ring before it knows which family it is; `calibrating` adds d, V
//  and the pitch. `fit` is exactly `calibrating(fitRing(...), d:...)`. Since 0.1.0 the fit runs with
//  `requireRingEvidence: true`: the report carries the ring evidence and the annulus's radii, so
//  a ring outside its annulus, a profile that is not a peak, an inconsistent geometry and too few sectors holding the ring
//  are refused by name, and one ring sigma at the 0.30 px floor is kept when the evidence passes.
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
    /// Exactly `calibrating(fitRing(...), d:highTensionVolts:pixelPitchMM:)`, with d, V and the pitch checked first.
    /// Two policies act: `policy` (the calibration's refusal bounds, with the evidence bounds) and `ellipsePolicy` (the
    /// ellipse fit's own rules); see `fitRing` for how they combine.
    /// Synchronous and CPU-bound (0.1 s at R ≈ 100 px; tens of seconds for a band of ±90 px at R 600 px on a 2048² pattern):
    /// a caller on the main actor runs it from a `@concurrent` function, and may pass `progress` (0 to 1) and cancel the
    /// task or pass `isCancelled` (it then throws `CancellationError`).
    public static func fit(pattern: DiffractionPattern, centerQX: Double, centerQY: Double,
                            innerRadius: Double, outerRadius: Double,
                            d: Double, highTensionVolts: Double?, pixelPitchMM: Double?,
                            policy: Policy = .default, ellipsePolicy: EllipseCalibration.Policy = .default,
                            progress: (@Sendable (Double) -> Void)? = nil,
                            isCancelled: (@Sendable () -> Bool)? = nil) throws -> RingFitCalibration {
        try validate(d: d, highTensionVolts: highTensionVolts, pixelPitchMM: pixelPitchMM)
        let ring = try fitRing(pattern: pattern, centerQX: centerQX, centerQY: centerQY, innerRadius: innerRadius,
                               outerRadius: outerRadius, policy: policy, ellipsePolicy: ellipsePolicy,
                               progress: progress, isCancelled: isCancelled)
        return try calibrating(ring, d: d, highTensionVolts: highTensionVolts, pixelPitchMM: pixelPitchMM)
    }

    /// Fits the ring in the annulus about the start centre and decides the fit's own refusals, without a spacing: the
    /// radius of any ring, before its family is known (the reference implementation's ask G1). Throws the fit's own error when no ring was
    /// fitted (and `CancellationError` when cancelled); a ring refused for coverage before fitting comes back with nil
    /// fit and report.
    ///
    /// **Two policies, one evidence.** `policy.evidence` REPLACES `ellipsePolicy.evidence`: the evidence is measured and
    /// judged by the calibration policy's bounds, whatever the ellipse policy holds there (also when it holds an invalid
    /// one). Every other field of `ellipsePolicy` is checked: an invalid one throws `EllipseCalibration.FitError.invalidInput`
    /// (and an invalid `policy`, `Refusal.invalidInput`). The coverage rule is applied twice, once by each policy, so
    /// loosening it through this call needs both: the ellipse fit refuses a ring before fitting it when fewer than
    /// `ellipsePolicy.degeneracyBoundBins` bins are occupied (the refusal `insufficientAngularCoverage` with no fit),
    /// and the calibration then refuses a fitted ring that holds fewer than `policy.minimumOccupiedAngularBins` bins or
    /// whose evidence holds fewer than `policy.evidence.minimumOccupiedBins` sectors. Lowering only the calibration
    /// policy's bounds changes nothing for a ring the ellipse fit refuses; lowering only the ellipse policy's lets the
    /// fit run and the calibration's bounds refuse it.
    public static func fitRing(pattern: DiffractionPattern, centerQX: Double, centerQY: Double,
                                innerRadius: Double, outerRadius: Double,
                                policy: Policy = .default, ellipsePolicy: EllipseCalibration.Policy = .default,
                                progress: (@Sendable (Double) -> Void)? = nil,
                                isCancelled: (@Sendable () -> Bool)? = nil) throws -> RingRadiusFit {
        guard policy.isValid else { throw Refusal.invalidInput("the calibration policy has a field outside its domain") }
        // One evidence policy: the ring evidence is measured and judged by the calibration policy's.
        var fitPolicy = ellipsePolicy
        fitPolicy.evidence = policy.evidence
        let fit: EllipseCalibrationFit
        do {
            fit = try EllipseCalibration.fitBestAvailable(
                pattern: pattern, centerQX: centerQX, centerQY: centerQY,
                innerRadius: innerRadius, outerRadius: outerRadius, acceptSparseCoverage: false,
                requireRingEvidence: true, policy: fitPolicy, progress: progress, isCancelled: isCancelled)
        } catch EllipseCalibration.FitError.insufficientAngularCoverage(let bins) {
            return RingRadiusFit(fit: nil, report: nil, refusals: [.insufficientAngularCoverage(occupiedBins: bins)])
        }
        let report = FitReport(
            model: fit.model, semiAxisA: fit.a, semiAxisB: fit.b,
            occupiedAngularBins: fit.occupiedAngularBins, sparseCoverage: fit.sparseCoverage,
            normalizedResidual: fit.normalizedResidual,
            innerSigma: fit.profile?.innerSigma, outerSigma: fit.profile?.outerSigma,
            annulusWidth: outerRadius - innerRadius,
            annulusInnerRadius: innerRadius, annulusOuterRadius: outerRadius, ringEvidence: fit.ringEvidence)
        return RingRadiusFit(fit: fit, report: report, refusals: fitRefusals(report, policy: policy))
    }

    /// The calibration of a fitted ring with d (Å), V (volts; nil when not stated) and the pixel pitch (mm; nil when not
    /// stated): the ring's refusals, `.noHighTension` when V is nil, and Δq and L unless a refusal withholds Δq. Throws
    /// `Refusal.invalidInput` for an impossible d, V or pitch.
    public static func calibrating(_ ring: RingRadiusFit, d: Double, highTensionVolts: Double?,
                                   pixelPitchMM: Double?) throws -> RingFitCalibration {
        try validate(d: d, highTensionVolts: highTensionVolts, pixelPitchMM: pixelPitchMM)
        var refused = ring.refusals
        let invalid = refused.first?.name == "invalidInput"
        if highTensionVolts == nil, !invalid { refused.append(.noHighTension) }
        guard let report = ring.report, !refused.contains(where: \.withholdsDeltaQ) else {
            return RingFitCalibration(fit: ring.fit, report: ring.report, refusals: refused, calibration: nil)
        }
        let calibration = try calibrate(radius: report.radius, d: d, highTensionVolts: highTensionVolts,
                                        pixelPitchMM: pixelPitchMM)
        return RingFitCalibration(fit: ring.fit, report: report, refusals: refused, calibration: calibration)
    }

    static func validate(d: Double, highTensionVolts: Double?, pixelPitchMM: Double?) throws {
        guard d.isFinite, d > 0 else { throw Refusal.invalidInput("the ring spacing must be a positive number of Å") }
        if let volts = highTensionVolts, !(volts.isFinite && volts > 0) {
            throw Refusal.invalidInput("the high tension must be a positive number of volts")
        }
        if let pitch = pixelPitchMM, !(pitch.isFinite && pitch > 0) {
            throw Refusal.invalidInput("the pixel pitch must be a positive number of millimetres")
        }
    }
}

/// A ring fitted without a spacing (`RingCalibration.fitRing`): the fit (nil when the ring was refused before fitting),
/// the numbers its refusals were decided on, and the fit's own refusals (never `.noHighTension`). `radius` is the ring's
/// (a + b)/2 in pixels when a fit exists, whether or not a refusal applies: read it beside `refusals`.
public struct RingRadiusFit: Sendable, Equatable {
    public let fit: EllipseCalibrationFit?
    public let report: RingCalibration.FitReport?
    public let refusals: [RingCalibration.Refusal]

    public init(fit: EllipseCalibrationFit?, report: RingCalibration.FitReport?, refusals: [RingCalibration.Refusal]) {
        self.fit = fit
        self.report = report
        self.refusals = refusals
    }

    public var radius: Double? { report?.radius }
    /// Whether the radius may be used as a measurement (no refusal withholds Δq).
    public var isUsable: Bool { report != nil && !refusals.contains(where: \.withholdsDeltaQ) }
}
