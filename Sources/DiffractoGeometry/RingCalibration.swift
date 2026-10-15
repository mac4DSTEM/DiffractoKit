//
//  RingCalibration.swift
//  Role: the calibration a diffraction ring of known d gives (Δq in 1/Å per pixel, and the camera length) as closed forms,
//        and the decision of when a ring is NOT to be calibrated from. Slice 3a of DiffractoKit
//        (docs/registration-slice-3.md); ported from macTEM 9482831 `RingCalibration.swift` with every formula, constant,
//        NaN rule, refusal order and message kept. The fit that fills a `FitReport` is slice 3b, not here.
//
//  THE FORMULAS (exact; nothing here is a small-angle form):
//    λ     = hc / √(eV (2 m₀c² + eV))        the relativistic electron wavelength, Å, V in volts (`ElectronWavelength`)
//    2θ    = 2 arcsin(λ / 2d)                Bragg's law, d in Å
//    L_px  = R / tan 2θ                      the camera length in pixels, R the ring's radius in pixels (a flat detector
//                                            perpendicular to the beam: R·p = L·tan 2θ, so L_px = L/p = R / tan 2θ)
//    Δq    = (1/d) / R                       Å⁻¹ per pixel; q = 1/d (no 2π), needs no wavelength
//    L_mm  = p · L_px                        only where the detector's pixel pitch p (mm) is stated
//  The small-angle form λL = Rd is off by 4.3e-5 at 200 kV Au {111} and 9.4e-4 at 80 kV Au {422}.
//  L is the EFFECTIVE camera length, never the file's indicated one. L ∝ 1/λ moves 0.58 % per 1 % of V at 200 kV and a
//  stored V is nominal, so `calibrate` returns the V it used beside L, and Δq (which needs no V) leads.
//
//  The closed forms are total: outside their domain they return NaN, never a plausible number. `calibrate` and
//  `refusals` check the domain first and say why, by name.
//
//  TWO REFUSALS beyond the six of the original registration (model, coverage, sigma, residual, ellipticity, high tension):
//  `.braggLimit` refuses λ/2d ≥ 1/√2, i.e. 2θ ≥ 90° (2θ ≥ 90° makes tan 2θ negative, and a flat detector perpendicular
//  to the beam shows no ring beyond 90° anyway), and `.invalidInput` refuses a number that is not finite and positive
//  (a NaN passes every comparison).
//

import Foundation

/// The model a ring fit used: the conic is a fallback whose radius is biased, the asymmetric-Gaussian ring is the profile fit.
public enum EllipseCalibrationModel: String, Sendable, Equatable {
    case conic = "Intensity-weighted conic"
    case asymmetricGaussian = "Asymmetric Gaussian ring"
}

public struct RingCalibration: Sendable, Equatable {
    /// The ring's spacing (Å) and radius ((a + b)/2 of the fitted ellipse, pixels), as the caller gave them.
    public let d: Double
    public let radius: Double
    /// Δq = (1/d)/R, Å⁻¹ per pixel. Always given once the ring is: it needs no wavelength.
    public let deltaQ: Double
    /// The high tension (volts) λ and L were computed at, shown beside L. nil when none was stated.
    public let highTensionVolts: Double?
    public let wavelengthAngstroms: Double?
    /// 2θ in radians.
    public let twoTheta: Double?
    public let lengthPixels: Double?
    /// p · L_px; nil without a pixel pitch (there is no field for one in `AcquisitionMetadata` yet).
    public let lengthMillimetres: Double?
    /// Why L is absent (`.noHighTension` or `.braggLimit`); Δq is given either way.
    public let lengthRefusal: Refusal?

    // MARK: Refusals

    /// Why a ring is not calibrated from: the six registered refusals, and the two additions `.braggLimit` and
    /// `.invalidInput` (see the header). Every case but `.noHighTension` and `.braggLimit` withholds Δq as well; those
    /// two withhold only L (Δq needs no wavelength).
    public enum Refusal: Error, Equatable, Sendable, LocalizedError {
        /// The fit fell back to the conic: its radius is not the profile's peak.
        case modelNotAsymmetricGaussian
        /// Fewer than 30 of 36 azimuthal bins hold ring signal, or the fit was accepted only on the caller's assertion
        /// of one ring in a sparse annulus (`sparseCoverage`). With ≥ 30 bins the payload means the assertion.
        case insufficientAngularCoverage(occupiedBins: Int)
        /// A ring sigma (the larger of the inner and outer) at or above the annulus width: the ring is not resolved
        /// inside the annulus, and the fit's peak can read a percent off at a small residual.
        case sigmaAtAnnulusWidth(sigma: Double, annulusWidth: Double)
        case excessiveResidual(Double)
        /// a/b − 1 above `maximumEllipticity`: (a + b)/2 and √(ab) would differ by 0.1 %.
        case excessiveEllipticity(Double)
        case noHighTension
        /// λ/2d at or above `maximumLambdaOver2d`: no ring of this spacing reaches a flat detector at this wavelength.
        case braggLimit(lambdaOver2d: Double)
        /// A number that is not finite and positive (or a count outside its range): the caller's mistake, not data.
        case invalidInput(String)
        /// The fitted radius (a + b)/2 lies outside the annulus the fit was given: the ring term fitted something other
        /// than the ring the caller drew (DiffractoKit 0.0.4, registration-0.0.4 §1.4; slice 3b's F2, a central blob).
        case ringOutsideAnnulus(radius: Double, inner: Double, outer: Double)
        /// The data's profile about the fit does not fall by a quarter of the ring's amplitude on both sides: an edge or a
        /// slope, not a ring (registration-0.0.4 §1.3; slice 3b's F2, a disk edge). nil where a side had no pixels.
        case notARingPeak(innerDrop: Double?, outerDrop: Double?)
        /// The sectors put the ring off the fitted centre or ellipse by more than 0.25 px and 4 standard errors: the fit
        /// was pulled by where the intensity is, not by where the ring is (registration-0.0.4 §1.3; slice 3b's F1 gap
        /// noise, F2 texture).
        case ringGeometryInconsistent(centreShift: Double, ellipseShift: Double)
        /// The annulus holds a second ring: the pooled profile rises again beyond the fitted ring by `height` of its
        /// amplitude (registration-0.0.4 §8, after Gate B: a neighbouring ring pulled the fit onto itself or between the
        /// two). Narrow the annulus to one ring.
        case secondRingInAnnulus(height: Double)

        /// The case's name, without its numbers: a stable key for logs and tools/ring-calibration-parity.
        public var name: String {
            switch self {
            case .modelNotAsymmetricGaussian: return "modelNotAsymmetricGaussian"
            case .insufficientAngularCoverage: return "insufficientAngularCoverage"
            case .sigmaAtAnnulusWidth: return "sigmaAtAnnulusWidth"
            case .excessiveResidual: return "excessiveResidual"
            case .excessiveEllipticity: return "excessiveEllipticity"
            case .noHighTension: return "noHighTension"
            case .braggLimit: return "braggLimit"
            case .invalidInput: return "invalidInput"
            case .ringOutsideAnnulus: return "ringOutsideAnnulus"
            case .notARingPeak: return "notARingPeak"
            case .ringGeometryInconsistent: return "ringGeometryInconsistent"
            case .secondRingInAnnulus: return "secondRingInAnnulus"
            }
        }

        /// Whether Δq is withheld too. False only where the missing or impossible thing is the wavelength.
        public var withholdsDeltaQ: Bool {
            switch self {
            case .noHighTension, .braggLimit: return false
            default: return true
            }
        }

        /// One sentence, for the person reading the result.
        public var reason: String {
            switch self {
            case .modelNotAsymmetricGaussian:
                return "The fit fell back to the conic, whose radius is biased high by about 1.75 (s/R)² and is not offered; only the asymmetric-Gaussian ring's peak is."
            case .insufficientAngularCoverage(let bins):
                if bins < RingCalibration.minimumOccupiedAngularBins {
                    return "The ring signal covers only \(bins) of \(RingCalibration.angularBinCount) azimuthal bins; at least \(RingCalibration.minimumOccupiedAngularBins) are needed."
                }
                return "The fit was accepted only because the caller asserted one ring in a sparse annulus, which a calibration does not accept."
            case .sigmaAtAnnulusWidth(let sigma, let width):
                return String(format: "A ring sigma of %.1f px is not less than the annulus width of %.1f px, so the ring is not resolved inside the annulus and its peak is not trusted.", sigma, width)
            case .excessiveResidual(let residual):
                return String(format: "The fit residual is %.3f, above the %.2f cap, so the model does not describe this ring.", residual, RingCalibration.maximumResidual)
            case .excessiveEllipticity(let ellipticity):
                return String(format: "The ring's a/b − 1 is %.3f, above %.3f, where the mean of the semi-axes and their geometric mean differ by 0.1 %%, so the radius would depend on the convention.", ellipticity, RingCalibration.maximumEllipticity)
            case .noHighTension:
                return "No high tension is stated, so the wavelength and the camera length are not computed; Δq in Å⁻¹ per pixel is still given."
            case .braggLimit(let ratio):
                return String(format: "At this wavelength no ring of this spacing reaches a flat detector: λ/2d is %.4f and must stay below %.4f (2θ below 90°).", ratio, RingCalibration.maximumLambdaOver2d)
            case .invalidInput(let detail):
                return "Cannot calibrate: \(detail)."
            case .ringOutsideAnnulus(let radius, let inner, let outer):
                return String(format: "The fitted ring radius of %.1f px lies outside the annulus from %.1f to %.1f px, so the fit describes something other than the ring in the annulus.", radius, inner, outer)
            case .notARingPeak(let innerDrop, let outerDrop):
                func percent(_ drop: Double?) -> String { drop.map { String(format: "%.0f %%", 100 * $0) } ?? "nothing" }
                return "The intensity about the fitted ring falls by \(percent(innerDrop)) of the ring's amplitude inside it and by \(percent(outerDrop)) outside it; a ring falls by at least \(Int(100 * RingEvidence.minimumDrop)) % on both sides, so this is an edge or a slope, not a ring."
            case .ringGeometryInconsistent(let centre, let ellipse):
                return String(format: "Read sector by sector, the ring lies %.2f px off the fitted centre and %.2f px off the fitted ellipse (the bound is %.2f px), so the fit follows where the intensity is rather than where the ring is.", centre, ellipse, RingEvidence.maximumShift)
            case .secondRingInAnnulus(let height):
                return String(format: "The annulus holds a second ring (a peak of %.0f %% of the fitted ring's amplitude beyond it), so the fit may sit on either ring or between them; narrow the annulus to one ring.", 100 * height)
            }
        }

        public var errorDescription: String? { reason }
    }

    /// The numbers a ring fit reports, as plain values (slice 3b fills them from a ring fit and the annulus
    /// it was given; this file has no fit in it).
    public struct FitReport: Sendable, Equatable {
        public var model: EllipseCalibrationModel
        /// The ellipse's semi-axes, pixels, in either order.
        public var semiAxisA: Double
        public var semiAxisB: Double
        public var occupiedAngularBins: Int
        public var sparseCoverage: Bool
        public var normalizedResidual: Double
        /// The ring's inner and outer Gaussian sigmas, pixels; nil for a conic fit, which has none.
        public var innerSigma: Double?
        public var outerSigma: Double?
        /// outerRadius − innerRadius of the annulus the fit was given, pixels.
        public var annulusWidth: Double
        /// The annulus's radii (DiffractoKit 0.0.4): with them, a radius outside the annulus is refused. nil skips that.
        public var annulusInnerRadius: Double?
        public var annulusOuterRadius: Double?
        /// The fit's ring evidence (DiffractoKit 0.0.4, `EllipseCalibration.ringEvidence`): with it, too few sectors
        /// holding the ring, no peak and an inconsistent geometry are refused. nil skips those.
        public var ringEvidence: RingEvidence?

        // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
        public init(model: EllipseCalibrationModel, semiAxisA: Double, semiAxisB: Double, occupiedAngularBins: Int,
                     sparseCoverage: Bool, normalizedResidual: Double, innerSigma: Double?, outerSigma: Double?,
                     annulusWidth: Double, annulusInnerRadius: Double? = nil, annulusOuterRadius: Double? = nil,
                     ringEvidence: RingEvidence? = nil) {
            self.annulusInnerRadius = annulusInnerRadius
            self.annulusOuterRadius = annulusOuterRadius
            self.ringEvidence = ringEvidence
            self.model = model
            self.semiAxisA = semiAxisA
            self.semiAxisB = semiAxisB
            self.occupiedAngularBins = occupiedAngularBins
            self.sparseCoverage = sparseCoverage
            self.normalizedResidual = normalizedResidual
            self.innerSigma = innerSigma
            self.outerSigma = outerSigma
            self.annulusWidth = annulusWidth
        }

        /// R = (a + b)/2, py4DSTEM's convention (`R = np.mean((a, b))`); a/b is printed beside it, not folded in.
        public var radius: Double { (semiAxisA + semiAxisB) / 2 }
        /// a/b − 1 with a the larger axis.
        public var ellipticity: Double { max(semiAxisA, semiAxisB) / min(semiAxisA, semiAxisB) - 1 }
    }

    // MARK: Registered bounds

    /// The azimuthal bins the fit divides the annulus into, and the fewest that must hold signal (a gap above about
    /// 60°). The same numbers as the ring fit's own (slice 3b holds them equal in a test).
    public static let angularBinCount = 36
    public static let minimumOccupiedAngularBins = 30
    /// The profile fit's cap on the normalised residual.
    public static let maximumResidual = 0.22
    /// Past a/b − 1 = 0.093, (a + b)/2 and √(ab) differ by more than 0.1 %. Derived, not measured: the root of
    /// (2 + ε) / (2√(1 + ε)) − 1 = 1e-3 is ε = 0.0936 (1.0009e-3 there; 9.886e-4 at 0.093), registered rounded down.
    public static let maximumEllipticity = 0.093
    /// sin 45° = 1/√2: λ/2d at or above it puts 2θ at or beyond 90° (an addition beyond the registration: see the header).
    public static let maximumLambdaOver2d = 0.5.squareRoot()

    // MARK: Constants (CODATA 2018, stated in the header)

    public static let planckConstant = ElectronWavelength.planckConstant                 // J s
    public static let speedOfLight = ElectronWavelength.speedOfLight                     // m/s
    public static let elementaryCharge = ElectronWavelength.elementaryCharge             // C
    public static let electronMass = ElectronWavelength.electronMass                     // kg
    /// hc/e in eV·Å: 12398.41984332003.
    public static let planckTimesLightEVAngstroms = ElectronWavelength.planckTimesLightEVAngstroms
    /// m₀c²/e in eV: 510998.9499961642.
    public static let electronRestEnergyEV = ElectronWavelength.electronRestEnergyEV

    // MARK: The closed forms

    /// λ in Å at an accelerating voltage in VOLTS (200 kV is 200_000). NaN unless volts is finite and positive.
    public static func wavelengthAngstroms(volts: Double) -> Double {
        ElectronWavelength.angstroms(volts: volts)
    }

    /// 2θ in radians for a ring of spacing d (Å) at wavelength λ (Å). NaN unless both are finite and positive and
    /// λ/2d is below `maximumLambdaOver2d`.
    public static func twoTheta(wavelengthAngstroms: Double, d: Double) -> Double {
        guard wavelengthAngstroms.isFinite, wavelengthAngstroms > 0, d.isFinite, d > 0 else { return .nan }
        let ratio = wavelengthAngstroms / (2 * d)
        guard ratio < maximumLambdaOver2d else { return .nan }
        return 2 * asin(ratio)
    }

    /// L_px = R / tan 2θ, pixels. NaN unless R is finite and positive and 2θ lies in (0, π/2).
    /// VALIDITY: the form is exact, but 2θ is a double and tan 2θ diverges at 90°. The roundings of λ, of λ/2d and of asin
    /// (about 5e-16 rad in all at λ/2d near 0.7) reach L_px multiplied by 1/(sin 2θ cos 2θ), so against 60-digit
    /// arithmetic its relative error is at most about 5.5e-16/|cos 2θ|: a sweep of 9000 random (V, d) draws at λ/2d from
    /// 0.706 to 0.70699 found 5.5e-16/|cos 2θ| at worst (9.4e-13 at 0.7069, 1.6e-12 at 0.70699), and one probe at
    /// λ/2d = 0.707 gave 3.5e-13. The 1e-12 band of tools/ring-calibration-parity therefore holds for 2θ below about
    /// 89.97°, i.e. d above about 0.018 Å at 200 kV (0.030 Å at 80 kV, 0.014 Å at 300 kV), a spacing no lattice has.
    /// Closer to 90° the number is right only to that error; `.braggLimit` is the only guard.
    public static func cameraLengthPixels(radius: Double, twoTheta: Double) -> Double {
        guard radius.isFinite, radius > 0, twoTheta.isFinite, twoTheta > 0, twoTheta < .pi / 2 else { return .nan }
        let length = radius / tan(twoTheta)
        return length.isFinite ? length : .nan   // an overflow is NaN, never infinity (Gate B on slice 3a, F2)
    }

    /// Δq = (1/d)/R, Å⁻¹ per pixel, for d in Å and R in pixels. NaN unless both are finite and positive.
    public static func qPerPixel(d: Double, radius: Double) -> Double {
        guard d.isFinite, d > 0, radius.isFinite, radius > 0 else { return .nan }
        let q = (1 / d) / radius
        return q.isFinite ? q : .nan   // a subnormal d or R overflows: NaN, never infinity (Gate B on slice 3a, F2)
    }

    /// L = p · L_px, millimetres, for a pixel pitch p in MILLIMETRES (a 14 µm pixel is 0.014). NaN unless p is finite
    /// and positive.
    public static func cameraLengthMillimetres(pixelPitchMM: Double, cameraLengthPixels: Double) -> Double {
        guard pixelPitchMM.isFinite, pixelPitchMM > 0 else { return .nan }
        return pixelPitchMM * cameraLengthPixels
    }

    // MARK: The calibration

    /// Δq, and L where V (volts) is stated and the ring can exist at it. A V that is nil is "not stated": it is never
    /// defaulted, and the result carries `.noHighTension` with Δq. Throws `.invalidInput` for a radius, spacing, V or
    /// pitch that is not finite and positive. The refusals that come from the FIT (`refusals`) are the caller's to
    /// apply before it calls this.
    public static func calibrate(radius: Double, d: Double, highTensionVolts: Double?,
                                  pixelPitchMM: Double? = nil) throws -> RingCalibration {
        guard radius.isFinite, radius > 0 else { throw Refusal.invalidInput("the ring radius must be a positive number of pixels") }
        guard d.isFinite, d > 0 else { throw Refusal.invalidInput("the ring spacing must be a positive number of Å") }
        if let volts = highTensionVolts, !(volts.isFinite && volts > 0) {
            throw Refusal.invalidInput("the high tension must be a positive number of volts")
        }
        if let pitch = pixelPitchMM, !(pitch.isFinite && pitch > 0) {
            throw Refusal.invalidInput("the pixel pitch must be a positive number of millimetres")
        }
        let deltaQ = qPerPixel(d: d, radius: radius)
        guard let volts = highTensionVolts else {
            return RingCalibration(d: d, radius: radius, deltaQ: deltaQ, highTensionVolts: nil, wavelengthAngstroms: nil,
                                   twoTheta: nil, lengthPixels: nil, lengthMillimetres: nil, lengthRefusal: .noHighTension)
        }
        let lambda = wavelengthAngstroms(volts: volts)
        let angle = twoTheta(wavelengthAngstroms: lambda, d: d)
        guard angle.isFinite else {
            return RingCalibration(d: d, radius: radius, deltaQ: deltaQ, highTensionVolts: volts, wavelengthAngstroms: lambda,
                                   twoTheta: nil, lengthPixels: nil, lengthMillimetres: nil,
                                   lengthRefusal: .braggLimit(lambdaOver2d: lambda / (2 * d)))
        }
        let lengthPixels = cameraLengthPixels(radius: radius, twoTheta: angle)
        return RingCalibration(d: d, radius: radius, deltaQ: deltaQ, highTensionVolts: volts, wavelengthAngstroms: lambda,
                               twoTheta: angle, lengthPixels: lengthPixels,
                               lengthMillimetres: pixelPitchMM.map { cameraLengthMillimetres(pixelPitchMM: $0, cameraLengthPixels: lengthPixels) },
                               lengthRefusal: nil)
    }

    /// Every registered refusal that applies to a fit's numbers and the stated V (volts, nil when none), in the order
    /// model, coverage, sigma, residual, ellipticity, outside the annulus, not a peak, geometry (the last three only
    /// when the report carries the annulus's radii and the ring evidence, DiffractoKit 0.0.4), high tension; empty when
    /// none applies. A number that is not finite (or a count out of range) returns `[.invalidInput]` alone: a NaN would
    /// otherwise pass every comparison. A caller withholds Δq when any refusal's `withholdsDeltaQ` is true.
    public static func refusals(_ fit: FitReport, highTensionVolts: Double?) -> [Refusal] {
        guard highTensionVolts.map({ $0.isFinite && $0 > 0 }) ?? true else {
            return [.invalidInput("a fit number is not finite and positive, or the bin count is out of range")]
        }
        var refused = fitRefusals(fit)
        if refused.first?.name == "invalidInput" { return refused }
        if highTensionVolts == nil { refused.append(.noHighTension) }
        return refused
    }

    /// The refusals that come from the fit alone (everything but `.noHighTension`), in `refusals`' order.
    static func fitRefusals(_ fit: FitReport) -> [Refusal] {
        let sigmas = [fit.innerSigma, fit.outerSigma].compactMap { $0 }
        let valid = fit.semiAxisA.isFinite && fit.semiAxisA > 0 && fit.semiAxisB.isFinite && fit.semiAxisB > 0
            && (0...angularBinCount).contains(fit.occupiedAngularBins)
            && fit.normalizedResidual.isFinite && fit.normalizedResidual >= 0
            && sigmas.allSatisfy { $0.isFinite && $0 > 0 }
            && fit.annulusWidth.isFinite && fit.annulusWidth > 0
            && [fit.annulusInnerRadius, fit.annulusOuterRadius].allSatisfy { $0.map { $0.isFinite && $0 >= 0 } ?? true }
        guard valid else { return [.invalidInput("a fit number is not finite and positive, or the bin count is out of range")] }

        var refused: [Refusal] = []
        if fit.model != .asymmetricGaussian { refused.append(.modelNotAsymmetricGaussian) }
        // The coverage is the fit's pre-fit count and, with the evidence, the sectors that hold the ring at the fit
        // (DiffractoKit 0.0.4: noise in a gap makes pixels bright, not a ring at the fitted radius).
        let evidenceBins = fit.ringEvidence?.occupiedBins ?? angularBinCount
        if fit.occupiedAngularBins < minimumOccupiedAngularBins || fit.sparseCoverage
            || evidenceBins < RingEvidence.minimumOccupiedBins {
            refused.append(.insufficientAngularCoverage(occupiedBins: min(fit.occupiedAngularBins, evidenceBins)))
        }
        if let widest = sigmas.max(), widest >= fit.annulusWidth {
            refused.append(.sigmaAtAnnulusWidth(sigma: widest, annulusWidth: fit.annulusWidth))
        }
        if fit.normalizedResidual > maximumResidual { refused.append(.excessiveResidual(fit.normalizedResidual)) }
        if fit.ellipticity > maximumEllipticity { refused.append(.excessiveEllipticity(fit.ellipticity)) }
        if let inner = fit.annulusInnerRadius, let outer = fit.annulusOuterRadius,
           !(fit.radius > inner && fit.radius < outer) {
            refused.append(.ringOutsideAnnulus(radius: fit.radius, inner: inner, outer: outer))
        }
        if let evidence = fit.ringEvidence {
            if !evidence.isPeak { refused.append(.notARingPeak(innerDrop: evidence.innerDrop, outerDrop: evidence.outerDrop)) }
            if !evidence.isSingleRing, let height = evidence.secondaryPeak { refused.append(.secondRingInAnnulus(height: height)) }
            // Undefined shifts come with fewer than 6 sectors holding the ring, which the coverage refusal names.
            if !evidence.geometryHolds, let centre = evidence.centreShift, let ellipse = evidence.ellipseShift {
                refused.append(.ringGeometryInconsistent(centreShift: centre, ellipseShift: ellipse))
            }
        }
        return refused
    }
}
