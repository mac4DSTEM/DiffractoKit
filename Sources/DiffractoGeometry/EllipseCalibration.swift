//
//  EllipseCalibration.swift
//  Role: detector-ellipse calibration of a ring pattern, in detector PIXELS. The conic path matches py4DSTEM
//        `fit_ellipse_1D`; the profile path ports `fit_ellipse_amorphous_ring`'s 11-parameter central Gaussian +
//        asymmetric ("Janus") Gaussian ring model; `fitBestAvailable` chains them. Slice 3b of DiffractoKit
//        (docs/registration-slice-3.md section 2).
//
//  SEED: mac4DSTEM 7d6fdbd1 `mac4DSTEM/Core/Analysis/EllipseCalibration.swift` (853 lines), ported verbatim in arithmetic
//  (every constant, every boundary, operation order): coverage 36 / 30 / 12 azimuthal bins, residual caps 0.2 / 0.4 / 0.22,
//  the one-ring ratio 1.10, the MAD cut 4 * 1.4826 and 70 % inliers, axis ratio 0.2 and 5, and ADR 069's two starts and its
//  0.30 px ring-width floor guard. NOT ported, because they never landed there: the annulus guard (`ringOutsideAnnulus`)
//  and the seed check. The only changes are the access (`public`), `DiffractionPattern`'s field names (`pattern.height` is
//  the seed's `pattern.qy`, the ROW count; `pattern.width` its `pattern.qx`, the COLUMN count) and the removal of
//  `EllipseCalibrationModel`, which `RingCalibration.swift` holds. The tests prove the arithmetic by bit-identity against
//  the seed compiled as it was (Tests/.../EllipseAppsParityTests.swift, tools/slice-3b/apps-parity).
//
//  COORDINATES (py4DSTEM's order, as the seed): `qx` is the ROW index and `qy` the COLUMN index everywhere, so the loops
//  below read `for qx in 0..<pattern.height` and `pixels[qx * pattern.width + qy]`. `theta` lies in [0, pi) and is measured
//  from the row axis towards the column axis; the major axis `a` lies along (cos theta, sin theta) in (row, column). The
//  conic path returns a >= b; the profile path keeps a >= b by swapping and adding pi/2 to theta.
//
//  FitError messages name no application; the one that reads "select an annulus" speaks of the annulus the caller gave.
//

import Foundation

public struct EllipseProfileParameters: Sendable, Equatable {
    public let centralIntensity: Double
    public let ringIntensity: Double
    public let centralSigma: Double
    public let innerSigma: Double
    public let outerSigma: Double
    public let background: Double

    // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
    public init(centralIntensity: Double, ringIntensity: Double, centralSigma: Double, innerSigma: Double, outerSigma: Double, background: Double) {
        self.centralIntensity = centralIntensity
        self.ringIntensity = ringIntensity
        self.centralSigma = centralSigma
        self.innerSigma = innerSigma
        self.outerSigma = outerSigma
        self.background = background
    }
}

public struct EllipseCalibrationFit: Sendable, Equatable {
    /// Fitted center in py4DSTEM detector order: qx is row, qy is column.
    public let centerQX: Double
    public let centerQY: Double
    public let a: Double
    public let b: Double
    public let theta: Double
    /// Dimensionless residual for the model that was accepted.
    public let normalizedResidual: Double
    public let conicResidual: Double
    public let sampleCount: Int
    public let occupiedAngularBins: Int
    public let model: EllipseCalibrationModel
    public let profile: EllipseProfileParameters?
    /// Populated when profile refinement was attempted but the conic fallback
    /// was retained. This is diagnostic provenance, not a scientific failure.
    public let profileFallbackReason: String?
    /// True only when this fit was accepted below the angular-coverage
    /// degeneracy bound because the caller asserted the annulus holds one
    /// ring (`acceptSparseCoverage`). Inert (false) whenever coverage was
    /// already full. 2026-09-15.
    public let sparseCoverage: Bool

    // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
    public init(centerQX: Double, centerQY: Double, a: Double, b: Double, theta: Double, normalizedResidual: Double, conicResidual: Double, sampleCount: Int, occupiedAngularBins: Int, model: EllipseCalibrationModel, profile: EllipseProfileParameters?, profileFallbackReason: String?, sparseCoverage: Bool = false) {
        self.centerQX = centerQX
        self.centerQY = centerQY
        self.a = a
        self.b = b
        self.theta = theta
        self.normalizedResidual = normalizedResidual
        self.conicResidual = conicResidual
        self.sampleCount = sampleCount
        self.occupiedAngularBins = occupiedAngularBins
        self.model = model
        self.profile = profile
        self.profileFallbackReason = profileFallbackReason
        self.sparseCoverage = sparseCoverage
    }
}

public enum EllipseCalibration {
    /// Azimuthal bins the coverage checks below divide the annulus into, and
    /// the two bounds read against them: full coverage needed by default, and
    /// the restored pre-2026-09-14 floor below which even "fit anyway" refuses.
    public static let angularBinCount = 36
    public static let degeneracyBoundBins = angularBinCount * 5 / 6
    public static let sparseFloorBins = angularBinCount / 3

    public enum FitError: LocalizedError, Equatable {
        case invalidInput(String)
        case insufficientSignal
        case insufficientAngularCoverage(Int)
        case moreThanOneRing(minRadius: Double, maxRadius: Double)
        case didNotConverge
        case invalidEllipse
        case excessiveResidual(Double)
        case ringWidthAtFloor

        public var errorDescription: String? {
            switch self {
            case .invalidInput(let detail):
                return "Cannot fit detector ellipse: \(detail)."
            case .insufficientSignal:
                return "Cannot fit detector ellipse: the selected annulus has no resolved ring signal."
            case .insufficientAngularCoverage(let bins):
                return "Cannot fit detector ellipse: ring signal covers only \(bins) of "
                    + "\(EllipseCalibration.angularBinCount) angular bins. An ellipse fitted to spots is decided by where the "
                    + "grains happen to be as much as by the detector — select an annulus "
                    + "holding one continuous ring, or use an amorphous standard."
            case .moreThanOneRing(let minRadius, let maxRadius):
                return String(format: "Cannot fit anyway: the strong signal in this annulus "
                    + "lies at radii from %.1f to %.1f px, more than one ring apart — no "
                    + "detector is that elliptical. Narrow the annulus to one ring.",
                    minRadius, maxRadius)
            case .didNotConverge:
                return "Detector ellipse fitting did not converge."
            case .invalidEllipse:
                return "Detector ellipse fitting produced a non-physical conic."
            case .excessiveResidual(let residual):
                return String(format: "Detector ellipse fit residual is too large (%.3f).", residual)
            case .ringWidthAtFloor:
                return "The ring profile fit lost the ring: its fitted width ended within 0.05 px of the 0.25 px floor."
            }
        }
    }

    /// Port of py4DSTEM `fit_ellipse_1D`'s operative model:
    /// `(A dx² + B dx dy + C dy² - 1) * intensity`.
    /// Coordinates deliberately stay in py4DSTEM order (qx=row, qy=column).
    ///
    /// RADIUS BIAS, a property of the port, not a deviation (2026-10-12,
    /// macTEM's report): the objective's minimum lies OUTSIDE a ring's peak.
    /// For a Gaussian ring of width s at radius R, no background, an annulus
    /// many s wide, (a + b)/2 = R (1 + 1.75 (s/R)²) — the r dr measure, the
    /// conic's r² and the squared intensity (R_f² = ∫r⁵I² dr / ∫r³I² dr).
    /// `tools/ellipse-calibration-test` prints it (+0.2456 % at R 80 px, s 3 px).
    /// That is a floor, not a correction: background and the annulus width
    /// move it, and the coefficient is the Gaussian's. A caller applies only
    /// b/a and θ (`EllipseTransform`), so a and b are not a
    /// calibrated ring radius (`archive/v5/ellipse-stop-rule-gateD-2026-10-12.md`).
    /// width × height equals the pixel count, without an Int overflow trap (Gate B on slice 3b, F5: a width of
    /// Int.max/2 + 1 trapped in the product). DEVIATION from the apps: they multiply unchecked.
    static func patternSizeIsConsistent(_ pattern: DiffractionPattern) -> Bool {
        let (count, overflow) = pattern.height.multipliedReportingOverflow(by: pattern.width)
        return !overflow && pattern.pixels.count == count
    }

    /// The outer radius is finite and at most twice the pattern's diagonal: beyond that no pixel is added, and an infinite
    /// or huge radius trapped in `Int(...)` or allocated gigabytes of bins (Gate B on slice 3b, F4). DEVIATION from the
    /// apps, which check only outer > inner; no annulus inside the pattern is affected.
    static func outerRadiusIsSane(_ outerRadius: Double, _ pattern: DiffractionPattern) -> Bool {
        let diagonal = (Double(pattern.width) * Double(pattern.width) + Double(pattern.height) * Double(pattern.height)).squareRoot()
        return outerRadius.isFinite && outerRadius <= 2 * diagonal
    }

    public static func fit1D(
        pattern: DiffractionPattern,
        centerQX: Double,
        centerQY: Double,
        innerRadius: Double,
        outerRadius: Double,
        maximumResidual: Double = 0.2,
        acceptSparseCoverage: Bool = false
    ) throws -> EllipseCalibrationFit {
        guard pattern.height > 2, pattern.width > 2,
              patternSizeIsConsistent(pattern),
              centerQX.isFinite, centerQY.isFinite,
              innerRadius >= 0, outerRadius > innerRadius, outerRadiusIsSane(outerRadius, pattern),
              maximumResidual > 0 else {
            throw FitError.invalidInput("dimensions, center, or fit radii are invalid")
        }

        struct Sample {
            let x: Double
            let y: Double
            let value: Double
        }
        var samples = [Sample]()
        samples.reserveCapacity(pattern.pixels.count / 2)
        var minimum = Double.greatestFiniteMagnitude
        var maximum = -Double.greatestFiniteMagnitude
        for qx in 0..<pattern.height {
            let dx = Double(qx) - centerQX
            for qy in 0..<pattern.width {
                let dy = Double(qy) - centerQY
                let radius = hypot(dx, dy)
                guard radius > innerRadius, radius <= outerRadius else { continue }
                let value = Double(pattern.pixels[qx * pattern.width + qy])
                guard value.isFinite else { continue }
                samples.append(Sample(x: Double(qx), y: Double(qy), value: value))
                minimum = min(minimum, value)
                maximum = max(maximum, value)
            }
        }
        guard samples.count >= 20 else {
            throw FitError.invalidInput("the fitting annulus contains too few pixels")
        }
        let dynamicRange = maximum - minimum
        guard dynamicRange.isFinite,
              dynamicRange > max(1, abs(maximum)) * 1e-6 else {
            throw FitError.insufficientSignal
        }

        // Reject spot-only/partial rings before optimization.
        //
        // DEVIATION from py4DSTEM, deliberate and now the stricter of two.
        // `fit_ellipse_1D` (`process/calibration/ellipse.py`) has no guard at
        // all: it least-squares whatever it is handed and returns five numbers.
        // Its answer to the degenerate case is `constrain_degenerate_ellipse`,
        // which adds a constraint rather than refusing. This port refuses,
        // because a caller has no place to carry "these numbers came from a
        // constraint you did not choose" and the numbers go straight into
        // Bragg coordinates and strain.
        //
        // WHY FIVE SIXTHS AND NOT A THIRD (2026-09-14). A third was written to
        // stop "a plausible conic from four peaks", and it does. It does not
        // stop one from several GRAINS: on the demo cube the owner's fit
        // reported a = 43.68, b = 39.72 — 10 % distortion — on a detector that
        // is isotropic by construction, because grains at different radii in
        // different azimuthal sectors let one ellipse thread through them.
        //
        // Three statistics were measured against a fixture sweep and all three
        // were refuted (`docs/open-items.md`): azimuthal contrast, the fit's
        // own residual, and radial multiplicity. The reason none of them works
        // is visible in the sweep: a three-grain annulus and a LEGITIMATE
        // six-azimuth ring on a distorted detector occupy the same 12 bins,
        // carry residuals of 0.082 and 0.118, and differ in nothing a
        // statistic can read — only in the answer. They are the same
        // measurement. An ellipse has five free parameters, and spots at a
        // dozen azimuths determine it no better than the three radii they
        // happen to lie on.
        //
        // So this is a DEGENERACY bound, not a separation: the fit is refused
        // wherever the data cannot decide, whether or not the answer would
        // have been right. It costs the sparse legitimate case, and the
        // message says what to do instead. What it buys is that a 10 %
        // distortion measured from the arrangement of three grains can no
        // longer reach a strain map.
        //
        // THE ANYWAY PATH (2026-09-15). The owner decided a caller may assert,
        // per annulus, "this holds exactly one ring" (`acceptSparseCoverage`)
        // and get the fit anyway, MARKED `sparseCoverage`. That re-admits the
        // legitimate sparse case the bound above costs, down to a restored
        // pre-2026-09-14 floor of one third of the bins: fewer than that still
        // refuses outright, because no per-bin check can rescue coverage this
        // thin. Above the floor, an accepted fit is checked for exactly what
        // the coverage bound cannot tell apart — several grains at different
        // radii in one annulus — by the per-azimuthal-bin radius check below
        // the fit. Full coverage (>= five sixths of the bins) makes the flag
        // and that check inert: the default path is unchanged either way.
        let angularBinCount = EllipseCalibration.angularBinCount
        let strongThreshold = minimum + 0.2 * dynamicRange
        var occupied = [Bool](repeating: false, count: angularBinCount)
        for sample in samples where sample.value >= strongThreshold {
            var angle = atan2(sample.y - centerQY, sample.x - centerQX)
            if angle < 0 { angle += 2 * .pi }
            let index = min(angularBinCount - 1, Int(angle / (2 * .pi) * Double(angularBinCount)))
            occupied[index] = true
        }
        let occupiedCount = occupied.filter { $0 }.count
        let degeneracyBound = EllipseCalibration.degeneracyBoundBins
        let sparseCoverage: Bool
        if acceptSparseCoverage {
            // Pre-2026-09-14 rule, restored ONLY for the anyway path.
            guard occupiedCount >= EllipseCalibration.sparseFloorBins else {
                throw FitError.insufficientAngularCoverage(occupiedCount)
            }
            sparseCoverage = occupiedCount < degeneracyBound
        } else {
            guard occupiedCount >= degeneracyBound else {
                throw FitError.insufficientAngularCoverage(occupiedCount)
            }
            sparseCoverage = false
        }

        let initialRadius = (innerRadius + outerRadius) / 2
        var parameters = [centerQX, centerQY,
                          1 / (initialRadius * initialRadius), 0,
                          1 / (initialRadius * initialRadius)]

        func costAndNormal(_ p: [Double], buildNormal: Bool)
            -> (cost: Double, weight: Double, normal: [[Double]], gradient: [Double]) {
            var cost = 0.0
            var weight = 0.0
            var normal = Array(repeating: Array(repeating: 0.0, count: 5), count: 5)
            var gradient = [Double](repeating: 0, count: 5)
            for sample in samples {
                let dx = sample.x - p[0], dy = sample.y - p[1]
                let value = sample.value
                let residual = (p[2] * dx * dx + p[3] * dx * dy
                                + p[4] * dy * dy - 1) * value
                cost += residual * residual
                weight += value * value
                guard buildNormal else { continue }
                let jacobian = [
                    (-2 * p[2] * dx - p[3] * dy) * value,
                    (-p[3] * dx - 2 * p[4] * dy) * value,
                    dx * dx * value,
                    dx * dy * value,
                    dy * dy * value,
                ]
                for row in 0..<5 {
                    gradient[row] += jacobian[row] * residual
                    for column in row..<5 {
                        normal[row][column] += jacobian[row] * jacobian[column]
                    }
                }
            }
            for row in 0..<5 {
                for column in 0..<row { normal[row][column] = normal[column][row] }
            }
            return (cost, weight, normal, gradient)
        }

        var current = costAndNormal(parameters, buildNormal: true)
        guard current.weight > 0, current.cost.isFinite else { throw FitError.insufficientSignal }
        var damping = 1e-3
        var converged = false
        for _ in 0..<100 {
            var system = current.normal
            for diagonal in 0..<5 {
                system[diagonal][diagonal] += damping * max(system[diagonal][diagonal], 1e-12)
            }
            guard let delta = solve(system, rhs: current.gradient.map { -$0 }) else {
                damping *= 10
                continue
            }
            let candidate = zip(parameters, delta).map(+)
            guard isPhysicalConic(candidate) else {
                damping *= 10
                continue
            }
            let trial = costAndNormal(candidate, buildNormal: false)
            let relativeStep = sqrt(delta.reduce(0) { $0 + $1 * $1 })
                / max(1, sqrt(candidate.reduce(0) { $0 + $1 * $1 }))
            if trial.cost < current.cost {
                parameters = candidate
                current = costAndNormal(parameters, buildNormal: true)
                damping = max(1e-12, damping / 3)
                if relativeStep < 1e-10 {
                    converged = true
                    break
                }
            } else {
                // STOP RULE ON A REFUSED STEP TOO (2026-10-12). At the optimum
                // the cost change of a proposed step is roundoff, so the test
                // above refuses it, damping grows x10 per pass, the step shrinks
                // until the candidate equals the point, and with `converged`
                // set only on an accepted step the 100-pass cap threw
                // didNotConverge AT the optimum: 12 of 84 clean planted rings
                // (macTEM's report, `docs/archive/v5/ellipse-stop-rule-gateD-2026-10-12.md`).
                // A step this small that roundoff will not let lower the cost is
                // converged; the fit returns the last accepted point. py4DSTEM's
                // `leastsq` (MINPACK) tests its xtol whether or not the step is
                // accepted.
                if relativeStep < 1e-10 {
                    converged = true
                    break
                }
                damping *= 10
            }
        }
        guard converged else { throw FitError.didNotConverge }
        guard let ellipse = userParameters(A: parameters[2], B: parameters[3], C: parameters[4])
        else { throw FitError.invalidEllipse }
        let centerDistance = hypot(parameters[0] - centerQX, parameters[1] - centerQY)
        guard centerDistance <= max(outerRadius - innerRadius, 2),
              ellipse.b / ellipse.a >= 0.2,
              ellipse.a <= outerRadius * 1.5,
              ellipse.b >= max(1, innerRadius * 0.25) else {
            throw FitError.invalidEllipse
        }
        let residual = sqrt(current.cost / current.weight)
        guard residual.isFinite, residual <= maximumResidual else {
            throw FitError.excessiveResidual(residual)
        }

        // ONE-RING CHECK (2026-09-15), sparse path only, after every guard
        // above has passed. The coverage bound exists because it cannot tell
        // a single distorted ring from several grains at different radii in
        // one annulus; once the caller has asserted "one ring" and bypassed
        // it, this catches the one thing that assertion does not make true.
        // About the FITTED centre, take the intensity-weighted mean radius of
        // the strong samples in each of the 36 azimuthal bins, and compare the
        // largest such mean to the smallest over the bins that hold any. A
        // ratio past 1.10 is more ellipticity than any STEM camera distortion
        // this calibration targets — py4DSTEM's own worked examples run a few
        // percent — and less than the closest common ring pair, fcc
        // {111}/{200} at 15.5%; hcp {100}/{002} at 6% is inside the bound and
        // is the caller's responsibility under "fit anyway".
        // DEVIATION: py4DSTEM `fit_ellipse_1D` has no such path; its answer to
        // degeneracy is `constrain_degenerate_ellipse`, not a refusal.
        if sparseCoverage {
            var weightedRadius = [Double](repeating: 0, count: angularBinCount)
            var weight = [Double](repeating: 0, count: angularBinCount)
            for sample in samples where sample.value >= strongThreshold {
                let dx = sample.x - parameters[0], dy = sample.y - parameters[1]
                var angle = atan2(dy, dx)
                if angle < 0 { angle += 2 * .pi }
                let index = min(angularBinCount - 1, Int(angle / (2 * .pi) * Double(angularBinCount)))
                weightedRadius[index] += hypot(dx, dy) * sample.value
                weight[index] += sample.value
            }
            var minMeanRadius = Double.greatestFiniteMagnitude
            var maxMeanRadius = -Double.greatestFiniteMagnitude
            for index in 0..<angularBinCount where weight[index] > 0 {
                let meanRadius = weightedRadius[index] / weight[index]
                minMeanRadius = min(minMeanRadius, meanRadius)
                maxMeanRadius = max(maxMeanRadius, meanRadius)
            }
            if minMeanRadius > 0, maxMeanRadius / minMeanRadius > 1.10 {
                throw FitError.moreThanOneRing(minRadius: minMeanRadius, maxRadius: maxMeanRadius)
            }
        }

        return EllipseCalibrationFit(
            centerQX: parameters[0], centerQY: parameters[1],
            a: ellipse.a, b: ellipse.b, theta: ellipse.theta,
            normalizedResidual: residual, conicResidual: residual,
            sampleCount: samples.count, occupiedAngularBins: occupiedCount,
            model: .conic, profile: nil, profileFallbackReason: nil,
            sparseCoverage: sparseCoverage
        )
    }

    /// Prefer the physically richer amorphous-ring profile, but never discard
    /// a valid conic calibration merely because the profile is absent,
    /// underconstrained, or incompatible with crystalline/overlapping signal.
    public static func fitBestAvailable(
        pattern: DiffractionPattern,
        centerQX: Double,
        centerQY: Double,
        innerRadius: Double,
        outerRadius: Double,
        maximumConicResidual: Double = 0.2,
        maximumProfileResidual: Double = 0.22,
        acceptSparseCoverage: Bool = false
    ) throws -> EllipseCalibrationFit {
        // A broad/asymmetric ring can legitimately exceed the strict conic
        // publication residual while still providing a useful profile seed.
        let conic = try fit1D(
            pattern: pattern,
            centerQX: centerQX, centerQY: centerQY,
            innerRadius: innerRadius, outerRadius: outerRadius,
            maximumResidual: max(0.4, maximumConicResidual),
            acceptSparseCoverage: acceptSparseCoverage
        )
        do {
            return try fitAmorphousRing(
                pattern: pattern, initial: conic,
                centerQX: centerQX, centerQY: centerQY,
                innerRadius: innerRadius, outerRadius: outerRadius,
                maximumResidual: maximumProfileResidual
            )
        } catch {
            guard conic.conicResidual <= maximumConicResidual else { throw error }
            return EllipseCalibrationFit(
                centerQX: conic.centerQX, centerQY: conic.centerQY,
                a: conic.a, b: conic.b, theta: conic.theta,
                normalizedResidual: conic.normalizedResidual,
                conicResidual: conic.conicResidual,
                sampleCount: conic.sampleCount,
                occupiedAngularBins: conic.occupiedAngularBins,
                model: .conic, profile: nil,
                profileFallbackReason: error.localizedDescription,
                sparseCoverage: conic.sparseCoverage
            )
        }
    }

    /// Bounded Levenberg–Marquardt refinement of py4DSTEM's calibration-module
    /// 11-parameter model. We optimize equivalent user-facing `(a,b,theta)`
    /// variables instead of unconstrained `(A,B,C)`, then evaluate the exact
    /// canonical model. This is an intentional numerical-stability deviation:
    /// every candidate remains a physical ellipse.
    public static func fitAmorphousRing(
        pattern: DiffractionPattern,
        initial: EllipseCalibrationFit,
        centerQX: Double,
        centerQY: Double,
        innerRadius: Double,
        outerRadius: Double,
        maximumResidual: Double = 0.22
    ) throws -> EllipseCalibrationFit {
        struct Sample {
            let x: Double
            let y: Double
            let value: Double
        }
        guard pattern.height > 2, pattern.width > 2,
              patternSizeIsConsistent(pattern),
              innerRadius >= 0, outerRadius > innerRadius, outerRadiusIsSane(outerRadius, pattern),
              maximumResidual > 0 else {
            throw FitError.invalidInput("dimensions or profile bounds are invalid")
        }
        var samples: [Sample] = []
        var minimum = Double.greatestFiniteMagnitude
        var maximum = -Double.greatestFiniteMagnitude
        var globalMaximum = maximum
        for value in pattern.pixels where value.isFinite {
            globalMaximum = max(globalMaximum, Double(value))
        }
        for qx in 0..<pattern.height {
            let dx = Double(qx) - centerQX
            for qy in 0..<pattern.width {
                let dy = Double(qy) - centerQY
                let radius = hypot(dx, dy)
                guard radius > innerRadius, radius < outerRadius else { continue }
                let value = Double(pattern.pixels[qx * pattern.width + qy])
                guard value.isFinite else { continue }
                samples.append(Sample(x: Double(qx), y: Double(qy), value: value))
                minimum = min(minimum, value)
                maximum = max(maximum, value)
            }
        }
        let dynamicRange = maximum - minimum
        guard samples.count >= 32, dynamicRange.isFinite,
              dynamicRange > max(1, abs(maximum)) * 1e-6 else {
            throw FitError.insufficientSignal
        }

        let ringWidth = outerRadius - innerRadius
        let intensityLimit = 5 * max(1, abs(globalMaximum), dynamicRange)
        let sigmaLimit = max(2, outerRadius * 2)
        let centerLimit = max(2, ringWidth)
        let axisMinimum = max(1, innerRadius * 0.5)
        let axisMaximum = max(axisMinimum + 1, outerRadius * 1.5)
        let backgroundLower = minimum - 2 * dynamicRange
        let backgroundUpper = maximum + 2 * dynamicRange

        // DEVIATION: py4DSTEM `fit_ellipse_amorphous_ring`
        // (process/calibration/ellipse.py:189-199) seeds absolute intensities
        // I0=max(data), I1=max(data*mask), c=min(data), and derives its start
        // geometry from a radial integral (a circle a=b=R, theta=0). This port
        // runs *after* the conic fit, so it seeds (a,b,theta) from `initial`
        // rather than a radial-integral R, and uses baseline-subtracted,
        // annulus-restricted amplitudes: I0=globalMax-annularMin,
        // I1=annularMax-annularMin, c=annularMin. Same model and residual, only
        // the LM start point differs; the effect on convergence is unmeasured.
        // [I0,I1,sigma0,sigma1,sigma2,c,x0,y0,a,b,theta]
        let parameters = [
            max(0, globalMaximum - minimum), dynamicRange,
            max(1, innerRadius / 2), max(1, ringWidth / 4), max(1, ringWidth / 4),
            minimum, initial.centerQX, initial.centerQY,
            initial.a, initial.b, initial.theta,
        ]
        let scales = [
            max(1, dynamicRange), max(1, dynamicRange),
            max(1, innerRadius / 2), max(1, ringWidth / 4), max(1, ringWidth / 4),
            max(1, dynamicRange), 1, 1,
            max(1, initial.a), max(1, initial.b), 1,
        ]

        func bounded(_ input: [Double]) -> [Double] {
            var p = input
            p[0] = min(intensityLimit, max(0, p[0]))
            p[1] = min(intensityLimit, max(0, p[1]))
            for index in 2...4 { p[index] = min(sigmaLimit, max(0.25, p[index])) }
            p[5] = min(backgroundUpper, max(backgroundLower, p[5]))
            p[6] = min(Double(pattern.height - 1), max(0, p[6]))
            p[6] = min(centerQX + centerLimit, max(centerQX - centerLimit, p[6]))
            p[7] = min(Double(pattern.width - 1), max(0, p[7]))
            p[7] = min(centerQY + centerLimit, max(centerQY - centerLimit, p[7]))
            p[8] = min(axisMaximum, max(axisMinimum, p[8]))
            p[9] = min(axisMaximum, max(axisMinimum, p[9]))
            if p[8] < p[9] {
                p.swapAt(8, 9)
                p[10] += .pi / 2
            }
            p[10].formTruncatingRemainder(dividingBy: .pi)
            if p[10] < 0 { p[10] += .pi }
            return p
        }

        func predictions(_ p: [Double], indices: [Int]) -> [Double] {
            let a = p[8], b = p[9], theta = p[10]
            let sinTheta = sin(theta), cosTheta = cos(theta)
            let a2 = a * a, b2 = b * b
            let sin2 = sinTheta * sinTheta, cos2 = cosTheta * cosTheta
            var A = sin2 / b2 + cos2 / a2
            var B = 2 * (b2 - a2) * sinTheta * cosTheta / (a2 * b2)
            var C = cos2 / b2 + sin2 / a2
            let radius = (a + b) / 2
            let radius2 = radius * radius
            A *= radius2; B *= radius2; C *= radius2
            let twoSigma0 = 2 * p[2] * p[2]
            let twoSigma1 = 2 * p[3] * p[3]
            let twoSigma2 = 2 * p[4] * p[4]
            return indices.map { index in
                let sample = samples[index]
                let dx = sample.x - p[6], dy = sample.y - p[7]
                let ellipticalRadius2 = max(0, A * dx * dx + B * dx * dy + C * dy * dy)
                let delta = sqrt(ellipticalRadius2) - radius
                let ringSigma = delta < 0 ? twoSigma1 : twoSigma2
                return p[0] * exp(-ellipticalRadius2 / twoSigma0)
                    + p[1] * exp(-(delta * delta) / ringSigma) + p[5]
            }
        }

        func optimize(_ start: [Double], indices: [Int]) throws -> [Double] {
            var p = bounded(start)
            var predicted = predictions(p, indices: indices)
            var residual = zip(predicted, indices).map { value, index in
                value - samples[index].value
            }
            var cost = residual.reduce(0) { $0 + $1 * $1 }
            guard cost.isFinite else { throw FitError.didNotConverge }
            var damping = 1e-3
            var acceptedSteps = 0
            let epsilon = 1e-4
            for _ in 0..<80 {
                var columns = [[Double]]()
                columns.reserveCapacity(parameters.count)
                for parameter in parameters.indices {
                    var plus = p, minus = p
                    plus[parameter] += scales[parameter] * epsilon
                    minus[parameter] -= scales[parameter] * epsilon
                    plus = bounded(plus); minus = bounded(minus)
                    let high = predictions(plus, indices: indices)
                    let low = predictions(minus, indices: indices)
                    columns.append(zip(high, low).map { ($0 - $1) / (2 * epsilon) })
                }
                var normal = Array(
                    repeating: [Double](repeating: 0, count: parameters.count),
                    count: parameters.count
                )
                var gradient = [Double](repeating: 0, count: parameters.count)
                for row in parameters.indices {
                    for sampleIndex in residual.indices {
                        gradient[row] += columns[row][sampleIndex] * residual[sampleIndex]
                    }
                    for column in row..<parameters.count {
                        var value = 0.0
                        for sampleIndex in residual.indices {
                            value += columns[row][sampleIndex] * columns[column][sampleIndex]
                        }
                        normal[row][column] = value
                        normal[column][row] = value
                    }
                }
                for diagonal in parameters.indices {
                    normal[diagonal][diagonal] += damping
                        * max(normal[diagonal][diagonal], 1e-12)
                }
                guard let step = solve(normal, rhs: gradient.map { -$0 }) else {
                    damping *= 10
                    continue
                }
                var candidate = p
                for index in parameters.indices { candidate[index] += step[index] * scales[index] }
                candidate = bounded(candidate)
                let trialPrediction = predictions(candidate, indices: indices)
                let trialResidual = zip(trialPrediction, indices).map { value, index in
                    value - samples[index].value
                }
                let trialCost = trialResidual.reduce(0) { $0 + $1 * $1 }
                if trialCost < cost {
                    let relativeImprovement = (cost - trialCost) / max(cost, 1e-20)
                    p = candidate; predicted = trialPrediction; residual = trialResidual
                    cost = trialCost; acceptedSteps += 1
                    damping = max(1e-10, damping / 3)
                    let stepNorm = sqrt(step.reduce(0) { $0 + $1 * $1 })
                    if stepNorm < 1e-7 || relativeImprovement < 1e-10 { break }
                } else {
                    damping *= 10
                }
            }
            guard acceptedSteps > 0, p.allSatisfy(\.isFinite) else {
                throw FitError.didNotConverge
            }
            return p
        }

        let allIndices = Array(samples.indices)

        /// One refinement from `start`: the fit, py4DSTEM's optional robust refit, and the acceptance checks.
        func refine(_ start: [Double]) throws -> EllipseCalibrationFit {
            var parameters = try optimize(start, indices: allIndices)

            // One robust refit mirrors py4DSTEM's optional outlier refinement while
            // keeping the final acceptance residual defined over every selected pixel.
            let firstPrediction = predictions(parameters, indices: allIndices)
            let firstResiduals = zip(firstPrediction, samples).map { value, sample in
                value - sample.value
            }
            let centerResidual = median(firstResiduals) ?? 0
            let mad = median(firstResiduals.map { abs($0 - centerResidual) }) ?? 0
            if mad > 0 {
                let threshold = 4 * 1.4826 * mad
                let inliers = allIndices.filter { abs(firstResiduals[$0] - centerResidual) <= threshold }
                if inliers.count * 10 >= samples.count * 7, inliers.count < samples.count {
                    parameters = try optimize(parameters, indices: inliers)
                }
            }

            let finalPrediction = predictions(parameters, indices: allIndices)
            let finalCost = zip(finalPrediction, samples).reduce(0.0) { partial, pair in
                let error = pair.0 - pair.1.value
                return partial + error * error
            }
            let residual = sqrt(finalCost / Double(samples.count)) / dynamicRange
            guard residual.isFinite, residual <= maximumResidual else {
                throw FitError.excessiveResidual(residual)
            }
            guard hypot(parameters[6] - centerQX, parameters[7] - centerQY) <= centerLimit,
                  parameters[8] / parameters[9] <= 5 else {
                throw FitError.invalidEllipse
            }
            // A ring sigma within 0.05 px of its 0.25 px floor is a fit that drained the ring into
            // the background and the central Gaussian: a, b and the centre then
            // carry no ring information (`archive/v5/ring-fit-two-start-gateD-2026-10-14.md`).
            guard min(parameters[3], parameters[4]) > 0.25 + 0.05 else {
                throw FitError.ringWidthAtFloor
            }
            return EllipseCalibrationFit(
                centerQX: parameters[6], centerQY: parameters[7],
                a: parameters[8], b: parameters[9], theta: parameters[10],
                normalizedResidual: residual, conicResidual: initial.conicResidual,
                sampleCount: samples.count,
                occupiedAngularBins: initial.occupiedAngularBins,
                model: .asymmetricGaussian,
                profile: EllipseProfileParameters(
                    centralIntensity: parameters[0], ringIntensity: parameters[1],
                    centralSigma: parameters[2], innerSigma: parameters[3],
                    outerSigma: parameters[4], background: parameters[5]
                ),
                profileFallbackReason: nil,
                sparseCoverage: initial.sparseCoverage
            )
        }

        var starts = [parameters]
        // DEVIATION: py4DSTEM refines from one start, both ring sigmas at
        // (ro - ri)/4. On a noisy ring that start (10-17 px for a 3 px ring)
        // can end in a minimum with no ring information, and a start at the
        // ring's own half-maximum width can end worse on a sloping background
        // (`archive/v5/ring-fit-two-start-gateD-2026-10-14.md`). Both are
        // refined and the fit with the lower cost on the same pixels is kept.
        if let ringSigma = halfMaximumRingSigma(
            values: samples.map(\.value), xs: samples.map(\.x), ys: samples.map(\.y),
            seed: initial, innerRadius: innerRadius, outerRadius: outerRadius
        ), abs(ringSigma - parameters[3]) > 1e-9 {
            var second = parameters
            second[3] = ringSigma
            second[4] = ringSigma
            starts.append(second)
        }
        var best: EllipseCalibrationFit?
        var firstError: Error?
        for start in starts {
            do {
                let fit = try refine(start)
                if fit.normalizedResidual < best?.normalizedResidual ?? .infinity { best = fit }
            } catch {
                if firstError == nil { firstError = error }
            }
        }
        guard let best else { throw firstError ?? FitError.didNotConverge }
        return best
    }

    /// The ring sigma whose Gaussian has the half-maximum width of the
    /// annulus's radial profile: 1 px bins of the seed ellipse's elliptical
    /// radius (scaled to (a + b)/2, as the profile model measures it), a
    /// 3-bin mean, the half maximum between the profile's minimum and its
    /// peak, crossings interpolated linearly, FWHM / 2.3548, at least 1 px.
    /// Nil when the peak has no crossing on both sides inside the annulus
    /// (a ring wider than the annulus, or no ring), or fewer than 5 bins.
    private static func halfMaximumRingSigma(
        values: [Double], xs: [Double], ys: [Double],
        seed: EllipseCalibrationFit, innerRadius: Double, outerRadius: Double
    ) -> Double? {
        let a = seed.a, b = seed.b
        let binCount = Int((outerRadius - innerRadius).rounded(.up))
        guard a > 0, b > 0, binCount >= 5,
              values.count == xs.count, values.count == ys.count else { return nil }
        let sinTheta = sin(seed.theta), cosTheta = cos(seed.theta)
        let a2 = a * a, b2 = b * b
        let radius2 = (a + b) * (a + b) / 4
        let A = (sinTheta * sinTheta / b2 + cosTheta * cosTheta / a2) * radius2
        let B = 2 * (b2 - a2) * sinTheta * cosTheta / (a2 * b2) * radius2
        let C = (cosTheta * cosTheta / b2 + sinTheta * sinTheta / a2) * radius2
        var sums = [Double](repeating: 0, count: binCount)
        var counts = [Double](repeating: 0, count: binCount)
        for index in values.indices {
            let dx = xs[index] - seed.centerQX, dy = ys[index] - seed.centerQY
            let radius = max(0, A * dx * dx + B * dx * dy + C * dy * dy).squareRoot()
            let bin = Int((radius - innerRadius).rounded(.down))
            guard bin >= 0, bin < binCount else { continue }
            sums[bin] += values[index]
            counts[bin] += 1
        }
        let profile = (0..<binCount).map { counts[$0] > 0 ? sums[$0] / counts[$0] : Double.nan }
        let smooth: [Double] = (0..<binCount).map { bin in
            let window = (max(0, bin - 1)...min(binCount - 1, bin + 1)).map { profile[$0] }.filter(\.isFinite)
            return window.isEmpty ? Double.nan : window.reduce(0, +) / Double(window.count)
        }
        let finite = smooth.indices.filter { smooth[$0].isFinite }
        guard let peak = finite.max(by: { smooth[$0] < smooth[$1] }),
              let floor = finite.map({ smooth[$0] }).min(),
              smooth[peak] > floor else { return nil }
        let half = floor + (smooth[peak] - floor) / 2
        var left = peak
        while left > 0, smooth[left - 1].isFinite, smooth[left - 1] >= half { left -= 1 }
        var right = peak
        while right < binCount - 1, smooth[right + 1].isFinite, smooth[right + 1] >= half { right += 1 }
        guard left > 0, smooth[left - 1].isFinite,
              right < binCount - 1, smooth[right + 1].isFinite else { return nil }
        let leftCrossing = Double(left - 1)
            + (half - smooth[left - 1]) / (smooth[left] - smooth[left - 1])
        let rightCrossing = Double(right)
            + (smooth[right] - half) / (smooth[right] - smooth[right + 1])
        let width = rightCrossing - leftCrossing
        guard width.isFinite, width > 0 else { return nil }
        return max(1, width / (2 * (2 * log(2.0)).squareRoot()))
    }

    private static func median(_ values: [Double]) -> Double? {
        guard !values.isEmpty else { return nil }
        let sorted = values.sorted()
        let middle = sorted.count / 2
        return sorted.count.isMultiple(of: 2)
            ? (sorted[middle - 1] + sorted[middle]) / 2 : sorted[middle]
    }

    private static func isPhysicalConic(_ p: [Double]) -> Bool {
        guard p.count == 5, p.allSatisfy(\.isFinite) else { return false }
        return p[2] > 0 && p[4] > 0 && 4 * p[2] * p[4] - p[3] * p[3] > 0
    }

    /// Exact `convert_ellipse_params` convention from py4DSTEM.
    private static func userParameters(A: Double, B: Double, C: Double)
        -> (a: Double, b: Double, theta: Double)? {
        let separation = hypot(A - C, B)
        let smallEigenvalue = (A + C - separation) / 2
        let largeEigenvalue = (A + C + separation) / 2
        guard smallEigenvalue > 0, largeEigenvalue > 0 else { return nil }
        let a = 1 / sqrt(smallEigenvalue)
        let b = 1 / sqrt(largeEigenvalue)
        var theta: Double
        if abs(B) < 1e-15 {
            theta = A < C ? 0 : .pi / 2
        } else {
            theta = atan2(C - A - separation, B)
        }
        theta.formTruncatingRemainder(dividingBy: .pi)
        if theta < 0 { theta += .pi }
        return (a, b, theta)
    }

    private static func solve(_ matrix: [[Double]], rhs: [Double]) -> [Double]? {
        let count = rhs.count
        guard matrix.count == count, matrix.allSatisfy({ $0.count == count }) else { return nil }
        var augmented = matrix.enumerated().map { index, row in row + [rhs[index]] }
        for column in 0..<count {
            var pivot = column
            for row in (column + 1)..<count
                where abs(augmented[row][column]) > abs(augmented[pivot][column]) {
                pivot = row
            }
            guard abs(augmented[pivot][column]) > 1e-20 else { return nil }
            if pivot != column { augmented.swapAt(pivot, column) }
            let divisor = augmented[column][column]
            for index in column...count { augmented[column][index] /= divisor }
            for row in 0..<count where row != column {
                let factor = augmented[row][column]
                for index in column...count {
                    augmented[row][index] -= factor * augmented[column][index]
                }
            }
        }
        return augmented.map { $0[count] }
    }
}
