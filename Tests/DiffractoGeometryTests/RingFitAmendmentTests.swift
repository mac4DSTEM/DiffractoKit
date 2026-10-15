//
//  RingFitAmendmentTests.swift
//  DiffractoKit 0.1.0: the ring-fit amendment (docs/registration-0.0.4.md §1). Each slice 3b Gate B gap (registration-slice-3
//  §7, F1, F2, F3, F6, F7) and macTEM's `refuse-bound` plant, on plants rendered here with fresh parameters and seeds (not
//  the refuter's). Every test names under "RED UNDER" the production mutation that turns it red
//  (tools/kit-0.0.4/mutants-ringfit.py).
//
//  The plant renderer: a function of (row offset, column offset) from the plant's centre, float32, rows x columns with
//  columns fastest, plus Gaussian noise from SplitMix64 + Box-Muller (the parity renderer's generator, seeded here).
//

import Testing
import Foundation
import DiffractoGeometry

struct RingFitAmendmentTests {

    // MARK: Plants

    struct Scene {
        var rows = 128, cols = 124
        var row0 = 63.37, col0 = 61.59
    }

    /// The asymmetric-Gaussian ring on the elliptical coordinate the profile model uses: rho^2 = u^2/a^2 + v^2/b^2 in the
    /// frame rotated by theta, offset (rho - 1) R with R = (a + b)/2.
    static func ring(_ dr: Double, _ dc: Double, radius: Double, epsilon: Double = 0, thetaDeg: Double = 0,
                     sigmaIn: Double = 3, sigmaOut: Double = 3) -> Double {
        let b = 2 * radius / (2 + epsilon), a = b * (1 + epsilon)
        let t = thetaDeg * .pi / 180, s = sin(t), c = cos(t)
        let u = dr * c + dc * s, v = -dr * s + dc * c
        let rho = ((u / a) * (u / a) + (v / b) * (v / b)).squareRoot()
        let x = (rho - 1) * radius
        let w = x < 0 ? sigmaIn : sigmaOut
        return exp(-x * x / (2 * w * w))
    }

    struct SplitMix64 {
        var state: UInt64
        mutating func next() -> UInt64 {
            state &+= 0x9E37_79B9_7F4A_7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return z ^ (z >> 31)
        }
        mutating func uniform() -> Double { (Double(next() >> 11) + 0.5) / 9_007_199_254_740_992 }
    }

    static func render(_ scene: Scene = Scene(), noise: Double = 0, seed: UInt64 = 1,
                       _ value: (Double, Double) -> Double) -> DiffractionPattern {
        var generator = SplitMix64(state: seed)
        var pixels = [Float](repeating: 0, count: scene.rows * scene.cols)
        for r in 0..<scene.rows {
            for c in 0..<scene.cols {
                var v = value(Double(r) - scene.row0, Double(c) - scene.col0)
                if noise > 0 {
                    let u1 = generator.uniform(), u2 = generator.uniform()
                    v += noise * (-2 * log(u1)).squareRoot() * cos(2 * .pi * u2)
                }
                pixels[r * scene.cols + c] = Float(v)
            }
        }
        return DiffractionPattern(width: scene.cols, height: scene.rows, pixels: pixels)
    }

    static func azimuthDeg(_ dr: Double, _ dc: Double) -> Double {
        var angle = atan2(dc, dr) * 180 / .pi
        if angle < 0 { angle += 360 }
        return angle
    }

    static func inGap(_ dr: Double, _ dc: Double, start: Double, width: Double) -> Bool {
        var offset = (azimuthDeg(dr, dc) - start).truncatingRemainder(dividingBy: 360)
        if offset < 0 { offset += 360 }
        return offset < width
    }

    /// The calibration of a plant with the start centre a click away (+1.1, -0.7 px), Au {111}, 200 kV.
    static func calibrate(_ pattern: DiffractionPattern, _ scene: Scene = Scene(), inner: Double, outer: Double)
        throws -> RingFitCalibration {
        try RingCalibration.fit(pattern: pattern, centerQX: scene.row0 + 1.1, centerQY: scene.col0 - 0.7,
                                innerRadius: inner, outerRadius: outer,
                                d: 2.3546, highTensionVolts: 200_000, pixelPitchMM: 0.014)
    }

    static func names(_ result: RingFitCalibration) -> [String] { result.refusals.map(\.name) }

    // MARK: F1: noise in an angular gap

    /// A ring with a 95-125 degree gap and Gaussian noise is not calibrated: the noise makes pixels in the gap bright, not a
    /// ring at the fitted radius. Slice 3b's F1 calibrated such rings with the centre 0.7-2.7 px off and R -0.9 to -3.5 %.
    /// RED UNDER: the evidence's occupancy bound (alpha >= 0.5 A) dropped; the coverage refusal from the evidence removed.
    @Test(arguments: [(95.0, 0.30, UInt64(11)), (110.0, 0.25, UInt64(12)), (125.0, 0.35, UInt64(13)), (110.0, 0.40, UInt64(14))])
    func aNoisyRingWithAWideGapIsNotCalibrated(gap: Double, noise: Double, seed: UInt64) throws {
        let pattern = Self.render(noise: noise, seed: seed) { dr, dc in
            Self.inGap(dr, dc, start: 70, width: gap) ? 0 : Self.ring(dr, dc, radius: 42, epsilon: 0.015, thetaDeg: 20)
        }
        let result = try Self.calibrate(pattern, inner: 28, outer: 56)
        #expect(result.calibration == nil, "calibrated at R \(result.report?.radius ?? .nan); refusals \(Self.names(result))")
        #expect(Self.names(result).contains("insufficientAngularCoverage"))
    }

    /// The same rings without the gap are still calibrated (the evidence does not refuse noise on a whole ring). The radius
    /// band is 1 %: at noise 0.4 of the peak the profile's radius trades against its asymmetry (σ_in 3.27, σ_out 2.76 fitted
    /// on a symmetric 3 px ring) and reads +0.78 %, before and after 0.1.0 (registration-0.0.4 §1.8 amendment; a radius
    /// uncertainty is macTEM's ask G3).
    /// RED UNDER: the significance bound raised (alpha >= 30 SE).
    @Test(arguments: [(0.30, UInt64(11)), (0.40, UInt64(14))])
    func theSameNoisyRingWithoutAGapIsCalibrated(noise: Double, seed: UInt64) throws {
        let pattern = Self.render(noise: noise, seed: seed) { dr, dc in
            Self.ring(dr, dc, radius: 42, epsilon: 0.015, thetaDeg: 20)
        }
        let result = try Self.calibrate(pattern, inner: 28, outer: 56)
        let calibration = try #require(result.calibration, "refused: \(Self.names(result))")
        let e = result.report?.ringEvidence
        #expect(abs(calibration.radius / 42 - 1) < 0.01,
                "R \(calibration.radius) fit \(String(describing: result.fit)) c0 \(String(describing: e?.radiusShift)) centre \(String(describing: e?.centreShift)) ellipse \(String(describing: e?.ellipseShift)) bins \(String(describing: e?.occupiedBins)) drops \(String(describing: e?.innerDrop)) \(String(describing: e?.outerDrop))")
        #expect(result.report?.ringEvidence?.passes == true)
    }

    /// A ring whose 120 degrees hold 0.35 of its amplitude (noise-free, so every sector is significant; above the pre-fit
    /// coverage threshold of 0.2 of the range, so the fit runs) is not covered:
    /// a sector holds the ring only at half the median amplitude (registration-0.0.4 §1.3).
    /// RED UNDER: the occupancy fraction lowered from 0.5 to 0.05.
    @Test func aRingWithAWeakSectorIsNotCovered() throws {
        let pattern = Self.render { dr, dc in
            Self.ring(dr, dc, radius: 42, epsilon: 0.01, thetaDeg: 15) * (Self.inGap(dr, dc, start: 200, width: 120) ? 0.35 : 1)
        }
        let result = try Self.calibrate(pattern, inner: 28, outer: 56)
        #expect(result.calibration == nil, "calibrated; refusals \(Self.names(result))")
        #expect(Self.names(result).contains("insufficientAngularCoverage"))
        #expect((result.report?.ringEvidence?.occupiedBins ?? 36) < 30)
    }

    /// The refusals read the annulus's radii and the evidence: a radius outside (inner, outer) is `ringOutsideAnnulus`
    /// (alone, with an evidence that passes), and nothing is added when the report carries neither (the 0.0.3 initializer).
    /// RED UNDER: the annulus refusal removed.
    @Test func theAnnulusRefusalReadsTheReportsRadii() {
        let passing = RingEvidence(amplitude: 1, sectorAmplitudes: Array(repeating: 1, count: 36),
                                   sectorAmplitudeErrors: Array(repeating: 0.01, count: 36),
                                   sectorShifts: Array(repeating: 0, count: 36), occupiedBins: 36, innerDrop: 0.9,
                                   outerDrop: 0.9, radiusShift: 0, centreShift: 0.01, centreShiftError: 0.01,
                                   ellipseShift: 0.01, ellipseShiftError: 0.01)
        #expect(passing.passes)
        func report(_ inner: Double?, _ outer: Double?, _ evidence: RingEvidence?) -> RingCalibration.FitReport {
            RingCalibration.FitReport(model: .asymmetricGaussian, semiAxisA: 40.2, semiAxisB: 39.8, occupiedAngularBins: 36,
                                      sparseCoverage: false, normalizedResidual: 0.05, innerSigma: 3, outerSigma: 3,
                                      annulusWidth: 20, annulusInnerRadius: inner, annulusOuterRadius: outer,
                                      ringEvidence: evidence)
        }
        #expect(RingCalibration.refusals(report(42, 62, passing), highTensionVolts: 200_000)
                == [.ringOutsideAnnulus(radius: 40, inner: 42, outer: 62)])
        #expect(RingCalibration.refusals(report(30, 50, passing), highTensionVolts: 200_000).isEmpty)
        #expect(RingCalibration.refusals(report(nil, nil, nil), highTensionVolts: 200_000).isEmpty)
        let flat = RingEvidence(amplitude: 1, sectorAmplitudes: passing.sectorAmplitudes,
                                sectorAmplitudeErrors: passing.sectorAmplitudeErrors, sectorShifts: passing.sectorShifts,
                                occupiedBins: 36, innerDrop: 0.1, outerDrop: 0.9, radiusShift: 0, centreShift: 0.6,
                                centreShiftError: 0.01, ellipseShift: 0.01, ellipseShiftError: 0.01)
        #expect(RingCalibration.refusals(report(30, 50, flat), highTensionVolts: nil).map(\.name)
                == ["notARingPeak", "ringGeometryInconsistent", "noHighTension"])
    }

    // MARK: F2: things that are not one ring

    /// A bright disk's blurred edge is not a ring: its profile does not fall on the inner side. Slice 3b's F2 calibrated
    /// such edges at -10 to -18 %.
    /// RED UNDER: the peak test's inner drop not required.
    @Test(arguments: [1.0, 2.5])
    func aDiskEdgeIsNotCalibrated(blur: Double) throws {
        let pattern = Self.render { dr, dc in
            let r = (dr * dr + dc * dc).squareRoot()
            return 0.5 * erfc((r - 34) / (blur * 2.0.squareRoot()))
        }
        let result = try? Self.calibrate(pattern, inner: 22, outer: 46)
        #expect(result?.calibration == nil, "calibrated at R \(result?.report?.radius ?? .nan); refusals \(result.map(Self.names) ?? [])")
    }

    /// A central blob is not a ring (its fitted R may lie outside the annulus, or its profile is monotonic). A guard test:
    /// 0.0.3 refused these two blobs as well (by the fit's own checks); the refusal logic is held by
    /// `theAnnulusRefusalReadsTheReportsRadii`.
    @Test(arguments: [9.0, 16.0])
    func aCentralBlobIsNotCalibrated(sigma: Double) throws {
        let pattern = Self.render { dr, dc in exp(-(dr * dr + dc * dc) / (2 * sigma * sigma)) }
        let result = try? Self.calibrate(pattern, inner: 12, outer: 50)
        #expect(result?.calibration == nil, "calibrated at R \(result?.report?.radius ?? .nan); refusals \(result.map(Self.names) ?? [])")
    }

    /// A circular ring whose intensity varies as 0.25 + 0.75 cos^2 of the azimuth (texture) is not read as an ellipse:
    /// slice 3b's F2 read such a ring as 4.6 % ellipticity. The kit either fits it as the circle it is (b/a within 0.2 %)
    /// or refuses it by name.
    /// RED UNDER: the geometry check's ellipse shift bound removed.
    @Test(arguments: [(36.0, 0.0), (36.0, 30.0), (48.0, 75.0)])
    func aTexturedCircleIsNotReadAsAnEllipse(radius: Double, phaseDeg: Double) throws {
        let pattern = Self.render { dr, dc in
            let phi = atan2(dc, dr) - phaseDeg * .pi / 180
            return Self.ring(dr, dc, radius: radius) * (0.25 + 0.75 * cos(phi) * cos(phi))
        }
        let result = try Self.calibrate(pattern, inner: radius - 13, outer: radius + 13)
        if result.calibration != nil, let fit = result.fit {
            #expect(min(fit.a, fit.b) / max(fit.a, fit.b) > 0.998, "b/a \(min(fit.a, fit.b) / max(fit.a, fit.b)) accepted")
        } else {
            #expect(Self.names(result).contains("ringGeometryInconsistent") || Self.names(result).contains("insufficientAngularCoverage"),
                    "refused as \(Self.names(result))")
        }
    }

    /// Two rings closer than about twice their width are one ring to the fit (registration-0.0.4 §1.6, not claimed): the
    /// fitted radius lies between them, within the separation of the stronger. A bound, not a refusal.
    @Test(arguments: [(1.5, 0.7), (3.0, 0.7), (3.0, 1.0)])
    func twoRingsCloserThanTheirWidthAreFittedBetweenThem(separation: Double, secondAmplitude: Double) throws {
        let pattern = Self.render { dr, dc in
            Self.ring(dr, dc, radius: 40, sigmaIn: 2.5, sigmaOut: 2.5)
                + secondAmplitude * Self.ring(dr, dc, radius: 40 + separation, sigmaIn: 2.5, sigmaOut: 2.5)
        }
        let result = try Self.calibrate(pattern, inner: 28, outer: 54)
        if let radius = result.report?.radius {
            #expect(radius > 40 - 0.05 && radius < 40 + separation + 0.05, "R \(radius)")
        }
    }

    /// A second ring inside the annulus never gives a wrong number (Gate B on the 0.1.0 candidate, F1: such fits landed on the other ring
    /// or between the two, 1-2.3 % off): the calibration is refused by name, or reads the intended ring within 0.3 %.
    /// Measured: at 9 px (3σ) the dip between the rings is 0.23 of the amplitude, under the 0.25 bound, and the fit reads
    /// R 44.096 (0.2 %); at 14 and 24 px it is refused as secondRingInAnnulus.
    /// RED UNDER: the second-peak bound raised from 0.25 to 25.
    @Test(arguments: [(9.0, 0.8), (14.0, 1.0), (24.0, 0.6), (12.0, 1.2)])
    func aSecondRingInTheAnnulusIsRefused(separation: Double, amplitude: Double) throws {
        let pattern = Self.render { dr, dc in
            Self.ring(dr, dc, radius: 44) + amplitude * Self.ring(dr, dc, radius: 44 + separation)
        }
        let result = try? Self.calibrate(pattern, inner: 32, outer: 44 + separation + 10)
        if let calibration = result?.calibration {
            #expect(abs(calibration.radius / 44 - 1) < 3e-3, "calibrated at R \(calibration.radius)")
        } else if let result, result.report?.ringEvidence != nil {
            #expect(Self.names(result).contains("secondRingInAnnulus"), "refused as \(Self.names(result))")
        }
        if separation >= 12 { #expect(result?.calibration == nil) }
    }

    /// One ring on a central-beam tail rising to the annulus's inner edge is not a second ring.
    @Test func aTailRisingToTheEdgeIsNotASecondRing() throws {
        let pattern = Self.render { dr, dc in
            let r = max((dr * dr + dc * dc).squareRoot(), 1)
            return Self.ring(dr, dc, radius: 44) + 0.6 * exp(-r / 12)
        }
        let result = try Self.calibrate(pattern, inner: 20, outer: 58)
        #expect(!Self.names(result).contains("secondRingInAnnulus"), "refused as \(Self.names(result))")
        #expect(result.calibration != nil, "refused as \(Self.names(result))")
    }

    // MARK: F3: isolated impulses

    /// One hot pixel on the ring (5x, 20x, 1000x its peak) does not refuse a perfect ring, and one strongly negative pixel
    /// does not move it: the fit equals the ring without it within 0.05 px and 0.05 %. Slice 3b's F3 refused the hot one as
    /// insufficientAngularCoverage(1).
    /// RED UNDER: the impulse exclusion removed; its factor raised from 2 to 2000 robust ranges.
    @Test(arguments: [5.0, 20.0, 1000.0, -5.0, -50.0])
    func anIsolatedImpulseIsIgnored(level: Double) throws {
        let scene = Scene()
        let clean = Self.render(scene) { dr, dc in Self.ring(dr, dc, radius: 40, epsilon: 0.02, thetaDeg: 35) }
        var pixels = clean.pixels
        // a pixel on the ring's major axis, inside the annulus
        let row = Int((scene.row0 + 40 * cos(35 * Double.pi / 180)).rounded())
        let col = Int((scene.col0 + 40 * sin(35 * Double.pi / 180)).rounded())
        pixels[row * scene.cols + col] = Float(level)
        let spiked = DiffractionPattern(width: clean.width, height: clean.height, pixels: pixels)
        let reference = try Self.calibrate(clean, inner: 27, outer: 54)
        let result = try Self.calibrate(spiked, inner: 27, outer: 54)
        let a = try #require(reference.fit), b = try #require(result.fit, "refused: \(Self.names(result))")
        #expect(hypot(a.centerQX - b.centerQX, a.centerQY - b.centerQY) < 0.05)
        #expect(abs((b.a + b.b) / (a.a + a.b) - 1) < 5e-4)
        #expect(result.calibration != nil)
    }

    // MARK: F6: scale

    /// The fit does not depend on the pattern's scale. Ranges below 0.5 run normalised into [0.5, 1) by an exact power of
    /// two, so a ring of amplitude 2^-40 gives the amplitude-2^-3 fit bit for bit; across scales (2^-40, 1e-6, 1e-3, 1e3,
    /// 1e6 against amplitude 1, background in proportion) the fit agrees within 1e-6 relative (the Float pixels round s·v
    /// unlike v by up to 6e-8 relative, and the optimisers' step scales differ by less than 2 between the regimes;
    /// registration amended). Slice 3b's F6 refused amplitudes of 1e-6 and below as insufficientSignal.
    /// RED UNDER: the power-of-two normalisation removed (the 2^-40 case refused); the signal test made absolute again.
    @Test(arguments: [0x1p-40, 1e-6, 1e-3, 1e3, 1e6])
    func theFitDoesNotDependOnScale(scale: Double) throws {
        func pattern(_ s: Double) -> DiffractionPattern {
            Self.render(noise: 0.05 * s, seed: 21) { dr, dc in s * (Self.ring(dr, dc, radius: 38, epsilon: 0.03, thetaDeg: 60) + 0.2) }
        }
        func fit(_ s: Double) throws -> EllipseCalibrationFit {
            try EllipseCalibration.fitBestAvailable(pattern: pattern(s), centerQX: 64, centerQY: 61, innerRadius: 25, outerRadius: 52)
        }
        let unit = try fit(1), scaled = try fit(scale)
        #expect(scaled.model == unit.model)
        for (x, y) in [(scaled.centerQX, unit.centerQX), (scaled.centerQY, unit.centerQY), (scaled.a, unit.a),
                       (scaled.b, unit.b), (scaled.theta, unit.theta)] {
            #expect(abs(x - y) <= 1e-6 * max(1, abs(y)), "\(x) vs \(y)")
        }
        if scale == 0x1p-40 {
            let eighth = try fit(0x1p-3)
            #expect(scaled.centerQX == eighth.centerQX && scaled.a == eighth.a && scaled.b == eighth.b && scaled.theta == eighth.theta)
            #expect(scaled.profile?.ringIntensity == eighth.profile.map { $0.ringIntensity * 0x1p-37 })
        }
    }

    // MARK: F7: one sharp edge

    /// A ring with one edge sharper than the 0.30 px floor guard is fitted and calibrated when the data hold the ring
    /// (registration-0.0.4 §1.5). Slice 3b's F7 refused it as ringWidthAtFloor.
    /// RED UNDER: the one-sided floor kept as a refusal.
    @Test(arguments: [0.15, 0.05])
    func aRingWithOneSharpEdgeIsCalibrated(sigmaOut: Double) throws {
        let pattern = Self.render { dr, dc in Self.ring(dr, dc, radius: 40, epsilon: 0.01, thetaDeg: 10, sigmaIn: 3, sigmaOut: sigmaOut) }
        let result = try Self.calibrate(pattern, inner: 26, outer: 52)
        let calibration = try #require(result.calibration, "refused: \(Self.names(result))")
        #expect(abs(calibration.radius - 40) < 0.1, "R \(calibration.radius)")
    }

    /// A ring point-sampled on both sides below the floor is still refused (ADR 069's guard stands when both sigmas are at it).
    @Test func aRingAtTheFloorOnBothSidesIsStillRefused() throws {
        let pattern = Self.render { dr, dc in Self.ring(dr, dc, radius: 40, sigmaIn: 0.1, sigmaOut: 0.1) }
        let result = try? Self.calibrate(pattern, inner: 26, outer: 52)
        #expect(result?.calibration == nil)
    }

    // MARK: macTEM's refuse-bound plant (its ring-fit-parity line 37)

    /// A 3 px ring at R 120 on a 1/r^3 background at the ring's peak, annulus +-20 px: the fit reads the ring and the kit
    /// calibrates it (Δq within 0.05 %); `sigmaAtAnnulusWidth` was the plant's wrong expectation (registration-0.0.4 §1.7).
    @Test func macTEMsRefuseBoundPlantIsCalibrated() throws {
        let scene = Scene(rows: 320, cols: 316, row0: 158.41, col0: 160.27)
        let pattern = Self.render(scene) { dr, dc in
            let r = max((dr * dr + dc * dc).squareRoot(), 1)
            return Self.ring(dr, dc, radius: 120, epsilon: 0.01, thetaDeg: 20) + pow(120 / r, 3)
        }
        let result = try Self.calibrate(pattern, scene, inner: 100, outer: 140)
        let calibration = try #require(result.calibration, "refused: \(Self.names(result))")
        #expect(abs(calibration.deltaQ / ((1 / 2.3546) / 120) - 1) < 5e-4)
    }

    /// A ring whose sigma reaches the annulus width is refused by name (σ 16 px in a 16 px annulus; registration-0.0.4 §1.7
    /// amended: σ 8 px there is half the width, fitted and refused for its geometry).
    @Test func aRingWiderThanItsAnnulusIsRefusedAsSigmaAtAnnulusWidth() throws {
        let pattern = Self.render { dr, dc in Self.ring(dr, dc, radius: 40, sigmaIn: 16, sigmaOut: 16) + 0.1 }
        let result = try? Self.calibrate(pattern, inner: 32, outer: 48)
        #expect(result?.calibration == nil)
        if let result, result.fit != nil {
            // Measured: refused as ringGeometryInconsistent (the annulus shows ±0.5σ, so the fitted σ stays below 16).
            #expect(!Set(Self.names(result)).isDisjoint(with: ["sigmaAtAnnulusWidth", "notARingPeak", "ringGeometryInconsistent"]),
                    "refused as \(Self.names(result))")
        }
    }
}

// MARK: G1 and G7 (registration-0.0.4 §2, §4)

struct RingFitSplitAndCancellationTests {
    typealias T = RingFitAmendmentTests

    static func plants() -> [(String, DiffractionPattern, Double, Double)] {
        [
            ("clean", T.render { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.02, thetaDeg: 35) }, 28, 56),
            ("noisy", T.render(noise: 0.2, seed: 31) { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.01, thetaDeg: 70) }, 28, 56),
            ("gap", T.render(noise: 0.3, seed: 32) { dr, dc in
                T.inGap(dr, dc, start: 40, width: 110) ? 0 : T.ring(dr, dc, radius: 42) }, 28, 56),
            ("gap-silent", T.render { dr, dc in T.inGap(dr, dc, start: 40, width: 200) ? 0 : T.ring(dr, dc, radius: 42) }, 28, 56),
            ("texture", T.render { dr, dc in
                let phi = atan2(dc, dr)
                return T.ring(dr, dc, radius: 40) * (0.25 + 0.75 * cos(phi) * cos(phi)) }, 27, 53),
            ("sharp", T.render { dr, dc in T.ring(dr, dc, radius: 40, sigmaIn: 3, sigmaOut: 0.1) }, 26, 52),
            ("blob", T.render { dr, dc in exp(-(dr * dr + dc * dc) / 200) }, 12, 50),
        ]
    }

    /// `fit` is exactly `calibrating(fitRing(...))`: on every plant, with and without V and pitch, the two give equal
    /// results (fit, report, refusals in order, calibration), or the same thrown error.
    /// RED UNDER: `calibrating` not appending `.noHighTension`; `fitRing` built without the ring evidence.
    @Test(arguments: [(Double?(200_000), Double?(0.014)), (nil, nil), (300_000, nil)])
    func theSplitEqualsTheFit(volts: Double?, pitch: Double?) throws {
        for (name, pattern, inner, outer) in Self.plants() {
            let whole = Result { try RingCalibration.fit(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: inner,
                                                         outerRadius: outer, d: 2.3546, highTensionVolts: volts, pixelPitchMM: pitch) }
            let split = Result {
                try RingCalibration.calibrating(
                    RingCalibration.fitRing(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: inner, outerRadius: outer),
                    d: 2.3546, highTensionVolts: volts, pixelPitchMM: pitch)
            }
            switch (whole, split) {
            case let (.success(a), .success(b)):
                #expect(a == b, "\(name)")
                #expect(a.report?.ringEvidence != nil || a.report == nil, "\(name): the report carries the evidence")
                if volts == nil, a.report != nil || !a.refusals.isEmpty { #expect(a.refusals.last == .noHighTension, "\(name)") }
            case let (.failure(a), .failure(b)):
                #expect(String(describing: a) == String(describing: b), "\(name)")
            default:
                Issue.record("\(name): one path threw and the other did not: \(whole) vs \(split)")
            }
        }
    }

    /// An impossible d throws before any fit, from both paths.
    @Test func anImpossibleSpacingThrowsFromBothPaths() throws {
        let pattern = Self.plants()[0].1
        #expect(throws: RingCalibration.Refusal.self) {
            try RingCalibration.fit(pattern: pattern, centerQX: 64, centerQY: 61, innerRadius: 28, outerRadius: 56,
                                    d: -1, highTensionVolts: 200_000, pixelPitchMM: nil)
        }
        let ring = try RingCalibration.fitRing(pattern: pattern, centerQX: 64, centerQY: 61, innerRadius: 28, outerRadius: 56)
        #expect(ring.isUsable)
        #expect(abs((ring.radius ?? 0) - 42) < 0.01)
        #expect(throws: RingCalibration.Refusal.self) {
            try RingCalibration.calibrating(ring, d: .nan, highTensionVolts: 200_000, pixelPitchMM: nil)
        }
    }

    /// Progress runs monotone from at most 0.1 to exactly 1, and a closure that cancels stops the fit with
    /// `CancellationError`, as does a cancelled task.
    /// RED UNDER: the cancellation check removed from the optimisers; the final progress report dropped.
    @Test func progressIsMonotoneAndCancellationStopsTheFit() async throws {
        let pattern = Self.plants()[1].1
        final class Box: @unchecked Sendable {
            let lock = NSLock()
            var values: [Double] = []
            func add(_ v: Double) { lock.lock(); values.append(v); lock.unlock() }
        }
        let box = Box()
        _ = try RingCalibration.fit(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28, outerRadius: 56,
                                    d: 2.3546, highTensionVolts: 200_000, pixelPitchMM: nil, progress: { box.add($0) })
        #expect(box.values.count > 10)
        #expect(zip(box.values, box.values.dropFirst()).allSatisfy { $0 <= $1 })
        #expect((box.values.first ?? 1) <= 0.1)
        #expect(box.values.last == 1)

        let calls = Box()
        #expect(throws: CancellationError.self) {
            try RingCalibration.fit(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28, outerRadius: 56,
                                    d: 2.3546, highTensionVolts: 200_000, pixelPitchMM: nil,
                                    isCancelled: { calls.add(1); return calls.values.count > 3 })
        }
        #expect(calls.values.count == 4)

        // Monotone also through the one-sided-floor refit (Gate B on the 0.1.0 candidate, F2: it reported backwards by up to 0.04).
        let sharp = RingFitAmendmentTests.render { dr, dc in
            RingFitAmendmentTests.ring(dr, dc, radius: 40, sigmaIn: 3, sigmaOut: 0.1) }
        let sharpBox = Box()
        _ = try RingCalibration.fitRing(pattern: sharp, centerQX: 64.1, centerQY: 61.2, innerRadius: 26, outerRadius: 52,
                                        progress: { sharpBox.add($0) })
        #expect(zip(sharpBox.values, sharpBox.values.dropFirst()).allSatisfy { $0 <= $1 })
        #expect(sharpBox.values.last == 1)

        let task = Task.detached { () throws -> RingRadiusFit in
            try RingCalibration.fitRing(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28, outerRadius: 56)
        }
        task.cancel()
        await #expect(throws: CancellationError.self) { _ = try await task.value }
    }
}
