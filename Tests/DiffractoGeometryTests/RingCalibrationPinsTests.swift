//
//  RingCalibrationPinsTests.swift
//  DiffractoKit slice 3a (docs/registration-slice-3.md section 3): `ElectronWavelength` and `RingCalibration` against the
//  pins tools/slice-3/ring-pins.py writes (Fixtures/slice-3-ring-pins.json): the closed forms evaluated in 60-digit decimal
//  arithmetic, 320 cases, the refusal decision from the registered thresholds (68 plants), the Bragg limit (18 plants) and the
//  radius (10 plants). Band 1e-7 relative (D17, 0.9.1; it was 1e-12); a difference is a failure of the port, never a band to widen.
//  Ported from reference commit 9482831 tools/ring-calibration-parity/main.swift. The fixture is read by #filePath (the Fixtures folder
//  is excluded from the target).
//

import Testing
import Foundation
import DiffractoGeometry

private struct RingPins: Decodable {
    struct Constants: Decodable {
        let h_Js: Double, c_ms: Double, e_C: Double, m0_kg: Double, hc_over_e_eV_A: Double, m0c2_over_e_eV: Double
    }
    struct Wavelength: Decodable { let kv: Int, volts: Double, lambda_A: Double }
    struct Case: Decodable {
        let name: String
        let d_A: Double, kv: Int, volts: Double, r_px: Double, pitch_mm: Double
        let lambda_A: Double, two_theta_rad: Double, l_px: Double, dq_per_A_per_px: Double, l_mm: Double
    }
    struct Refusal: Decodable {
        let name: String, check: String, model: String
        let a_px: Double, b_px: Double, occupied_bins: Int, sparse_coverage: Bool, residual: Double
        let inner_sigma: Double?, outer_sigma: Double?, annulus_width_px: Double, volts: Double?
        let expect_names: [String], expect_values: [Double?], withholds_dq: Bool
    }
    struct Axes: Decodable { let name: String, a_px: Double, b_px: Double, radius_px: Double }
    struct Domain: Decodable {
        let name: String, volts: Double, d_A: Double, r_px: Double, pitch_mm: Double, inside: Bool, lambda_over_2d: Double
    }
    struct Truth: Decodable {
        let lambda_A: [String: Double]
        let au111_200kv_r150_p014: Au111
        let d_au111: Double
        let small_angle_excess_au111_200kv: Double, small_angle_excess_au422_80kv: Double
        struct Au111: Decodable { let two_theta_rad: Double, l_px: Double, dq_per_A_per_px: Double, l_mm: Double }
        enum CodingKeys: String, CodingKey {
            case lambda_A, d_au111, small_angle_excess_au111_200kv, small_angle_excess_au422_80kv
            case au111_200kv_r150_p014 = "au111_200kv_r150_p0.014"
        }
    }
    let numpy: String, tolerance: Double
    let constants: Constants
    let wavelengths: [Wavelength]
    let cases: [Case]
    let refusals: [Refusal]
    let radius: [Axes]
    let domain: [Domain]
    let truth: Truth

    static let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
        .appendingPathComponent("Fixtures/slice-3-ring-pins.json")

    static let shared: RingPins = {
        let decoder = JSONDecoder()
        decoder.nonConformingFloatDecodingStrategy = .convertFromString(positiveInfinity: "inf", negativeInfinity: "-inf", nan: "nan")
        do { return try decoder.decode(RingPins.self, from: Data(contentsOf: url)) } catch { fatalError("cannot read \(url.path): \(error)") }
    }()
}

private func relative(_ x: Double, _ truth: Double) -> Double { abs(x - truth) / abs(truth) }

/// Equal, or within the relative tolerance (a reported 0 must equal 0, where the relative error is 0/0).
private func agrees(_ x: Double, _ truth: Double, _ tol: Double) -> Bool { x == truth || relative(x, truth) <= tol }

/// The number a refusal reports, the same one ring-pins.py's `decide` lists (nil where it reports none).
private func reported(_ refusal: RingCalibration.Refusal) -> Double? {
    switch refusal {
    case .insufficientAngularCoverage(let bins): return Double(bins)
    case .sigmaAtAnnulusWidth(let sigma, _): return sigma
    case .excessiveResidual(let residual): return residual
    case .excessiveEllipticity(let ellipticity): return ellipticity
    case .braggLimit(let ratio): return ratio
    case .modelNotAsymmetricGaussian, .noHighTension, .invalidInput: return nil
    // DiffractoKit 0.1.0's refusals come from a fit's ring evidence; ring-pins.py's closed-form cases never meet them.
    case .ringOutsideAnnulus(let radius, _, _): return radius
    case .notARingPeak, .ringGeometryInconsistent, .secondRingInAnnulus: return nil
    }
}

struct RingCalibrationPinsTests {
    private let pins = RingPins.shared
    /// The registration's band; the fixture's own tolerance must not be looser.
    private let band = 1e-7  // D17: assert_allclose default (rtol 1e-7); was 1e-12, measured worst is far below

    @Test func testThePinsAreTheRegisteredSetAndTheirBandIsNotLoosened() {
        #expect(pins.tolerance <= band)
        #expect(pins.cases.count == 320)
        #expect(pins.refusals.count == 68)
        #expect(pins.domain.count == 18)
        #expect(pins.radius.count == 10)
        #expect(pins.wavelengths.map(\.kv) == [60, 80, 120, 200, 300])
    }

    @Test func testTheCODATAConstantsEqualThePinsAndTheDerivedEnergiesFollow() {
        let k = pins.constants
        #expect(ElectronWavelength.planckConstant == k.h_Js)
        #expect(ElectronWavelength.speedOfLight == k.c_ms)
        #expect(ElectronWavelength.elementaryCharge == k.e_C)
        #expect(ElectronWavelength.electronMass == k.m0_kg)
        #expect(relative(ElectronWavelength.planckTimesLightEVAngstroms, k.hc_over_e_eV_A) <= 1e-14)
        #expect(relative(ElectronWavelength.electronRestEnergyEV, k.m0c2_over_e_eV) <= 1e-14)
        // RingCalibration's constants are the same numbers.
        #expect(RingCalibration.planckConstant == ElectronWavelength.planckConstant)
        #expect(RingCalibration.electronRestEnergyEV == ElectronWavelength.electronRestEnergyEV)
    }

    @Test func testTheWavelengthAtEachPinnedVoltageIsWithinTheBand() {
        for w in pins.wavelengths {
            let lambda = ElectronWavelength.angstroms(volts: w.volts)
            #expect(relative(lambda, w.lambda_A) <= band, "\(w.kv) kV: \(lambda) vs \(w.lambda_A)")
            // RingCalibration delegates: the same bits.
            #expect(RingCalibration.wavelengthAngstroms(volts: w.volts) == lambda, "\(w.kv) kV delegates")
        }
    }

    @Test func testEachClosedFormGivenThePinnedInputsIsWithinTheBand() {
        for c in pins.cases {
            let twoTheta = RingCalibration.twoTheta(wavelengthAngstroms: c.lambda_A, d: c.d_A)
            let length = RingCalibration.cameraLengthPixels(radius: c.r_px, twoTheta: c.two_theta_rad)
            let q = RingCalibration.qPerPixel(d: c.d_A, radius: c.r_px)
            let mm = RingCalibration.cameraLengthMillimetres(pixelPitchMM: c.pitch_mm, cameraLengthPixels: c.l_px)
            #expect(relative(twoTheta, c.two_theta_rad) <= band, "\(c.name) 2θ \(twoTheta) vs \(c.two_theta_rad)")
            #expect(relative(length, c.l_px) <= band, "\(c.name) L_px \(length) vs \(c.l_px)")
            #expect(relative(q, c.dq_per_A_per_px) <= band, "\(c.name) Δq \(q) vs \(c.dq_per_A_per_px)")
            #expect(relative(mm, c.l_mm) <= band, "\(c.name) L_mm \(mm) vs \(c.l_mm)")
        }
    }

    @Test func testCalibrateEndToEndAtEveryCaseWithAndWithoutPitchAndHighTension() throws {
        for c in pins.cases {
            let full = try RingCalibration.calibrate(radius: c.r_px, d: c.d_A, highTensionVolts: c.volts, pixelPitchMM: c.pitch_mm)
            func close(_ what: String, _ x: Double?, _ truth: Double) {
                guard let x else { Issue.record("\(c.name): \(what) absent"); return }
                #expect(relative(x, truth) <= band, "\(c.name) \(what): \(x) vs \(truth)")
            }
            close("λ", full.wavelengthAngstroms, c.lambda_A)
            close("2θ", full.twoTheta, c.two_theta_rad)
            close("L_px", full.lengthPixels, c.l_px)
            close("Δq", full.deltaQ, c.dq_per_A_per_px)
            close("L_mm", full.lengthMillimetres, c.l_mm)
            #expect(full.highTensionVolts == c.volts, "\(c.name): V is returned")
            #expect(full.lengthRefusal == nil, "\(c.name): L not refused")
            #expect(full.d == c.d_A && full.radius == c.r_px, "\(c.name): d and R echoed")

            let noPitch = try RingCalibration.calibrate(radius: c.r_px, d: c.d_A, highTensionVolts: c.volts, pixelPitchMM: nil)
            close("L_px without a pitch", noPitch.lengthPixels, c.l_px)
            #expect(noPitch.lengthMillimetres == nil, "\(c.name): no pitch, no L_mm")

            let noV = try RingCalibration.calibrate(radius: c.r_px, d: c.d_A, highTensionVolts: nil, pixelPitchMM: c.pitch_mm)
            close("Δq without V", noV.deltaQ, c.dq_per_A_per_px)
            #expect(noV.lengthRefusal == .noHighTension, "\(c.name): V is never defaulted")
            #expect(noV.lengthPixels == nil && noV.lengthMillimetres == nil && noV.wavelengthAngstroms == nil
                    && noV.twoTheta == nil && noV.highTensionVolts == nil, "\(c.name): nothing of L without V")
        }
    }

    @Test func testEveryRefusalPlantGivesTheRegisteredNamesInOrderWithTheirNumbers() {
        for r in pins.refusals {
            let fit = RingCalibration.FitReport(
                model: r.model == "asymmetricGaussian" ? .asymmetricGaussian : .conic, semiAxisA: r.a_px, semiAxisB: r.b_px,
                occupiedAngularBins: r.occupied_bins, sparseCoverage: r.sparse_coverage, normalizedResidual: r.residual,
                innerSigma: r.inner_sigma, outerSigma: r.outer_sigma, annulusWidth: r.annulus_width_px)
            let got = RingCalibration.refusals(fit, highTensionVolts: r.volts)
            #expect(got.map(\.name) == r.expect_names,
                    "\(r.name) \(r.check): refused as \(got.map(\.name)), the registered thresholds say \(r.expect_names)")
            if got.map(\.name) == r.expect_names {
                for (refusal, want) in zip(got, r.expect_values) {
                    let have = reported(refusal)
                    let ok = (have == nil) == (want == nil) && (have == nil || agrees(have!, want!, band))
                    #expect(ok, "\(r.name): \(refusal.name) reports \(String(describing: have)), expected \(String(describing: want))")
                    if case .sigmaAtAnnulusWidth(_, let width) = refusal {
                        #expect(width == r.annulus_width_px, "\(r.name): reports the annulus width \(width)")
                    }
                }
            }
            #expect(got.contains(where: \.withholdsDeltaQ) == r.withholds_dq, "\(r.name): withholdsDeltaQ")
        }
    }

    @Test func testTheRadiusIsNumpysMeanOfTheTwoSemiAxesOnEveryPlant() {
        for p in pins.radius {
            let fit = RingCalibration.FitReport(model: .asymmetricGaussian, semiAxisA: p.a_px, semiAxisB: p.b_px,
                                                occupiedAngularBins: 36, sparseCoverage: false, normalizedResidual: 0.05,
                                                innerSigma: 3, outerSigma: 3, annulusWidth: 40)
            #expect(relative(fit.radius, p.radius_px) <= band, "\(p.name): R \(fit.radius) vs \(p.radius_px)")
        }
    }

    @Test func testTheBraggLimitIsRefusedByNameOutsideAndNotInsideOnBothSidesOfTheEdge() throws {
        for b in pins.domain {
            let lambda = ElectronWavelength.angstroms(volts: b.volts)
            let angle = RingCalibration.twoTheta(wavelengthAngstroms: lambda, d: b.d_A)
            #expect(angle.isFinite == b.inside, "\(b.name): 2θ \(angle), the ring is \(b.inside ? "inside" : "outside") the limit")
            let result = try RingCalibration.calibrate(radius: b.r_px, d: b.d_A, highTensionVolts: b.volts, pixelPitchMM: b.pitch_mm)
            #expect(result.deltaQ == RingCalibration.qPerPixel(d: b.d_A, radius: b.r_px) && result.deltaQ.isFinite,
                    "\(b.name): Δq is given beside a Bragg refusal")
            if b.inside {
                #expect(result.lengthRefusal == nil, "\(b.name): refused \(String(describing: result.lengthRefusal))")
                #expect(result.lengthPixels.map { $0.isFinite && $0 > 0 } ?? false, "\(b.name): L_px \(String(describing: result.lengthPixels))")
            } else if case .braggLimit(let ratio)? = result.lengthRefusal {
                #expect(relative(ratio, b.lambda_over_2d) <= band, "\(b.name): λ/2d \(ratio) vs \(b.lambda_over_2d)")
                #expect(result.lengthPixels == nil && result.lengthMillimetres == nil && result.twoTheta == nil, "\(b.name): no L or 2θ")
                #expect(result.highTensionVolts == b.volts && result.wavelengthAngstroms != nil, "\(b.name): V and λ beside the refusal")
            } else {
                Issue.record("\(b.name): not refused as braggLimit: \(String(describing: result.lengthRefusal))")
            }
        }
    }

    @Test func testTheRegisteredTruthLiteralsAreWhatThePinsRecomputeAndTheCodeGives() throws {
        let t = pins.truth
        for (kv, want) in t.lambda_A {
            let volts = (Double(kv) ?? 0) * 1000
            #expect(relative(ElectronWavelength.angstroms(volts: volts), want) <= band, "λ at \(kv) kV")
        }
        let r = try RingCalibration.calibrate(radius: 150, d: t.d_au111, highTensionVolts: 200_000, pixelPitchMM: 0.014)
        let au = t.au111_200kv_r150_p014
        #expect(relative(try #require(r.twoTheta), au.two_theta_rad) <= band)
        #expect(relative(try #require(r.lengthPixels), au.l_px) <= band)
        #expect(relative(r.deltaQ, au.dq_per_A_per_px) <= band)
        #expect(relative(try #require(r.lengthMillimetres), au.l_mm) <= band)
        let excess = 150 * t.d_au111 / ElectronWavelength.angstroms(volts: 200_000) / (try #require(r.lengthPixels)) - 1
        #expect(abs(excess - t.small_angle_excess_au111_200kv) <= 1e-11)
    }
}
