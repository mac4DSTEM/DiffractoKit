//
//  RingCalibrationQTests.swift
//  DiffractoKit 0.1.0 (docs/registration-0.0.4.md section 3, G2): `RingCalibration.q(radius:cameraLengthPixels:wavelength
//  Angstroms:)` and its inverse `radius(q:cameraLengthPixels:wavelengthAngstroms:)` against the pins tools/kit-0.0.4/q-pins.py
//  writes (Fixtures/q-pins.json): the closed forms q = 2 sin(atan(R/L)/2)/λ and R = L tan(2 asin(qλ/2)) in 60-digit decimal
//  arithmetic (two independent evaluations that agree to 1e-50), on 100 forward and 120 inverse cases; the round trip; the
//  identity q(R, L_px, λ) = 1/d on slice 3a's 320 cases; and the NaN plants. Band 1e-12 relative, the registration's; a
//  difference is a failure of the port, never a band to widen. The fixture is read by #filePath.
//

import Testing
import Foundation
import DiffractoGeometry

private struct QPins: Decodable {
    struct Forward: Decodable {
        let name: String
        let kv: Int, volts: Double, lambda_A: Double, l_px: Double, r_px: Double, q: Double, numpy_q: Double
    }
    struct Inverse: Decodable {
        let name: String, kind: String
        let kv: Int, volts: Double, lambda_A: Double, l_px: Double, q: Double, r_px: Double, numpy_r: Double
        let start_r_px: Double?
    }
    struct Identity: Decodable {
        let name: String
        let d_A: Double, r_px: Double, l_px: Double, lambda_A: Double, inv_d: Double, q: Double
    }
    struct Plant: Decodable {
        let name: String, fn: String, note: String
        let a: Double, l_px: Double, lambda_A: Double
    }
    struct Inside: Decodable { let name: String, q: Double, l_px: Double, lambda_A: Double }
    let numpy: String, tolerance: Double, limit_float: Double
    let forward: [Forward]
    let inverse: [Inverse]
    let identity: [Identity]
    let nan: [Plant]
    let just_inside: Inside

    static let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
        .appendingPathComponent("Fixtures/q-pins.json")

    static let shared: QPins = {
        let decoder = JSONDecoder()
        decoder.nonConformingFloatDecodingStrategy = .convertFromString(positiveInfinity: "inf", negativeInfinity: "-inf", nan: "nan")
        do { return try decoder.decode(QPins.self, from: Data(contentsOf: url)) } catch { fatalError("cannot read \(url.path): \(error)") }
    }()
}

private func relative(_ x: Double, _ truth: Double) -> Double { abs(x - truth) / abs(truth) }

struct RingCalibrationQTests {
    private let pins = QPins.shared
    /// The registration's band; the fixture's own tolerance must not be looser.
    private let band = 1e-12

    @Test func testThePinsAreTheRegisteredSetAndTheirBandIsNotLoosened() {
        #expect(pins.tolerance <= band)
        #expect(pins.forward.count == 100)           // 5 V x 4 L x 5 R
        #expect(pins.inverse.count == 120)           // the q those give, and 20 at 0.99 of the Bragg limit
        #expect(pins.inverse.filter { $0.kind == "bragg99" }.count == 20)
        #expect(pins.identity.count == 320)
        #expect(pins.nan.count == 51)
        #expect(Set(pins.forward.map(\.kv)) == [60, 80, 120, 200, 300])
        #expect(Set(pins.forward.map(\.l_px)) == [500, 2000, 14082.307097295722, 1e5])
        #expect(Set(pins.forward.map(\.r_px)) == [1, 40, 150, 600, 2000])
        // The Swift limit is the float64 the pins plant around.
        #expect(RingCalibration.maximumLambdaOver2d == pins.limit_float)
    }

    @Test func testQFromARadiusIsWithinTheBandOfTheDecimalTruth() {
        for c in pins.forward {
            let q = RingCalibration.q(radius: c.r_px, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
            #expect(relative(q, c.q) <= band, "\(c.name): q \(q) vs \(c.q)")
        }
    }

    @Test func testTheRadiusFromAQIsWithinTheBandOfTheDecimalTruth() {
        for c in pins.inverse {
            let r = RingCalibration.radius(q: c.q, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
            #expect(relative(r, c.r_px) <= band, "\(c.name) (\(c.kind)): R \(r) vs \(c.r_px)")
        }
    }

    @Test func testTheRoundTripRadiusOfQOfRIsTheRadiusToTheBand() {
        for c in pins.forward {
            let q = RingCalibration.q(radius: c.r_px, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
            let r = RingCalibration.radius(q: q, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
            #expect(relative(r, c.r_px) <= band, "\(c.name): R \(r) after the round trip")
        }
    }

    @Test func testQOfTheCalibratingRingsOwnRadiusLengthAndWavelengthIsOneOverD() {
        for c in pins.identity {
            let q = RingCalibration.q(radius: c.r_px, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
            #expect(relative(q, c.inv_d) <= band, "\(c.name): q \(q) vs 1/d \(c.inv_d)")
            #expect(relative(q, c.q) <= band, "\(c.name): q \(q) vs the decimal q \(c.q)")
        }
    }

    @Test func testTheIdentityHoldsThroughCalibrateEndToEnd() throws {
        // The kit's own L_px, from calibrate, closes the identity as well (not only the pinned L_px).
        for c in pins.identity where c.name.hasSuffix("-p0.014") {
            let kv = Double(c.name.split(separator: "-")[1].dropLast(2))!
            let cal = try RingCalibration.calibrate(radius: c.r_px, d: c.d_A, highTensionVolts: kv * 1000)
            let length = try #require(cal.lengthPixels)
            let lambda = try #require(cal.wavelengthAngstroms)
            let q = RingCalibration.q(radius: c.r_px, cameraLengthPixels: length, wavelengthAngstroms: lambda)
            #expect(relative(q, 1 / c.d_A) <= band, "\(c.name): q \(q) vs 1/d \(1 / c.d_A)")
        }
    }

    @Test func testEveryNaNPlantIsNaNNeverAnInfinityOrANumber() {
        for p in pins.nan {
            let x: Double
            switch p.fn {
            case "q": x = RingCalibration.q(radius: p.a, cameraLengthPixels: p.l_px, wavelengthAngstroms: p.lambda_A)
            case "radius": x = RingCalibration.radius(q: p.a, cameraLengthPixels: p.l_px, wavelengthAngstroms: p.lambda_A)
            default: Issue.record("\(p.name): unknown function \(p.fn)"); continue
            }
            #expect(x.isNaN, "\(p.name) (\(p.note)): \(x)")
        }
        // Both functions are exercised, and the Bragg limit and the overflow plants are among them.
        #expect(pins.nan.contains { $0.fn == "q" } && pins.nan.contains { $0.fn == "radius" })
        #expect(pins.nan.contains { $0.name == "radius-bragg-at-limit" } && pins.nan.contains { $0.name.contains("overflow") })
    }

    @Test func testJustInsideTheBraggLimitIsAFinitePositiveNumber() {
        let c = pins.just_inside
        let r = RingCalibration.radius(q: c.q, cameraLengthPixels: c.l_px, wavelengthAngstroms: c.lambda_A)
        #expect(r.isFinite && r > 0, "\(c.name): \(r)")
    }

    @Test func testBothFunctionsAreMonotoneInTheirFirstArgument() {
        // A cheap structural check beside the pins: q grows with R, and R grows with q, at 200 kV and L = 2000 px.
        let lambda = 0.025079340450548003
        let radii: [Double] = [1, 10, 100, 1000, 10_000]
        let qs = radii.map { RingCalibration.q(radius: $0, cameraLengthPixels: 2000, wavelengthAngstroms: lambda) }
        #expect(zip(qs, qs.dropFirst()).allSatisfy { $0 < $1 })
        let back = qs.map { RingCalibration.radius(q: $0, cameraLengthPixels: 2000, wavelengthAngstroms: lambda) }
        #expect(zip(back, back.dropFirst()).allSatisfy { $0 < $1 })
    }
}
