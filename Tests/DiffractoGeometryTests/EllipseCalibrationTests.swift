//
//  EllipseCalibrationTests.swift
//  `EllipseCalibration` on small clean patterns planted here, in
//  Swift, with their truth known: the constants, the input guards, the annulus boundaries, the coverage bound and the "fit
//  anyway" path, the one-ring check, the 0.30 px ring-width floor guard (ADR 069), the profile residual cap, a clean circle's
//  radius on both paths, and the row / column convention on a non-square pattern. The noisy plants of the reference implementations are
//  EllipseNoisyRingTests, the stop-rule rings EllipseStopRuleTests.
//
//  Written for the kit (the seed has these checks only as tools/ellipse-calibration-test, a scientific harness); the numbers
//  marked "measured" were measured on reference commit 7d6fdbd1's own EllipseCalibration compiled -O (and equalled the kit's in
//  every bit when they were measured).
//
//  Every test names under "RED UNDER" the production mutation that turns it red.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import DiffractoGeometry

struct EllipseCalibrationTests {

    // MARK: Planted patterns

    /// A clean Gaussian ring on a rows x columns pattern: I = exp(-((rho - 1) R)^2 / 2 s^2), rho the conic's elliptical radius in
    /// py4DSTEM's (A, B, C) convention with qx the row and a along theta, R = (a + b)/2.
    private func ring(rows: Int, columns: Int, centerQX: Double, centerQY: Double, a: Double, b: Double, theta: Double,
                      sigma: Double) -> DiffractionPattern {
        let sin2 = sin(theta) * sin(theta), cos2 = cos(theta) * cos(theta)
        let A = sin2 / (b * b) + cos2 / (a * a)
        let C = cos2 / (b * b) + sin2 / (a * a)
        let B = 2 * (b * b - a * a) * sin(theta) * cos(theta) / (a * a * b * b)
        let mean = (a + b) / 2
        var pixels = [Float](repeating: 0, count: rows * columns)
        for qx in 0..<rows {
            for qy in 0..<columns {
                let dx = Double(qx) - centerQX, dy = Double(qy) - centerQY
                let rho = max(0, A * dx * dx + B * dx * dy + C * dy * dy).squareRoot()
                let distance = (rho - 1) * mean
                pixels[qx * columns + qy] = Float(exp(-(distance * distance) / (2 * sigma * sigma)))
            }
        }
        return DiffractionPattern(width: columns, height: rows, pixels: pixels)
    }

    /// A circle of radius R in each azimuth sector (degrees from the row axis towards the column axis, [from, to)): single arcs
    /// or "grains" at different radii about one centre. Zero outside the sectors.
    private func arcs(size: Int, center: Double, sigma: Double, sectors: [(from: Double, to: Double, radius: Double)]) -> DiffractionPattern {
        var pixels = [Float](repeating: 0, count: size * size)
        for qx in 0..<size {
            for qy in 0..<size {
                let dx = Double(qx) - center, dy = Double(qy) - center
                var angle = atan2(dy, dx) * 180 / .pi
                if angle < 0 { angle += 360 }
                let r = hypot(dx, dy)
                var value = 0.0
                for sector in sectors where angle >= sector.from && angle < sector.to {
                    value = exp(-(r - sector.radius) * (r - sector.radius) / (2 * sigma * sigma))
                }
                pixels[qx * size + qy] = Float(value)
            }
        }
        return DiffractionPattern(width: size, height: size, pixels: pixels)
    }

    // MARK: The constants

    /// 36 azimuthal bins, 30 needed by default (five sixths), 12 on the anyway path (one third); the calibration's own bounds are
    /// the same numbers (RingCalibrationTests holds them to the literals; this is the equality its header promised).
    /// RED UNDER: any of the three constants changed; the ring fit's bin count different from `RingCalibration`'s.
    @Test func testTheBinConstantsAreTheDocumentedOnesAndEqualTheCalibrationsOwn() {
        #expect(EllipseCalibration.angularBinCount == 36)
        #expect(EllipseCalibration.degeneracyBoundBins == 30)
        #expect(EllipseCalibration.sparseFloorBins == 12)
        #expect(EllipseCalibration.angularBinCount == RingCalibration.angularBinCount)
        #expect(EllipseCalibration.degeneracyBoundBins == RingCalibration.minimumOccupiedAngularBins)
    }

    // MARK: Input guards

    /// Every guard of `fit1D` and `fitAmorphousRing` throws `invalidInput` before any work; none of the messages names an
    /// application.
    /// RED UNDER: a guard removed (a NaN centre or an inner radius below zero would then be fitted or crash).
    @Test func testImpossibleInputsAreRefusedByName() throws {
        let good = ring(rows: 40, columns: 40, centerQX: 20.3, centerQY: 19.6, a: 10, b: 10, theta: 0, sigma: 1.5)
        let tiny = DiffractionPattern(width: 2, height: 40, pixels: [Float](repeating: 1, count: 80))
        let torn = DiffractionPattern(width: 40, height: 40, pixels: [Float](repeating: 1, count: 1599))
        func refused(_ body: () throws -> EllipseCalibrationFit) -> Bool {
            do { _ = try body() } catch let caught as EllipseCalibration.FitError { if case .invalidInput = caught { return true }; return false } catch { return false }
            return false
        }
        #expect(refused { try EllipseCalibration.fit1D(pattern: tiny, centerQX: 20, centerQY: 1, innerRadius: 1, outerRadius: 5) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: torn, centerQX: 20, centerQY: 20, innerRadius: 5, outerRadius: 9) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: .nan, centerQY: 20, innerRadius: 5, outerRadius: 9) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: 20, centerQY: .infinity, innerRadius: 5, outerRadius: 9) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: 20, centerQY: 20, innerRadius: -1, outerRadius: 9) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: 20, centerQY: 20, innerRadius: 9, outerRadius: 9) })
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: 20, centerQY: 20, innerRadius: 5, outerRadius: 9,
                                                       maximumResidual: 0) })
        // fewer than 20 pixels in the annulus
        #expect(refused { try EllipseCalibration.fit1D(pattern: good, centerQX: 20, centerQY: 20, innerRadius: 1, outerRadius: 2.1) })
        let seed = EllipseCalibrationFit(centerQX: 20, centerQY: 20, a: 10, b: 10, theta: 0, normalizedResidual: 0,
                                         conicResidual: 0, sampleCount: 1, occupiedAngularBins: 36, model: .conic, profile: nil,
                                         profileFallbackReason: nil)
        #expect(refused { try EllipseCalibration.fitAmorphousRing(pattern: tiny, initial: seed, centerQX: 20, centerQY: 1,
                                                                  innerRadius: 1, outerRadius: 5) })
        #expect(refused { try EllipseCalibration.fitAmorphousRing(pattern: good, initial: seed, centerQX: 20, centerQY: 20,
                                                                  innerRadius: 9, outerRadius: 5) })
        #expect(refused { try EllipseCalibration.fitAmorphousRing(pattern: good, initial: seed, centerQX: 20, centerQY: 20,
                                                                  innerRadius: 5, outerRadius: 9, maximumResidual: -1) })
        let errors: [EllipseCalibration.FitError] = [
            .invalidInput("x"), .insufficientSignal, .insufficientAngularCoverage(5), .moreThanOneRing(minRadius: 1, maxRadius: 2),
            .didNotConverge, .invalidEllipse, .excessiveResidual(0.3), .ringWidthAtFloor,
        ]
        for error in errors {
            let text = try #require(error.errorDescription)
            for word in ["mac4DSTEM", "macTEM", "ACOM", "SpectroKit"] { #expect(!text.contains(word), "\(text)") }
        }
    }

    /// A flat annulus has no ring: `insufficientSignal` from the conic, from the profile and from the chain.
    /// RED UNDER: the dynamic-range guard removed.
    @Test func testAFlatAnnulusHasNoSignal() throws {
        let flat = DiffractionPattern(width: 48, height: 48, pixels: [Float](repeating: 0, count: 48 * 48))
        #expect(throws: EllipseCalibration.FitError.insufficientSignal) {
            try EllipseCalibration.fit1D(pattern: flat, centerQX: 24, centerQY: 24, innerRadius: 8, outerRadius: 18)
        }
        #expect(throws: EllipseCalibration.FitError.insufficientSignal) {
            try EllipseCalibration.fitBestAvailable(pattern: flat, centerQX: 24, centerQY: 24, innerRadius: 8, outerRadius: 18)
        }
        let seed = EllipseCalibrationFit(centerQX: 24, centerQY: 24, a: 12, b: 12, theta: 0, normalizedResidual: 0,
                                         conicResidual: 0, sampleCount: 1, occupiedAngularBins: 36, model: .conic, profile: nil,
                                         profileFallbackReason: nil)
        #expect(throws: EllipseCalibration.FitError.insufficientSignal) {
            try EllipseCalibration.fitAmorphousRing(pattern: flat, initial: seed, centerQX: 24, centerQY: 24, innerRadius: 8,
                                                    outerRadius: 18)
        }
    }

    // MARK: The annulus boundaries

    /// The conic's annulus is `inner < r <= outer` and the profile's `inner < r < outer` (kept for bit identity with the reference implementations).
    /// On a lattice centred on a pixel the radii 10 and 20 are hit exactly by pixels such as (6, 8) and
    /// (12, 16), so the counts differ by exactly the pixels on the boundaries; the expected counts are worked in integers.
    /// RED UNDER: `<=` for `<` or `<` for `<=` at either outer boundary; `>=` for `>` at the inner one.
    @Test func testTheAnnulusBoundariesAreInnerExcludedAndOuterIncludedForTheConicOnly() throws {
        let size = 65, centre = 32
        // the plant sits a fraction off the lattice centre; the annulus is cut about the START centre, which is the lattice point
        let pattern = ring(rows: size, columns: size, centerQX: 32.3, centerQY: 31.8, a: 15, b: 15, theta: 0, sigma: 1.5)
        var halfOpen = 0, open = 0, onInner = 0, onOuter = 0
        for qx in 0..<size {
            for qy in 0..<size {
                let d2 = (qx - centre) * (qx - centre) + (qy - centre) * (qy - centre)
                if d2 > 100 && d2 <= 400 { halfOpen += 1 }
                if d2 > 100 && d2 < 400 { open += 1 }
                if d2 == 100 { onInner += 1 }
                if d2 == 400 { onOuter += 1 }
            }
        }
        #expect(onInner == 12 && onOuter == 12, "the lattice must put pixels on both boundaries: \(onInner), \(onOuter)")
        #expect(halfOpen - open == onOuter)
        let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 32, centerQY: 32, innerRadius: 10, outerRadius: 20)
        #expect(conic.sampleCount == halfOpen, "conic: \(conic.sampleCount) against inner < r <= outer = \(halfOpen)")
        let best = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 32, centerQY: 32, innerRadius: 10,
                                                           outerRadius: 20)
        #expect(best.model == .asymmetricGaussian)
        #expect(best.sampleCount == open, "profile: \(best.sampleCount) against inner < r < outer = \(open)")
        #expect(abs((best.a + best.b) / 2 - 15) < 1e-3)
    }

    // MARK: Coverage

    /// A ring with a 120 degree arc missing (a plant of the former apps-parity set: 128 px, R 40, eps 2 %, theta 30 degrees) holds 25 of
    /// 36 bins (measured): refused by default with that count, fitted on the "anyway" path and MARKED `sparseCoverage`, with its centre
    /// 0.89 px from the plant's (measured; a gapped conic carries a centre shift, which is why the one-ring check and the caller's
    /// assertion exist) and b/a within 0.0012 of the plant's.
    /// RED UNDER: the 30-bin bound lowered to 24 or lower (the default path then fits it); `sparseCoverage` not set below the
    /// bound; the one-third floor changed (see the next test).
    @Test func testAGappedRingIsRefusedByDefaultAndFittedMarkedOnTheAnywayPath() throws {
        let plant = EllipsePlant(size: 128, radius: 40, eps: 0.02, thetaDegrees: 30, gap: 120, noise: 0, seed: 201)
        let pattern = plant.render()
        let start = (plant.centerQX + 1.3, plant.centerQY - 0.9)
        #expect(throws: EllipseCalibration.FitError.insufficientAngularCoverage(25)) {
            try EllipseCalibration.fit1D(pattern: pattern, centerQX: start.0, centerQY: start.1, innerRadius: plant.b - 20,
                                         outerRadius: plant.a + 20)
        }
        let fit = try EllipseCalibration.fit1D(pattern: pattern, centerQX: start.0, centerQY: start.1, innerRadius: plant.b - 20,
                                               outerRadius: plant.a + 20, maximumResidual: 0.4, acceptSparseCoverage: true)
        #expect(fit.sparseCoverage)
        #expect(fit.occupiedAngularBins == 25)
        #expect(fit.model == .conic)
        #expect(hypot(fit.centerQX - plant.centerQX, fit.centerQY - plant.centerQY) <= 1.0)
        #expect(abs(fit.b / fit.a - plant.b / plant.a) <= 0.01)
        // full coverage makes the flag inert
        let full = EllipsePlant(size: 128, radius: 40, eps: 0.02, thetaDegrees: 30, gap: 0, noise: 0, seed: 201)
        let whole = try EllipseCalibration.fit1D(pattern: full.render(), centerQX: full.centerQX + 1.3, centerQY: full.centerQY - 0.9,
                                                 innerRadius: full.b - 20, outerRadius: full.a + 20, acceptSparseCoverage: true)
        #expect(!whole.sparseCoverage)
        #expect(whole.occupiedAngularBins >= EllipseCalibration.degeneracyBoundBins)
    }

    /// Under one third of the bins (a 270 degree arc missing: 9 of 36, measured) even the anyway path refuses, with the count.
    /// RED UNDER: the one-third floor lowered (the anyway path would try to fit nine bins).
    @Test func testTheAnywayPathStillRefusesBelowOneThirdOfTheBins() throws {
        let plant = EllipsePlant(size: 128, radius: 40, eps: 0.02, thetaDegrees: 30, gap: 270, noise: 0, seed: 205)
        let pattern = plant.render()
        var count = -1
        do {
            _ = try EllipseCalibration.fit1D(pattern: pattern, centerQX: plant.centerQX + 1.3, centerQY: plant.centerQY - 0.9,
                                             innerRadius: plant.b - 20, outerRadius: plant.a + 20, maximumResidual: 0.4,
                                             acceptSparseCoverage: true)
        } catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let bins) = caught else { throw caught }
            count = bins
        }
        #expect(count >= 0 && count < EllipseCalibration.sparseFloorBins, "occupied bins \(count)")
    }

    // MARK: The one-ring check

    /// Three arcs at radii 24, 33 and 42 about one centre (20 occupied bins, measured) are refused for coverage by default;
    /// the anyway path catches them with `moreThanOneRing` (mean radii about the fitted centre 28.8 to 37.6 px, measured), the
    /// one thing the assertion "this holds one ring" does not make true; two arcs of ONE ring pass and are marked sparse.
    /// RED UNDER: the one-ring check removed or its 1.10 ratio raised to 1.5.
    @Test func testThreeGrainsAtDifferentRadiiAreNotOneRing() throws {
        let grains = arcs(size: 128, center: 64, sigma: 2, sectors: [(0, 60, 24), (120, 180, 33), (240, 300, 42)])
        #expect(throws: EllipseCalibration.FitError.insufficientAngularCoverage(20)) {
            try EllipseCalibration.fit1D(pattern: grains, centerQX: 64.4, centerQY: 63.6, innerRadius: 12, outerRadius: 54)
        }
        do {
            _ = try EllipseCalibration.fit1D(pattern: grains, centerQX: 64.4, centerQY: 63.6, innerRadius: 12, outerRadius: 54,
                                             maximumResidual: 0.4, acceptSparseCoverage: true)
            Issue.record("three grains were fitted as one ring")
        } catch let caught as EllipseCalibration.FitError { guard case .moreThanOneRing(let low, let high) = caught else { throw caught }
            #expect(high / low > 1.10, "ratio \(high / low)")
        }
        let oneRing = arcs(size: 128, center: 64, sigma: 2, sectors: [(0, 90, 30), (180, 270, 30)])
        let fit = try EllipseCalibration.fit1D(pattern: oneRing, centerQX: 64.4, centerQY: 63.6, innerRadius: 20, outerRadius: 48,
                                               maximumResidual: 0.4, acceptSparseCoverage: true)
        #expect(fit.sparseCoverage && fit.occupiedAngularBins == 19)
        #expect(abs((fit.a + fit.b) / 2 - 30) < 0.5)
    }

    // MARK: ADR 069: the 0.30 px floor guard

    /// A ring 0.27 px wide (a point-sampled ring: the profile fit ends at the ring's own 0.27 px) lies within 0.30 px of the
    /// 0.25 px floor: `fitAmorphousRing` refuses it with `ringWidthAtFloor`, and `fitBestAvailable` falls back to the conic and
    /// says so. A ring 0.35 px wide is above the floor and is accepted with its own width (measured: 0.3499999986 and
    /// 0.3500000016 px). The guard is `> 0.25 + 0.05`.
    /// RED UNDER: the guard removed, or its margin dropped (`> 0.25`) or shrunk to 0.26: the 0.27 ring is then accepted.
    @Test func testARingFittedWithinPointThreePxOfTheFloorIsRefusedAndFallsBackToTheConic() throws {
        func thin(_ sigma: Double) -> DiffractionPattern {
            ring(rows: 64, columns: 64, centerQX: 31.4, centerQY: 32.3, a: 20, b: 20, theta: 0, sigma: sigma)
        }
        let narrow = thin(0.27)
        let conic = try EllipseCalibration.fit1D(pattern: narrow, centerQX: 32.5, centerQY: 31.5, innerRadius: 10, outerRadius: 30,
                                                 maximumResidual: 0.4)
        #expect(throws: EllipseCalibration.FitError.ringWidthAtFloor) {
            try EllipseCalibration.fitAmorphousRing(pattern: narrow, initial: conic, centerQX: 32.5, centerQY: 31.5,
                                                    innerRadius: 10, outerRadius: 30)
        }
        let fallback = try EllipseCalibration.fitBestAvailable(pattern: narrow, centerQX: 32.5, centerQY: 31.5, innerRadius: 10,
                                                               outerRadius: 30)
        #expect(fallback.model == .conic)
        #expect(fallback.profile == nil)
        #expect(fallback.profileFallbackReason?.contains("0.25 px floor") ?? false, "\(fallback.profileFallbackReason ?? "nil")")

        let wide = thin(0.35)
        let accepted = try EllipseCalibration.fitBestAvailable(pattern: wide, centerQX: 32.5, centerQY: 31.5, innerRadius: 10,
                                                               outerRadius: 30)
        #expect(accepted.model == .asymmetricGaussian)
        let profile = try #require(accepted.profile)
        #expect(abs(profile.innerSigma - 0.35) < 1e-5 && abs(profile.outerSigma - 0.35) < 1e-5)
        #expect(abs((accepted.a + accepted.b) / 2 - 20) < 1e-4)
    }

    // MARK: The profile residual cap

    /// Two equal rings 8 px apart, then 9 px apart, inside one annulus (s = 2 px, R 24, centre (39.6, 40.3), annulus 10..40
    /// about (40.5, 39.5)): the one-ring profile leaves a normalised residual of 0.2031 at 8 px, which the 0.22 cap accepts, and
    /// 0.2459 at 9 px, which it refuses with `excessiveResidual` (measured; the conic seed is 0.27 and 0.28 of a dynamic range,
    /// above the 0.2 conic cap, so no conic fallback hides the refusal).
    /// RED UNDER: the profile cap 0.22 raised to 0.30 (the 9 px pair is accepted) or lowered to 0.20 (the 8 px pair is refused).
    @Test func testTheProfileResidualCapAcceptsPointTwoZeroThreeAndRefusesPointTwoFourSix() throws {
        func pair(_ separation: Double) -> DiffractionPattern {
            var pixels = [Float](repeating: 0, count: 80 * 80)
            for qx in 0..<80 {
                for qy in 0..<80 {
                    let r = hypot(Double(qx) - 39.6, Double(qy) - 40.3)
                    let value = exp(-(r - 24) * (r - 24) / 8) + exp(-(r - 24 - separation) * (r - 24 - separation) / 8)
                    pixels[qx * 80 + qy] = Float(value)
                }
            }
            return DiffractionPattern(width: 80, height: 80, pixels: pixels)
        }
        let near = try EllipseCalibration.fitBestAvailable(pattern: pair(8), centerQX: 40.5, centerQY: 39.5, innerRadius: 10, outerRadius: 40)
        #expect(near.model == .asymmetricGaussian)
        #expect(near.normalizedResidual > 0.20 && near.normalizedResidual < 0.22, "\(near.normalizedResidual)")
        do {
            _ = try EllipseCalibration.fitBestAvailable(pattern: pair(9), centerQX: 40.5, centerQY: 39.5, innerRadius: 10, outerRadius: 40)
            Issue.record("a residual above the cap was accepted")
        } catch let caught as EllipseCalibration.FitError { guard case .excessiveResidual(let residual) = caught else { throw caught }
            #expect(residual > 0.22 && residual < 0.26, "\(residual)")
        }
    }

    // MARK: Truth: a clean circle

    /// A clean circular ring of width s about a known centre, on both paths: the conic's objective puts its radius at
    /// R (1 + 1.75 (s/R)^2) (a property of the objective, in the header: +0.4375 % at R 30, s 3), the profile fit
    /// recovers R and its width s; a/b is 1 on both. Plant: the reference implementations' stop-rule ring 1 (96 px, R 30, s 3).
    /// RED UNDER: the profile model's radius taken from the conic; the profile's sigmas swapped with the central one.
    @Test func testACleanCircleGivesTheBiasedRadiusOnTheConicAndTheTrueRadiusOnTheProfile() throws {
        let plant = EllipsePlant(size: 96, radius: 30, eps: 0, thetaDegrees: 0, gap: 0, noise: 0, seed: 1)
        let pattern = plant.render()
        let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: plant.centerQX + 1.3, centerQY: plant.centerQY - 0.9,
                                                 innerRadius: 10, outerRadius: 50)
        #expect(abs((conic.a + conic.b) / 2 - 30 * (1 + 1.75 * 0.01)) <= 0.02, "\((conic.a + conic.b) / 2)")
        #expect(abs(conic.a / conic.b - 1) < 1e-8)
        let best = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: plant.centerQX + 1.3,
                                                           centerQY: plant.centerQY - 0.9, innerRadius: 10, outerRadius: 50)
        #expect(best.model == .asymmetricGaussian)
        #expect(abs((best.a + best.b) / 2 - 30) < 1e-4, "\((best.a + best.b) / 2)")
        #expect(abs(best.a / best.b - 1) < 1e-4)
        let profile = try #require(best.profile)
        #expect(abs(profile.innerSigma - 3) < 1e-3 && abs(profile.outerSigma - 3) < 1e-3, "\(profile)")
        #expect(abs(best.centerQX - plant.centerQX) < 1e-4 && abs(best.centerQY - plant.centerQY) < 1e-4)
    }

    // MARK: Row and column

    /// qx is the ROW and qy the COLUMN, on a pattern that is not square. The ring of the reference implementations' reference harness
    /// ("non_square_rotated": 72 rows by 110 columns, centre (34.5, 58.25), a 24, b 17, theta 0.42 rad, s 1.1 px) is found with
    /// its centre in (row, column) order; the same pattern TRANSPOSED (110 rows by 72 columns, centre (58.25, 34.5)) gives the axes
    /// unchanged and theta' = pi/2 - theta (the major axis direction (cos theta, sin theta) in (row, column) becomes
    /// (sin theta, cos theta)).
    /// RED UNDER: the pattern's width and height exchanged in the loops; qx and qy exchanged in the sample coordinates;
    /// theta measured from the column axis.
    @Test func testRowIsQxAndColumnIsQyOnANonSquarePattern() throws {
        let theta = 0.42
        let wide = ring(rows: 72, columns: 110, centerQX: 34.5, centerQY: 58.25, a: 24, b: 17, theta: theta, sigma: 1.1)
        let fit = try EllipseCalibration.fitBestAvailable(pattern: wide, centerQX: 35.3, centerQY: 57.65, innerRadius: 9, outerRadius: 33)
        #expect(fit.model == .asymmetricGaussian)
        #expect(abs(fit.centerQX - 34.5) < 1e-3 && abs(fit.centerQY - 58.25) < 1e-3, "(\(fit.centerQX), \(fit.centerQY))")
        #expect(abs(fit.a - 24) < 1e-2 && abs(fit.b - 17) < 1e-2, "a \(fit.a) b \(fit.b)")
        #expect(abs(fit.theta - theta) < 1e-3, "theta \(fit.theta)")

        var transposed = [Float](repeating: 0, count: 110 * 72)
        for row in 0..<72 { for column in 0..<110 { transposed[column * 72 + row] = wide.pixels[row * 110 + column] } }
        let tall = DiffractionPattern(width: 72, height: 110, pixels: transposed)
        let other = try EllipseCalibration.fitBestAvailable(pattern: tall, centerQX: 57.65, centerQY: 35.3, innerRadius: 9, outerRadius: 33)
        #expect(abs(other.centerQX - 58.25) < 1e-3 && abs(other.centerQY - 34.5) < 1e-3, "(\(other.centerQX), \(other.centerQY))")
        #expect(abs(other.a - 24) < 1e-2 && abs(other.b - 17) < 1e-2)
        #expect(abs(other.theta - (.pi / 2 - theta)) < 1e-3, "theta \(other.theta)")
    }

    // MARK: The clean sweep against SciPy

    private struct Sweep: Decodable {
        struct Ring: Decodable {
            let name: String, size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, seed: Int
            let startQX: Double, startQY: Double, innerRadius: Double, outerRadius: Double
            let pixelHash: String, sampleCount: Int, occupiedBins: Int
            let optimumCenterQX: Double, optimumCenterQY: Double, optimumA: Double, optimumB: Double, optimumTheta: Double
        }
        let scipy: String
        let rings: [Ring]
        static let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .appendingPathComponent("Fixtures/scipy-ellipse-sweep-pins.json")
        static let shared: Sweep = {
            do { return try JSONDecoder().decode(Sweep.self, from: Data(contentsOf: url)) } catch { fatalError("scipy-ellipse-sweep-pins.json: \(error)") }
        }()
    }

    private static func angleError(_ x: Double, _ y: Double) -> Double {
        let raw = abs(x - y).truncatingRemainder(dividingBy: .pi)
        return min(raw, .pi - raw)
    }

    /// The 20 clean sweep rings (square 112-144 px, R 32-52, eps 0 / 1 / 3 %, unbroken or a 60 degree arc missing; plant and start
    /// read from the sweep fixture) fitted by `fit1D` from the click's error equal SciPy 1.18.1's `least_squares` optimum of the same
    /// weighted-conic objective, tools/pins/ellipse/sweep.py: centre 1e-5 px, a and b 1e-7 relative, theta 1e-6 rad (not on a circle,
    /// where it is undefined), the same pixel and bin counts, and the centre within 0.05 px of the plant (0.1 px for the rings with an arc missing, see below). These are the bands of
    /// reference tools/ellipse-calibration-test (measured there: 1.1e-7 px, 6.4e-10, 1.7e-8 rad). This pins the OBJECTIVE: an
    /// unweighted or sqrt-weighted variant lands on a different minimum and fails here.
    /// RED UNDER: the residual weighted by sqrt(value) or not at all; the Jacobian wrong; a different start; the annulus boundary.
    @Test func testCleanSweepEqualsSciPysOptimum() throws {
        let sweep = Sweep.shared
        #expect(sweep.rings.count == 20)
        var worstPlant = 0.0, worstCentre = 0.0, worstAxes = 0.0, worstTheta = 0.0
        for ring in sweep.rings {
            let plant = EllipsePlant(size: ring.size, radius: ring.radius, eps: ring.eps, thetaDegrees: ring.thetaDegrees, gap: ring.gap,
                                     noise: 0, seed: ring.seed)
            let fit = try EllipseCalibration.fit1D(pattern: plant.render(), centerQX: ring.startQX, centerQY: ring.startQY,
                                                   innerRadius: ring.innerRadius, outerRadius: ring.outerRadius)
            #expect(fit.sampleCount == ring.sampleCount && fit.occupiedAngularBins == ring.occupiedBins, "\(ring.name)")
            let toPlant = max(abs(fit.centerQX - plant.centerQX), abs(fit.centerQY - plant.centerQY))
            let centre = max(abs(fit.centerQX - ring.optimumCenterQX), abs(fit.centerQY - ring.optimumCenterQY))
            let axes = max(abs(fit.a - ring.optimumA) / ring.optimumA, abs(fit.b - ring.optimumB) / ring.optimumB)
            let theta = ring.eps == 0 ? 0 : Self.angleError(fit.theta, ring.optimumTheta)
            worstPlant = max(worstPlant, toPlant); worstCentre = max(worstCentre, centre)
            worstAxes = max(worstAxes, axes); worstTheta = max(worstTheta, theta)
            // The reference implementations' clause is 0.05 px on rings of R 80-200. A conic fitted to a ring with a 60 degree arc missing sits off the
            // plant by an amount that grows as R shrinks (the objective, not the port: SciPy's optimum agrees to 7e-8 px); on these
            // R 32-52 rings it is 0.055-0.094 px (measured), so a gapped ring is held to 0.1 px and an unbroken one to 0.05 px.
            #expect(toPlant < (ring.gap == 0 ? 0.05 : 0.1), "\(ring.name): centre \(toPlant) px from the plant")
            #expect(centre < 1e-5, "\(ring.name): centre \(centre) px from scipy's")
            #expect(axes < 1e-7, "\(ring.name): axes \(axes) relative")
            #expect(theta < 1e-6, "\(ring.name): theta \(theta) rad")
        }
        print("clean sweep vs scipy \(sweep.scipy): plant \(worstPlant) px, centre \(worstCentre) px, axes \(worstAxes), theta \(worstTheta) rad")
    }
}
