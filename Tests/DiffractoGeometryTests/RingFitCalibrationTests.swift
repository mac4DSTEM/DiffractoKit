//
//  RingFitCalibrationTests.swift
//  The fit-driven calibration (`RingCalibration.fit`) on rings planted here, in Swift, with their truth
//  known. Ported from reference commit 9482831 Tests/RingFitCalibrationTests.swift (8 tests, same names, same plant, same bands),
//  RUN ON THE KIT'S FIT, which is reference commit 7d6fdbd1's `EllipseCalibration` (ADR 069: two starts and the 0.30 px floor guard);
//  the reference implementation's own copy there is the pre-ADR-069 file. No number of these eight moved: every assertion holds as the reference implementation wrote it
//  (said they might, the measured values are in tools/slice-3b/README.md).
//
//  The plant is the one tools/ring-fit-parity's plants.py draws: a Gaussian ring on an elliptical coordinate,
//  I = exp(-((rho - 1) R)^2 / 2s^2), rho^2 = u^2/a^2 + v^2/b^2 in the frame rotated by theta, R = (a + b)/2, inside the
//  asymmetric-Gaussian model's family (sigma_in = sigma_out = s, no central disc, no background). Rows are the fit's "qx",
//  columns its "qy"; the centre is off the grid and its row differs from its column.
//
//  Every test names under "RED UNDER" the production mutation that turns it red.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import DiffractoGeometry

struct RingFitCalibrationTests {

    // MARK: The plant

    private struct Plant {
        let rows: Int, cols: Int
        let row0: Double, col0: Double
        let a: Double, b: Double, theta: Double
        let sigma: Double
        var radius: Double { (a + b) / 2 }
    }

    /// The ring, float32, rows x columns with columns fastest; `gapDegrees` cuts the signal over that azimuth (measured
    /// from the row axis towards the column axis) starting at `gapStart`.
    private static func render(_ p: Plant, gapDegrees: Double = 0, gapStart: Double = 0) -> DiffractionPattern {
        var pixels = [Float](repeating: 0, count: p.rows * p.cols)
        let s = sin(p.theta), c = cos(p.theta)
        for r in 0..<p.rows {
            for col in 0..<p.cols {
                let dr = Double(r) - p.row0, dc = Double(col) - p.col0
                if gapDegrees > 0 {
                    var angle = atan2(dc, dr) * 180 / .pi
                    if angle < 0 { angle += 360 }
                    var offset = (angle - gapStart).truncatingRemainder(dividingBy: 360)
                    if offset < 0 { offset += 360 }
                    if offset < gapDegrees { continue }
                }
                let u = dr * c + dc * s, v = -dr * s + dc * c
                let rho = ((u / p.a) * (u / p.a) + (v / p.b) * (v / p.b)).squareRoot()
                let x = (rho - 1) * p.radius
                pixels[r * p.cols + col] = Float(exp(-x * x / (2 * p.sigma * p.sigma)))
            }
        }
        return DiffractionPattern(width: p.cols, height: p.rows, pixels: pixels)
    }

    /// 2 % ellipticity (a/b = 1.02), the major axis at 35 degrees from the rows, R = 60 px, s = 2.5 px; centre (81.31, 77.73).
    private static let plant: Plant = {
        let radius = 60.0, epsilon = 0.02
        let b = 2 * radius / (2 + epsilon)
        return Plant(rows: 160, cols: 156, row0: 81.31, col0: 77.73, a: b * (1 + epsilon), b: b,
                     theta: 35 * .pi / 180, sigma: 2.5)
    }()

    /// The annulus (b - 15, a + 15) about a start centre a click away from the plant's, (+1.3, -0.9) px.
    private static func fit(_ pattern: DiffractionPattern, _ p: Plant = plant, d: Double = auD220,
                            volts: Double? = 200_000, pitchMM: Double? = 0.1) throws -> RingFitCalibration {
        try RingCalibration.fit(pattern: pattern, centerQX: p.row0 + 1.3, centerQY: p.col0 - 0.9,
                                innerRadius: p.b - 15, outerRadius: p.a + 15,
                                d: d, highTensionVolts: volts, pixelPitchMM: pitchMM)
    }

    /// Au {220}: a/sqrt(8) with a = 4.0782 A.
    private static let auD220 = 4.0782 / 8.0.squareRoot()

    // MARK: A clean ring

    /// The planted ring calibrates to its own Δq and L within the documented 0.1 %, its centre within 0.05 px, and the
    /// whole chain is the closed forms of RingCalibration applied to the fitted radius.
    /// RED UNDER: the radius handed to `calibrate` taken as `fit.a` (1.0 % off at 2 % ellipticity); the d argument
    /// ignored for a fixed one. (Swapping the START centre's row and column is not caught: the fit converges from that
    /// 5 px-off start to the same ring; measured by the reference implementation 2026-10-13.)
    @Test func testAPlantedRingCalibratesToItsTruth() throws {
        let p = Self.plant
        let result = try Self.fit(Self.render(p))
        #expect(result.refusals == [])
        let fit = try #require(result.fit)
        #expect(fit.model == .asymmetricGaussian)
        #expect(abs(fit.centerQX - p.row0) <= 0.05, "the row")
        #expect(abs(fit.centerQY - p.col0) <= 0.05, "the column")
        let calibration = try #require(result.calibration)
        let truth = try RingCalibration.calibrate(radius: p.radius, d: Self.auD220, highTensionVolts: 200_000,
                                                  pixelPitchMM: 0.1)
        #expect(abs(calibration.deltaQ / truth.deltaQ - 1) <= 1e-3, "Δq")
        let calibratedPixels: Double = try #require(calibration.lengthPixels)
        let truthPixels: Double = try #require(truth.lengthPixels)
        let pixelsRatio: Double = calibratedPixels / truthPixels - 1
        #expect(abs(pixelsRatio) <= 1e-3)
        let calibratedMillimetres: Double = try #require(calibration.lengthMillimetres)
        let truthMillimetres: Double = try #require(truth.lengthMillimetres)
        let millimetresRatio: Double = calibratedMillimetres / truthMillimetres - 1
        #expect(abs(millimetresRatio) <= 1e-3)
        #expect(calibration.d == Self.auD220)
        #expect(calibration.highTensionVolts == 200_000)
    }

    /// The report the refusals are decided on carries the fit's own numbers and the annulus it was given, and the
    /// calibration's radius is the report's (a + b)/2 to the last bit.
    /// RED UNDER: the annulus width taken as the outer radius; the sigmas not copied from the profile (nil); the radius
    /// taken as one axis.
    @Test func testTheReportCarriesTheFitAndTheAnnulus() throws {
        let p = Self.plant
        let result = try Self.fit(Self.render(p))
        let fit = try #require(result.fit)
        let report = try #require(result.report)
        #expect(report.model == fit.model)
        #expect(report.semiAxisA == fit.a)
        #expect(report.semiAxisB == fit.b)
        #expect(report.occupiedAngularBins == fit.occupiedAngularBins)
        #expect(report.sparseCoverage == fit.sparseCoverage)
        #expect(report.normalizedResidual == fit.normalizedResidual)
        #expect(report.innerSigma == fit.profile?.innerSigma)
        #expect(report.outerSigma == fit.profile?.outerSigma)
        #expect(report.innerSigma != nil)
        #expect(report.annulusWidth == (p.a + 15) - (p.b - 15))
        #expect(try #require(result.calibration).radius == (fit.a + fit.b) / 2)
    }

    // MARK: Refusals through the fit

    /// No high tension: Δq is given, L is not, and the one refusal is `noHighTension`.
    /// RED UNDER: a missing V defaulted to 200 kV.
    @Test func testWithoutAHighTensionOnlyDeltaQIsGiven() throws {
        let result = try Self.fit(Self.render(Self.plant), volts: nil)
        #expect(result.refusals.map(\.name) == ["noHighTension"])
        let calibration = try #require(result.calibration)
        #expect(calibration.highTensionVolts == nil)
        #expect(calibration.lengthPixels == nil)
        #expect(calibration.lengthMillimetres == nil)
        #expect(calibration.lengthRefusal == .noHighTension)
        #expect(abs(calibration.deltaQ * Self.plant.radius * Self.auD220 - 1) <= 1e-3)
    }

    /// A ring with a 120 degree gap (24 of 36 bins) is refused by its coverage, by name and with its bin count, and nothing is
    /// calibrated from it; the fit's own refusal is not thrown at the caller as an error.
    /// RED UNDER: the coverage error rethrown instead of reported; the fit run with `acceptSparseCoverage: true` (Δq is
    /// still withheld, by `refusals`' own bin bound, but a fit comes back).
    @Test func testAGappedRingIsRefusedByItsCoverage() throws {
        let result = try Self.fit(Self.render(Self.plant, gapDegrees: 120, gapStart: 10))
        #expect(result.refusals.map(\.name) == ["insufficientAngularCoverage"])
        guard case .insufficientAngularCoverage(let bins) = try #require(result.refusals.first) else {
            Issue.record("not a coverage refusal")
            return
        }
        #expect(bins < RingCalibration.minimumOccupiedAngularBins)
        #expect(bins >= 22)
        #expect(result.fit == nil)
        #expect(result.report == nil)
        #expect(result.calibration == nil)
    }

    /// The refusals are listed whole on the coverage path too: a gapped ring with no high tension is refused for both,
    /// coverage first.
    /// RED UNDER: `noHighTension` not appended when the fit refuses the coverage.
    @Test func testTheCoveragePathListsEveryRefusal() throws {
        let result = try Self.fit(Self.render(Self.plant, gapDegrees: 120, gapStart: 10), volts: nil)
        #expect(result.refusals.map(\.name) == ["insufficientAngularCoverage", "noHighTension"])
        #expect(result.calibration == nil)
    }

    /// A ring the fit finds but whose a/b − 1 is past 0.093 (12 % planted) is refused by name and nothing is calibrated
    /// from it, though Δq could be computed.
    /// RED UNDER: the calibration computed whatever the refusals say (the guard on `withholdsDeltaQ` removed).
    @Test func testAnEllipticRingIsRefusedAndNotCalibrated() throws {
        let radius = 60.0, epsilon = 0.12
        let b = 2 * radius / (2 + epsilon)
        let p = Plant(rows: 172, cols: 168, row0: 86.41, col0: 83.62, a: b * (1 + epsilon), b: b,
                      theta: 70 * .pi / 180, sigma: 2.5)
        let result = try Self.fit(Self.render(p), p)
        #expect(result.fit != nil)
        #expect(result.refusals.map(\.name).contains("excessiveEllipticity"), "\(result.refusals.map(\.name))")
        #expect(result.calibration == nil)
    }

    /// An empty pattern is not a refusal of a ring but no ring at all: the fit's error reaches the caller.
    /// RED UNDER: every fit error swallowed into an empty result.
    @Test func testAnEmptyPatternThrowsTheFitsError() {
        let p = Self.plant
        let empty = DiffractionPattern(width: p.cols, height: p.rows, pixels: [Float](repeating: 0, count: p.rows * p.cols))
        #expect(throws: EllipseCalibration.FitError.insufficientSignal) { try Self.fit(empty) }
    }

    /// The d, V and pitch are checked before any fit runs: an impossible one is the caller's mistake, by name.
    /// RED UNDER: the inputs checked only after the fit (an empty pattern would then throw the fit's error instead).
    @Test func testAnInvalidInputIsRefusedBeforeTheFit() {
        let p = Self.plant
        let empty = DiffractionPattern(width: p.cols, height: p.rows, pixels: [Float](repeating: 0, count: p.rows * p.cols))
        for (d, volts, pitch) in [(-1.0, 200_000.0, 0.1), (1.44, -5.0, 0.1), (1.44, 200_000.0, .nan)] as [(Double, Double?, Double?)] {
            do {
                _ = try Self.fit(empty, d: d, volts: volts, pitchMM: pitch)
                Issue.record("d \(d) V \(String(describing: volts)) was not refused")
            } catch {
                #expect((error as? RingCalibration.Refusal)?.name == "invalidInput", "d \(d) V \(String(describing: volts))")
            }
        }
    }
}
