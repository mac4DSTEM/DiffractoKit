//
//  Geometry004Checks.swift
//  Role: the public API DiffractoKit 0.1.0 adds to DiffractoGeometry, used from outside the package with a real call and a
//        result a reader verifies by hand: q(radius:...) and its inverse, the ring evidence, `fitRing` / `calibrating`
//        and their equality with `fit`, `requireRingEvidence`, the three new refusals, `FitReport`'s new fields, and
//        `progress:` / `isCancelled:` on every fit. `runGeometry004Checks` is called from main.swift.
//

import Foundation
import DiffractoCore
import DiffractoGeometry

/// A progress sink a `@Sendable` closure may write to.
private final class ProgressLog: @unchecked Sendable {
    private let lock = NSLock()
    private var values: [Double] = []
    func add(_ value: Double) { lock.lock(); values.append(value); lock.unlock() }
    var all: [Double] { lock.lock(); defer { lock.unlock() }; return values }
    var isMonotoneWithinUnit: Bool {
        let v = all
        return !v.isEmpty && v.allSatisfy { $0 >= 0 && $0 <= 1 } && zip(v, v.dropFirst()).allSatisfy { $0 <= $1 }
    }
}

/// Whether `body` throws a `CancellationError`.
private func throwsCancellation(_ body: () throws -> Void) -> Bool {
    do { try body() } catch is CancellationError { return true } catch {}
    return false
}

func runGeometry004Checks(_ report: inout Report) {
    let size = 128, R = 40.0, s = 2.0, cr = 63.4, cc = 64.2
    let pattern = plantedRing(size: size, centerRow: cr, centerColumn: cc, radius: R, sigma: s)
    let half = plantedRing(size: size, centerRow: cr, centerColumn: cc, radius: R, sigma: s, keepAngle: { $0 > 0 })
    let d = 2.3546, volts = 200_000.0
    let lambda = ElectronWavelength.angstroms(volts: volts)

    report.section("0.1.0: RingCalibration.q and radius") { r in
        let twoTheta = RingCalibration.twoTheta(wavelengthAngstroms: lambda, d: d)
        let lengthPx = RingCalibration.cameraLengthPixels(radius: 150, twoTheta: twoTheta)
        // For the calibrating ring's own (R, L_px, lambda), q = 2 sin(theta) / lambda = 1 / d (Bragg's law).
        r.near("q(150 px, L_px, lambda) = 1/d for Au {111}", RingCalibration.q(radius: 150, cameraLengthPixels: lengthPx,
                                                                              wavelengthAngstroms: lambda), 1 / d, 1e-13)
        r.near("radius(1/d, L_px, lambda) = 150 px", RingCalibration.radius(q: 1 / d, cameraLengthPixels: lengthPx,
                                                                           wavelengthAngstroms: lambda), 150, 1e-9)
        // q = 2 sin(atan2(R, L) / 2) / lambda, by hand for R = L: atan 1 = pi/4, q = 2 sin(pi/8) / lambda.
        r.near("R = L: q = 2 sin(pi/8) / lambda", RingCalibration.q(radius: 1000, cameraLengthPixels: 1000, wavelengthAngstroms: lambda),
               2 * sin(.pi / 8) / lambda, 1e-12)
        let q = RingCalibration.q(radius: 77, cameraLengthPixels: 9000, wavelengthAngstroms: lambda)
        r.near("radius inverts q", RingCalibration.radius(q: q, cameraLengthPixels: 9000, wavelengthAngstroms: lambda), 77, 1e-10)
        r.equal("q of a zero radius is NaN", RingCalibration.q(radius: 0, cameraLengthPixels: 9000, wavelengthAngstroms: lambda).isNaN, true)
        r.equal("q of a negative camera length is NaN", RingCalibration.q(radius: 10, cameraLengthPixels: -1, wavelengthAngstroms: lambda).isNaN, true)
        r.equal("q of a NaN wavelength is NaN", RingCalibration.q(radius: 10, cameraLengthPixels: 100, wavelengthAngstroms: .nan).isNaN, true)
        // q lambda / 2 at or above 1/sqrt 2 (2 theta at 90 degrees) shows no ring on a flat detector.
        r.equal("radius at the Bragg limit is NaN", RingCalibration.radius(q: 2 * RingCalibration.maximumLambdaOver2d / lambda,
                                                                           cameraLengthPixels: 100, wavelengthAngstroms: lambda).isNaN, true)
        r.equal("radius of a zero q is NaN", RingCalibration.radius(q: 0, cameraLengthPixels: 100, wavelengthAngstroms: lambda).isNaN, true)
    }

    report.section("0.1.0: RingEvidence by hand and the new refusals") { r in
        r.equal("registered bounds", RingEvidence.sectorCount == 36 && RingEvidence.minimumSectorPixels == 8
                && RingEvidence.minimumAmplitudeFraction == 0.5 && RingEvidence.minimumSignificance == 3
                && RingEvidence.minimumOccupiedBins == 30 && RingEvidence.minimumDrop == 0.25
                && RingEvidence.maximumShift == 0.25 && RingEvidence.shiftSignificance == 4, true)
        func evidence(bins: Int = 36, inner: Double? = 0.9, outer: Double? = 0.9, centre: Double? = 0.05, centreError: Double? = 0.02,
                      ellipse: Double? = 0.05, ellipseError: Double? = 0.02, secondary: Double? = nil) -> RingEvidence {
            RingEvidence(amplitude: 100, sectorAmplitudes: Array(repeating: 100, count: 36),
                         sectorAmplitudeErrors: Array(repeating: 1, count: 36), sectorShifts: Array(repeating: 0, count: 36),
                         occupiedBins: bins, innerDrop: inner, outerDrop: outer, radiusShift: 0, centreShift: centre,
                         centreShiftError: centreError, ellipseShift: ellipse, ellipseShiftError: ellipseError,
                         secondaryPeak: secondary)
        }
        let good = evidence()
        r.equal("the fields are the ones given", good.amplitude == 100 && good.sectorAmplitudes.count == 36 && good.sectorAmplitudeErrors.count == 36
                && good.sectorShifts.count == 36 && good.occupiedBins == 36 && good.innerDrop == 0.9 && good.outerDrop == 0.9
                && good.radiusShift == 0 && good.centreShift == 0.05 && good.centreShiftError == 0.02
                && good.ellipseShift == 0.05 && good.ellipseShiftError == 0.02, true)
        r.equal("a clean evidence passes", good.isPeak && good.geometryHolds && good.passes, true)
        r.equal("drops 0.25 and 0.25 are a peak, 0.2499 on one side is not",
                evidence(inner: 0.25, outer: 0.25).isPeak && !evidence(inner: 0.2499).isPeak && !evidence(outer: 0.2499).isPeak, true)
        r.equal("a side with no bin is not a peak", !evidence(inner: nil).isPeak && !evidence(outer: nil).isPeak, true)
        r.equal("29 of 36 sectors do not pass, 30 do", !evidence(bins: 29).passes && evidence(bins: 30).passes, true)
        // A shift fails only beyond 0.25 px AND 4 standard errors.
        r.equal("a 0.3 px centre shift at SE 0.05 fails", !evidence(centre: 0.3, centreError: 0.05).geometryHolds, true)
        r.equal("a 0.3 px centre shift at SE 0.1 holds (within 4 SE)", evidence(centre: 0.3, centreError: 0.1).geometryHolds, true)
        r.equal("a 0.3 px ellipse shift at SE 0.05 fails", !evidence(ellipse: 0.3, ellipseError: 0.05).geometryHolds, true)
        r.equal("an undefined shift does not hold", !evidence(centre: nil).geometryHolds && !evidence(ellipseError: nil).geometryHolds, true)
        r.equal("an evidence is Equatable", good == evidence() && good != evidence(bins: 35), true)

        // The three new refusals, decided by `refusals` on a FitReport carrying the annulus and the evidence.
        func fitReport(a: Double = 40, b: Double = 40, inner: Double? = 28, outer: Double? = 52,
                       evidence: RingEvidence? = good) -> RingCalibration.FitReport {
            RingCalibration.FitReport(model: .asymmetricGaussian, semiAxisA: a, semiAxisB: b, occupiedAngularBins: 36,
                                      sparseCoverage: false, normalizedResidual: 0.05, innerSigma: 2, outerSigma: 2,
                                      annulusWidth: 24, annulusInnerRadius: inner, annulusOuterRadius: outer, ringEvidence: evidence)
        }
        let clean = fitReport()
        r.equal("FitReport carries the annulus and the evidence", clean.annulusInnerRadius == 28 && clean.annulusOuterRadius == 52
                && clean.ringEvidence == good, true)
        r.equal("FitReport's new parameters default to nil", RingCalibration.FitReport(
            model: .asymmetricGaussian, semiAxisA: 40, semiAxisB: 40, occupiedAngularBins: 36, sparseCoverage: false,
            normalizedResidual: 0.05, innerSigma: 2, outerSigma: 2, annulusWidth: 24).ringEvidence == nil, true)
        r.equal("a clean report with the evidence has no refusal", RingCalibration.refusals(clean, highTensionVolts: volts), [])

        let outside = RingCalibration.refusals(fitReport(a: 60, b: 60), highTensionVolts: volts)
        r.equal("a radius of 60 px in the annulus (28, 52) is ringOutsideAnnulus", outside, [.ringOutsideAnnulus(radius: 60, inner: 28, outer: 52)])
        if case .ringOutsideAnnulus(let radius, let inner, let outer)? = outside.first {
            r.equal("ringOutsideAnnulus carries its numbers", radius == 60 && inner == 28 && outer == 52, true)
        } else { r.check("ringOutsideAnnulus carries its numbers", false) }
        r.equal("a radius on the annulus's edge is outside it", RingCalibration.refusals(fitReport(a: 52, b: 52), highTensionVolts: volts).map(\.name),
                ["ringOutsideAnnulus"])
        r.equal("without the annulus's radii the check is skipped", RingCalibration.refusals(fitReport(a: 60, b: 60, inner: nil, outer: nil),
                                                                                           highTensionVolts: volts), [])

        let edge = RingCalibration.refusals(fitReport(evidence: evidence(inner: 0.1, outer: nil)), highTensionVolts: volts)
        r.equal("drops (0.1, none) are notARingPeak", edge, [.notARingPeak(innerDrop: 0.1, outerDrop: nil)])
        if case .notARingPeak(let inner, let outer)? = edge.first {
            r.equal("notARingPeak carries its drops", inner == 0.1 && outer == nil, true)
        } else { r.check("notARingPeak carries its drops", false) }

        let pulled = RingCalibration.refusals(fitReport(evidence: evidence(centre: 0.6, centreError: 0.05)), highTensionVolts: volts)
        r.equal("a 0.6 px centre shift is ringGeometryInconsistent", pulled, [.ringGeometryInconsistent(centreShift: 0.6, ellipseShift: 0.05)])
        if case .ringGeometryInconsistent(let centre, let ellipse)? = pulled.first {
            r.equal("ringGeometryInconsistent carries its shifts", centre == 0.6 && ellipse == 0.05, true)
        } else { r.check("ringGeometryInconsistent carries its shifts", false) }

        let twoRings = RingCalibration.refusals(fitReport(evidence: evidence(secondary: 0.6)), highTensionVolts: volts)
        r.equal("a second peak of 0.6 A is secondRingInAnnulus", twoRings, [.secondRingInAnnulus(height: 0.6)])
        let secondary: Double? = evidence(secondary: 0.6).secondaryPeak
        r.equal("secondaryPeak is the one given", secondary == 0.6, true)
        r.equal("a second peak below the bound is one ring", evidence(secondary: 0.2).isSingleRing, true)
        r.equal("the second-peak bound is 0.25", RingEvidence.maximumSecondaryPeak, 0.25)
        let second = RingCalibration.Refusal.secondRingInAnnulus(height: 0.6)
        r.equal("secondRingInAnnulus is named, withholds Δq and says why", second.name == "secondRingInAnnulus", true)
        r.equal("secondRingInAnnulus withholds Δq", second.withholdsDeltaQ, true)
        r.equal("secondRingInAnnulus says why", second.reason.contains("a peak of 60 % of the fitted ring's amplitude"), true)

        let fewSectors = RingCalibration.refusals(fitReport(evidence: evidence(bins: 20)), highTensionVolts: volts)
        r.equal("20 sectors holding the ring is a coverage refusal at 20", fewSectors, [.insufficientAngularCoverage(occupiedBins: 20)])

        let all = [RingCalibration.Refusal.ringOutsideAnnulus(radius: 60, inner: 28, outer: 52),
                   .notARingPeak(innerDrop: 0.1, outerDrop: nil), .ringGeometryInconsistent(centreShift: 0.6, ellipseShift: 0.05)]
        r.equal("the new refusals' names", all.map(\.name), ["ringOutsideAnnulus", "notARingPeak", "ringGeometryInconsistent"])
        r.equal("the new refusals withhold Δq", all.allSatisfy(\.withholdsDeltaQ), true)
        r.equal("ringOutsideAnnulus says why", all[0].reason.hasPrefix("The fitted ring radius of 60.0 px lies outside the annulus from 28.0 to 52.0 px"), true)
        r.equal("notARingPeak says why (nothing for a missing side)", all[1].reason.contains("falls by 10 % of the ring's amplitude inside it and by nothing outside it"), true)
        r.equal("ringGeometryInconsistent says why", all[2].reason.contains("0.60 px off the fitted centre and 0.05 px off the fitted ellipse"), true)
        r.equal("the new refusals are LocalizedErrors that speak their reason", all.allSatisfy { $0.errorDescription == $0.reason && !$0.reason.isEmpty }, true)
        r.equal("Refusal's localizedDescription is its reason", all[0].localizedDescription == all[0].reason, true)
    }

    report.section("0.1.0: ring evidence on a planted ring, requireRingEvidence") { r in
        let best = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                                           requireRingEvidence: true)
        let plain = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("without requireRingEvidence the fit carries none", plain.ringEvidence, nil)
        let evidence = try unwrap(best.ringEvidence)
        r.equal("fitBestAvailable(requireRingEvidence:) carries the evidence", best.model, .asymmetricGaussian)
        r.near("the radius is the plain fit's within 0.5 px", (best.a + best.b) / 2, (plain.a + plain.b) / 2, 0.5)
        r.equal("a clean ring: 36 of 36 sectors hold it", evidence.occupiedBins, 36)
        r.equal("a clean ring: peak, geometry and passes", evidence.isPeak && evidence.geometryHolds && evidence.passes, true)
        r.near("the median sector amplitude is the planted peak 100", evidence.amplitude, 100, 5)
        r.equal("36 sectors in each array", evidence.sectorAmplitudes.count == 36 && evidence.sectorAmplitudeErrors.count == 36
                && evidence.sectorShifts.count == 36, true)
        r.equal("every sector reads and holds the ring", evidence.sectorAmplitudes.allSatisfy { ($0 ?? 0) > 50 }
                && evidence.sectorAmplitudeErrors.allSatisfy { ($0 ?? .infinity) < 10 } && evidence.sectorShifts.allSatisfy { $0 != nil }, true)
        r.equal("both drops are a fraction of the amplitude, above 0.25", (evidence.innerDrop ?? 0) > 0.25 && (evidence.outerDrop ?? 0) > 0.25
                && (evidence.innerDrop ?? 2) <= 1.1 && (evidence.outerDrop ?? 2) <= 1.1, true)
        r.near("the consensus radial offset is zero", try unwrap(evidence.radiusShift), 0, 0.1)
        r.near("the centre shift is zero", try unwrap(evidence.centreShift), 0, 0.1)
        r.near("the ellipse shift is zero", try unwrap(evidence.ellipseShift), 0, 0.1)
        r.equal("the shifts have standard errors", (evidence.centreShiftError ?? -1) >= 0 && (evidence.ellipseShiftError ?? -1) >= 0, true)

        let direct = EllipseCalibration.ringEvidence(pattern: pattern, fit: best, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("EllipseCalibration.ringEvidence is the fit's own", direct, evidence)
        r.equal("EllipseCalibrationFit.ringEvidence is public and settable through init", EllipseCalibrationFit(
            centerQX: 1, centerQY: 2, a: 3, b: 4, theta: 0, normalizedResidual: 0, conicResidual: 0, sampleCount: 5, occupiedAngularBins: 36,
            model: .conic, profile: nil, profileFallbackReason: nil, ringEvidence: evidence).ringEvidence, evidence)
        r.equal("the init's ringEvidence defaults to nil", EllipseCalibrationFit(
            centerQX: 1, centerQY: 2, a: 3, b: 4, theta: 0, normalizedResidual: 0, conicResidual: 0, sampleCount: 5, occupiedAngularBins: 36,
            model: .conic, profile: nil, profileFallbackReason: nil).ringEvidence, nil)

        // A measurement: a degenerate annulus gives an empty evidence and never throws.
        let empty = EllipseCalibration.ringEvidence(pattern: pattern, fit: best, centerQX: 60, centerQY: 60, innerRadius: 30, outerRadius: 30)
        r.equal("an empty annulus gives an empty evidence that does not pass", empty.occupiedBins == 0 && empty.amplitude == 0 && !empty.passes
                && empty.sectorAmplitudes.count == 36, true)

        // The conic (fit1D) has no evidence of its own; the sectors read the same ring.
        let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("fit1D carries no evidence", conic.ringEvidence, nil)
        let conicEvidence = EllipseCalibration.ringEvidence(pattern: pattern, fit: conic, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("the conic's ring is read by the sectors too: 36 of 36", conicEvidence.occupiedBins, 36)

        let refined = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 60, centerQY: 60,
                                                              innerRadius: 28, outerRadius: 52, requireRingEvidence: true)
        // (fitAmorphousRing(requireRingEvidence:) only widens the floor rule; fitBestAvailable attaches the evidence. The
        // evidence of its fit is measured by the call below: the doc comment of `ringEvidence` says the fit carries it.)
        r.equal("fitAmorphousRing(requireRingEvidence:) is the profile fit", refined.model, .asymmetricGaussian)
        r.near("fitAmorphousRing(requireRingEvidence:) recovers R within 1 %", (refined.a + refined.b) / 2, R, 0.01 * R)
        r.equal("its ring evidence, measured, passes", EllipseCalibration.ringEvidence(
            pattern: pattern, fit: refined, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52).passes, true)

        // A central blob: no ring in the annulus. It is refused or not usable, never calibrated from.
        var blob = [Float](repeating: 0, count: size * size)
        for row in 0..<size { for column in 0..<size {
            let dr = Double(row) - cr, dc = Double(column) - cc
            blob[row * size + column] = Float(100 * exp(-(dr * dr + dc * dc) / (2 * 30 * 30)))
        } }
        let blobPattern = DiffractionPattern(width: size, height: size, pixels: blob)
        let blobRing = try? RingCalibration.fitRing(pattern: blobPattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("a central blob is not calibrated from", blobRing == nil || blobRing?.isUsable == false, true)
    }

    report.section("0.1.0: fitRing, calibrating, FitReport on a planted ring") { r in
        let ring = try RingCalibration.fitRing(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        let fit = try unwrap(ring.fit)
        let rep = try unwrap(ring.report)
        r.equal("a clean ring has no refusal and is usable", ring.refusals.isEmpty && ring.isUsable, true)
        r.near("radius is the planted R within 1 %", try unwrap(ring.radius), R, 0.01 * R)
        r.equal("radius is the report's (a + b)/2", ring.radius, rep.radius)
        r.equal("the report names the annulus it was given", rep.annulusInnerRadius == 28 && rep.annulusOuterRadius == 52
                && rep.annulusWidth == 24, true)
        r.equal("the report carries the fit's evidence, which passes", rep.ringEvidence == fit.ringEvidence && rep.ringEvidence?.passes == true, true)
        r.equal("the report is the fit's numbers", rep.semiAxisA == fit.a && rep.semiAxisB == fit.b && rep.model == fit.model
                && rep.occupiedAngularBins == fit.occupiedAngularBins, true)

        let calibrated = try RingCalibration.calibrating(ring, d: d, highTensionVolts: volts, pixelPitchMM: 0.014)
        let whole = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                            d: d, highTensionVolts: volts, pixelPitchMM: 0.014)
        r.equal("fit equals calibrating(fitRing(...))", calibrated, whole)
        r.near("calibrating: delta q = (1/d) / fitted R", try unwrap(calibrated.calibration).deltaQ, (1 / d) / rep.radius, 1e-15)
        r.equal("calibrating a clean ring has no refusal", calibrated.refusals, [])
        let noVolts = try RingCalibration.calibrating(ring, d: d, highTensionVolts: nil, pixelPitchMM: nil)
        r.equal("calibrating without V: only .noHighTension, delta q kept", noVolts.refusals == [.noHighTension] && noVolts.calibration != nil, true)
        var thrown: RingCalibration.Refusal?
        do { _ = try RingCalibration.calibrating(ring, d: -1, highTensionVolts: volts, pixelPitchMM: nil) } catch let e as RingCalibration.Refusal { thrown = e } catch {}
        r.equal("calibrating with an impossible d throws invalidInput", thrown?.name, "invalidInput")
        thrown = nil
        do { _ = try RingCalibration.calibrating(ring, d: d, highTensionVolts: volts, pixelPitchMM: 0) } catch let e as RingCalibration.Refusal { thrown = e } catch {}
        r.equal("calibrating with a zero pitch throws invalidInput", thrown?.name, "invalidInput")

        let sparse = try RingCalibration.fitRing(pattern: half, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("a half ring is refused for coverage before any fit", sparse.fit == nil && sparse.report == nil
                && sparse.refusals.map(\.name) == ["insufficientAngularCoverage"], true)
        r.equal("a refused ring has no radius and is not usable", sparse.radius == nil && !sparse.isUsable, true)
        let sparseCalibrated = try RingCalibration.calibrating(sparse, d: d, highTensionVolts: volts, pixelPitchMM: nil)
        r.equal("calibrating a refused ring gives no calibration", sparseCalibrated.calibration == nil && sparseCalibrated.fit == nil
                && sparseCalibrated.refusals == sparse.refusals, true)

        let constructed = RingRadiusFit(fit: nil, report: rep, refusals: [.excessiveResidual(0.5)])
        r.equal("RingRadiusFit is constructible; a report with a withholding refusal is not usable", constructed.radius == rep.radius
                && !constructed.isUsable && constructed.refusals == [.excessiveResidual(0.5)] && constructed.fit == nil, true)
        r.equal("RingRadiusFit with only a noHighTension refusal stays usable",
                RingRadiusFit(fit: nil, report: rep, refusals: [.noHighTension]).isUsable, true)
        r.equal("RingRadiusFit is Equatable", ring == ring && ring != sparse, true)
    }

    report.section("0.1.0: progress and isCancelled") { r in
        let always: @Sendable () -> Bool = { true }
        let never: @Sendable () -> Bool = { false }
        // progress: monotone within [0, 1], ending at 1.
        let p1 = ProgressLog()
        _ = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                         progress: { p1.add($0) })
        r.equal("fit1D reports progress, monotone in [0, 1], ending at 1", p1.isMonotoneWithinUnit && p1.all.last == 1, true)
        let p2 = ProgressLog()
        _ = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                                    requireRingEvidence: true, progress: { p2.add($0) })
        r.equal("fitBestAvailable reports progress, monotone in [0, 1], ending at 1", p2.isMonotoneWithinUnit && p2.all.last == 1, true)
        let p3 = ProgressLog()
        _ = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                    d: d, highTensionVolts: volts, pixelPitchMM: nil, progress: { p3.add($0) })
        r.equal("RingCalibration.fit reports progress, monotone in [0, 1], ending at 1", p3.isMonotoneWithinUnit && p3.all.last == 1, true)
        let p4 = ProgressLog()
        _ = try RingCalibration.fitRing(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                        progress: { p4.add($0) })
        r.equal("fitRing reports progress, monotone in [0, 1], ending at 1", p4.isMonotoneWithinUnit && p4.all.last == 1, true)
        let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        let p5 = ProgressLog()
        _ = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 60, centerQY: 60, innerRadius: 28,
                                                    outerRadius: 52, progress: { p5.add($0) })
        r.equal("fitAmorphousRing reports progress, monotone in [0, 1], ending at 1", p5.isMonotoneWithinUnit && p5.all.last == 1, true)

        // isCancelled: true throws CancellationError, from every entry point.
        r.equal("fit1D(isCancelled:) throws CancellationError", throwsCancellation {
            _ = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52, isCancelled: always) }, true)
        r.equal("fitBestAvailable(isCancelled:) throws CancellationError", throwsCancellation {
            _ = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                                        requireRingEvidence: true, isCancelled: always) }, true)
        r.equal("fitAmorphousRing(isCancelled:) throws CancellationError", throwsCancellation {
            _ = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 60, centerQY: 60, innerRadius: 28,
                                                        outerRadius: 52, isCancelled: always) }, true)
        r.equal("RingCalibration.fit(isCancelled:) throws CancellationError", throwsCancellation {
            _ = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                        d: d, highTensionVolts: volts, pixelPitchMM: nil, isCancelled: always) }, true)
        r.equal("fitRing(isCancelled:) throws CancellationError", throwsCancellation {
            _ = try RingCalibration.fitRing(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52, isCancelled: always) }, true)
        r.equal("a false isCancelled changes nothing", try RingCalibration.fitRing(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28,
                                                                                outerRadius: 52, isCancelled: never),
                try RingCalibration.fitRing(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52))
    }
}
