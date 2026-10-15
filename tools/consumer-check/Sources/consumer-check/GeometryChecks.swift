//
//  GeometryChecks.swift
//  Role: DiffractoGeometry's public API (slice 3, DiffractoKit 0.0.3) used from outside the package, against values a
//        reader verifies by hand: lambda at 200 kV, one Au {111} ring at R 150 px, indexing a q at Au {111}'s g, the beta''
//        [0 0 1] axis, a clean ring rendered here and recovered, and the 2:1 ellipse's closed-form correction.
//        `runGeometryChecks` is called from main.swift.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

/// A clean Gaussian ring of radius `radius` and width `sigma` (pixels) about (row, column) in a size x size pattern.
func plantedRing(size: Int, centerRow: Double, centerColumn: Double, radius: Double, sigma: Double,
                         keepAngle: ((Double) -> Bool)? = nil) -> DiffractionPattern {
    var pixels = [Float](repeating: 0, count: size * size)
    for row in 0..<size {
        for column in 0..<size {
            let dr = Double(row) - centerRow, dc = Double(column) - centerColumn
            let r = hypot(dr, dc)
            if let keep = keepAngle, !keep(atan2(dr, dc)) { continue }
            pixels[row * size + column] = Float(100 * exp(-(r - radius) * (r - radius) / (2 * sigma * sigma)))
        }
    }
    return DiffractionPattern(width: size, height: size, pixels: pixels)
}

/// The refusal a typed-throwing indexing call gives, nil when it gives an assignment.
private func refusal(_ body: () throws(RingIndexing.Refusal) -> Void) -> RingIndexing.Refusal? {
    do throws(RingIndexing.Refusal) { try body() } catch { return error }
    return nil
}

func runGeometryChecks(_ report: inout Report) {
    let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: fcc.map { AtomSite(z: 79, fractional: $0) })

    report.section("ElectronWavelength") { r in
        let lambda200 = ElectronWavelength.angstroms(volts: 200_000)
        r.near("lambda at 200 kV = 0.025079340450548003 A", lambda200, 0.025079340450548003, 1e-17)
        r.near("lambda at 80 kV", ElectronWavelength.angstroms(volts: 80_000), 0.041757160772834212, 1e-16)
        r.near("lambda at 120 kV", ElectronWavelength.angstroms(volts: 120_000), 0.033492152726723738, 1e-16)
        r.near("lambda at 300 kV", ElectronWavelength.angstroms(volts: 300_000), 0.019687489006848796, 1e-16)
        // lambda = lambda_classical / sqrt(1 + eV / (2 m0 c^2)): the closed form a reader checks by hand.
        let eV = 200_000.0
        let classical = ElectronWavelength.planckTimesLightEVAngstroms / (2 * ElectronWavelength.electronRestEnergyEV * eV).squareRoot()
        r.near("lambda = classical / sqrt(1 + eV/2m0c^2)", lambda200,
               classical / (1 + eV / (2 * ElectronWavelength.electronRestEnergyEV)).squareRoot(), 1e-16)
        r.near("m0 c^2 / e = 510998.95 eV", ElectronWavelength.electronRestEnergyEV, 510_998.95, 1)
        r.near("hc/e = 12398.4198 eV A", ElectronWavelength.planckTimesLightEVAngstroms, 12_398.4198, 1e-4)
        r.equal("CODATA constants are exact (h, c, e)",
                ElectronWavelength.planckConstant == 6.626_070_15e-34 && ElectronWavelength.speedOfLight == 299_792_458
                    && ElectronWavelength.elementaryCharge == 1.602_176_634e-19, true)
        r.near("electron mass 9.1093837015e-31 kg", ElectronWavelength.electronMass, 9.109_383_7015e-31, 1e-40)
        r.equal("zero volts is NaN", ElectronWavelength.angstroms(volts: 0).isNaN, true)
        r.equal("negative volts is NaN", ElectronWavelength.angstroms(volts: -1).isNaN, true)
        r.equal("NaN volts is NaN", ElectronWavelength.angstroms(volts: .nan).isNaN, true)
        r.equal("overflowing volts is NaN, not 0", ElectronWavelength.angstroms(volts: 1e200).isNaN, true)
        r.equal("RingCalibration.wavelengthAngstroms is the same function",
                RingCalibration.wavelengthAngstroms(volts: 200_000) == lambda200, true)
    }

    report.section("RingCalibration closed forms and calibrate") { r in
        let d = 2.3546, radius = 150.0, volts = 200_000.0
        let lambda = ElectronWavelength.angstroms(volts: volts)
        let twoTheta = RingCalibration.twoTheta(wavelengthAngstroms: lambda, d: d)
        r.near("Au {111} 2 theta = 0.010651260937060768", twoTheta, 0.010651260937060768, 1e-15)
        r.near("2 theta = 2 arcsin(lambda / 2d)", twoTheta, 2 * asin(lambda / (2 * d)), 1e-18)
        let lengthPx = RingCalibration.cameraLengthPixels(radius: radius, twoTheta: twoTheta)
        r.near("L_px = 14082.307097295722 within 1e-9", lengthPx, 14082.307097295722, 1e-9)
        r.near("L_px = R / tan 2 theta", lengthPx, radius / tan(twoTheta), 1e-9)
        // Bragg's small-angle form is 4.25e-5 relative above the exact length at 200 kV Au {111}.
        r.near("exact L is below R d / lambda by 1.5 (lambda/2d)^2", 1 - lengthPx / (radius * d / lambda),
               1.5 * (lambda / (2 * d)) * (lambda / (2 * d)), 5e-8)
        let dq = RingCalibration.qPerPixel(d: d, radius: radius)
        r.near("delta q = (1/d)/R = 0.0028313372405787253", dq, 0.0028313372405787253, 1e-17)
        r.near("L_mm at 14 um = 197.15229936214011", RingCalibration.cameraLengthMillimetres(pixelPitchMM: 0.014, cameraLengthPixels: lengthPx),
               197.15229936214011, 1e-9)

        let cal = try RingCalibration.calibrate(radius: radius, d: d, highTensionVolts: volts, pixelPitchMM: 0.014)
        r.near("calibrate: delta q", cal.deltaQ, dq, 0)
        r.near("calibrate: L_px", try unwrap(cal.lengthPixels), 14082.307097295722, 1e-9)
        r.near("calibrate: L_mm", try unwrap(cal.lengthMillimetres), 197.15229936214011, 1e-9)
        r.near("calibrate: 2 theta", try unwrap(cal.twoTheta), 0.010651260937060768, 1e-15)
        r.near("calibrate: lambda", try unwrap(cal.wavelengthAngstroms), 0.025079340450548003, 1e-17)
        r.equal("calibrate carries the high tension it used", cal.highTensionVolts, volts)
        r.equal("calibrate echoes d and R", cal.d == d && cal.radius == radius, true)
        r.equal("no refusal on a good ring", cal.lengthRefusal, nil)

        let noVolts = try RingCalibration.calibrate(radius: radius, d: d, highTensionVolts: nil)
        r.equal("no high tension: delta q still given", noVolts.deltaQ, dq)
        r.equal("no high tension: L absent, named", noVolts.lengthRefusal, RingCalibration.Refusal.noHighTension)
        r.equal("no high tension: nothing defaulted", noVolts.highTensionVolts == nil && noVolts.lengthPixels == nil, true)
        let noPitch = try RingCalibration.calibrate(radius: radius, d: d, highTensionVolts: volts)
        r.equal("no pixel pitch: L_px given, L_mm nil", noPitch.lengthPixels != nil && noPitch.lengthMillimetres == nil, true)

        let bragg = try RingCalibration.calibrate(radius: radius, d: 0.01, highTensionVolts: volts)
        var braggLimit = false
        if case .braggLimit? = bragg.lengthRefusal { braggLimit = true }
        r.equal("Bragg limit: d = 0.01 A at 200 kV refuses L", braggLimit, true)
        r.equal("Bragg limit: delta q still given", bragg.deltaQ, (1 / 0.01) / radius)
        r.near("lambda/2d limit is 1/sqrt 2", RingCalibration.maximumLambdaOver2d, 0.5.squareRoot(), 0)
        r.equal("twoTheta at the limit is NaN", RingCalibration.twoTheta(wavelengthAngstroms: 0.7071068, d: 0.5).isNaN, true)
        r.equal("cameraLengthPixels outside (0, pi/2) is NaN", RingCalibration.cameraLengthPixels(radius: 10, twoTheta: 1.6).isNaN, true)
        r.equal("qPerPixel of a zero radius is NaN", RingCalibration.qPerPixel(d: 2, radius: 0).isNaN, true)
        r.equal("cameraLengthMillimetres of a negative pitch is NaN",
                RingCalibration.cameraLengthMillimetres(pixelPitchMM: -1, cameraLengthPixels: 1).isNaN, true)

        var thrown: RingCalibration.Refusal?
        do { _ = try RingCalibration.calibrate(radius: 0, d: d, highTensionVolts: volts) } catch let e as RingCalibration.Refusal { thrown = e } catch {}
        r.equal("a zero radius throws invalidInput", thrown?.name, "invalidInput")
        thrown = nil
        do { _ = try RingCalibration.calibrate(radius: 1, d: d, highTensionVolts: -5) } catch let e as RingCalibration.Refusal { thrown = e } catch {}
        r.equal("a negative high tension throws invalidInput", thrown?.name, "invalidInput")
        r.equal("registered bounds", RingCalibration.angularBinCount == 36 && RingCalibration.minimumOccupiedAngularBins == 30
                && RingCalibration.maximumResidual == 0.22 && RingCalibration.maximumEllipticity == 0.093, true)
        r.equal("the model enum", EllipseCalibrationModel.asymmetricGaussian.rawValue, "Asymmetric Gaussian ring")
        r.equal("the conic fallback", EllipseCalibrationModel.conic.rawValue, "Intensity-weighted conic")
    }

    report.section("RingCalibration.FitReport and refusals") { r in
        func report(model: EllipseCalibrationModel = .asymmetricGaussian, a: Double = 100, b: Double = 100, bins: Int = 36,
                    sparse: Bool = false, residual: Double = 0.05, sigma: Double? = 2, width: Double = 24) -> RingCalibration.FitReport {
            RingCalibration.FitReport(model: model, semiAxisA: a, semiAxisB: b, occupiedAngularBins: bins, sparseCoverage: sparse,
                                      normalizedResidual: residual, innerSigma: sigma, outerSigma: sigma, annulusWidth: width)
        }
        let ok = report()
        r.near("report radius = (a + b)/2", report(a: 110, b: 90).radius, 100, 0)
        r.near("report ellipticity = a/b - 1 (a the larger)", report(a: 90, b: 110).ellipticity, 110.0 / 90 - 1, 1e-15)
        r.equal("a clean fit with a high tension has no refusal", RingCalibration.refusals(ok, highTensionVolts: 200_000), [])
        r.equal("no high tension is the only refusal", RingCalibration.refusals(ok, highTensionVolts: nil), [.noHighTension])
        r.equal("conic model refused", RingCalibration.refusals(report(model: .conic, sigma: nil), highTensionVolts: 200_000),
                [.modelNotAsymmetricGaussian])
        r.equal("29 of 36 bins refused, 30 accepted",
                RingCalibration.refusals(report(bins: 29), highTensionVolts: 1).map(\.name) == ["insufficientAngularCoverage"]
                    && RingCalibration.refusals(report(bins: 30), highTensionVolts: 1).isEmpty, true)
        r.equal("a sparse-coverage assertion is refused even at 36 bins",
                RingCalibration.refusals(report(sparse: true), highTensionVolts: 1).map(\.name), ["insufficientAngularCoverage"])
        r.equal("sigma at the annulus width refused",
                RingCalibration.refusals(report(sigma: 24, width: 24), highTensionVolts: 1).map(\.name), ["sigmaAtAnnulusWidth"])
        r.equal("residual 0.23 refused, 0.22 accepted",
                RingCalibration.refusals(report(residual: 0.23), highTensionVolts: 1).map(\.name) == ["excessiveResidual"]
                    && RingCalibration.refusals(report(residual: 0.22), highTensionVolts: 1).isEmpty, true)
        r.equal("ellipticity 0.094 refused, 0.093 accepted",
                RingCalibration.refusals(report(a: 109.4, b: 100), highTensionVolts: 1).map(\.name) == ["excessiveEllipticity"]
                    && RingCalibration.refusals(report(a: 109.3, b: 100), highTensionVolts: 1).isEmpty, true)
        r.equal("refusals come in the registered order",
                RingCalibration.refusals(report(model: .conic, a: 120, bins: 10, residual: 0.5, sigma: 30), highTensionVolts: nil).map(\.name),
                ["modelNotAsymmetricGaussian", "insufficientAngularCoverage", "sigmaAtAnnulusWidth", "excessiveResidual",
                 "excessiveEllipticity", "noHighTension"])
        r.equal("a NaN number gives invalidInput alone",
                RingCalibration.refusals(report(residual: .nan), highTensionVolts: nil).map(\.name), ["invalidInput"])
        r.equal("withholdsDeltaQ is false only for no high tension and the Bragg limit",
                [RingCalibration.Refusal.noHighTension, .braggLimit(lambdaOver2d: 1)].allSatisfy { !$0.withholdsDeltaQ }
                    && [RingCalibration.Refusal.modelNotAsymmetricGaussian, .excessiveResidual(1), .invalidInput("x"),
                        .insufficientAngularCoverage(occupiedBins: 1), .excessiveEllipticity(1),
                        .sigmaAtAnnulusWidth(sigma: 1, annulusWidth: 1)].allSatisfy(\.withholdsDeltaQ), true)
        r.equal("a refusal says why", RingCalibration.Refusal.noHighTension.reason.contains("Δq"), true)
        r.equal("a refusal is a LocalizedError", RingCalibration.Refusal.noHighTension.errorDescription
                == RingCalibration.Refusal.noHighTension.reason, true)
    }

    report.section("RingIndexing") { r in
        let list = try gold.kinematicReflections(kMax: 1.2)
        let g111 = 3.0.squareRoot() / 4.0782
        let sigma = 1e-3 * g111
        let hit = try RingIndexing.index(q: g111, sigma: sigma, in: list)
        r.equal("q = g(111), sigma = 1e-3 q: accepted as {1 1 1}", hit.best.family.name, "{1 1 1}")
        r.near("z = 0 at the ring", hit.best.z, 0, 1e-12)
        r.near("deviation 0 at the ring", hit.best.deviation, 0, 1e-12)
        r.equal("the next family is {2 0 0}", hit.next?.family.name, "{2 0 0}")
        // z of the next = |q - g200| / sigma = (2/a - sqrt3/a) / sigma; margin = z^2 / 2
        let zNext = (2 / 4.0782 - g111) / sigma
        r.near("next z by hand", try unwrap(hit.next).z, zNext, 1e-9 * zNext)
        r.near("margin = z_next^2 / 2", try unwrap(hit.margin), zNext * zNext / 2, 1e-9 * zNext * zNext)
        r.equal("margin exceeds ln 19", try unwrap(hit.margin) > RingIndexing.minimumMargin, true)
        r.near("ln 19", RingIndexing.minimumMargin, log(19.0), 0)
        r.equal("default coverage is 3", RingIndexing.defaultCoverage == 3 && hit.coverage == 3, true)
        r.near("reach = sqrt(k^2 + 2 ln 19)", RingIndexing.reach(coverage: 3), (9 + 2 * log(19.0)).squareRoot(), 1e-15)
        r.equal("assignment echoes q, sigma", hit.q == g111 && hit.sigma == sigma, true)
        r.near("measured d = 1/q", hit.dMeasured, 1 / g111, 1e-12)
        r.equal("validation says none", hit.validation == RingIndexing.validation && RingIndexing.validation.hasPrefix("none"), true)
        r.near("relative deviation of a q 0.5 sigma high", try RingIndexing.index(q: g111 * (1 + 5e-4), sigma: sigma, in: list)
            .best.relativeDeviation, 5e-4 / (1 + 5e-4), 1e-9)
        r.equal("the list overload equals the families overload",
                try RingIndexing.index(q: g111, sigma: sigma, in: list.families, kMax: list.kMax) == hit, true)
        let rings = list.families.map(\.family)
        r.equal("the ReflectionFamily overload equals it",
                try RingIndexing.index(q: g111, sigma: sigma, in: rings, kMax: 1.2) == hit, true)
        r.equal("a Reflection list's own families index too",
                try RingIndexing.index(q: g111, sigma: sigma, in: try gold.reflections(kMax: 1.2).families(), kMax: 1.2).best.family.name, "{1 1 1}")
        // The apps' macTEM T1 plant: a ring at 3 sigma exactly is accepted, at 3.01 sigma refused as none-within-tolerance.
        func verdict(_ q: Double, _ s: Double) -> String {
            refusal { () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: q, sigma: s, in: list) }?.name ?? "accepted"
        }
        r.equal("a q midway between {1 1 1} and {2 0 0} is ambiguous at a wide sigma", verdict((g111 + 2 / 4.0782) / 2, 0.05 * g111), "ambiguous")
        r.equal("a q between rings at a tight sigma is none within tolerance", verdict(0.45, 1e-4), "noneWithinTolerance")
        r.equal("a q at k_max is beyond k_max", verdict(1.2, 1e-3), "beyondKMax")
        r.equal("a negative q is invalid", verdict(-1, 1e-3), "invalidInput")
        r.equal("a zero sigma is invalid", verdict(g111, 0), "invalidInput")
        let nanCoverage = refusal { () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: g111, sigma: sigma, in: list, coverage: Double.nan) }
        r.equal("a NaN coverage is invalid", nanCoverage?.name, "invalidInput")
        let ambiguous = refusal { () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: (g111 + 2 / 4.0782) / 2, sigma: 0.05 * g111, in: list) }
        if case .ambiguous(let candidates)? = ambiguous {
            r.equal("ambiguous names two candidates by (z, g)", candidates.count >= 2 && candidates[0].z <= candidates[1].z, true)
        } else { r.check("ambiguous names its candidates", false, "got \(String(describing: ambiguous))") }
        r.equal("a refusal says why", ambiguous?.reason.hasPrefix("More than one family") ?? false, true)
        let none = refusal { () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: 0.45, sigma: 1e-4, in: list) }
        if case .noneWithinTolerance(let nearest)? = none {
            r.equal("none-within-tolerance names the nearest", nearest?.family.name, "{1 1 1}")
        } else { r.check("none-within-tolerance names the nearest", false, "got \(String(describing: none))") }
        let empty = refusal { () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: 0.45, sigma: 1e-3, in: [KinematicFamily](), kMax: 1.2) }
        if case .noneWithinTolerance(let nearest)? = empty { r.equal("an empty list has no nearest", nearest == nil, true) }
        else { r.check("an empty list has no nearest", false, "got \(String(describing: empty))") }
        let manual = RingIndexing.Candidate(family: rings[0], z: 1, deviation: -0.02)
        r.equal("Candidate is constructible, |deviation| is relativeDeviation", manual.relativeDeviation, 0.02)
        if case .beyondKMax(let reach, let kMax)? = refusal({ () throws(RingIndexing.Refusal) in _ = try RingIndexing.index(q: 1.2, sigma: 1e-3, in: list) }) {
            r.equal("beyond k_max carries its reach", reach > kMax && kMax == 1.2, true)
        } else { r.check("beyond k_max carries its reach", false) }
    }

    report.section("ZoneAxis") { r in
        r.equal("49 low-index axes", ZoneAxis.lowIndex.count, 49)
        r.equal("the first and last axes (mac4DSTEM's order)", ZoneAxis.lowIndex.first == ZoneAxis(-1, -1, -1) && ZoneAxis.lowIndex.last == ZoneAxis(0, 0, -1), true)
        r.equal("axes are distinct", Set(ZoneAxis.lowIndex).count, 49)
        r.equal("description", ZoneAxis(1, 1, 0).description, "[1 1 0]")
        let beta = 105.3
        let cell = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: beta)
        let n001 = try unwrap(ZoneAxis(0, 0, 1).cartesian(in: cell))
        r.near("beta'' [0 0 1] x = cos beta", n001.x, cos(beta * .pi / 180), 1e-12)
        r.near("beta'' [0 0 1] y = 0", n001.y, 0, 1e-12)
        r.near("beta'' [0 0 1] z = sin beta", n001.z, sin(beta * .pi / 180), 1e-12)
        r.equal("a cartesian axis is a unit direction", ZoneAxis.isUnitDirection(n001), true)
        let n110 = try unwrap(ZoneAxis(1, 1, 0).cartesian(in: UnitCell.cubic(a: 4.0782)))
        r.near("cubic [1 1 0] = (1, 1, 0)/sqrt 2", n110.x, 0.5.squareRoot(), 1e-15)
        r.equal("[0 0 0] has no direction", ZoneAxis(0, 0, 0).cartesian(in: cell) == nil, true)
        r.equal("isUnitDirection refuses a long, a zero and a NaN vector",
                !ZoneAxis.isUnitDirection([0, 0, 2]) && !ZoneAxis.isUnitDirection([0, 0, 0]) && !ZoneAxis.isUnitDirection([.nan, 0, 1]), true)

        // [0 0 1] of a cubic cell: reference x, e1 = n x x = y, e2 = n x e1 = -x.
        let basis = ZoneAxis.detectorBasis([0, 0, 1])
        r.equal("detector basis of [0 0 1]: e1 = y", basis.columns.0, SIMD3(0, 1, 0))
        r.equal("detector basis of [0 0 1]: e2 = -x", basis.columns.1, SIMD3(-1, 0, 0))
        r.equal("detector basis of [0 0 1]: third column is n", basis.columns.2, SIMD3(0, 0, 1))
        r.near("detector basis is right-handed", basis.determinant, 1, 1e-15)
        let tilted = ZoneAxis.detectorBasis(simd_normalize(SIMD3(1, 0.3, -0.2)))
        r.near("a beam along x uses y as the reference: orthonormal", simd_length(tilted.columns.0) + simd_length(tilted.columns.1), 2, 1e-12)
        r.near("e1 . e2 = 0", simd_dot(tilted.columns.0, tilted.columns.1), 0, 1e-12)
        r.near("determinant +1 for every low-index axis of beta''", ZoneAxis.lowIndex.map { axis -> Double in
            (axis.cartesian(in: cell).map { ZoneAxis.detectorBasis($0).determinant } ?? 0) }.map { abs($0 - 1) }.max() ?? 1, 0, 1e-12)

        let list = try gold.kinematicReflections(kMax: 1.2)
        let flat = ZoneAxis.project(list.kept, onto: [0, 0, 1], sgMax: 1e-9)
        r.equal("flat sphere: the zero-order Laue zone is the l = 0 plane", !flat.isEmpty && flat.allSatisfy { $0.reflection.reflection.hkl.l == 0 && abs($0.sg) < 1e-12 }, true)
        r.equal("flat sphere: nothing off the zone law", flat.count, list.kept.filter { $0.reflection.hkl.l == 0 }.count)
        let f200 = try unwrap(flat.first { $0.reflection.reflection.hkl == MillerIndex(2, 0, 0) })
        r.near("(2 0 0) lands at x = g . e1 = 0", f200.x, 0, 1e-15)
        r.near("(2 0 0) lands at y = g . e2 = -2/a", f200.y, -2 / 4.0782, 1e-15)
        let lambda = ElectronWavelength.angstroms(volts: 200_000)
        let curved = ZoneAxis.project(list.kept, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: 0.3)
        let c200 = try unwrap(curved.first { $0.reflection.reflection.hkl == MillerIndex(2, 0, 0) })
        r.near("curved sphere: s_g = lambda |g|^2 / 2 for g . n = 0", c200.sg, lambda * (2 / 4.0782) * (2 / 4.0782) / 2, 1e-15)
        let c111 = try unwrap(curved.first { $0.reflection.reflection.hkl == MillerIndex(1, 1, 1) })
        r.near("curved sphere: s_g = g . n + lambda |g|^2 / 2 for (1 1 1)", c111.sg,
               1 / 4.0782 + lambda * (3 / (4.0782 * 4.0782)) / 2, 1e-15)
        r.equal("a non-unit n keeps nothing", ZoneAxis.project(list.kept, onto: [0, 0, 2], sgMax: 1).isEmpty, true)
        r.equal("a negative wavelength keeps nothing", ZoneAxis.project(list.kept, onto: [0, 0, 1], wavelengthAngstroms: -1, sgMax: 1).isEmpty, true)
        r.equal("sgMax bounds the slab", ZoneAxis.project(list.kept, onto: [0, 0, 1], sgMax: 0.3).allSatisfy { abs($0.sg) <= 0.3 }, true)
    }

    report.section("DiffractionPattern and EllipseTransform") { r in
        let p = DiffractionPattern(width: 3, height: 2, pixels: [1, 2, 3, 4, 5, 6])
        r.equal("pattern fields", p.width == 3 && p.height == 2 && p.pixels.count == 6, true)
        r.equal("pixel at row 1, column 2 is pixels[1 * width + 2]", p.pixels[1 * p.width + 2], 6)
        r.equal("equatable", p == DiffractionPattern(width: 3, height: 2, pixels: [1, 2, 3, 4, 5, 6]), true)

        // a = 2, b = 1, theta = pi/2: s = 0, c = 1, t00 = 1, t01 = 0, t11 = b/a = 1/2; (x, y) = (t11 dx, t00 dy).
        let t = try unwrap(EllipseTransform(a: 2, b: 1, theta: .pi / 2))
        let corrected = t.correctedOffset(dx: 4, dy: 3)
        r.near("a 2:1 ellipse halves dx", Double(corrected.x), 2, 1e-6)
        r.near("and leaves dy", Double(corrected.y), 3, 1e-6)
        let back = t.uncorrectedOffset(dx: corrected.x, dy: corrected.y)
        r.near("uncorrected x returns dx", Double(back.x), 4, 1e-5)
        r.near("uncorrected y returns dy", Double(back.y), 3, 1e-5)
        r.equal("the transform keeps its triple", t.a == 2 && t.b == 1 && t.theta == .pi / 2, true)
        let circle = try unwrap(EllipseTransform(a: 5, b: 5, theta: 0.7))
        let same = circle.correctedOffset(dx: 3, dy: -2)
        r.near("a circle corrects nothing (x)", Double(same.x), 3, 1e-6)
        r.near("a circle corrects nothing (y)", Double(same.y), -2, 1e-6)
        r.equal("a NaN theta is not an ellipse", EllipseTransform(a: 2, b: 1, theta: .nan) == nil, true)
        r.equal("a missing axis is not an ellipse", EllipseTransform(a: nil, b: 1, theta: 0) == nil, true)
        r.equal("a negative axis is not an ellipse", EllipseTransform.accepted(a: -1, b: 1, theta: 0) == nil, true)
        r.equal("an ellipse is accepted as given", EllipseTransform.accepted(a: 2, b: 1, theta: 0.5)?.theta, 0.5)
        let identity = EllipseTransform.correctedOffset(dx: 4, dy: 3, a: nil, b: nil, theta: nil)
        r.equal("no ellipse corrects nothing", identity.x == 4 && identity.y == 3, true)
        let viaStatic = EllipseTransform.correctedOffset(dx: 4, dy: 3, a: 2, b: 1, theta: .pi / 2)
        r.equal("the static form equals the instance", viaStatic.x == corrected.x && viaStatic.y == corrected.y, true)
    }

    report.section("EllipseCalibration and RingCalibration.fit on a planted ring") { r in
        let size = 128, R = 40.0, s = 2.0, cr = 63.4, cc = 64.2
        let pattern = plantedRing(size: size, centerRow: cr, centerColumn: cc, radius: R, sigma: s)
        // The registered band is 1 % on the radius, 0.01 on b/a, 1 px on the centre (mac4DSTEM EllipseNoisyRingTests).
        let best = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("the profile model is used on a clean ring", best.model, .asymmetricGaussian)
        r.near("the clean ring's radius is recovered within 1 %", (best.a + best.b) / 2, R, 0.01 * R)
        r.near("the clean ring is circular: b/a within 0.01", best.a / best.b, 1, 0.01)
        r.near("centre row within 1 px", best.centerQX, cr, 1)
        r.near("centre column within 1 px", best.centerQY, cc, 1)
        r.equal("36 of 36 angular bins hold the ring", best.occupiedAngularBins, EllipseCalibration.angularBinCount)
        r.equal("the profile is reported, with a fitted width near the planted one",
                abs((best.profile?.innerSigma ?? 0) - s) < 0.5 && abs((best.profile?.outerSigma ?? 0) - s) < 0.5, true)
        r.equal("no fallback reason", best.profileFallbackReason, nil)
        r.equal("not an asserted sparse ring", best.sparseCoverage, false)

        let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52)
        r.equal("fit1D is the conic", conic.model, .conic)
        // The conic's radius is biased high by R * 1.75 (s/R)^2 = 0.175 px here (a property of the objective, documented).
        r.near("fit1D radius sits R * 1.75 (s/R)^2 above the plant", (conic.a + conic.b) / 2, R * (1 + 1.75 * (s / R) * (s / R)), 0.05)
        let refined = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 60, centerQY: 60,
                                                              innerRadius: 28, outerRadius: 52)
        r.equal("fitAmorphousRing is the profile fit", refined.model, .asymmetricGaussian)
        r.near("fitAmorphousRing recovers R", (refined.a + refined.b) / 2, R, 0.01 * R)
        r.equal("constants", EllipseCalibration.angularBinCount == 36 && EllipseCalibration.degeneracyBoundBins == 30
                && EllipseCalibration.sparseFloorBins == 12, true)

        // Refusals of the fit itself, by name.
        func fitError(_ body: () throws -> Void) -> EllipseCalibration.FitError? {
            do { try body() } catch let e as EllipseCalibration.FitError { return e } catch {}
            return nil
        }
        r.equal("an empty annulus is invalid input", {
            if case .invalidInput? = fitError({ _ = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 30, outerRadius: 30) }) { return true }
            return false }(), true)
        let blank = DiffractionPattern(width: size, height: size, pixels: [Float](repeating: 0, count: size * size))
        r.equal("a blank pattern has no ring signal", fitError({
            _ = try EllipseCalibration.fit1D(pattern: blank, centerQX: 64, centerQY: 64, innerRadius: 28, outerRadius: 52) }) != nil, true)
        let half = plantedRing(size: size, centerRow: cr, centerColumn: cc, radius: R, sigma: s, keepAngle: { $0 > 0 })
        if case .insufficientAngularCoverage(let bins)? = fitError({
            _ = try EllipseCalibration.fit1D(pattern: half, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52) }) {
            r.equal("a half ring covers 18 of 36 bins and is refused", bins <= 20, true)
        } else { r.check("a half ring is refused for coverage", false) }
        r.equal("the fit errors speak", EllipseCalibration.FitError.insufficientSignal.errorDescription?.hasPrefix("Cannot fit detector ellipse") ?? false, true)
        r.equal("the profile parameters are constructible", EllipseProfileParameters(centralIntensity: 1, ringIntensity: 2, centralSigma: 3,
                                                                                    innerSigma: 4, outerSigma: 5, background: 6).outerSigma, 5)

        // End to end: the fit -> report -> refusals -> calibrate chain, with the ring's own d and a 200 kV beam.
        let d = 2.3546
        let result = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                             d: d, highTensionVolts: 200_000, pixelPitchMM: 0.014)
        r.equal("a clean ring has no refusal", result.refusals, [])
        let report = try unwrap(result.report)
        r.near("the report's radius is the fit's (a + b)/2", report.radius, (try unwrap(result.fit).a + (try unwrap(result.fit)).b) / 2, 1e-12)
        r.equal("the report's annulus width is outer - inner", report.annulusWidth, 24)
        let calibration = try unwrap(result.calibration)
        r.near("delta q = (1/d) / fitted R", calibration.deltaQ, (1 / d) / report.radius, 1e-15)
        r.near("delta q is within 1 % of the planted ring's (1/d)/R", calibration.deltaQ, (1 / d) / R, 0.01 * (1 / d) / R)
        r.near("L_px = R / tan 2 theta at the fitted R", try unwrap(calibration.lengthPixels),
               report.radius / tan(2 * asin(ElectronWavelength.angstroms(volts: 200_000) / (2 * d))), 1e-9)
        let noVolts = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                              d: d, highTensionVolts: nil, pixelPitchMM: nil)
        r.equal("without a high tension: only .noHighTension, delta q kept", noVolts.refusals == [.noHighTension] && noVolts.calibration != nil, true)
        let sparse = try RingCalibration.fit(pattern: half, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                             d: d, highTensionVolts: 200_000, pixelPitchMM: nil)
        r.equal("a half ring is refused for coverage before any fit, no calibration",
                sparse.refusals.map(\.name) == ["insufficientAngularCoverage"] && sparse.calibration == nil && sparse.fit == nil, true)
        var thrown: RingCalibration.Refusal?
        do { _ = try RingCalibration.fit(pattern: pattern, centerQX: 60, centerQY: 60, innerRadius: 28, outerRadius: 52,
                                         d: 0, highTensionVolts: 200_000, pixelPitchMM: nil) } catch let e as RingCalibration.Refusal { thrown = e } catch {}
        r.equal("an impossible d throws invalidInput before any fit", thrown?.name, "invalidInput")
        let value = RingFitCalibration(fit: nil, report: nil, refusals: [.noHighTension], calibration: nil)
        r.equal("RingFitCalibration is constructible", value.refusals, [.noHighTension])
    }
}
