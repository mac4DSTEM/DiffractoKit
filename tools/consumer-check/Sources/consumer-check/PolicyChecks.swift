//
//  PolicyChecks.swift
//  Role: the public API DiffractoKit 0.2.0 adds, used from outside the package: every `Policy` type with its `.default`, one
//        non-default use that moves a result a reader can verify by hand, and the refusals of an invalid value; the position
//        tolerance's domain; the metric tolerances; `RingIndexing`'s `minimumMargin:`; `retainingAll:`. `runPolicyChecks` is
//        called from main.swift. (Negative controls: controls.py.)
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

private let goldSites: [AtomSite] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]].map {
    AtomSite(z: 79, fractional: SIMD3<Double>($0[0], $0[1], $0[2]))
}

private let goldCIF = """
data_gold
_cell_length_a 4.0782
_cell_length_b 4.0782
_cell_length_c 4.0782
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Au1 Au 0 0 0
Au2 Au 0.5 0.5 0
Au3 Au 0.5 0 0.5
Au4 Au 0 0.5 0.5
"""

func runPolicyChecks(_ report: inout Report) {
    let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: goldSites)

    report.section("0.2.0: ReflectionList.Policy") { r in
        let d = ReflectionList.Policy.default
        r.equal("default k_max cap", d.maximumKMax, 4)
        r.equal("default reflection cap", d.maximumReflections, 10_000_000)
        r.equal("default enumeration cap", d.maximumEnumeration, 100_000_000)
        r.equal("default family tolerance", d.familyTolerance, 1e-9)
        r.equal("default absence tolerance", d.absenceTolerance, 1e-9)
        r.equal("an empty initialiser is the default", ReflectionList.Policy(), d)
        var refused = false
        do { _ = try gold.reflections(kMax: 5) } catch let e as ReflectionListError { if case .invalidKMax = e { refused = true } } catch {}
        r.equal("k_max 5 is refused at the default cap", refused, true)
        let raised = try gold.reflections(kMax: 5, policy: ReflectionList.Policy(maximumKMax: 6))
        r.check("k_max 5 is built at a cap of 6", raised.reflections.count > 0)
        r.equal("the list keeps its policy", raised.policy.maximumKMax, 6)
        let big = CrystalStructure.detectingSymmetry(cell: .cubic(a: 40), sites: goldSites)
        var tooLarge: ReflectionListError?
        do { _ = try big.reflections(kMax: 3.9, policy: ReflectionList.Policy(maximumReflections: 1000)) }
        catch let e as ReflectionListError { tooLarge = e } catch {}
        r.check("a reflection cap of 1000 refuses a 40 A cell at k_max 3.9", tooLarge != nil)
        r.check("the refusal names the count, the default limit and the policy field",
                tooLarge.map { "\($0)".contains("1e7") && "\($0)".contains("maximumReflections") } ?? false)
        let list = try gold.kinematicReflections(kMax: 1.5, policy: ReflectionList.Policy(maximumKMax: 2))
        r.check("the kinematic list takes the policy", list.kept.count > 0)
    }

    report.section("0.2.0: CrystalStructure.ValidationPolicy") { r in
        let d = CrystalStructure.ValidationPolicy.default
        r.equal("default contact limit", d.closeContactLimit, 0.5)
        r.equal("default split-site tolerance", d.splitSiteOccupancyTolerance, 0.02)
        r.equal("an empty initialiser is the default", CrystalStructure.ValidationPolicy(), d)
        let close = CrystalStructure(cell: .cubic(a: 4),
                                     sites: [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.075, 0, 0])],
                                     operations: [.identity])
        r.equal("0.3 A apart is not usable at the default", close.isUsable(under: .default), false)
        r.equal("... and is at a limit of 0.2", close.isUsable(under: CrystalStructure.ValidationPolicy(closeContactLimit: 0.2)), true)
        r.near("the shortest contact", close.shortestCloseContact(using: .default) ?? .nan, 0.3, 1e-9)
        r.equal("no contact under a limit of 0.2", close.shortestCloseContact(using: CrystalStructure.ValidationPolicy(closeContactLimit: 0.2)), nil)
        r.equal("validationIssues(using:) names the contact", close.validationIssues(using: .default).isEmpty, false)
        var refused = false
        do { _ = try close.kinematicReflections(kMax: 1) } catch let e as KinematicError { if case .invalidStructure = e { refused = true } } catch {}
        r.equal("the kinematic list refuses it", refused, true)
        let relaxed = try close.kinematicReflections(kMax: 1, validation: CrystalStructure.ValidationPolicy(closeContactLimit: 0.2))
        r.check("... and builds under the relaxed policy", relaxed.all.count > 0)
    }

    report.section("0.2.0: SymmetryDetection.Policy, Tolerance and SymmetryOperation(xyz:maximumMagnitude:)") { r in
        let d = SymmetryDetection.Policy.default
        r.equal("default occupancy tolerance", d.occupancyTolerance, 0.02)
        r.equal("default translation grid", d.translationGrid, 1.0 / 24)
        r.equal("default snap factor", d.gridSnapFactor, 3)
        r.equal("default precision limit factor", d.precisionLimitFactor, 10)
        // two sites on one position whose occupancies differ by 0.01: one species at the default tolerance, two at 0
        let sites = [AtomSite(z: 13, fractional: [0, 0, 0], occupancy: 0.50), AtomSite(z: 13, fractional: [0.5, 0.5, 0.5], occupancy: 0.51)]
        let loose = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4), sites: sites, conventional: false)
        let strict = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4), sites: sites, conventional: false,
                                                        policy: SymmetryDetection.Policy(occupancyTolerance: 0))
        r.check("a stricter occupancy tolerance finds fewer operations", strict.operations.count < loose.operations.count)
        let result = SymmetryDetection.detect(sites: sites, cell: .cubic(a: 4), policy: .default)
        r.equal("detect(policy:) agrees with detectingSymmetry", result.operations.count, loose.operations.count)

        let cell = UnitCell.cubic(a: 4)
        for bad in [SymmetryDetection.Tolerance.fractional(-1), .fractional(.nan), .fractional(0), .fractional(1e-300),
                    .fractional(0.5), .cartesian(.infinity)] {
            r.check("\(bad) is outside the domain", bad.refusal(for: cell) != nil && !bad.isValid(for: cell))
            let s = CrystalStructure.detectingSymmetry(cell: cell, sites: goldSites, tolerance: bad)
            r.check("... and the structure is the identity alone, unusable, with the reason",
                    s.operations == [.identity] && s.symmetryRefusal != nil && !s.isUsable)
        }
        r.check("the default is in the domain", SymmetryDetection.Tolerance.default.isValid(for: cell))
        r.equal("... and detects gold's 192 operations", CrystalStructure.detectingSymmetry(cell: cell, sites: goldSites, tolerance: .default).operations.count, 192)

        r.equal("default magnitude", SymmetryOperation.defaultMaximumMagnitude, 1000)
        var refused = false
        do { _ = try SymmetryOperation(xyz: "x+2000,y,z") } catch SymmetryOperationError.unparsable { refused = true } catch {}
        r.equal("a number of 2000 is unparsable at the default", refused, true)
        let wide = try SymmetryOperation(xyz: "x+2000,y,z", maximumMagnitude: 5000)
        r.equal("... and parsed at 5000, as the identity mod 1", wide.isIdentity, true)
    }

    report.section("0.2.0: CIFImport.Policy") { r in
        let d = CIFImport.Policy.default
        r.equal("default merge distance", d.maximumMergeDistance, 0.10)
        r.equal("default minimum merge tolerance", d.minimumMergeTolerance, 1e-6)
        r.equal("default half step cap", d.maximumCoordinateHalfStep, 1e-2)
        r.equal("default right-angle tolerance", d.rightAngleToleranceDeg, 0.05)
        r.equal("default operation magnitude", d.maximumOperationMagnitude, 1000)
        r.equal("default position tolerance", d.positionTolerance, SymmetryDetection.Tolerance.default)
        r.equal("default metric length tolerance", d.metricLengthTolerance, 1e-9)
        r.equal("default metric angle tolerance", d.metricAngleToleranceDeg, 1e-9)
        r.equal("default nested policies", d.symmetry == .default && d.validation == .default, true)
        let plain = try CIFImport.read(goldCIF, fileBaseName: "gold")
        let withDefault = try CIFImport.read(goldCIF, fileBaseName: "gold", policy: .default)
        r.equal("a read given the default equals a read given none", plain, withDefault)
        r.equal("gold has 192 operations", plain.structure.operations.count, 192)
        var tight = CIFImport.Policy.default
        tight.positionTolerance = .fractional(-1)
        var refusal: CIFImportError?
        do { _ = try CIFImport.read(goldCIF, policy: tight) } catch let e as CIFImportError { refusal = e } catch {}
        if case .invalidStructure(let messages)? = refusal {
            r.check("an unusable position tolerance is refused with its name", messages.contains { $0.contains("position tolerance") })
        } else {
            r.check("an unusable position tolerance is refused as invalidStructure", false, "\(String(describing: refusal))")
        }
        var nan = CIFImport.Policy.default
        nan.metricLengthTolerance = .nan
        let zero = CIFImport.Policy(metricLengthTolerance: 0)
        r.equal("a NaN metric tolerance reads as 0", try CIFImport.read(goldCIF, policy: nan), try CIFImport.read(goldCIF, policy: zero))
        var strict = CIFImport.Policy.default
        strict.validation = CrystalStructure.ValidationPolicy(closeContactLimit: 5)
        var contact = false
        do { _ = try CIFImport.read(goldCIF, policy: strict) } catch let e as CIFImportError { if case .invalidStructure = e { contact = true } } catch {}
        r.equal("a contact limit of 5 A refuses gold (2.88 A)", contact, true)
    }

    report.section("0.2.0: EllipseCalibration.Policy, RingEvidence.Policy, RingCalibration.Policy") { r in
        let e = EllipseCalibration.Policy.default
        r.equal("ellipse default bins", [e.angularBinCount, e.degeneracyBoundBins, e.sparseFloorBins], [36, 30, 12])
        r.equal("ellipse default residual caps", [e.maximumConicResidual, e.maximumProfileResidual, e.profileSeedResidual], [0.2, 0.22, 0.4])
        r.equal("ellipse default axis ratio, one-ring ratio", [e.maximumAxisRatio, e.oneRingRatio], [5, 1.10])
        r.equal("ellipse default floors", [e.ringWidthFloor, e.ringWidthFloorTolerance, e.loweredRingWidthFloor], [0.25, 0.05, 0.01])
        r.equal("ellipse default is valid", e.isValid, true)
        r.equal("ellipse default evidence", e.evidence, RingEvidence.Policy.default)
        let derived = EllipseCalibration.Policy(angularBinCount: 72)
        r.equal("the initialiser derives the other two bin counts", [derived.degeneracyBoundBins, derived.sparseFloorBins], [60, 24])
        var mutated = EllipseCalibration.Policy.default
        mutated.angularBinCount = 72
        r.equal("a mutated bin count keeps 30 and 12", [mutated.degeneracyBoundBins, mutated.sparseFloorBins], [30, 12])
        r.equal("any Int is accepted without a trap, and is invalid outside 6...3600",
                [Int.max, Int.min, 0, 5, 3601].map { EllipseCalibration.Policy(angularBinCount: $0).isValid }, [false, false, false, false, false])

        let p = RingEvidence.Policy.default
        r.equal("evidence defaults (pixels, fraction, significance, bins)", [Double(p.minimumSectorPixels), p.minimumAmplitudeFraction, p.minimumSignificance, Double(p.minimumOccupiedBins)], [8, 0.5, 3, 30])
        r.equal("evidence defaults (drop, shift, shift significance, secondary peak)", [p.minimumDrop, p.maximumShift, p.shiftSignificance, p.maximumSecondaryPeak], [0.25, 0.25, 4, 0.25])
        r.equal("evidence default is valid", p.isValid, true)

        let c = RingCalibration.Policy.default
        r.equal("calibration defaults", [Double(c.angularBinCount), Double(c.minimumOccupiedAngularBins), c.maximumResidual, c.maximumEllipticity], [36, 30, 0.22, 0.093])
        r.equal("calibration default evidence", c.evidence, p)
        r.equal("calibration default is valid", c.isValid, true)
        r.equal("an invalid calibration policy is not valid", RingCalibration.Policy(maximumResidual: .nan).isValid, false)

        // a ring with a 120 degree gap: refused at the default coverage bound, fitted at a bound of 20
        let size = 128
        let gapped = plantedRing(size: size, centerRow: 63.4, centerColumn: 64.2, radius: 40, sigma: 2.5,
                                 keepAngle: { !($0 > 0 && $0 < 2 * .pi / 3) })
        var refused = false
        do { _ = try EllipseCalibration.fit1D(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55) }
        catch let error as EllipseCalibration.FitError { if case .insufficientAngularCoverage = error { refused = true } } catch {}
        r.equal("the 120 degree gap is refused at the default", refused, true)
        let loose = EllipseCalibration.Policy(degeneracyBoundBins: 20)
        let fit = try EllipseCalibration.fit1D(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55, policy: loose)
        r.near("... and fitted under a bound of 20 bins (the radius of the arc it covers)", (fit.a + fit.b) / 2, 40, 1.5)
        var invalid = false
        do { _ = try EllipseCalibration.fit1D(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55, policy: EllipseCalibration.Policy(angularBinCount: Int.max)) }
        catch let error as EllipseCalibration.FitError { if case .invalidInput = error { invalid = true } } catch {}
        r.equal("an invalid ellipse policy is refused with invalidInput", invalid, true)

        // the ring evidence under a policy; fitRing with the two policies
        let whole = plantedRing(size: size, centerRow: 63.4, centerColumn: 64.2, radius: 40, sigma: 2.5)
        let ring = try RingCalibration.fitRing(pattern: whole, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55,
                                               policy: .default, ellipsePolicy: .default)
        r.equal("fitRing with the defaults is usable", ring.isUsable, true)
        let evidence = try unwrap(ring.report?.ringEvidence)
        r.equal("the evidence passes the default policy", evidence.passes(.default), true)
        r.equal("... and fails a policy demanding a drop of 1.5 amplitudes on both sides", evidence.isPeak(RingEvidence.Policy(minimumDrop: 1.5)), false)
        let measured = EllipseCalibration.ringEvidence(pattern: whole, fit: try unwrap(ring.fit), centerQX: 63, centerQY: 64,
                                                       innerRadius: 25, outerRadius: 55, policy: RingEvidence.Policy(minimumSectorPixels: 100_000))
        r.equal("a sector needing 100000 pixels reads no sector", measured.occupiedBins, 0)
        r.equal("a refusal under a strict calibration policy names the residual",
                RingCalibration.refusals(try unwrap(ring.report), highTensionVolts: 200_000, policy: RingCalibration.Policy(maximumResidual: 1e-9)).map(\.name),
                ["excessiveResidual"])
        r.check("the reason quotes the policy's bound",
                RingCalibration.Refusal.excessiveResidual(0.3).reason(policy: RingCalibration.Policy(maximumResidual: 0.5)).contains("0.50"))
        let gappedRing = try RingCalibration.fitRing(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55,
                                                     policy: .default, ellipsePolicy: .default)
        r.equal("a gapped ring is refused for coverage by fitRing", gappedRing.refusals.map(\.name), ["insufficientAngularCoverage"])
        var calibration = RingCalibration.Policy.default
        calibration.minimumOccupiedAngularBins = 20
        calibration.evidence.minimumOccupiedBins = 20
        let calibrationOnly = try RingCalibration.fitRing(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55,
                                                          policy: calibration, ellipsePolicy: .default)
        r.equal("loosening the calibration policy alone changes nothing", calibrationOnly.refusals.map(\.name), ["insufficientAngularCoverage"])
        let both = try RingCalibration.fitRing(pattern: gapped, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55,
                                               policy: calibration, ellipsePolicy: loose)
        r.check("loosening both lets the ring through the coverage rule", !both.refusals.map(\.name).contains("insufficientAngularCoverage"))
        var badCalibration = false
        do { _ = try RingCalibration.fitRing(pattern: whole, centerQX: 63, centerQY: 64, innerRadius: 25, outerRadius: 55,
                                             policy: RingCalibration.Policy(maximumResidual: .nan)) }
        catch let error as RingCalibration.Refusal { if case .invalidInput = error { badCalibration = true } } catch {}
        r.equal("an invalid calibration policy throws invalidInput", badCalibration, true)
    }

    report.section("0.2.0: RingIndexing.minimumMargin") { r in
        r.near("default margin is ln 19", RingIndexing.minimumMargin, log(19), 1e-15)
        r.near("reach at k = 3", RingIndexing.reach(coverage: 3), (9 + 2 * log(19.0)).squareRoot(), 1e-12)
        r.equal("reach with the default margin passed", RingIndexing.reach(coverage: 3, minimumMargin: RingIndexing.minimumMargin), RingIndexing.reach(coverage: 3))
        r.check("a larger margin reaches further", RingIndexing.reach(coverage: 3, minimumMargin: 10) > RingIndexing.reach(coverage: 3))
        r.check("reach is NaN where index refuses (negative margin)", RingIndexing.reach(coverage: 3, minimumMargin: -1).isNaN)
        let families = [KinematicFamily](gold_families())
        let sigma = 0.001
        let g = families[0].family.g
        let accepted = try? RingIndexing.index(q: g, sigma: sigma, in: families, kMax: 1.5, minimumMargin: 0.5)
        r.check("a ring on gold's first family is indexed under a margin of 0.5", accepted != nil)
        var refused = false
        do throws(RingIndexing.Refusal) { _ = try RingIndexing.index(q: g, sigma: sigma, in: families, kMax: 1.5, minimumMargin: -1) }
        catch { if case .invalidInput = error { refused = true } }
        r.equal("a negative margin is invalidInput", refused, true)
    }

    report.section("0.2.0: kinematicReflections(retainingAll:)") { r in
        let all = try gold.kinematicReflections(kMax: 1.5)
        let none = try gold.kinematicReflections(kMax: 1.5, retainingAll: false)
        r.check("the default holds every reflection before the cut", all.all.count > all.kept.count)
        r.equal("retainingAll: false holds none of them", none.all.count, 0)
        r.equal("kept is identical", none.kept, all.kept)
        r.equal("families are identical", none.families, all.families)
        let explicit = try gold.kinematicReflections(kMax: 1.5, retainingAll: true)
        r.equal("retainingAll: true is the default", explicit, all)
        var refused: KinematicError?
        do { _ = try gold.kinematicReflections(kMax: 1.5, cut: -1, retainingAll: false) } catch let e as KinematicError { refused = e } catch {}
        r.equal("the refusals are the same", refused, .invalidCut(-1))
    }
}

/// Gold's kinematic families (the first rings), for the indexing check.
private func gold_families() -> [KinematicFamily] {
    let gold = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: goldSites)
    return (try? gold.kinematicReflections(kMax: 1.5, retainingAll: false).families) ?? []
}
