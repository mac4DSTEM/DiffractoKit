//
//  PolicyTests.swift
//  DiffractoKit 0.2.0 (docs/registration-0.2.0.md §1, §2): the policies of DiffractoGeometry are parameters with the suite's
//  values as defaults. Tier 2 (the default path is held bit for bit by the existing suite, EllipseAppsParityTests included):
//  here a non-default value moves the outcome in the stated direction, and an invalid one is refused as documented.
//  Every test names under "RED UNDER" the production mutation that turns it red (tools/kit-0.2.0/mutants-policy.py).
//

import Testing
import Foundation
import DiffractoGeometry
import DiffractoCore
import DiffractoScattering

// MARK: - EllipseCalibration.Policy

struct EllipsePolicyTests {
    typealias T = RingFitAmendmentTests
    typealias Policy = EllipseCalibration.Policy

    static func clean(epsilon: Double = 0.02, thetaDeg: Double = 35) -> DiffractionPattern {
        T.render { dr, dc in T.ring(dr, dc, radius: 42, epsilon: epsilon, thetaDeg: thetaDeg) }
    }

    static func conic(_ pattern: DiffractionPattern, policy: Policy = .default, sparse: Bool = false,
                      residual: Double? = nil) throws -> EllipseCalibrationFit {
        try EllipseCalibration.fit1D(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28, outerRadius: 56,
                                     maximumResidual: residual, acceptSparseCoverage: sparse, policy: policy)
    }

    static func gapped(_ gap: Double) -> EllipsePlant {
        EllipsePlant(size: 128, radius: 40, eps: 0.02, thetaDegrees: 30, gap: gap, noise: 0, seed: 201)
    }

    static func conic(_ plant: EllipsePlant, policy: Policy, sparse: Bool = false) throws -> EllipseCalibrationFit {
        try EllipseCalibration.fit1D(pattern: plant.render(), centerQX: plant.centerQX + 1.3, centerQY: plant.centerQY - 0.9,
                                     innerRadius: plant.b - 20, outerRadius: plant.a + 20, maximumResidual: 0.4,
                                     acceptSparseCoverage: sparse, policy: policy)
    }

    /// The documented defaults are the constants the fit used before 0.2.0, and the derived bounds follow the bin count.
    /// RED UNDER: any default of `Policy.init` changed (tools/kit-0.2.0/mutants-policy.py, ellipse-*).
    @Test func theDefaultsAreTheDocumentedValues() {
        let p: Policy = .default
        #expect(p.angularBinCount == EllipseCalibration.angularBinCount)
        #expect(p.degeneracyBoundBins == EllipseCalibration.degeneracyBoundBins)
        #expect(p.sparseFloorBins == EllipseCalibration.sparseFloorBins)
        #expect(p.strongSignalFraction == 0.2)
        #expect(p.oneRingRatio == 1.10)
        #expect(p.minimumConicSamples == 20)
        #expect(p.minimumProfileSamples == 32)
        #expect(p.minimumRelativeSignal == 1e-6)
        #expect(p.impulseRangeFactor == 2)
        #expect(p.maximumConicResidual == 0.2)
        #expect(p.maximumProfileResidual == 0.22)
        #expect(p.profileSeedResidual == 0.4)
        #expect(p.maximumAxisRatio == 5)
        #expect(p.robustCutMultiplier == 4 * 1.4826)
        #expect(p.minimumInlierPercent == 70)
        #expect(p.ringWidthFloor == 0.25)
        #expect(p.ringWidthFloorTolerance == 0.05)
        #expect(p.loweredRingWidthFloor == 0.01)
        #expect(p.evidence == RingEvidence.Policy.default)
        #expect(p.isValid)
        let custom = Policy(angularBinCount: 72)
        #expect(custom.degeneracyBoundBins == 60)
        #expect(custom.sparseFloorBins == 24)
        #expect(custom.isValid)
    }

    /// A fit given `.default` equals one given no policy, on every entry point.
    @Test func aFitGivenTheDefaultEqualsAFitGivenNone() throws {
        let pattern = Self.clean()
        let a = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28,
                                                        outerRadius: 56, requireRingEvidence: true)
        let b = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28,
                                                        outerRadius: 56, requireRingEvidence: true, policy: .default)
        #expect(a == b)
        let c = try Self.conic(pattern)
        let d = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28, outerRadius: 56)
        #expect(c == d)
    }

    static func invalidPolicies() -> [Policy] {
        var list: [Policy] = []
        list.append(Policy(angularBinCount: 5))
        list.append(Policy(angularBinCount: 36, degeneracyBoundBins: 40))
        list.append(Policy(angularBinCount: 36, degeneracyBoundBins: 10, sparseFloorBins: 12))
        list.append(Policy(angularBinCount: 36, sparseFloorBins: 0))
        list.append(Policy(strongSignalFraction: .nan))
        list.append(Policy(strongSignalFraction: 1))
        list.append(Policy(oneRingRatio: 0.9))
        list.append(Policy(minimumConicSamples: 2))
        list.append(Policy(minimumProfileSamples: 5))
        list.append(Policy(minimumRelativeSignal: -1))
        list.append(Policy(impulseRangeFactor: -2))
        list.append(Policy(maximumConicResidual: 0))
        list.append(Policy(maximumProfileResidual: .nan))
        list.append(Policy(profileSeedResidual: -0.4))
        list.append(Policy(maximumAxisRatio: 0.5))
        list.append(Policy(robustCutMultiplier: 0))
        list.append(Policy(minimumInlierPercent: 101))
        list.append(Policy(ringWidthFloor: 0))
        list.append(Policy(ringWidthFloorTolerance: -0.05))
        list.append(Policy(loweredRingWidthFloor: 0.5))
        list.append(Policy(evidence: RingEvidence.Policy(minimumSectorPixels: 2)))
        return list
    }

    /// An invalid policy is refused, never clamped: `isValid` is false, every fit throws `invalidInput`, and the ring
    /// evidence is empty.
    /// RED UNDER: the validity guard removed from `fit1D` (the first case then fits).
    @Test func anInvalidPolicyIsRefused() throws {
        let pattern = Self.clean()
        let valid = try Self.conic(pattern)
        for policy in Self.invalidPolicies() {
            #expect(!policy.isValid)
            var refusedConic = false
            do { _ = try Self.conic(pattern, policy: policy) } catch let caught as EllipseCalibration.FitError { guard case .invalidInput = caught else { throw caught }; refusedConic = true }
            #expect(refusedConic, "\(policy)")
            var refusedBest = false
            do {
                _ = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28,
                                                            outerRadius: 56, policy: policy)
            } catch let caught as EllipseCalibration.FitError { guard case .invalidInput = caught else { throw caught }; refusedBest = true }
            #expect(refusedBest)
            var refusedProfile = false
            do {
                _ = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: valid, centerQX: 64.1, centerQY: 61.2,
                                                            innerRadius: 28, outerRadius: 56, policy: policy)
            } catch let caught as EllipseCalibration.FitError { guard case .invalidInput = caught else { throw caught }; refusedProfile = true }
            #expect(refusedProfile)
        }
        let bad = RingEvidence.Policy(minimumSectorPixels: 2)
        let evidence = EllipseCalibration.ringEvidence(pattern: pattern, fit: valid, centerQX: 64.1, centerQY: 61.2,
                                                       innerRadius: 28, outerRadius: 56, policy: bad)
        #expect(evidence.occupiedBins == 0)
        #expect(!evidence.passes(bad))
    }

    /// The coverage bound: a ring with a 120 degree gap (25 of 36 bins, measured) is refused at the default 30 and fitted,
    /// unmarked, at a bound of 20; a ring with a 40 degree gap is fitted at 30 and refused at a bound of 36.
    /// RED UNDER: the degeneracy bound ignored (`policy.degeneracyBoundBins` read as 30); the default bin count changed.
    @Test func theCoverageBoundMovesWhatIsAccepted() throws {
        let wide = Self.gapped(120)
        var refused = -1
        do { _ = try Self.conic(wide, policy: .default) } catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; refused = n }
        #expect(refused == 25)
        let loose = Policy(degeneracyBoundBins: 20)
        let fit = try Self.conic(wide, policy: loose)
        #expect(fit.occupiedAngularBins == 25)
        #expect(!fit.sparseCoverage)

        let narrow = Self.gapped(40)
        let accepted = try Self.conic(narrow, policy: .default)
        #expect(accepted.occupiedAngularBins >= 30)
        let strict = Policy(degeneracyBoundBins: 36)
        var strictRefused = -1
        do { _ = try Self.conic(narrow, policy: strict) } catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; strictRefused = n }
        #expect(strictRefused == accepted.occupiedAngularBins)
    }

    /// The anyway path's floor: 25 bins pass the default floor of 12 (marked sparse) and are refused at a floor of 26; 9 bins
    /// are refused at the default and fitted (marked) at a floor of 6.
    /// RED UNDER: the sparse floor ignored.
    @Test func theSparseFloorMovesWhatTheAnywayPathAccepts() throws {
        let plant = Self.gapped(120)
        let marked = try Self.conic(plant, policy: .default, sparse: true)
        #expect(marked.sparseCoverage)
        let raised = Policy(angularBinCount: 36, degeneracyBoundBins: 30, sparseFloorBins: 26)
        var raisedCount = -1
        do { _ = try Self.conic(plant, policy: raised, sparse: true) }
        catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; raisedCount = n }
        #expect(raisedCount == 25)

        let sparse = Self.gapped(270)
        var defaultCount = -1
        do { _ = try Self.conic(sparse, policy: .default, sparse: true) }
        catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; defaultCount = n }
        #expect(defaultCount >= 0 && defaultCount < 12)
        let lowered = Policy(angularBinCount: 36, degeneracyBoundBins: 30, sparseFloorBins: 6)
        var loweredAnswer = "fitted"
        do { _ = try Self.conic(sparse, policy: lowered, sparse: true) } catch { loweredAnswer = "\(error)" }
        #expect(!loweredAnswer.contains("insufficientAngularCoverage"), "\(loweredAnswer)")
    }

    /// The strong-signal fraction decides which pixels occupy a bin: a ring whose 120 degrees hold 0.35 of its amplitude
    /// occupies all bins at 0.2 of the range and 24 or fewer at 0.5.
    /// RED UNDER: the strong-signal fraction ignored.
    @Test func theStrongSignalFractionDecidesWhichBinsAreOccupied() throws {
        let pattern = T.render { dr, dc in
            T.ring(dr, dc, radius: 42, epsilon: 0.01, thetaDeg: 15) * (T.inGap(dr, dc, start: 200, width: 120) ? 0.35 : 1)
        }
        let low = try Self.conic(pattern, residual: 0.4)
        #expect(low.occupiedAngularBins >= 30)
        let high = Policy(strongSignalFraction: 0.5)
        var count = -1
        do { _ = try Self.conic(pattern, policy: high, residual: 0.4) }
        catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; count = n }
        #expect(count >= 0 && count < 30, "\(count)")
    }

    /// The one-ring ratio of the anyway path: three grains are more than one ring at 1.10, and pass at 10; two arcs of one
    /// ring pass at 1.10 and are called more than one ring at the least ratio, 1.
    /// RED UNDER: the one-ring ratio ignored.
    @Test func theOneRingRatioDecidesWhatTheAnywayPathCallsOneRing() throws {
        func arcs(_ sectors: [(from: Double, to: Double, radius: Double)]) -> DiffractionPattern {
            var pixels = [Float](repeating: 0, count: 128 * 128)
            for qx in 0..<128 {
                for qy in 0..<128 {
                    let dx = Double(qx) - 64, dy = Double(qy) - 64
                    var angle = atan2(dy, dx) * 180 / .pi
                    if angle < 0 { angle += 360 }
                    let r = hypot(dx, dy)
                    var value = 0.0
                    for s in sectors where angle >= s.from && angle < s.to { value = exp(-(r - s.radius) * (r - s.radius) / 8) }
                    pixels[qx * 128 + qy] = Float(value)
                }
            }
            return DiffractionPattern(width: 128, height: 128, pixels: pixels)
        }
        let grains = arcs([(0, 60, 24), (120, 180, 33), (240, 300, 42)])
        var strictAnswer = "fitted"
        do { _ = try EllipseCalibration.fit1D(pattern: grains, centerQX: 64.4, centerQY: 63.6, innerRadius: 12, outerRadius: 54,
                                              maximumResidual: 0.4, acceptSparseCoverage: true) }
        catch let caught as EllipseCalibration.FitError { guard case .moreThanOneRing = caught else { throw caught }; strictAnswer = "moreThanOneRing" } catch { strictAnswer = "\(error)" }
        #expect(strictAnswer == "moreThanOneRing")
        var looseAnswer = "fitted"
        do { _ = try EllipseCalibration.fit1D(pattern: grains, centerQX: 64.4, centerQY: 63.6, innerRadius: 12, outerRadius: 54,
                                              maximumResidual: 0.4, acceptSparseCoverage: true, policy: Policy(oneRingRatio: 10)) }
        catch let caught as EllipseCalibration.FitError { guard case .moreThanOneRing = caught else { throw caught }; looseAnswer = "moreThanOneRing" } catch { looseAnswer = "\(error)" }
        #expect(looseAnswer != "moreThanOneRing", "\(looseAnswer)")

        let oneRing = arcs([(0, 90, 30), (180, 270, 30)])
        let passed = try EllipseCalibration.fit1D(pattern: oneRing, centerQX: 64.4, centerQY: 63.6, innerRadius: 20, outerRadius: 48,
                                                  maximumResidual: 0.4, acceptSparseCoverage: true)
        #expect(passed.sparseCoverage)
        var tightAnswer = "fitted"
        do { _ = try EllipseCalibration.fit1D(pattern: oneRing, centerQX: 64.4, centerQY: 63.6, innerRadius: 20, outerRadius: 48,
                                              maximumResidual: 0.4, acceptSparseCoverage: true, policy: Policy(oneRingRatio: 1)) }
        catch let caught as EllipseCalibration.FitError { guard case .moreThanOneRing = caught else { throw caught }; tightAnswer = "moreThanOneRing" } catch { tightAnswer = "\(error)" }
        #expect(tightAnswer == "moreThanOneRing")
    }

    /// The residual caps: a caller's explicit cap wins over the policy's; without one the policy's decides. The conic's cap
    /// is read by `fit1D`, the profile's by `fitBestAvailable` (two rings 8 px apart leave 0.2031, measured), and the
    /// conic that seeds the profile is held to the seed residual (0.27 there).
    /// RED UNDER: the policy's residual caps ignored; the seed residual ignored.
    @Test func theResidualCapsOfThePolicyAreTheDefaults() throws {
        let pattern = Self.clean()
        var tinyRefused = false
        do { _ = try Self.conic(pattern, policy: Policy(maximumConicResidual: 1e-9)) }
        catch let caught as EllipseCalibration.FitError { guard case .excessiveResidual = caught else { throw caught }; tinyRefused = true }
        #expect(tinyRefused)
        let explicit = try Self.conic(pattern, policy: Policy(maximumConicResidual: 1e-9), residual: 0.2)
        let baseline = try Self.conic(pattern)
        #expect(explicit == baseline)

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
        func best(_ separation: Double, _ policy: Policy) throws -> EllipseCalibrationFit {
            try EllipseCalibration.fitBestAvailable(pattern: pair(separation), centerQX: 40.5, centerQY: 39.5, innerRadius: 10,
                                                    outerRadius: 40, policy: policy)
        }
        let near = try best(8, .default)
        #expect(near.model == .asymmetricGaussian)
        var capped = false
        do { _ = try best(8, Policy(maximumProfileResidual: 0.2)) } catch let caught as EllipseCalibration.FitError { guard case .excessiveResidual = caught else { throw caught }; capped = true }
        #expect(capped)
        let farDefault = try? best(9, .default)
        #expect(farDefault == nil)
        let farLoose = try best(9, Policy(maximumProfileResidual: 0.30))
        #expect(farLoose.model == .asymmetricGaussian)
        var seedRefused = false
        do { _ = try best(8, Policy(profileSeedResidual: 0.25)) } catch let caught as EllipseCalibration.FitError { guard case .excessiveResidual = caught else { throw caught }; seedRefused = true }
        #expect(seedRefused)
    }

    /// The axis ratio: an ellipse with a/b = 1.3 is fitted at the default ratio 5 and refused as `invalidEllipse` at 1.2.
    /// RED UNDER: the axis ratio ignored.
    @Test func theAxisRatioBoundsTheEllipse() throws {
        let pattern = T.render { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.3, thetaDeg: 20) }
        let fit = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 22, outerRadius: 64)
        #expect(fit.a / fit.b > 1.25)
        var refused = false
        do {
            _ = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 22, outerRadius: 64,
                                             policy: Policy(maximumAxisRatio: 1.2))
        } catch let caught as EllipseCalibration.FitError { guard case .invalidEllipse = caught else { throw caught }; refused = true }
        #expect(refused)
    }

    /// The pre-fit thresholds: too few pixels, no relative signal, and the impulse exclusion.
    /// RED UNDER: the sample minimums, the relative signal or the impulse factor ignored.
    @Test func theSampleSignalAndImpulseThresholdsAreParameters() throws {
        let pattern = Self.clean()
        var fewConic = false
        do { _ = try Self.conic(pattern, policy: Policy(minimumConicSamples: 100_000)) }
        catch let caught as EllipseCalibration.FitError { guard case .invalidInput = caught else { throw caught }; fewConic = true }
        #expect(fewConic)
        let seed = try Self.conic(pattern)
        var fewProfile = false
        do {
            _ = try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: seed, centerQX: 64.1, centerQY: 61.2,
                                                        innerRadius: 28, outerRadius: 56, policy: Policy(minimumProfileSamples: 100_000))
        } catch let caught as EllipseCalibration.FitError { guard case .insufficientSignal = caught else { throw caught }; fewProfile = true }
        #expect(fewProfile)

        // a ring of amplitude 1 on a background of 1e5: range / maximum is 1e-5
        let lifted = T.render { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.02, thetaDeg: 35) + 1e5 }
        _ = try Self.conic(lifted, residual: 0.4)
        var noSignal = false
        do { _ = try Self.conic(lifted, policy: Policy(minimumRelativeSignal: 1e-2), residual: 0.4) }
        catch let caught as EllipseCalibration.FitError { guard case .insufficientSignal = caught else { throw caught }; noSignal = true }
        #expect(noSignal)

        // one hot pixel on the ring: ignored by default, a coverage of one bin when the exclusion is off
        var pixels = pattern.pixels
        let row = Int((63.37 + 42 * cos(35 * Double.pi / 180)).rounded())
        let col = Int((61.59 + 42 * sin(35 * Double.pi / 180)).rounded())
        pixels[row * pattern.width + col] = 1000
        let spiked = DiffractionPattern(width: pattern.width, height: pattern.height, pixels: pixels)
        let ignored = try Self.conic(spiked)
        #expect(ignored.occupiedAngularBins >= 30)
        var coverage = -1
        do { _ = try Self.conic(spiked, policy: Policy(impulseRangeFactor: 0)) }
        catch let caught as EllipseCalibration.FitError { guard case .insufficientAngularCoverage(let n) = caught else { throw caught }; coverage = n }
        #expect(coverage >= 0 && coverage < 5, "\(coverage)")
    }

    /// The ring-width floor: a 0.27 px ring is refused as `ringWidthAtFloor` by default and accepted when the tolerance
    /// around the floor is 0; a 0.35 px ring is accepted by default and refused when the tolerance is 0.15.
    /// RED UNDER: the floor, its tolerance or the lowered floor ignored.
    @Test func theRingWidthFloorAndItsToleranceAreParameters() throws {
        func thin(_ sigma: Double) -> DiffractionPattern {
            let scene = RingFitAmendmentTests.Scene(rows: 64, cols: 64, row0: 31.4, col0: 32.3)
            return T.render(scene) { dr, dc in exp(-(hypot(dr, dc) - 20) * (hypot(dr, dc) - 20) / (2 * sigma * sigma)) }
        }
        func profile(_ pattern: DiffractionPattern, _ policy: Policy) throws -> EllipseCalibrationFit {
            let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 32.5, centerQY: 31.5, innerRadius: 10,
                                                     outerRadius: 30, maximumResidual: 0.4)
            return try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 32.5, centerQY: 31.5,
                                                           innerRadius: 10, outerRadius: 30, policy: policy)
        }
        var narrowRefused = false
        do { _ = try profile(thin(0.27), .default) } catch let caught as EllipseCalibration.FitError { guard case .ringWidthAtFloor = caught else { throw caught }; narrowRefused = true }
        #expect(narrowRefused)
        let narrow = try profile(thin(0.27), Policy(ringWidthFloorTolerance: 0))
        #expect(narrow.model == .asymmetricGaussian)

        let wide = try profile(thin(0.35), .default)
        #expect(wide.model == .asymmetricGaussian)
        var wideRefused = false
        do { _ = try profile(thin(0.35), Policy(ringWidthFloorTolerance: 0.15)) }
        catch let caught as EllipseCalibration.FitError { guard case .ringWidthAtFloor = caught else { throw caught }; wideRefused = true }
        #expect(wideRefused)
        var floorRaised = false
        do { _ = try profile(thin(0.35), Policy(ringWidthFloor: 0.5)) } catch let caught as EllipseCalibration.FitError { guard case .ringWidthAtFloor = caught else { throw caught }; floorRaised = true }
        #expect(floorRaised)
    }

    /// The lowered floor of the one-sided refit: a ring with one edge of 0.05 px is read within 0.1 px of its radius at the
    /// lowered floor of 0.01, and not when the lowered floor is the ordinary one (0.25).
    /// RED UNDER: the lowered floor ignored.
    @Test func theLoweredFloorKeepsTheRadiusOfAOneSidedRing() throws {
        let pattern = T.render { dr, dc in T.ring(dr, dc, radius: 40, epsilon: 0.01, thetaDeg: 10, sigmaIn: 3, sigmaOut: 0.05) }
        func radius(_ policy: RingCalibration.Policy, _ ellipse: EllipseCalibration.Policy) throws -> Double {
            let ring = try RingCalibration.fitRing(pattern: pattern, centerQX: 64.47, centerQY: 60.89, innerRadius: 26,
                                                   outerRadius: 52, policy: policy, ellipsePolicy: ellipse)
            return try #require(ring.radius)
        }
        let atDefault = try radius(.default, .default)
        #expect(abs(atDefault - 40) < 0.1, "\(atDefault)")
        let atFloor = try radius(.default, Policy(loweredRingWidthFloor: 0.25))
        #expect(abs(atFloor - 40) > 0.1, "\(atFloor)")
    }

    /// The robust refit: with outliers (bright 3 x 3 blocks off the ring) the refit changes the fit; with the cut very wide
    /// or the inlier share at 100 % no refit runs, and the two agree.
    /// RED UNDER: the robust cut or the inlier percentage ignored.
    @Test func theRobustRefitIsParameterised() throws {
        let scene = RingFitAmendmentTests.Scene()
        let base = T.render(scene, noise: 0.05, seed: 5) { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.03, thetaDeg: 40) }
        var pixels = base.pixels
        let blocks: [(Int, Int)] = [(40, 63), (92, 70), (66, 100), (62, 26), (45, 90), (85, 40)]
        for (r0, c0) in blocks {
            for dr in 0..<3 {
                for dc in 0..<3 { pixels[(r0 + dr) * base.width + (c0 + dc)] = 0.9 }
            }
        }
        let pattern = DiffractionPattern(width: base.width, height: base.height, pixels: pixels)
        func fit(_ policy: Policy) throws -> EllipseCalibrationFit {
            let conic = try EllipseCalibration.fit1D(pattern: pattern, centerQX: 64.1, centerQY: 61.2, innerRadius: 28,
                                                     outerRadius: 56, maximumResidual: 0.4)
            return try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: conic, centerQX: 64.1, centerQY: 61.2,
                                                           innerRadius: 28, outerRadius: 56, maximumResidual: 0.5, policy: policy)
        }
        let refit = try fit(.default)
        let wide = try fit(Policy(robustCutMultiplier: 1e6))
        let all = try fit(Policy(minimumInlierPercent: 100))
        #expect(refit != wide)
        #expect(wide == all)
    }
}

// MARK: - RingCalibration.Policy and RingEvidence.Policy

struct RingCalibrationPolicyTests {
    typealias P = RingCalibration.Policy

    static func evidence(occupied: Int = 36, innerDrop: Double = 0.9, outerDrop: Double = 0.9, centreShift: Double = 0.01,
                         centreError: Double = 0.01, secondary: Double? = nil) -> RingEvidence {
        RingEvidence(amplitude: 1, sectorAmplitudes: Array(repeating: 1, count: 36),
                     sectorAmplitudeErrors: Array(repeating: 0.01, count: 36), sectorShifts: Array(repeating: 0, count: 36),
                     occupiedBins: occupied, innerDrop: innerDrop, outerDrop: outerDrop, radiusShift: 0,
                     centreShift: centreShift, centreShiftError: centreError, ellipseShift: 0.01, ellipseShiftError: 0.01,
                     secondaryPeak: secondary)
    }

    static func report(a: Double = 40.2, b: Double = 39.8, bins: Int = 36, residual: Double = 0.05,
                       evidence: RingEvidence? = nil) -> RingCalibration.FitReport {
        RingCalibration.FitReport(model: .asymmetricGaussian, semiAxisA: a, semiAxisB: b, occupiedAngularBins: bins,
                                  sparseCoverage: false, normalizedResidual: residual, innerSigma: 3, outerSigma: 3,
                                  annulusWidth: 20, annulusInnerRadius: 30, annulusOuterRadius: 50, ringEvidence: evidence)
    }

    static func names(_ report: RingCalibration.FitReport, _ policy: P = .default) -> [String] {
        RingCalibration.refusals(report, highTensionVolts: 200_000, policy: policy).map(\.name)
    }

    /// The documented defaults equal the registered constants.
    /// RED UNDER: any default of `RingCalibration.Policy` or `RingEvidence.Policy` changed.
    @Test func theDefaultsAreTheRegisteredBounds() {
        let p: P = .default
        #expect(p.angularBinCount == 36)
        #expect(p.minimumOccupiedAngularBins == 30)
        #expect(p.maximumResidual == 0.22)
        #expect(p.maximumEllipticity == 0.093)
        let e: RingEvidence.Policy = p.evidence
        #expect(e.minimumSectorPixels == 8)
        #expect(e.minimumAmplitudeFraction == 0.5)
        #expect(e.minimumSignificance == 3.0)
        #expect(e.minimumOccupiedBins == 30)
        #expect(e.minimumDrop == 0.25)
        #expect(e.maximumShift == 0.25)
        #expect(e.shiftSignificance == 4.0)
        #expect(e.maximumSecondaryPeak == 0.25)
        #expect(p.isValid && e.isValid)
        #expect(Self.names(Self.report(evidence: Self.evidence())).isEmpty)
    }

    /// Invalid policies give `[.invalidInput]` alone and are not valid.
    /// RED UNDER: the validity guard of `fitRefusals` removed.
    @Test func anInvalidPolicyGivesInvalidInputAlone() {
        var bad: [P] = []
        bad.append(P(angularBinCount: 3))
        bad.append(P(minimumOccupiedAngularBins: 40))
        bad.append(P(maximumResidual: .nan))
        bad.append(P(maximumResidual: -0.1))
        bad.append(P(maximumEllipticity: -1))
        bad.append(P(evidence: RingEvidence.Policy(minimumDrop: -1)))
        bad.append(P(evidence: RingEvidence.Policy(minimumOccupiedBins: 37)))
        bad.append(P(evidence: RingEvidence.Policy(maximumShift: .nan)))
        for policy in bad {
            #expect(!policy.isValid)
            #expect(Self.names(Self.report(), policy) == ["invalidInput"])
        }
    }

    /// The residual and ellipticity caps, both directions.
    /// RED UNDER: the policy's residual cap or ellipticity ignored.
    @Test func theResidualAndEllipticityCapsMoveTheRefusals() {
        let high = Self.report(residual: 0.25)
        #expect(Self.names(high) == ["excessiveResidual"])
        #expect(Self.names(high, P(maximumResidual: 0.30)).isEmpty)
        #expect(Self.names(Self.report(residual: 0.15), P(maximumResidual: 0.10)) == ["excessiveResidual"])

        let elliptic = Self.report(a: 44, b: 40)
        #expect(Self.names(elliptic) == ["excessiveEllipticity"])
        #expect(Self.names(elliptic, P(maximumEllipticity: 0.15)).isEmpty)
        #expect(Self.names(Self.report(a: 41, b: 40), P(maximumEllipticity: 0.02)) == ["excessiveEllipticity"])
    }

    /// The coverage bound and the bin count.
    /// RED UNDER: the policy's occupied-bin bound or bin count ignored.
    @Test func theCoverageBoundsMoveTheRefusals() {
        let gapped = Self.report(bins: 28)
        #expect(Self.names(gapped) == ["insufficientAngularCoverage"])
        #expect(Self.names(gapped, P(minimumOccupiedAngularBins: 24)).isEmpty)
        #expect(Self.names(Self.report(bins: 32), P(minimumOccupiedAngularBins: 34)) == ["insufficientAngularCoverage"])
        // 40 occupied bins are out of range for 36 bins, in range for 72
        #expect(Self.names(Self.report(bins: 40)) == ["invalidInput"])
        #expect(Self.names(Self.report(bins: 40), P(angularBinCount: 72, minimumOccupiedAngularBins: 60)) == ["insufficientAngularCoverage"])
        #expect(Self.names(Self.report(bins: 66), P(angularBinCount: 72, minimumOccupiedAngularBins: 60)).isEmpty)
    }

    /// The evidence bounds, each both ways: sectors holding the ring, the peak drops, the shifts and their significance, the
    /// second peak. The instance judgements take the policy; the properties are the default's.
    /// RED UNDER: any of the evidence bounds ignored by the refusals or the judgements.
    @Test func theEvidenceBoundsMoveTheJudgements() {
        let sectors = Self.evidence(occupied: 28)
        #expect(Self.names(Self.report(evidence: sectors)) == ["insufficientAngularCoverage"])
        var loose = RingEvidence.Policy()
        loose.minimumOccupiedBins = 24
        #expect(Self.names(Self.report(evidence: sectors), P(evidence: loose)).isEmpty)
        #expect(!sectors.passes && sectors.passes(loose))

        let flat = Self.evidence(innerDrop: 0.2)
        #expect(Self.names(Self.report(evidence: flat)) == ["notARingPeak"])
        var shallow = RingEvidence.Policy()
        shallow.minimumDrop = 0.1
        #expect(Self.names(Self.report(evidence: flat), P(evidence: shallow)).isEmpty)
        #expect(!flat.isPeak && flat.isPeak(shallow))
        var steep = RingEvidence.Policy()
        steep.minimumDrop = 0.95
        #expect(Self.names(Self.report(evidence: Self.evidence()), P(evidence: steep)) == ["notARingPeak"])

        let shifted = Self.evidence(centreShift: 0.6)
        #expect(Self.names(Self.report(evidence: shifted)) == ["ringGeometryInconsistent"])
        var generous = RingEvidence.Policy()
        generous.maximumShift = 1.0
        #expect(Self.names(Self.report(evidence: shifted), P(evidence: generous)).isEmpty)
        #expect(!shifted.geometryHolds && shifted.geometryHolds(generous))

        let uncertain = Self.evidence(centreShift: 0.6, centreError: 0.2)
        #expect(uncertain.geometryHolds)
        var sharp = RingEvidence.Policy()
        sharp.shiftSignificance = 2
        #expect(!uncertain.geometryHolds(sharp))
        #expect(Self.names(Self.report(evidence: uncertain), P(evidence: sharp)) == ["ringGeometryInconsistent"])

        let second = Self.evidence(secondary: 0.3)
        #expect(Self.names(Self.report(evidence: second)) == ["secondRingInAnnulus"])
        var tolerant = RingEvidence.Policy()
        tolerant.maximumSecondaryPeak = 0.5
        #expect(Self.names(Self.report(evidence: second), P(evidence: tolerant)).isEmpty)
        #expect(!second.isSingleRing && second.isSingleRing(tolerant))
    }

    /// The reason names the bounds of the policy it is asked under.
    @Test func theReasonNamesThePoliciesBounds() {
        let refusal = RingCalibration.Refusal.excessiveResidual(0.25)
        #expect(refusal.reason.contains("0.22"))
        #expect(refusal.reason(policy: P(maximumResidual: 0.2)).contains("0.20"))
        let coverage = RingCalibration.Refusal.insufficientAngularCoverage(occupiedBins: 20)
        #expect(coverage.reason(policy: P(minimumOccupiedAngularBins: 25)).contains("at least 25"))
    }

    /// The measurement bounds act on the data: a ring whose 120 degrees hold 0.35 of its amplitude has fewer than 30 sectors
    /// at the default amplitude fraction and 30 or more at 0.05; sectors need pixels: none are read with 100000 demanded.
    /// RED UNDER: the amplitude fraction, the significance or the sector pixels ignored by `ringEvidence`.
    @Test func theMeasurementBoundsActOnTheSectors() throws {
        typealias T = RingFitAmendmentTests
        let pattern = T.render { dr, dc in
            T.ring(dr, dc, radius: 42, epsilon: 0.01, thetaDeg: 15) * (T.inGap(dr, dc, start: 200, width: 120) ? 0.35 : 1)
        }
        let fit = try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: 64.47, centerQY: 60.89, innerRadius: 28,
                                                          outerRadius: 56)
        func bins(_ policy: RingEvidence.Policy) -> Int {
            EllipseCalibration.ringEvidence(pattern: pattern, fit: fit, centerQX: 64.47, centerQY: 60.89, innerRadius: 28,
                                            outerRadius: 56, policy: policy).occupiedBins
        }
        let atDefault = bins(.default)
        #expect(atDefault < 30)
        var weak = RingEvidence.Policy()
        weak.minimumAmplitudeFraction = 0.05
        #expect(bins(weak) >= 30)
        var starved = RingEvidence.Policy()
        starved.minimumSectorPixels = 100_000
        #expect(bins(starved) == 0)
        var strict = RingEvidence.Policy()
        strict.minimumSignificance = 1e6
        #expect(bins(strict) == 0)
    }

    /// Through the fit: a residual cap of 1e-9 refuses a ring the default calibrates, and the same policy's evidence bound
    /// reaches the evidence the fit measures.
    /// RED UNDER: `fitRing` not passing the policy to the refusals; the evidence policy not reaching the fit.
    @Test func thePolicyReachesTheFit() throws {
        typealias T = RingFitAmendmentTests
        let pattern = T.render(noise: 0.05, seed: 3) { dr, dc in T.ring(dr, dc, radius: 42, epsilon: 0.02, thetaDeg: 30) }
        func ring(_ policy: P) throws -> RingRadiusFit {
            try RingCalibration.fitRing(pattern: pattern, centerQX: 64.47, centerQY: 60.89, innerRadius: 28, outerRadius: 56,
                                        policy: policy)
        }
        let ordinary = try ring(.default)
        #expect(ordinary.isUsable)
        let again = try ring(P())
        #expect(again == ordinary)
        let strict = try ring(P(maximumResidual: 1e-9))
        #expect(strict.refusals.map(\.name) == ["excessiveResidual"])
        let demanding = try ring(P(evidence: RingEvidence.Policy(minimumSectorPixels: 100_000)))
        #expect(demanding.report?.ringEvidence?.occupiedBins == 0)
        #expect(demanding.refusals.map(\.name).contains("insufficientAngularCoverage"))
        var refused = false
        do { _ = try ring(P(maximumResidual: .nan)) } catch let caught as RingCalibration.Refusal { guard case .invalidInput = caught else { throw caught }; refused = true }
        #expect(refused)
    }
}

// MARK: - RingIndexing's margin

struct RingIndexingPolicyTests {
    typealias Index = RingIndexing

    /// A measurement 2.9 sigma from one family and 3.1 from the next (z within k = 3 for one only): the margin is 0.6, under
    /// ln 19 (refused as ambiguous), accepted at a minimum margin of 0.5, refused again at 5. A negative or NaN margin is
    /// `invalidInput`.
    /// RED UNDER: the minimum margin ignored; its default changed from ln 19.
    @Test func theMinimumMarginDecidesAmbiguity() throws {
        let list = try IndexPresets.dyadic([1.0, 1.5])
        let sigma = 1.0 / 12.0
        let q = 1.0 + 2.9 * sigma
        #expect(Index.minimumMargin == log(19.0))
        var defaultAnswer = "accepted"
        do throws(Index.Refusal) { _ = try Index.index(q: q, sigma: sigma, in: list, kMax: 4) }
        catch { defaultAnswer = error.name }
        #expect(defaultAnswer == "ambiguous")
        let loose = try Index.index(q: q, sigma: sigma, in: list, kMax: 4, minimumMargin: 0.5)
        #expect(loose.best.family.g == 1.0)
        let margin = try #require(loose.margin)
        #expect(abs(margin - 0.6) < 1e-6, "\(margin)")
        var strictAnswer = "accepted"
        do throws(Index.Refusal) { _ = try Index.index(q: q, sigma: sigma, in: list, kMax: 4, minimumMargin: 5) }
        catch { strictAnswer = error.name }
        #expect(strictAnswer == "ambiguous")
        for bad in [-1.0, Double.nan, Double.infinity] {
            var answer = "accepted"
            do throws(Index.Refusal) { _ = try Index.index(q: q, sigma: sigma, in: list, kMax: 4, minimumMargin: bad) }
            catch { answer = error.name }
            #expect(answer == "invalidInput")
        }
        #expect(Index.reach(coverage: 3) == Index.reach(coverage: 3, minimumMargin: Index.minimumMargin))
        #expect(Index.reach(coverage: 3, minimumMargin: 10) > Index.reach(coverage: 3))
    }
}
