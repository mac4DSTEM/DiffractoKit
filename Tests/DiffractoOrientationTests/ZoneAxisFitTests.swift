//
//  ZoneAxisFitTests.swift
//  Role: asking the data which beam direction the specimen is on: the planted-family recovery, the two nulls, the total order
//        of the ranking.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct ZoneAxisFitTests {
    /// Planted truth: patterns built from aluminium's own ⟨110⟩ projection,
    /// at a sign-discriminating in-plane rotation. The fit must recover the
    /// FAMILY, and — because fcc is cubic — every sampled ⟨110⟩ equivalent
    /// must tie EXACTLY. The tie is the real assertion: it is what cubic
    /// symmetry requires, so a sweep that returns a winner without it is
    /// returning something rather than measuring something.
    ///
    /// Mutation: ranking by `meanDistance` alone instead of matched count
    /// first — the trap the source comment names, where one vector at
    /// 0.001 Å⁻¹ beats nine at 0.01.
    @Test func testTheFitRecoversThePlantedFamilyAndTiesAcrossIt() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let planted = SIMD3(1, 1, 0)
        let rotation = 37.2 * Double.pi / 180

        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, reference.kMaxInvAngstrom), crystal: al,
            zoneAxis: planted, settings: reference)
        #expect(base.count > 4)
        let turned = PhaseReferenceLibrary.rotate(base, by: rotation)

        let scale = 0.008
        let originX: Float = 128, originY: Float = 128
        let peaks: [[SIMD2<Double>]] = (0..<64).map { _ in
            [pixel(originX, originY)]
                + turned.map { v in
                    pixel(originX + Float(v.q.x / scale), originY + Float(v.q.y / scale))
                }
        }
        let scan = ScanFixture(scanWidth: 8, scanHeight: 8, peaks: peaks)

        let fits = PhaseVectorMatcher.fitZoneAxis(
            scan: scan, crystal: al, referenceSettings: reference,
            settings: PhaseVectorSettings(), originX: originX, originY: originY,
            invAngstromPerPixel: scale, inPlaneStepDeg: 2)
        let winner = try #require(fits.first)

        func isOneOneZeroFamily(_ a: SIMD3<Int>) -> Bool {
            let m = [abs(a.x), abs(a.y), abs(a.z)].sorted()
            return m == [0, 1, 1]
        }
        #expect(isOneOneZeroFamily(winner.zoneAxis),
            "fitted [\(winner.zoneAxis.x) \(winner.zoneAxis.y) \(winner.zoneAxis.z)], not a ⟨110⟩")
        #expect(abs(winner.explainedFraction - 1.0) <= 0.001)

        // THE TIE. Every ⟨110⟩ in the candidate list must score identically.
        let family = fits.filter { isOneOneZeroFamily($0.zoneAxis) }
        #expect(family.count >= 3,
            "too few ⟨110⟩ axes to test the tie")
        for fit in family {
            #expect(fit.matchedVectors == winner.matchedVectors,
            "[\(fit.zoneAxis.x) \(fit.zoneAxis.y) \(fit.zoneAxis.z)] broke the cubic tie")
        }
        // And a non-equivalent axis must NOT reach the same score, or the
        // sweep is not discriminating at all.
        let other = fits.first { !isOneOneZeroFamily($0.zoneAxis) }
        #expect(try #require(other).matchedVectors < winner.matchedVectors)
    }

    /// The chance floor (2026-09-14). A percentage alone cannot be read: at a
    /// tight tolerance every axis explains a few percent of ANY vectors, and
    /// a display presented a ⟨112⟩ family at 8 % the way it presents a real
    /// fit. `chanceMatchedVectors` is what separates them, and the two halves
    /// of this test are the two cases a reader cannot tell apart on screen.
    ///
    /// Mutations this names: computing the chance expectation at a fixed
    /// radius instead of each pattern's own reach; dropping the
    /// `chanceMatchMultiple` factor from `isAboveChance`; or returning zero
    /// chance, which would mark everything informative.
    @Test func testAFitOnRandomVectorsIsAtChanceAndAPlantedOneIsNot() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let settings = PhaseVectorSettings()
        let scale = 0.008
        let originX: Float = 128, originY: Float = 128

        // Half one: aluminium's own ⟨110⟩ projection, as the recovery test
        // plants it. This must read as information.
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, reference.kMaxInvAngstrom), crystal: al,
            zoneAxis: SIMD3(1, 1, 0), settings: reference)
        let turned = PhaseReferenceLibrary.rotate(base, by: 37.2 * .pi / 180)
        let plantedPeaks: [[SIMD2<Double>]] = (0..<64).map { _ in
            [pixel(originX, originY)]
                + turned.map { v in
                    pixel(originX + Float(v.q.x / scale), originY + Float(v.q.y / scale))
                }
        }
        let plantedFits = PhaseVectorMatcher.fitZoneAxis(
            scan: ScanFixture(scanWidth: 8, scanHeight: 8, peaks: plantedPeaks),
            crystal: al, referenceSettings: reference, settings: settings,
            originX: originX, originY: originY, invAngstromPerPixel: scale,
            inPlaneStepDeg: 2)
        let plantedWinner = try #require(plantedFits.first)
        #expect(plantedWinner.isAboveChance(multiple: settings.chanceMatchMultiple),
            "a perfectly planted ⟨110⟩ was marked as chance: \(plantedWinner.matchedVectors) matched against \(plantedWinner.chanceMatchedVectors) expected")
        #expect(plantedWinner.chanceMatchedVectors > 0,
            "the chance expectation is identically zero, so nothing could ever be marked at chance")

        // Half two: the same number of vectors pointing nowhere in
        // particular, drawn uniformly over the disc the data reaches — which
        // is the model `chanceMatchFraction` states for itself. Whatever axis
        // wins on these must NOT read as information.
        var seed: UInt64 = 0x9E3779B97F4A7C15
        func next() -> Double {                     // deterministic; no shared RNG state
            seed = seed &* 6364136223846793005 &+ 1442695040888963407
            return Double(seed >> 11) / Double(UInt64(1) << 53)
        }
        let reach = turned.map { simd_length($0.q) }.max() ?? 0.8
        func sweep(_ peaks: [[SIMD2<Double>]]) -> PhaseVectorMatcher.ZoneAxisFit? {
            PhaseVectorMatcher.fitZoneAxis(
                scan: ScanFixture(scanWidth: 8, scanHeight: 8, peaks: peaks),
                crystal: al, referenceSettings: reference, settings: settings,
                originX: originX, originY: originY, invAngstromPerPixel: scale,
                inPlaneStepDeg: 2).first
        }
        let uniformPeaks: [[SIMD2<Double>]] = (0..<64).map { _ in
            [pixel(originX, originY)]
                + turned.map { _ in
                    // sqrt keeps the draw uniform per unit AREA, not per radius.
                    let r = reach * next().squareRoot(), a = next() * 2 * .pi
                    return pixel(originX + Float(r * cos(a) / scale), originY + Float(r * sin(a) / scale))
                }
        }
        let uniformWinner = try #require(sweep(uniformPeaks))
        // The model's calibration, asserted rather than printed, because
        // the measurement record quotes it: 14 matched of 1 379 against 11.4
        // expected. A band, not the number, so a change to the sampling is a
        // failure and a change to the last digit is not.
        #expect(abs(Double(uniformWinner.matchedVectors) - uniformWinner.chanceMatchedVectors) <= 6.0,
            "the uniform-disc model no longer predicts what uniform vectors do: \(uniformWinner.matchedVectors) matched of \(uniformWinner.totalVectors), \(uniformWinner.chanceMatchedVectors) expected")
        #expect(!(uniformWinner.isAboveChance(multiple: settings.chanceMatchMultiple)),
            "vectors pointing nowhere produced a fit presented as information: \(uniformWinner.matchedVectors) matched, \(uniformWinner.chanceMatchedVectors) expected by chance")

        // MEASURED LIMIT, printed not asserted: vectors confined to the rings
        // the references occupy beat the disc model, because the model spreads
        // its area over the whole disc while the peaks sit where the
        // references are. `chanceMatchFraction`'s own comment says this
        // understatement is known; these two lines are what it costs here.
        let ringPeaks: [[SIMD2<Double>]] = (0..<64).map { _ in
            [pixel(originX, originY)]
                + turned.map { v in
                    let r = simd_length(v.q), a = next() * 2 * .pi
                    return pixel(originX + Float(r * cos(a) / scale), originY + Float(r * sin(a) / scale))
                }
        }
        let ringWinner = try #require(sweep(ringPeaks))
        // AND THE LIMIT, asserted too: 63 matched against 10.7 expected. The
        // model spreads its area over the whole disc while these vectors sit
        // where the references are, so it understates chance about sixfold —
        // and the guard therefore reads them as information. This assertion
        // exists to fail if that ever stops being true, in either direction.
        #expect(Double(ringWinner.matchedVectors) > (4 * ringWinner.chanceMatchedVectors),
            "the disc model no longer understates ring-concentrated vectors: \(ringWinner.matchedVectors) matched, \(ringWinner.chanceMatchedVectors) expected")
        #expect(ringWinner.isAboveChance(multiple: settings.chanceMatchMultiple),
            "ring-concentrated vectors are no longer read as information; the recorded measured limit needs rewriting")
    }

    /// THE MEASURED CASE, pinned (a diagnostic run, 2026-09-15 night): a wrong axis on a
    /// real crystal explains a large share of the vectors through SHARED
    /// reflections, clears the disc-chance rule by a wide margin, and used to
    /// be presented exactly like the true one. A ⟨110⟩ plant degraded the way
    /// data is (30 % of spots missing, 0.004 Å⁻¹ jitter, three spurious peaks
    /// per pattern): the ⟨112⟩ family explains ~16 % — five times the
    /// chance rule allows, so `isAboveChance` is TRUE, which is the defect
    /// and is asserted as the control — while it sits at ~1.2× the sweep's
    /// median, under the 2× bar. The true family sits at ~6×. Mutations this
    /// names: `informativeSweepRatio` lowered to 1 (every wrong axis passes),
    /// the median computed over the top few instead of the sweep, or
    /// `isInformative` dropping either rule.
    @Test func testAWrongAxisThatSharesReflectionsIsMarkedBelowTheSweep() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let settings = PhaseVectorSettings()
        let scale = 0.008
        let originX: Float = 128, originY: Float = 128
        var seed: UInt64 = 0x2545F4914F6CDD1D
        func next() -> Double {
            seed = seed &* 6364136223846793005 &+ 1442695040888963407
            return Double(seed >> 11) / Double(UInt64(1) << 53)
        }
        func gauss() -> Double {
            sqrt(-2 * log(max(1e-12, next()))) * cos(2 * .pi * next())
        }
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, reference.kMaxInvAngstrom), crystal: al,
            zoneAxis: SIMD3(1, 1, 0), settings: reference)
        let turned = PhaseReferenceLibrary.rotate(base, by: 37.2 * .pi / 180)
        let reach = turned.map { simd_length($0.q) }.max() ?? 0.8
        let peaks: [[SIMD2<Double>]] = (0..<64).map { _ in
            var out = [pixel(originX, originY)]
            for v in turned where next() >= 0.3 {
                out.append(pixel(originX + Float((v.q.x + gauss() * 0.004) / scale), originY + Float((v.q.y + gauss() * 0.004) / scale)))
            }
            for _ in 0..<3 {
                let r = reach * next().squareRoot(), a = next() * 2 * .pi
                out.append(pixel(originX + Float(r * cos(a) / scale), originY + Float(r * sin(a) / scale)))
            }
            return out
        }
        let fits = PhaseVectorMatcher.fitZoneAxis(
            scan: ScanFixture(scanWidth: 8, scanHeight: 8, peaks: peaks),
            crystal: al, referenceSettings: reference, settings: settings,
            originX: originX, originY: originY, invAngstromPerPixel: scale,
            inPlaneStepDeg: 2)
        let winner = try #require(fits.first)
        #expect(Set([abs(winner.zoneAxis.x), abs(winner.zoneAxis.y), abs(winner.zoneAxis.z)]) == Set([0, 1]),
            "the ⟨110⟩ plant was not won by a ⟨110⟩: \(winner.zoneAxis)")
        #expect(winner.isInformative(multiple: settings.chanceMatchMultiple),
            "the true axis was marked: \(winner.explainedFraction) against a sweep median of \(winner.sweepMedianFraction)")
        #expect(winner.sweepRatio > 4,
            "the true family no longer stands out")

        let wrong112 = try #require(fits.first {
            Set([abs($0.zoneAxis.x), abs($0.zoneAxis.y), abs($0.zoneAxis.z)]) == Set([1, 1, 2])
        })
        // The control: the disc rule alone reads it as information.
        #expect(wrong112.isAboveChance(multiple: settings.chanceMatchMultiple),
            "the defect this test pins no longer reproduces: ⟨112⟩ at \(wrong112.explainedFraction) is under 5× disc chance")
        #expect(wrong112.explainedFraction > 0.05,
            "⟨112⟩ explains too little for this to be the measured case")
        // The fix: against the sweep it is a wrong axis.
        #expect(!(wrong112.isAboveSweep),
            "⟨112⟩ at \(wrong112.explainedFraction) is \(wrong112.sweepRatio)× the sweep median of \(wrong112.sweepMedianFraction) and was not marked")
        #expect(!(wrong112.isInformative(multiple: settings.chanceMatchMultiple)))
        // And every axis under 30 % is a wrong one here; none may stand out.
        for fit in fits where fit.explainedFraction < 0.3 {
            #expect(!(fit.isAboveSweep),
            "\(fit.zoneAxis) at \(fit.explainedFraction) stands out at \(fit.sweepRatio)× the median")
        }
    }

    /// A review (2026-09-15 night): with two candidate axes the "median" was
    /// the top axis's own fraction, so its ratio saturated at exactly 1 and
    /// a perfect fit could never clear the bar. Below five axes the sweep
    /// rule is inert. Mutation this names: the `>= 5` guard dropped.
    @Test func testATinySweepHasNoMedianNull() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let scale = 0.008
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, reference.kMaxInvAngstrom), crystal: al,
            zoneAxis: SIMD3(1, 1, 0), settings: reference)
        let peaks: [[SIMD2<Double>]] = (0..<16).map { _ in
            [pixel(128, 128)] + base.map {
                pixel(128 + Float($0.q.x / scale), 128 + Float($0.q.y / scale))
            }
        }
        let fits = PhaseVectorMatcher.fitZoneAxis(
            scan: ScanFixture(scanWidth: 4, scanHeight: 4, peaks: peaks),
            crystal: al, referenceSettings: reference, settings: PhaseVectorSettings(),
            originX: 128, originY: 128, invAngstromPerPixel: scale,
            candidateAxes: [SIMD3(1, 1, 0), SIMD3(0, 0, 1)], inPlaneStepDeg: 2)
        let winner = try #require(fits.first)
        #expect(winner.sweepMedianFraction == 0,
            "two axes are not a null")
        #expect(winner.isAboveSweep,
            "a perfect ⟨110⟩ fit failed a sweep of two: ratio \(winner.sweepRatio)")
    }

    /// Mutation: the `pairs > 0` guard dropped, or an axis with no reachable
    /// reflection admitted — it would divide by zero or report a perfect fit
    /// for a direction that presents nothing.
    @Test func testAnAxisThatPresentsNothingIsNotRanked() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let scan = ScanFixture(scanWidth: 2, scanHeight: 2,
                                 peaks: Array(repeating: [], count: 4))
        #expect(PhaseVectorMatcher.fitZoneAxis(
            scan: scan, crystal: al, referenceSettings: reference,
            settings: PhaseVectorSettings(), originX: 0, originY: 0,
            invAngstromPerPixel: 0.008).isEmpty,
            "a scan with no peaks produced a zone-axis ranking")
        for fit in PhaseVectorMatcher.fitZoneAxis(
            scan: ScanFixture(scanWidth: 1, scanHeight: 1,
                                peaks: [[pixel(90, 0)]]),
            crystal: al, referenceSettings: reference, settings: PhaseVectorSettings(),
            originX: 0, originY: 0, invAngstromPerPixel: 0.008) {
            #expect(fit.meanDistance.isFinite,
            "an unranked axis leaked a NaN")
            #expect(fit.totalVectors > 0)
        }
    }

    // MARK: The ranking is a total order (S4, 2026-09-29)

    private func fit(_ axis: SIMD3<Int>, matched: Int,
                     mean: Double) -> PhaseVectorMatcher.ZoneAxisFit {
        PhaseVectorMatcher.ZoneAxisFit(
            zoneAxis: axis, inPlaneRotationRad: 0, matchedVectors: matched,
            totalVectors: 100, meanDistance: mean, chanceMatchedVectors: 1,
            sweepMedianFraction: 0)
    }

    private func ranked(_ fits: [PhaseVectorMatcher.ZoneAxisFit]) -> [SIMD3<Int>] {
        fits.sorted(by: PhaseVectorMatcher.ZoneAxisFit.ranksBefore).map(\.zoneAxis)
    }

    /// The axis indices decide ONLY what the score cannot. Each half puts the
    /// better score on the axis that sorts LATER by indices, so the index key
    /// would have to lose to it.
    ///
    /// Mutations: the axis key moved ahead of the matched count (red on the
    /// first pair); ahead of the mean distance (red on the second).
    @Test func testTheAxisTieBreakNeverOverrulesADifferentScore() throws {
        let early = SIMD3(-1, -1, 0), late = SIMD3(1, 1, 1)
        #expect(ranked([fit(early, matched: 9, mean: 0.010),
                               fit(late, matched: 10, mean: 0.020)]).first == late,
            "a fit with fewer matched vectors won on its axis indices")
        #expect(ranked([fit(early, matched: 10, mean: 0.0101),
                               fit(late, matched: 10, mean: 0.0100)]).first == late,
            "a fit 1e-4 Å⁻¹ further out won on its axis indices")
    }

    /// A tied family ranks in the same order whatever order it arrives in —
    /// including when the mean distances differ only in the last bit, which
    /// is what rotating two symmetry-equivalent reference sets leaves behind.
    ///
    /// Mutations: the axis tie-break removed (red: the input order leaks
    /// through the stable sort); the mean distance compared raw instead of on
    /// `meanDistanceTieResolution`'s grid (red: the one-ulp-larger [-1 -1 0]
    /// drops behind [0 -1 1]).
    @Test func testATiedFamilyRanksTheSameInAnyInputOrder() throws {
        let d = 0.012_345_678
        let family = [
            fit(SIMD3(0, -1, 1), matched: 39, mean: d),
            fit(SIMD3(-1, 0, 1), matched: 39, mean: d),
            fit(SIMD3(-1, -1, 0), matched: 39, mean: d.nextUp),
            fit(SIMD3(0, -1, -1), matched: 39, mean: d.nextDown),
            fit(SIMD3(1, 1, 1), matched: 12, mean: 0.001),
        ]
        let expected = [SIMD3(-1, -1, 0), SIMD3(-1, 0, 1), SIMD3(0, -1, -1),
                        SIMD3(0, -1, 1), SIMD3(1, 1, 1)]
        #expect(ranked(family) == expected)
        #expect(ranked(family.reversed()) == expected)
        #expect(ranked([family[3], family[0], family[4], family[2], family[1]]) == expected)
    }

    /// The sweep itself returns that order: the same planted ⟨110⟩ fitted
    /// against the same axes listed forwards and backwards must rank them
    /// identically. Mutation: `fitZoneAxis` back on its own count-then-
    /// distance closure (red wherever the family ties exactly, since the
    /// stable sort then keeps each input order).
    @Test func testTheSweepOrderDoesNotDependOnTheCandidateOrder() throws {
        let al = CrystalStructure.aluminum
        var reference = PhaseReferenceSettings()
        reference.kMaxInvAngstrom = 1.2
        let base = PhaseReferenceLibrary.projectedVectors(
            reflections: refl(al, reference.kMaxInvAngstrom), crystal: al,
            zoneAxis: SIMD3(1, 1, 0), settings: reference)
        let turned = PhaseReferenceLibrary.rotate(base, by: 37.2 * .pi / 180)
        let scale = 0.008
        let originX: Float = 128, originY: Float = 128
        let peaks: [[SIMD2<Double>]] = (0..<16).map { _ in
            [pixel(originX, originY)]
                + turned.map { v in
                    pixel(originX + Float(v.q.x / scale), originY + Float(v.q.y / scale))
                }
        }
        let scan = ScanFixture(scanWidth: 4, scanHeight: 4, peaks: peaks)
        let axes = PhaseReferenceLibrary.lowIndexZoneAxes.filter {
            let m = [abs($0.x), abs($0.y), abs($0.z)].sorted()
            return m == [0, 1, 1] || m == [0, 0, 1]
        }
        #expect(axes.count >= 6,
            "too few axes to test the order")
        func sweep(_ candidates: [SIMD3<Int>]) -> [SIMD3<Int>] {
            PhaseVectorMatcher.fitZoneAxis(
                scan: scan, crystal: al, referenceSettings: reference,
                settings: PhaseVectorSettings(), originX: originX, originY: originY,
                invAngstromPerPixel: scale, candidateAxes: candidates,
                inPlaneStepDeg: 2).map(\.zoneAxis)
        }
        let forwards = sweep(axes)
        #expect(forwards.count == axes.count)
        #expect(sweep(axes.reversed()) == forwards,
            "the ranking followed the order the axes were listed in")
    }
}
