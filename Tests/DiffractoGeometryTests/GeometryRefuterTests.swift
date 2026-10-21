//
//  GeometryRefuterTests.swift
//  The Gate B refuter's findings on slice 3a (a16cf84; registration slice 3 §7, amendment of 2026-10-15), each a test.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

struct GeometryRefuterTests {

    /// F1: above about 1.3e154 V the product overflows; λ is NaN, never 0. Mutant: the overflow guard removed.
    @Test func testWavelengthOverflowIsNaN() {
        for v in [1.4e154, 1e155, 1e300, Double.greatestFiniteMagnitude] {
            #expect(ElectronWavelength.angstroms(volts: v).isNaN, "\(v)")
        }
        #expect(ElectronWavelength.angstroms(volts: 2e5) == 0.025079340450548003)
    }

    /// F2: an overflowing camera length or Δq is NaN, never infinity. Mutant: the result guards removed.
    @Test func testOverflowingLengthsAreNaN() {
        #expect(RingCalibration.cameraLengthPixels(radius: 1e300, twoTheta: 1e-300).isNaN)
        #expect(RingCalibration.qPerPixel(d: 5e-324, radius: 1).isNaN)
        #expect(RingCalibration.qPerPixel(d: 1, radius: 5e-324).isNaN)
        #expect(RingCalibration.cameraLengthPixels(radius: 150, twoTheta: 0.010651260937060768).isFinite)
    }

    /// F4, F5: project keeps nothing for a non-unit, zero or non-finite n, or a negative or non-finite wavelength.
    @Test func testProjectRefusesBadDirectionsAndWavelengths() throws {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let au = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: fcc.map { AtomSite(z: 79, fractional: $0) })
        let kept = try au.kinematicReflections(kMax: 1.2).kept
        for n: SIMD3<Double> in [[0, 0, 0], [2, 0, 0], [0.5, 0, 0], [.nan, 0, 0], [.infinity, 0, 0]] {
            #expect(ZoneAxis.project(kept, onto: n, sgMax: 0.1).isEmpty, "\(n)")
            #expect(!ZoneAxis.isUnitDirection(n))
        }
        for lambda in [-0.025, Double.nan, .infinity] {
            #expect(ZoneAxis.project(kept, onto: [0, 0, 1], wavelengthAngstroms: lambda, sgMax: 0.1).isEmpty, "\(lambda)")
        }
        #expect(!ZoneAxis.project(kept, onto: [0, 0, 1], wavelengthAngstroms: 0.025, sgMax: 0.1).isEmpty)
        #expect(ZoneAxis.isUnitDirection([0, 0, 1]))
    }

    /// F6: with n the beam's travel direction, s_g = 0 is the Ewald condition |n/λ + g| = 1/λ (to first order exactly).
    @Test func testCurvedExcitationErrorIsTheEwaldConditionForTheTravelDirection() throws {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let au = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: fcc.map { AtomSite(z: 79, fractional: $0) })
        let lambda = 0.025
        let n = SIMD3<Double>(0, 0, 1)
        for p in ZoneAxis.project(try au.kinematicReflections(kMax: 1.2).kept, onto: n, wavelengthAngstroms: lambda, sgMax: 10) {
            let g = p.reflection.reflection.g
            let ewald = (simd_length_squared(n / lambda + g) - 1 / (lambda * lambda)) * lambda / 2   // = g·n + λ|g|²/2
            #expect(abs(p.sg - ewald) <= 1e-12, "\(p.reflection.reflection.hkl)")
        }
    }

    /// F7: the list overload uses the list's own k_max (the free-k_max overload with a larger k_max accepts what the
    /// list overload refuses as beyond k_max).
    @Test func testListOverloadUsesTheListsKMax() throws {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let au = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: fcc.map { AtomSite(z: 79, fractional: $0) })
        let list = try au.kinematicReflections(kMax: 0.5)
        let last = try #require(list.families.last)
        let q = last.family.g, sigma = q * 1e-2
        #expect(throws: RingIndexing.Refusal.self) { try RingIndexing.index(q: q, sigma: sigma, in: list) }
        #expect(throws: Never.self) { try RingIndexing.index(q: q, sigma: sigma, in: list.families, kMax: 4) }
    }
}

/// The Gate B refuter's two traps on slice 3b (704e5b0): an unbounded outer radius and an overflowing pattern size are
/// refused as invalid input, never trapped. Mutants: either guard removed.
struct EllipseRefuterTests {
    static func pattern() -> DiffractionPattern {
        let n = 64
        var pixels = [Float](repeating: 0, count: n * n)
        for r in 0..<n { for c in 0..<n {
            let d = (Double(r - 32) * Double(r - 32) + Double(c - 32) * Double(c - 32)).squareRoot()
            pixels[r * n + c] = Float(exp(-(d - 20) * (d - 20) / 8))
        } }
        return DiffractionPattern(width: n, height: n, pixels: pixels)
    }

    @Test func testUnboundedOuterRadiusIsInvalidInput() {
        let p = Self.pattern()
        for outer in [Double.infinity, 1e19, 1e9, 2 * (2.0 * 64 * 64).squareRoot() + 1] {
            #expect(throws: EllipseCalibration.FitError.self) {
                try EllipseCalibration.fitBestAvailable(pattern: p, centerQX: 32, centerQY: 32, innerRadius: 12, outerRadius: outer)
            }
            do { _ = try EllipseCalibration.fit1D(pattern: p, centerQX: 32, centerQY: 32, innerRadius: 12, outerRadius: outer) }
            catch { guard case EllipseCalibration.FitError.invalidInput = error else { Issue.record("\(outer): \(error)"); continue } }
        }
        #expect(throws: Never.self) {
            try EllipseCalibration.fitBestAvailable(pattern: p, centerQX: 32, centerQY: 32, innerRadius: 12, outerRadius: 28)
        }
    }

    @Test func testOverflowingPatternSizeIsInvalidInput() {
        let huge = DiffractionPattern(width: Int.max / 2 + 1, height: 3, pixels: [])
        do { _ = try EllipseCalibration.fit1D(pattern: huge, centerQX: 1, centerQY: 1, innerRadius: 1, outerRadius: 2); Issue.record("admitted") }
        catch { guard case EllipseCalibration.FitError.invalidInput = error else { Issue.record("\(error)"); return } }
    }
}
