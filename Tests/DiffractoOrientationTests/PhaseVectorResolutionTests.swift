//
//  PhaseVectorResolutionTests.swift
//  Role: what the tolerances in Å⁻¹ mean on a particular detector.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct PhaseVectorResolutionTests {
    /// A real cube: 0.045741 Å⁻¹ per detector pixel after 4× binning.
    /// The shipped 0.020 Å⁻¹ tolerances are 0.44 of one pixel there, matrix
    /// removal removed nothing at all, and 108 899 of 108 900 positions came
    /// back "not indexed".
    ///
    /// Mutation: `pairRadiusPixels` computed as `scale / setting` (inverted),
    /// which is the easiest way to get a unit conversion backwards and would
    /// read as a comfortable 2.29 px instead of a hopeless 0.44.
    @Test func testTolerancesAreReportedInDetectorPixels() throws {
        let resolution = PhaseVectorResolution(settings: PhaseVectorSettings(),
                                               invAngstromPerPixel: 0.045741477608680726)
        #expect(abs(resolution.pairRadiusPixels - 0.437) <= 0.001)
        #expect(abs(resolution.matrixRemovalPixels - 0.437) <= 0.001)
        #expect(abs(resolution.notIndexedAbovePixels - 0.328) <= 0.001)
        // A finer detector meets the same settings comfortably, which is what
        // makes the number worth showing rather than the setting worth banning.
        let fine = PhaseVectorResolution(settings: PhaseVectorSettings(),
                                         invAngstromPerPixel: 0.008)
        #expect(abs(fine.pairRadiusPixels - 2.5) <= 1e-9)
        #expect(fine.advice == nil,
            "a detector that can meet the tolerance was warned about")
    }

    /// Mutation: the `pairRadiusPixels < 1` guard removed or inverted. A
    /// warning that fires on every dataset is one nobody reads.
    @Test func testAdviceFiresOnlyWhenTheToleranceIsUnderOnePixel() throws {
        func advice(_ scale: Double) -> String? {
            PhaseVectorResolution(settings: PhaseVectorSettings(),
                                  invAngstromPerPixel: scale).advice
        }
        #expect(advice(0.045741477608680726) != nil)
        #expect(advice(0.0201) != nil)
        #expect(advice(0.0199) == nil,
            "0.02 Å⁻¹ over a 0.0199 Å⁻¹ pixel is 1.005 px")
        #expect(advice(0.008) == nil)
        #expect(advice(0) == nil,
            "an uncalibrated detector has no pixel to compare against")
        #expect(advice(0.045741477608680726)?.contains("0.44") ?? false,
            "the advice must carry the number: \(advice(0.045741477608680726) ?? "nil")")
    }

    /// Mutation: `scaledToDetector` setting the verdict distance equal to the
    /// pair radius rather than three quarters of it (half until 2026-09-15,
    /// when the half was measured to reject honest many-vector fits). The
    /// pair radius says what COULD be the same reflection; the verdict
    /// distance says what is close enough to call — equal values collapse
    /// the distinction and index everything the pair radius admits.
    @Test func testScalingToTheDetectorKeepsTheShippedRatio() throws {
        let scale = 0.045741477608680726
        let scaled = PhaseVectorResolution(settings: PhaseVectorSettings(),
                                           invAngstromPerPixel: scale)
            .scaledToDetector(PhaseVectorSettings())
        #expect(abs(scaled.pairRadiusInvAngstrom - scale) <= 1e-12)
        #expect(abs(scaled.matrixToleranceInvAngstrom - scale) <= 1e-12)
        #expect(abs(scaled.notIndexedAboveInvAngstrom - (scale * 0.75)) <= 1e-12)
        #expect(scaled.notIndexedAboveInvAngstrom < scaled.pairRadiusInvAngstrom)
        // The direct beam is only ever widened, never narrowed: a user who set
        // it wide for a big probe must not have it cut by a scaling action.
        var wide = PhaseVectorSettings()
        wide.directBeamRadiusInvAngstrom = 0.5
        let keptWide = PhaseVectorResolution(settings: wide, invAngstromPerPixel: scale)
            .scaledToDetector(wide)
        #expect(abs(keptWide.directBeamRadiusInvAngstrom - 0.5) <= 1e-12)
        // And a run with no Q calibration changes nothing at all.
        let uncalibrated = PhaseVectorResolution(settings: PhaseVectorSettings(),
                                                 invAngstromPerPixel: 0)
            .scaledToDetector(PhaseVectorSettings())
        #expect(uncalibrated == PhaseVectorSettings())
    }

}
