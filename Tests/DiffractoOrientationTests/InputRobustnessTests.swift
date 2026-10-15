//
//  InputRobustnessTests.swift
//  Role: inputs a caller can pass that must be refused or answered, never crash on or hang on: a tiny in-plane step, an
//        inconsistent library, a cancelled build. Each test fails on the code it was written against (a trap shows as a crashed
//        test process).
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct InPlaneStepRobustnessTests {

    private func phases() -> [PhaseDefinition] {
        [PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
         PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime, role: .candidate,
                         zoneAxes: [SIMD3(0, 1, 0)])]
    }

    private func refusal(_ step: Double) -> PhaseReferenceLibrary.Failure? {
        var s = PhaseReferenceSettings()
        s.kMaxInvAngstrom = 1.0
        s.inPlaneStepDeg = step
        do { _ = try PhaseReferenceLibrary.build(phases: phases(), settings: s) } catch { return error as? PhaseReferenceLibrary.Failure }
        return nil
    }

    /// 360 / 1e-300 is 3.6e302, which `Int(_:)` traps on; 360 / 1e-9 is 3.6e11 elements, which is an allocation the machine
    /// does not have. Both are refused as too large before anything is built.
    @Test(arguments: [1e-300, 1e-9, 1e-5])
    func buildRefusesATinyStepBeforeAllocating(_ step: Double) {
        guard case .libraryTooLarge(let requested, let limit)? = refusal(step) else {
            Issue.record("a step of \(step) degrees was not refused as too large"); return
        }
        #expect(limit == PhaseReferenceSettings().maximumEntries)
        #expect(requested > limit)
    }

    /// The array form never traps either: a count past the cap is held at the cap (documented on `inPlaneSteps`).
    @Test(arguments: [1e-300, 1e-9])
    func inPlaneStepsHoldsATinyStepAtTheCap(_ step: Double) {
        var s = PhaseReferenceSettings()
        s.inPlaneStepDeg = step
        let steps = PhaseReferenceLibrary.inPlaneSteps(s)
        #expect(steps.count == PhaseReferenceLibrary.maximumInPlaneStepCount)
        #expect(steps.first == 0)
        #expect(steps.last! < 2 * Double.pi)
    }

    /// A step at the cap is still exact: 360 / 0.001 = 360 000 steps.
    @Test func inPlaneStepsAtTheCapAreUnchanged() {
        var s = PhaseReferenceSettings()
        s.inPlaneStepDeg = 360.0 / Double(PhaseReferenceLibrary.maximumInPlaneStepCount)
        #expect(PhaseReferenceLibrary.inPlaneSteps(s).count == PhaseReferenceLibrary.maximumInPlaneStepCount)
    }

    /// `fitZoneAxis` refuses a step whose count exceeds the cap (it has no `maximumEntries`), with the count it was asked for.
    @Test(arguments: [1e-300, 1e-9])
    func fitZoneAxisRefusesATinyStep(_ step: Double) {
        let scan: [[SIMD2<Double>]] = [[SIMD2(0.2, 0.0), SIMD2(0.0, 0.2)]]
        do {
            _ = try PhaseVectorMatcher.fitZoneAxis(
                scan: scan, scanWidth: 1, scanHeight: 1, crystal: .aluminum,
                referenceSettings: PhaseReferenceSettings(), settings: PhaseVectorSettings(),
                coordinates: .calibrated, candidateAxes: [SIMD3(0, 0, 1)], inPlaneStepDeg: step)
            Issue.record("a step of \(step) degrees was not refused")
        } catch {
            guard case .libraryTooLarge(_, let limit) = error else { Issue.record("wrong refusal \(error)"); return }
            #expect(limit == PhaseReferenceLibrary.maximumInPlaneStepCount)
        }
    }

    /// An ordinary step is not refused and the library is as before (bit identity is held by the parity tests).
    @Test func anOrdinaryStepStillBuilds() throws {
        var s = PhaseReferenceSettings()
        s.kMaxInvAngstrom = 1.0
        s.inPlaneStepDeg = 30
        let lib = try PhaseReferenceLibrary.build(phases: phases(), settings: s)
        #expect(lib.entries.count == 24)
    }
}
