//
//  PhaseVectorReviewFindingTests.swift
//  Role: the two review findings (2026-09-12) that were fixed in the matcher rather than only recorded; both were invisible to
//        every check that existed at the time.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct PhaseVectorReviewFindingTests {
    /// Review finding 9. `gcd` floored itself at 1, so `gcd(0, 0)` returned 1
    /// and `[0 0 2]` never reduced — the list held it beside `[0 0 1]`, 50
    /// entries for 49 directions, and 180 redundant library entries per phase
    /// at the default step. Redundant entries are not merely wasted work: each
    /// is one more comparison in the best-of-N search every scan position runs.
    ///
    /// Mutation: `return a` in `gcd` put back to `return max(a, 1)`.
    @Test func testLowIndexZoneAxesAreReducedAndDistinct() throws {
        let axes = PhaseReferenceLibrary.lowIndexZoneAxes
        #expect(!(axes.isEmpty))
        func gcd(_ a: Int, _ b: Int) -> Int {
            var a = abs(a), b = abs(b)
            while b != 0 { (a, b) = (b, a % b) }
            return a
        }
        for axis in axes {
            let g = gcd(gcd(axis.x, axis.y), axis.z)
            #expect(g == 1,
            "[\(axis.x) \(axis.y) \(axis.z)] is reducible by \(g)")
        }
        // No two entries name the same direction, in either sense.
        for i in axes.indices {
            for j in (i + 1)..<axes.count {
                let a = axes[i], b = axes[j]
                let parallel = a.x * b.y == a.y * b.x
                    && a.y * b.z == a.z * b.y && a.x * b.z == a.z * b.x
                #expect(!(parallel),
            "[\(a.x) \(a.y) \(a.z)] and [\(b.x) \(b.y) \(b.z)] are the same axis")
            }
        }
    }

    /// Review finding 6. The chance expectation is computed over the area the
    /// DATA occupies, and that radius used to be `max |u|` — so one spurious
    /// maximum far from the others inflated the area, weakened the guard for
    /// every other vector in the same pattern, and turned refusals into
    /// labels. Measured by the reviewer at 10 false positives in 1024
    /// patterns; second-largest makes it 0.
    ///
    /// Mutation: `accessibleRadius` put back to the running maximum.
    @Test func testOneDistantSpuriousPeakDoesNotWeakenTheChanceGuard() throws {
        let library = try PhaseReferenceLibrary.build(phases: [
            PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum,
                            role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
            PhaseDefinition(id: "beta", displayName: "β″", crystal: .betaDoublePrime,
                            role: .candidate, zoneAxes: [SIMD3(0, 1, 0)]),
        ], settings: {
            var s = PhaseReferenceSettings()
            s.kMaxInvAngstrom = 1.2
            s.inPlaneStepDeg = 2
            return s
        }())
        let settings = PhaseVectorSettings()
        let scratch = PhaseVectorMatcher.Scratch(
            capacity: library.entries.map(\.vectors.count).max() ?? 1)

        // Deterministic: a check that is not reproducible is not a check.
        var state: UInt64 = 0xD1B5_4A32_D192_ED03
        func uniform() -> Double {
            state ^= state >> 12; state ^= state << 25; state ^= state >> 27
            return Double((state &* 2685821657736338717) >> 11) * (1.0 / 9007199254740992.0)
        }
        var indexedWithOutlier = 0, indexedWithout = 0
        for _ in 0..<1024 {
            var vectors: [SIMD2<Double>] = []
            for _ in 0..<8 {
                let r = 0.25 + 0.15 * uniform(), a = 2 * Double.pi * uniform()
                vectors.append(SIMD2(r * cos(a), r * sin(a)))
            }
            func verdict(_ v: [SIMD2<Double>]) -> PhaseVerdict {
                PhaseVectorMatcher.classify(
                    vectors: v, library: library, settings: settings, matrixEntry: nil,
                    candidateEntryIndices: library.candidateEntryIndices, scratch: scratch
                ).verdict
            }
            if verdict(vectors) == .indexed { indexedWithout += 1 }
            if verdict(vectors + [SIMD2(1.15, 0)]) == .indexed { indexedWithOutlier += 1 }
        }
        #expect(indexedWithout == 0,
            "random vectors were indexed even without an outlier")
        #expect(indexedWithOutlier == 0,
            "\(indexedWithOutlier) of 1024 random patterns became a phase because one distant peak widened the chance expectation")
    }
}
