//
//  PhaseReferenceLibrarySlabTests.swift
//  Role: the excitation slab is a property of a phase's shape: the per-phase override widens only its own phase's entries and
//        its absence changes nothing (registration slice 4, gate G5).
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct PhaseReferenceLibrarySlabTests {
    // MARK: Fixture
    //
    // A simple-cubic, single-atom-at-the-origin crystal, a = 10 Å, sampled at
    // its [001] zone with kMax 0.15 Å⁻¹. For this cell g_hkl = (h/10, k/10,
    // l/10) exactly (orthogonal, isotropic), so at zone axis [0 0 1]:
    //   • sg = g·n = l/10, independent of h and k;
    //   • the in-plane length |q| = √(h²+k²)/10, independent of l.
    // At kMax 0.15 Å⁻¹ the sphere |g| ≤ 0.15 admits only h²+k²+l² ∈ {1, 2}
    // (3 needs |g| = √3/10 = 0.173 > 0.15), enumerated once, by hand:
    //   l = 0, h²+k² ∈ {1,2}: (±1,0,0) (0,±1,0) (±1,±1,0)         →  8 vectors
    //   l = ±1, h²+k² ∈ {0,1}: (0,0,±1) has |q| = 0 (excluded,
    //          length > 1e-9 required); (±1,0,±1) (0,±1,±1)       → 4+4 = 8 more
    // `sg` for l = ±1 is 0.1 Å⁻¹: above the shipped default slab (0.05, so
    // excluded there) and within an override of 0.3 (so included). A single
    // atom at the origin has a real, positive structure factor at every
    // (h,k,l) (phase factor is 1), so nothing here is pruned by intensity —
    // `minimumIntensityFraction` is 0 and `maximumVectorsPerEntry` is well
    // above 16, so the ONLY gates in play are the slab and the length cut.
    // A wrong phase index (the override landing on the matrix instead of the
    // candidate, or vice versa) changes which of the two pinned counts below
    // moves, so it goes red rather than passing by coincidence.

    private func slabTestSettings() -> PhaseReferenceSettings {
        var settings = PhaseReferenceSettings()
        settings.kMaxInvAngstrom = 0.15
        settings.inPlaneStepDeg = 360     // one entry per zone axis; rotation is not this file's question
        settings.minimumIntensityFraction = 0
        settings.maximumVectorsPerEntry = 100
        return settings
    }

    private func slabMatrixPhase() -> PhaseDefinition {
        PhaseDefinition(id: "matrix", displayName: "Matrix", crystal: .sc(a: 10, z: 13),
                        role: .matrix, zoneAxes: [SIMD3(0, 0, 1)])
    }

    private func slabCandidatePhase(excitationSlabInvAngstrom: Double? = nil) -> PhaseDefinition {
        PhaseDefinition(id: "candidate", displayName: "Candidate", crystal: .sc(a: 10, z: 13),
                        role: .candidate, zoneAxes: [SIMD3(0, 0, 1)],
                        excitationSlabInvAngstrom: excitationSlabInvAngstrom)
    }

    private func vectorCount(_ library: PhaseReferenceLibrary, phaseIndex: Int) -> Int {
        library.entries.first(where: { $0.phaseIndex == phaseIndex })!.vectors.count
    }

    // MARK: (1) The override widens only its own phase's entry

    /// Mutation this catches: `build()` applying ANY phase's override to
    /// EVERY phase (M1 -- the matrix count would move from 8 to 16, failing
    /// the "untouched" assertion), or ignoring `PhaseDefinition
    /// .excitationSlabInvAngstrom` entirely (M2 -- the candidate count would
    /// stay 8 instead of reaching 16). Both were run and confirmed red
    /// (recorded when the test was written).
    @Test func testSlabOverrideAddsVectorsOnlyToItsOwnPhase() throws {
        let settings = slabTestSettings()
        let baseline = try PhaseReferenceLibrary.build(
            phases: [slabMatrixPhase(), slabCandidatePhase()], settings: settings)
        let overridden = try PhaseReferenceLibrary.build(
            phases: [slabMatrixPhase(), slabCandidatePhase(excitationSlabInvAngstrom: 0.3)],
            settings: settings)

        #expect(vectorCount(baseline, phaseIndex: 0) == 8,
            "matrix, baseline")
        #expect(vectorCount(baseline, phaseIndex: 1) == 8,
            "candidate, baseline (no override)")
        #expect(vectorCount(overridden, phaseIndex: 0) == 8,
            "the matrix carries no override and must be untouched by the candidate's")
        #expect(vectorCount(overridden, phaseIndex: 1) == 16,
            "the candidate's own override must admit its |l| = 1 reflections (sg = 0.1 Å⁻¹)")
        #expect(vectorCount(overridden, phaseIndex: 1) > vectorCount(baseline, phaseIndex: 1),
            "the override must ADD vectors relative to the same phase without it")
    }

    // MARK: (2) A nil override changes nothing

    /// A non-randomized fingerprint of a library's entries. Deliberately NOT
    /// Swift's `Hasher`/`hashValue`, which reseed every process launch by
    /// design and so cannot be pinned as an expected constant across runs.
    private func stableChecksum(of entries: [PhaseOrientationReference]) -> UInt64 {
        var hash: UInt64 = 0xcbf29ce484222325          // FNV-1a 64-bit offset basis
        let prime: UInt64 = 0x100000001b3
        func mix(_ s: String) {
            for byte in s.utf8 {
                hash ^= UInt64(byte)
                hash = hash &* prime
            }
        }
        for entry in entries {
            mix("phase\(entry.phaseIndex)|axis(\(entry.zoneAxis.x),\(entry.zoneAxis.y),\(entry.zoneAxis.z))"
                + "|rot\(String(format: "%.6f", entry.inPlaneRotationRad))")
            for v in entry.vectors {
                mix("h\(v.h)k\(v.k)l\(v.l)|q(\(String(format: "%.6f", v.q.x)),"
                    + "\(String(format: "%.6f", v.q.y)))|len\(String(format: "%.6f", v.length))"
                    + "|ri\(String(format: "%.6f", v.relativeIntensity))")
            }
        }
        return hash
    }

    /// Mutation this catches: any change to `build()` or `projectedVectors`
    /// that alters output on the path `phase.excitationSlabInvAngstrom ??
    /// settings.excitationSlabInvAngstrom` takes when the new field is left
    /// `nil` -- the path every pre-existing call site (every one before this
    /// slice) takes, unconditionally. The checksum was pinned by running this
    /// exact fixture once (recorded when the test was written) and was confirmed to change when
    /// the fixture's own `kMaxInvAngstrom` was perturbed, so it is not
    /// vacuously true.
    @Test func testNilOverrideReproducesThePreChangeEntriesExactly() throws {
        let settings = slabTestSettings()
        let library = try PhaseReferenceLibrary.build(
            phases: [slabMatrixPhase(), slabCandidatePhase()], settings: settings)
        // Pinned 2026-09-21 from this exact fixture (a one-off run of this fixture): 0x675b3e13f7e44cea.
        #expect(stableChecksum(of: library.entries) == (0x675b3e13f7e44cea as UInt64),
            "nil override must reproduce byte-identical entries to the pre-existing code path")
    }
}
