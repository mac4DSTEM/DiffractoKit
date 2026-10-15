//
//  ReflectionOrderTests.swift
//  Role: the order of the reflections a library is made from (registration slice 4, section 2.1). The density cap's sort and
//        the final length sort break ties on input order, so the input order is part of the definition: it must be the
//        order of the original reflection list (h, k, l enumerated ascending, filtered, then sorted by increasing |g| with
//        a stable sort), not the kit's own order of equal |g| (descending h, k, l). These tests re-derive the original
//        order from the kit's lattice by the original loop and compare it with `PhaseReferenceLibrary.reflections`.
//

import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

@Suite struct ReflectionOrderTests {

    /// The original loop: the per-axis index box, h then k then l ascending, 0 < |g| ≤ kMax, |F| above the cut, then a
    /// sort by |g| that keeps the enumeration order of equal |g| (the offset is the explicit tie key, so this does not
    /// depend on a sort being stable).
    private func originalOrder(_ s: CrystalStructure, _ kMax: Double) -> [MillerIndex] {
        let tile = s.cell.realLattice.map { max(0, Int((kMax * length($0)).rounded(.up))) }
        var out: [(hkl: MillerIndex, gLength: Double)] = []
        for h in -tile[0]...tile[0] {
            for k in -tile[1]...tile[1] {
                for l in -tile[2]...tile[2] {
                    if h == 0 && k == 0 && l == 0 { continue }
                    let hkl = MillerIndex(h, k, l)
                    let gLen = length(s.cell.gVector(hkl))
                    if gLen > kMax { continue }
                    guard let f = StructureFactors.compute(hkl, in: s) else { continue }
                    if f.intensity.squareRoot() <= 1e-4 { continue }
                    out.append((hkl, gLen))
                }
            }
        }
        return out.enumerated().sorted {
            $0.element.gLength != $1.element.gLength ? $0.element.gLength < $1.element.gLength : $0.offset < $1.offset
        }.map(\.element.hkl)
    }

    private var structures: [(name: String, structure: CrystalStructure)] {
        [
            ("Al fcc", .aluminum),
            ("Au fcc", .fcc(a: 4.0782, z: 79)),
            ("simple cubic 10 Å", .sc(a: 10, z: 13)),
            ("beta-double-prime", .betaDoublePrime),
            ("hexagonal Mg", makeCrystal(a: 3.2094, b: 3.2094, c: 5.2108, gammaDeg: 120, sites: [
                AtomSite(z: 12, fractional: [0, 0, 0]),
                AtomSite(z: 12, fractional: [2.0 / 3.0, 1.0 / 3.0, 0.5]),
            ])),
            ("triclinic", makeCrystal(a: 6.1, b: 7.3, c: 8.9, alphaDeg: 81, betaDeg: 97, gammaDeg: 104, sites: [
                AtomSite(z: 14, fractional: [0.1, 0.2, 0.3]),
                AtomSite(z: 8, fractional: [0.6, 0.45, 0.7]),
            ])),
        ]
    }

    /// Mutation: the kept list returned in the kit's own order (equal |g| by descending h, k, l), or the tie key dropped so
    /// the order of equal |g| is the sort's whim. Every structure, two k-max values.
    @Test func theLibraryReflectionsAreInTheOrderOfTheOriginalList() throws {
        for (name, structure) in structures {
            for kMax in [1.0, 1.6] {
                let got = try PhaseReferenceLibrary.reflections(of: structure, kMaxInvAngstrom: kMax).map(\.reflection.hkl)
                let expected = originalOrder(structure, kMax)
                #expect(!expected.isEmpty, "\(name) kMax \(kMax): no reflections, the test is empty")
                #expect(got == expected, "\(name) kMax \(kMax): the order differs from the original list")
            }
        }
    }

    /// The reordering is necessary: on a cubic cell the kit's own list orders equal |g| the other way round, so using it as
    /// it is would be a different library. If the two ever agree on every structure this test goes red and the helper's
    /// sort can be dropped.
    @Test func theKitsOwnOrderOfEqualLengthsIsDifferentSoTheReorderingIsNeeded() throws {
        let al = CrystalStructure.aluminum
        let kits = try al.kinematicReflections(kMax: 1.6).kept.map(\.reflection.hkl)
        let ours = try PhaseReferenceLibrary.reflections(of: al, kMaxInvAngstrom: 1.6).map(\.reflection.hkl)
        #expect(Set(kits) == Set(ours), "the same reflections, differently ordered")
        #expect(kits != ours)
        // {111}: the kit lists (1,1,1) before (-1,-1,-1); the original list has the opposite.
        #expect(kits.first == MillerIndex(1, 1, 1))
        #expect(ours.first == MillerIndex(-1, -1, -1))
    }

    /// Increasing |g| throughout, and each reflection's fields are the kit's: g Cartesian in Å⁻¹ from the reciprocal
    /// lattice, intensity |F|² above the cut.
    @Test func theFieldsAreTheKitsAndTheListIsSortedByLength() throws {
        let beta = CrystalStructure.betaDoublePrime
        let list = try PhaseReferenceLibrary.reflections(of: beta, kMaxInvAngstrom: 1.6)
        for pair in zip(list, list.dropFirst()) {
            #expect(pair.0.reflection.gLength <= pair.1.reflection.gLength)
        }
        for r in list {
            #expect(r.reflection.g == beta.cell.gVector(r.reflection.hkl))
            #expect(r.reflection.gLength <= 1.6)
            #expect(r.structureFactor.intensity.squareRoot() > StructureFactors.defaultCut)
        }
    }

    /// An element the scattering table does not cover is named, and a structure the kit refuses is refused with its
    /// reason, not silently scored as nothing.
    @Test func unsupportedElementsAndRefusedStructuresThrowWithTheirPhase() {
        let unknown = makeCrystal(a: 4, b: 4, c: 4, sites: [AtomSite(z: 104, fractional: [0, 0, 0])])
        #expect(PhaseReferenceLibrary.unsupportedElements(of: unknown) == [104])
        #expect(PhaseReferenceLibrary.unsupportedElements(of: .aluminum).isEmpty)
        #expect(throws: PhaseReferenceLibrary.Failure.unsupportedElements(phase: "X", z: [104])) {
            try PhaseReferenceLibrary.build(phases: [
                PhaseDefinition(id: "m", displayName: "M", crystal: .aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
                PhaseDefinition(id: "x", displayName: "X", crystal: unknown, role: .candidate, zoneAxes: [SIMD3(0, 0, 1)]),
            ])
        }

        // Two atoms on one site: a close contact, which the kit's structure validation refuses.
        let doubled = makeCrystal(a: 4, b: 4, c: 4, sites: [
            AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0, 0, 0]),
        ])
        do {
            _ = try PhaseReferenceLibrary.build(phases: [
                PhaseDefinition(id: "m", displayName: "M", crystal: .aluminum, role: .matrix, zoneAxes: [SIMD3(0, 0, 1)]),
                PhaseDefinition(id: "d", displayName: "Doubled", crystal: doubled, role: .candidate, zoneAxes: [SIMD3(0, 0, 1)]),
            ])
            Issue.record("a structure with two atoms on one site was accepted")
        } catch let error as PhaseReferenceLibrary.Failure {
            guard case .structureRefused(let phase, let reason) = error else {
                Issue.record("expected structureRefused, got \(error)")
                return
            }
            #expect(phase == "Doubled")
            #expect(!reason.isEmpty)
        } catch {
            Issue.record("unexpected error \(error)")
        }
    }
}
