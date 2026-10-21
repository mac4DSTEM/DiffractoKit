//
//  ScatteringRefuterTests.swift
//  The Gate B refuter's findings on slice 2 (be6b832; registration slice 2 §8, amendment of 2026-10-15), each a test.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

struct ScatteringRefuterTests {

    private static var rockSalt: CrystalStructure {
        let fcc = Presets.fcc
        return Presets.structure(.cubic(a: 5.64), fcc.map { AtomSite(z: 11, fractional: $0) }
                                 + fcc.map { AtomSite(z: 17, fractional: wrapFractional($0 + SIMD3(0.5, 0.5, 0.5))) })
    }

    /// F2: a NaN or infinite form factor (for every element or one) is refused at the first reflection, never scored as
    /// nothing kept. Mutant: the finiteness check removed.
    @Test func testNonFiniteFormFactorsAreRefused() {
        let s = Self.rockSalt
        for ff: FormFactor in [{ _, _ in .nan }, { _, _ in .infinity }, { z, g in z == 17 ? .nan : ElectronScatteringFactors.f(z: z, gSquared: g) }] {
            #expect(throws: KinematicError.self) { try s.kinematicReflections(kMax: 1.2, formFactor: ff) }
            do { _ = try s.kinematicReflections(kMax: 1.2, formFactor: ff); Issue.record("admitted") }
            catch { guard case .nonFiniteStructureFactor = error else { Issue.record("\(error)"); continue } }
        }
    }

    /// F2: an invalid structure (NaN occupancy, infinite coordinate, occupancy 1.5, no sites) is refused before any F.
    /// Mutant: the validation not called.
    @Test func testInvalidStructuresAreRefused() {
        let cell = UnitCell.cubic(a: 5)
        let cases: [[AtomSite]] = [
            [AtomSite(z: 26, fractional: [0, 0, 0], occupancy: .nan)],
            [AtomSite(z: 26, fractional: [.infinity, 0, 0])],
            [AtomSite(z: 26, fractional: [0, 0, 0], occupancy: 1.5)],
            [],
        ]
        for sites in cases {
            let s = CrystalStructure(cell: cell, sites: sites, operations: [.identity])
            do { _ = try s.kinematicReflections(kMax: 1.2); Issue.record("admitted \(sites)") }
            catch { guard case .invalidStructure = error else { Issue.record("\(error)"); continue } }
        }
    }

    /// F3: a form factor with no value at some g² names the element; one with no value for one element anywhere names it.
    @Test func testMissingFormFactorNamesTheElement() {
        let s = Self.rockSalt
        let noneAboveZero: FormFactor = { z, g in g > 0 ? nil : ElectronScatteringFactors.f(z: z, gSquared: g) }
        #expect(throws: KinematicError.unsupportedElements([11, 17])) { try s.kinematicReflections(kMax: 1.2, formFactor: noneAboveZero) }
        let noChlorine: FormFactor = { z, g in z == 17 ? nil : ElectronScatteringFactors.f(z: z, gSquared: g) }
        #expect(throws: KinematicError.unsupportedElements([17])) { try s.kinematicReflections(kMax: 1.2, formFactor: noChlorine) }
    }

    /// F4: magnesium with ⅓ written as 0.333 keeps 24 operation-absent reflections above the cut (|F| up to 3e-4), in
    /// extra rings (28 against 24, 306 kept against 270, the refuter's numbers); the list names them. With exact thirds there are none. Mutant: the property filtering nothing.
    @Test func testOperationAbsentReflectionsAboveTheCutAreNamed() throws {
        let coarse = Presets.structure(.hexagonal(a: 3.2094, c: 5.2108),
                                       [AtomSite(z: 12, fractional: [0.333, 0.667, 0.25]),
                                        AtomSite(z: 12, fractional: [0.667, 0.333, 0.75])])
        let list = try coarse.kinematicReflections(kMax: 1.2)
        #expect(coarse.operations.count == 24)
        #expect(list.keptButOperationAbsent.count == 24)
        #expect(list.keptButOperationAbsent.allSatisfy { $0.structureFactor.magnitude > 1e-4 && $0.structureFactor.magnitude < 3e-4 })
        let exact = try Presets.hcp(12, a: 3.2094, c: 5.2108).kinematicReflections(kMax: 1.2)
        #expect(exact.keptButOperationAbsent.isEmpty)
        #expect(list.kept.count == 306 && exact.kept.count == 270, "kept \(list.kept.count), \(exact.kept.count)")
        #expect(list.families.count == 28 && exact.families.count == 24, "rings \(list.families.count), \(exact.families.count)")
    }

    /// F5: a needle-shaped cell with a tiny volume but a long enumeration box is refused, not enumerated.
    @Test func testNeedleCellIsRefused() {
        let needle = CrystalStructure(cell: .orthorhombic(a: 1e7, b: 1e-3, c: 1e-3) , sites: [], operations: [.identity])
        #expect(throws: ReflectionListError.self) { try needle.reflections(kMax: 1.2) }
    }

    /// F6: a coplanar cell (α = β = γ = 120°) whose volume is rounding is invalid, so no structure factor of 1e6 comes out.
    @Test func testCoplanarCellIsInvalid() {
        let flat = UnitCell(a: 5, b: 5, c: 5, alphaDeg: 120, betaDeg: 120, gammaDeg: 120)
        #expect(flat.issues.contains(.degenerateVolume))
        let s = CrystalStructure(cell: flat, sites: [AtomSite(z: 26, fractional: [0, 0, 0])], operations: [.identity])
        #expect(throws: KinematicError.self) { try s.kinematicReflections(kMax: 1.2) }
    }

    /// F8: f_e is nil for a negative or non-finite g².
    @Test func testScatteringFactorDomain() {
        #expect(ElectronScatteringFactors.f(z: 13, gSquared: -0.1) == nil)
        #expect(ElectronScatteringFactors.f(z: 13, gSquared: .nan) == nil)
        #expect(ElectronScatteringFactors.f(z: 13, gSquared: .infinity) == nil)
        #expect(ElectronScatteringFactors.f(z: 13, gSquared: 0) != nil)
    }
}
