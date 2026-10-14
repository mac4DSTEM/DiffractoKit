//
//  SymmetryDetectionTests.swift
//  The operations of a structure given as a site list (registration 3.4b), on cells whose content is known by
//  construction: β″ expanded by hand (C2/m: 8 operations), gold P1 (Fm-3m: 192), iron (Im-3m: 96), magnesium
//  (P6₃/mmc: 24), the split-atom gold that must not be cubic, and a tetragonal cell with a = c.
//

import XCTest
import simd
import DiffractoCore

final class SymmetryDetectionTests: XCTestCase {

    static let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    /// β″'s 22 hand-expanded atoms give C2/m's eight operations and Laue class 2/m. Mutants: the C₂ about a; the
    /// target-distinct guard dropped (does not change this case; the split gold below).
    func testBetaDoublePrimeByHandIsC2OverM() {
        let structure = CrystalStructure.detectingSymmetry(cell: Fixtures.betaDoublePrimeCell, sites: Fixtures.betaDoublePrimeSites)
        XCTAssertEqual(structure.sites.count, 22)
        XCTAssertEqual(structure.formula, "Mg10 Si12")
        XCTAssertEqual(structure.operations.count, 8, "\(structure.operations.map(\.xyzString))")
        XCTAssertEqual(structure.laueClass, .twoOverM)
        XCTAssertEqual(structure.crystalSystem, .monoclinic)
        XCTAssertEqual(structure.symmetrySource, .detected)
        XCTAssertTrue(structure.unresolvedRotations.isEmpty)
        // The C-centring is among the operations, as a translation with the identity rotation.
        XCTAssertTrue(structure.operations.contains {
            $0.rotation.isIdentity && $0.translationMod1 == SIMD3(0.5, 0.5, 0)
        })
        XCTAssertTrue(structure.operations.contains { $0.isEquivalent(to: try! SymmetryOperation(xyz: "-x,y,-z")) })
        XCTAssertTrue(structure.isUsable)
    }

    func testGoldP1IsFm3m() {
        let sites = Self.fcc.map { AtomSite(z: 79, fractional: $0) }
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites)
        XCTAssertEqual(structure.operations.count, 192)
        XCTAssertEqual(structure.laueClass, .mThreeBarM)
        XCTAssertEqual(structure.laueGroup.order, 48)
        XCTAssertEqual(Set(structure.operations.map(\.rotation)).count, 48)
        XCTAssertEqual(structure.operations.filter { $0.rotation.isIdentity }.count, 4, "F-centring: four translations")
    }

    func testIronIsIm3m() {
        let sites = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 2.8665), sites: sites)
        XCTAssertEqual(structure.operations.count, 96)
        XCTAssertEqual(structure.laueClass, .mThreeBarM)
    }

    /// Magnesium's 6₃ screw carries a translation; the structure is still 6/mmm (the apps' reason to allow t ≠ 0).
    func testMagnesiumIsP63OverMmc() {
        let sites = [AtomSite(z: 12, fractional: [0, 0, 0]), AtomSite(z: 12, fractional: [2.0 / 3, 1.0 / 3, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites)
        XCTAssertEqual(structure.operations.count, 24)
        XCTAssertEqual(structure.laueClass, .sixOverMmm)
        // With an atom at the origin (the apps' cell, not ITA's origin on the 6₃ axis) the screw is
        // (x−y+⅔, x+⅓, z+½): a 6-fold rotation part with a half translation along c.
        let c6 = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        let screws = structure.operations.filter { $0.rotation == c6 }
        XCTAssertEqual(screws.count, 1)
        XCTAssertEqual(screws.first?.translationMod1.z ?? -1, 0.5, accuracy: 1e-9)
        XCTAssertEqual(screws.first?.translationMod1.x ?? -1, 2.0 / 3, accuracy: 1e-9)
    }

    /// Silicon's diamond cell: Fd-3m, 192 operations (origin at an atom: the 3-fold is a pure rotation, the 4-fold a
    /// screw with a glide-like translation).
    func testSiliconIsFd3m() {
        let shift = SIMD3<Double>(0.25, 0.25, 0.25)
        let sites = (Self.fcc + Self.fcc.map { $0 + shift }).map { AtomSite(z: 14, fractional: $0) }
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 5.4309), sites: sites)
        XCTAssertEqual(structure.operations.count, 192)
        XCTAssertEqual(structure.laueClass, .mThreeBarM)
    }

    /// An FCC cell whose face-centre atom is split into a dumbbell is at best tetragonal; the permutation test must
    /// not let two sites claim one target. Mutant: the target-distinct guard dropped.
    func testSplitSiteGoldIsNotCubic() {
        var sites = Self.fcc.map { AtomSite(z: 79, fractional: $0, occupancy: 1) }
        sites[1] = AtomSite(z: 79, fractional: [0.5, 0.5, 0.02], occupancy: 0.5)
        sites.append(AtomSite(z: 79, fractional: [0.5, 0.5, -0.02], occupancy: 0.5))
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites)
        XCTAssertNotEqual(structure.laueClass, .mThreeBarM)
        XCTAssertLessThan(structure.laueGroup.order, 48)
    }

    /// A metrically cubic cell (a = c) whose atoms carry only a 4-fold along c: 4/mmm, not m-3m (the apps refused it).
    func testTetragonalAtomsInACubicMetricAre4OverMmm() {
        let sites = [AtomSite(z: 22, fractional: [0, 0, 0]), AtomSite(z: 8, fractional: [0.5, 0.5, 0]),
                     AtomSite(z: 8, fractional: [0.5, 0.5, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4), sites: sites)
        XCTAssertEqual(structure.laueClass, .fourOverMmm)
        XCTAssertEqual(structure.crystalSystem, .tetragonal)
    }

    /// A triclinic cell with one atom: P-1, Laue -1.
    func testOneAtomInATriclinicCellIsPOneBar() {
        let cell = UnitCell.triclinic(a: 4, b: 5, c: 6, alphaDeg: 80, betaDeg: 85, gammaDeg: 95)
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 6, fractional: [0.1, 0.2, 0.3])])
        XCTAssertEqual(structure.operations.count, 2)
        XCTAssertEqual(structure.laueClass, .oneBar)
    }

    /// A rhombohedral cell with one atom at the origin has the full -3m of its axes (12 operations).
    func testOneAtomInRhombohedralAxesIsRThreeBarM() {
        let structure = CrystalStructure.detectingSymmetry(cell: .rhombohedral(a: 4, alphaDeg: 70),
                                                          sites: [AtomSite(z: 33, fractional: [0, 0, 0])])
        XCTAssertEqual(structure.operations.count, 12)
        XCTAssertEqual(structure.laueClass, .threeBarM)
    }

    /// Coordinates written too coarsely: the candidate fails at 5e-3, passes at 10 × the half-step, and is reported
    /// as unresolved, not taken.
    func testCoarseCoordinatesReportUnresolvedCandidates() {
        // ⅓ and ⅔ truncated to two decimals in magnesium: residual 3.3e-3 per coordinate, amplified by the rows.
        let sites = [AtomSite(z: 12, fractional: [0, 0, 0]), AtomSite(z: 12, fractional: [0.66, 0.33, 0.5])]
        let exact = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites)
        let withPrecision = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites,
                                                               coarsestHalfStep: 1e-2)
        XCTAssertEqual(exact.operations.count, withPrecision.operations.count)
        XCTAssertTrue(exact.unresolvedRotations.isEmpty)
        if exact.laueClass != .sixOverMmm {
            XCTAssertFalse(withPrecision.unresolvedRotations.isEmpty, "the failed candidates pass at 10 × 1e-2")
        }
    }

    func testCloseContactIsRefusedAndASplitSiteIsNot() {
        let close = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: 5), sites: [AtomSite(z: 14, fractional: [0.99, 0.5, 0.5]), AtomSite(z: 14, fractional: [0.01, 0.5, 0.5])])
        XCTAssertEqual(try XCTUnwrap(close.shortestCloseContact), 0.1, accuracy: 1e-6)
        XCTAssertEqual(close.validationIssues.map(\.code), ["atomic_sites_too_close"])
        let split = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: 5), sites: [AtomSite(z: 14, fractional: [0.49, 0.5, 0.5], occupancy: 0.5),
                                        AtomSite(z: 14, fractional: [0.51, 0.5, 0.5], occupancy: 0.5)])
        XCTAssertNil(split.shortestCloseContact)
        XCTAssertTrue(split.isUsable)
    }
}
