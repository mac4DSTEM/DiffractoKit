//
//  SymmetryDetectionTests.swift
//  The operations of a structure given as a site list, on cells whose content is known by
//  construction: β″ expanded by hand (C2/m: 8 operations), gold P1 (Fm-3m: 192), iron (Im-3m: 96), magnesium
//  (P6₃/mmc: 24), the split-atom gold that must not be cubic, and a tetragonal cell with a = c.
//

import DiffractoKitTestSupport
import Testing
import simd
import DiffractoCore

struct SymmetryDetectionTests {

    static let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    /// β″'s 22 hand-expanded atoms give C2/m's eight operations and Laue class 2/m. Mutants: the C₂ about a; the
    /// target-distinct guard dropped (does not change this case; the split gold below).
    @Test func testBetaDoublePrimeByHandIsC2OverM() {
        let structure = CrystalStructure.detectingSymmetry(cell: Fixtures.betaDoublePrimeCell, sites: Fixtures.betaDoublePrimeSites)
        #expect(structure.sites.count == 22)
        #expect(structure.formula == "Mg10 Si12")
        #expect(structure.operations.count == 8, "\(structure.operations.map(\.xyzString))")
        #expect(structure.laueClass == .twoOverM)
        #expect(structure.crystalSystem == .monoclinic)
        #expect(structure.symmetrySource == .detected)
        #expect(structure.unresolvedRotations.isEmpty)
        // The C-centring is among the operations, as a translation with the identity rotation.
        #expect(structure.operations.contains {
            $0.rotation.isIdentity && $0.translationMod1 == SIMD3(0.5, 0.5, 0)
        })
        #expect(structure.operations.contains { $0.isEquivalent(to: try! SymmetryOperation(xyz: "-x,y,-z")) })
        #expect(structure.isUsable)
    }

    @Test func testGoldP1IsFm3m() {
        let sites = Self.fcc.map { AtomSite(z: 79, fractional: $0) }
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites)
        #expect(structure.operations.count == 192)
        #expect(structure.laueClass == .mThreeBarM)
        #expect(structure.laueGroup.order == 48)
        #expect(Set(structure.operations.map(\.rotation)).count == 48)
        #expect(structure.operations.filter { $0.rotation.isIdentity }.count == 4, "F-centring: four translations")
    }

    @Test func testIronIsIm3m() {
        let sites = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 2.8665), sites: sites)
        #expect(structure.operations.count == 96)
        #expect(structure.laueClass == .mThreeBarM)
    }

    /// Magnesium's 6₃ screw carries a translation; the structure is still 6/mmm (the reference implementations' reason to allow t ≠ 0).
    @Test func testMagnesiumIsP63OverMmc() {
        let sites = [AtomSite(z: 12, fractional: [0, 0, 0]), AtomSite(z: 12, fractional: [2.0 / 3, 1.0 / 3, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites)
        #expect(structure.operations.count == 24)
        #expect(structure.laueClass == .sixOverMmm)
        // With an atom at the origin (the reference implementations' cell, not ITA's origin on the 6₃ axis) the screw is
        // (x−y+⅔, x+⅓, z+½): a 6-fold rotation part with a half translation along c.
        let c6 = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        let screws = structure.operations.filter { $0.rotation == c6 }
        #expect(screws.count == 1)
        #expect(abs((screws.first?.translationMod1.z ?? -1) - 0.5) <= 1e-9)
        #expect(abs((screws.first?.translationMod1.x ?? -1) - (2.0 / 3)) <= 1e-9)
    }

    /// Silicon's diamond cell: Fd-3m, 192 operations (origin at an atom: the 3-fold is a pure rotation, the 4-fold a
    /// screw with a glide-like translation).
    @Test func testSiliconIsFd3m() {
        let shift = SIMD3<Double>(0.25, 0.25, 0.25)
        let sites = (Self.fcc + Self.fcc.map { $0 + shift }).map { AtomSite(z: 14, fractional: $0) }
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 5.4309), sites: sites)
        #expect(structure.operations.count == 192)
        #expect(structure.laueClass == .mThreeBarM)
    }

    /// The reference implementations' bijection fixture (reference `CIFSymmetryBijectionTests`): FCC gold plus two extra Au at
    /// (0.004, 0, 0) and (0, 0.004, 0), inside the 5e-3 tolerance of the origin atom but 0.0057 from each other and
    /// from where the 3-fold sends them. No translation makes the 3-fold a permutation of these six sites, so the cell
    /// is not cubic; a matcher that lets two sites share one target says it is. Mutant: the target-distinct guard
    /// dropped. (A dumbbell split within the tolerance, ±0.002 about a face centre, IS cubic at that tolerance and is
    /// not a test of the guard; a ±0.02 split is refused by distance alone, second loop.)
    @Test func testSplitSiteGoldIsNotCubic() {
        var sites = Self.fcc.map { AtomSite(z: 79, fractional: $0, occupancy: 1) }
        sites.append(AtomSite(z: 79, fractional: [0.004, 0, 0], occupancy: 1))
        sites.append(AtomSite(z: 79, fractional: [0, 0.004, 0], occupancy: 1))
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: sites)
        #expect(structure.laueClass != .mThreeBarM)
        #expect(structure.laueGroup.order < 48)
        #expect(!structure.isUsable, "0.0163 Å apart is a close contact")

        var dumbbell = Self.fcc.map { AtomSite(z: 79, fractional: $0, occupancy: 0.5) }
        dumbbell[1] = AtomSite(z: 79, fractional: [0.5, 0.5, 0.02], occupancy: 0.5)
        dumbbell.append(AtomSite(z: 79, fractional: [0.5, 0.5, 0.98], occupancy: 0.5))
        let split = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: dumbbell)
        #expect(split.laueClass != .mThreeBarM)
    }

    /// A metrically cubic cell (a = c) whose atoms carry only a 4-fold along c: 4/mmm, not m-3m (the reference implementations refused it).
    @Test func testTetragonalAtomsInACubicMetricAre4OverMmm() {
        let sites = [AtomSite(z: 22, fractional: [0, 0, 0]), AtomSite(z: 8, fractional: [0.5, 0.5, 0]),
                     AtomSite(z: 8, fractional: [0.5, 0.5, 0.5])]
        let structure = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4), sites: sites)
        #expect(structure.laueClass == .fourOverMmm)
        #expect(structure.crystalSystem == .tetragonal)
    }

    /// A triclinic cell with one atom: P-1, Laue -1.
    @Test func testOneAtomInATriclinicCellIsPOneBar() {
        let cell = UnitCell.triclinic(a: 4, b: 5, c: 6, alphaDeg: 80, betaDeg: 85, gammaDeg: 95)
        let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 6, fractional: [0.1, 0.2, 0.3])])
        #expect(structure.operations.count == 2)
        #expect(structure.laueClass == .oneBar)
    }

    /// A rhombohedral cell with one atom at the origin has the full -3m of its axes (12 operations).
    @Test func testOneAtomInRhombohedralAxesIsRThreeBarM() {
        let structure = CrystalStructure.detectingSymmetry(cell: .rhombohedral(a: 4, alphaDeg: 70),
                                                          sites: [AtomSite(z: 33, fractional: [0, 0, 0])])
        #expect(structure.operations.count == 12)
        #expect(structure.laueClass == .threeBarM)
    }

    /// Coordinates written too coarsely: the candidate fails at 5e-3, passes at 10 × the half-step, and is reported
    /// as unresolved, not taken.
    @Test func testCoarseCoordinatesReportUnresolvedCandidates() {
        // ⅓ and ⅔ truncated to two decimals in magnesium: residual 3.3e-3 per coordinate, amplified by the rows.
        let sites = [AtomSite(z: 12, fractional: [0, 0, 0]), AtomSite(z: 12, fractional: [0.66, 0.33, 0.5])]
        let exact = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites)
        let withPrecision = CrystalStructure.detectingSymmetry(cell: .hexagonal(a: 3.2094, c: 5.2108), sites: sites,
                                                               coarsestHalfStep: 1e-2)
        #expect(exact.operations.count == withPrecision.operations.count)
        #expect(exact.unresolvedRotations.isEmpty)
        if exact.laueClass != .sixOverMmm {
            #expect(!withPrecision.unresolvedRotations.isEmpty, "the failed candidates pass at 10 × 1e-2")
        }
    }

    @Test func testCloseContactIsRefusedAndASplitSiteIsNot() throws {
        let close = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: 5), sites: [AtomSite(z: 14, fractional: [0.99, 0.5, 0.5]), AtomSite(z: 14, fractional: [0.01, 0.5, 0.5])])
        let contact = try #require(close.shortestCloseContact)
        #expect(abs(contact - 0.1) <= 1e-6)
        #expect(close.validationIssues.map(\.code) == ["atomic_sites_too_close"])
        let split = CrystalStructure.detectingSymmetry(
            cell: .cubic(a: 5), sites: [AtomSite(z: 14, fractional: [0.49, 0.5, 0.5], occupancy: 0.5),
                                        AtomSite(z: 14, fractional: [0.51, 0.5, 0.5], occupancy: 0.5)])
        #expect(split.shortestCloseContact == nil)
        #expect(split.isUsable)
    }
}
