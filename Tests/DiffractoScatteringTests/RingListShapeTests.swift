//
//  RingListShapeTests.swift
//  The shape of a crystal's ring list (`kinematicReflections(kMax:).families`) on cells near a higher symmetry, on lower-symmetry
//  cells and on a monoclinic cell with a hidden swap symmetry: the metric family each cell reads, its Laue class, how many rings
//  it has and how they are named. Moved from macTEM's app tests, which held them before the kit owned this code; each test
//  names the mutants that turn it red.
//

import Foundation
import Testing
import DiffractoCore
import DiffractoScattering

struct RingListShapeTests {

    private static let au = CrystalPresets.gold
    private static let auA = 4.0782

    private static func structure(_ cell: UnitCell, _ sites: [AtomSite]) -> CrystalStructure {
        CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
    }

    /// The rings of `s` after the |F| cut.
    private static func rings(_ s: CrystalStructure, kMax: Double = 1.2) throws -> [ReflectionFamily] {
        try s.kinematicReflections(kMax: kMax).families.map(\.family)
    }

    private static func family(_ l: [ReflectionFamily], _ m: MillerIndex) -> ReflectionFamily? { l.first { $0.members.contains(m) } }

    // MARK: The metric family and the Laue class

    /// The cell's family is read within 1e-9 (relative lengths, degrees), and the Laue class from the atoms: b = a(1 + 1e-11)
    /// stays m-3m (8 rings) and gamma = 120 deg + 1e-10 deg 6/mmm (24); b = a(1 + 1e-6) and 1e-8 read 4/mmm along b (15 rings),
    /// gamma + 1e-8 deg 2/m (24), alpha + 1e-8 deg -3m (8). Mutants: the metric tolerance 0 (the near cells fall to a lower
    /// class), any cell taken as cubic.
    @Test func metricToleranceKeepsNearCubicAndNearHexagonalCells() throws {
        let a = Self.auA
        let fcc = Self.au.sites
        func read(_ cell: UnitCell, _ sites: [AtomSite]) throws -> (LaueClass?, UnitCell.MetricFamily, Int) {
            let s = Self.structure(cell, sites)
            return (s.laueClass, s.cell.metricFamily(), try Self.rings(s).count)
        }
        let nearCubic = try read(UnitCell(a: a, b: a * (1 + 1e-11), c: a, alphaDeg: 90, betaDeg: 90, gammaDeg: 90), fcc)
        #expect(nearCubic.0 == .mThreeBarM)
        #expect(nearCubic.1 == .cubic)
        #expect(nearCubic.2 == 8)
        let mg = CrystalPresets.magnesium
        let nearHex = try read(UnitCell(a: mg.cell.a, b: mg.cell.b, c: mg.cell.c, alphaDeg: 90, betaDeg: 90, gammaDeg: 120 + 1e-10), mg.sites)
        #expect(nearHex.0 == .sixOverMmm)
        #expect(nearHex.1 == .hexagonal)
        #expect(nearHex.2 == 24)
        for eps in [1e-6, 1e-8] {
            // a = c with b different is tetragonal along b.
            let off = try read(UnitCell(a: a, b: a * (1 + eps), c: a, alphaDeg: 90, betaDeg: 90, gammaDeg: 90), fcc)
            #expect(off.0 == .fourOverMmm, "b = a(1 + \(eps))")
            #expect(off.1 == .tetragonal(uniqueAxis: .b), "b = a(1 + \(eps))")
            #expect(off.2 == 15, "b = a(1 + \(eps))")
        }
        let gamma = try read(UnitCell(a: mg.cell.a, b: mg.cell.b, c: mg.cell.c, alphaDeg: 90, betaDeg: 90, gammaDeg: 120 + 1e-8), mg.sites)
        #expect(gamma.0 == .twoOverM)
        #expect(gamma.1 == .monoclinic(uniqueAxis: .c))
        #expect(gamma.2 == 24)
        let alpha = try read(UnitCell(a: a, b: a, c: a, alphaDeg: 90 + 1e-8, betaDeg: 90, gammaDeg: 90), fcc)
        #expect(alpha.0 == .threeBarM)
        #expect(alpha.1 == .rhombohedral)
        #expect(alpha.2 == 8)
    }

    /// The cell's family from its six parameters (at 1e-9): unique axis b is monoclinic; beta within 1e-9 deg of 90 is
    /// orthorhombic, 1e-8 off is monoclinic; alpha or gamma 1e-8 off 90 with beta off 90 is triclinic; unique axes c and a are
    /// monoclinic about their own axis; cubic and hexagonal are unchanged. Mutants: the family tolerance 1e-7 (beta 90 deg + 1e-8
    /// reads orthorhombic), a unique axis taken as b, alpha / gamma ignored.
    @Test func latticeRecognition() {
        func f(_ a: Double, _ b: Double, _ c: Double, _ al: Double, _ be: Double, _ ga: Double) -> UnitCell.MetricFamily {
            UnitCell(a: a, b: b, c: c, alphaDeg: al, betaDeg: be, gammaDeg: ga).metricFamily()
        }
        let mP = Self.structure(UnitCell(a: 7.2, b: 4.9, c: 5.6, alphaDeg: 90, betaDeg: 104, gammaDeg: 90), [AtomSite(z: 13, fractional: [0, 0, 0])])
        #expect(mP.cell.metricFamily() == .monoclinic(uniqueAxis: .b))
        #expect(f(15.16, 4.05, 6.74, 90, 105.3, 90) == .monoclinic(uniqueAxis: .b))
        #expect(f(7.2, 4.9, 5.6, 90, 90 + 1e-8, 90) == .monoclinic(uniqueAxis: .b))
        #expect(f(7.2, 4.9, 5.6, 90, 96, 90) == .monoclinic(uniqueAxis: .b))
        #expect(f(7.2, 4.9, 5.6, 90, 90, 90) == .orthorhombic)
        #expect(f(7.2, 4.9, 5.6, 90, 90 + 1e-10, 90) == .orthorhombic, "within 1e-9 deg of 90")
        #expect(f(7.2, 4.9, 5.6, 90 + 1e-8, 104, 90) == .triclinic)
        #expect(f(7.2, 4.9, 5.6, 90 - 1e-8, 104, 90) == .triclinic)
        #expect(f(7.2, 4.9, 5.6, 90, 104, 90 + 1e-8) == .triclinic)
        #expect(f(7.2, 4.9, 5.6, 90, 90, 100) == .monoclinic(uniqueAxis: .c))
        #expect(f(7.2, 4.9, 5.6, 100, 90, 90) == .monoclinic(uniqueAxis: .a))
        #expect(f(7.2, 4.9, 5.6, 92, 95, 101) == .triclinic)
        #expect(CrystalPresets.gold.cell.metricFamily() == .cubic)
        #expect(CrystalPresets.magnesium.cell.metricFamily() == .hexagonal)
    }

    // MARK: Lower-symmetry cells

    /// Every crystal system has a ring list: a tetragonal cell (4, 4, 5; one Cu) reads 4/mmm with 58 rings and its {0 0 1} ring
    /// first at 1/c, an orthorhombic (3, 4, 5) mmm with 74, and b = a(1 + 1e-6) 4/mmm (tetragonal along b) with 15. An element
    /// beyond the table (Z = 104) is refused by the list as `unsupportedElements([104])`. Mutants: any cell taken as cubic,
    /// unsupported elements ignored.
    @Test func lowerSymmetryCellsHaveRingsAndAnUnsupportedElementIsRefused() throws {
        let cu = [AtomSite(z: 29, fractional: [0, 0, 0])]
        let tetra = Self.structure(.tetragonal(a: 4, c: 5), cu)
        #expect(tetra.laueClass == .fourOverMmm)
        let tr = try Self.rings(tetra)
        #expect(tr.count == 58)
        let first = try #require(tr.first)
        #expect(abs(first.g - 1.0 / 5) <= 1e-13)
        #expect(first.multiplicity == 2)
        let ortho = Self.structure(.orthorhombic(a: 3, b: 4, c: 5), cu)
        #expect(ortho.laueClass == .mmm)
        #expect(try Self.rings(ortho).count == 74)
        let off = Self.structure(UnitCell(a: Self.auA, b: Self.auA * (1 + 1e-6), c: Self.auA, alphaDeg: 90, betaDeg: 90, gammaDeg: 90),
                                 Self.au.sites)
        #expect(off.laueClass == .fourOverMmm)     // a = c, b different: tetragonal along b
        #expect(try Self.rings(off).count == 15)
        let withZ104 = Self.structure(.cubic(a: Self.auA), Self.au.sites + [AtomSite(z: 104, fractional: [0.5, 0.5, 0.5])])
        #expect(throws: KinematicError.unsupportedElements([104])) { try withZ104.kinematicReflections(kMax: 1.2) }
    }

    // MARK: A hidden swap symmetry

    /// With a = c, cos beta = -1/4 and one atom at the origin the a<->c swap (h,k,l) -> (l,k,h) is a true symmetry: the cell is
    /// mmm, so {1 0 0} is one orbit of 4 and {1 1 0} one of 8; what is still accidental is named by its orbits: {2 0 0} + {2 0 -1}
    /// (8, one |g|) and 27 more of the 83 families. 698 reflections. Mutants: grouping at 1e-4 or 0, the name listing the first
    /// orbit only.
    @Test func swapEquivalentReflectionsAreOneOrbitAndCoincidencesStayNamedByTheirOrbits() throws {
        let mSq = Self.structure(UnitCell(a: 5, b: 4, c: 5, alphaDeg: 90, betaDeg: acos(-0.25) * 180 / .pi, gammaDeg: 90),
                                 [AtomSite(z: 13, fractional: [0, 0, 0])])
        let l = try Self.rings(mSq)
        let a = try #require(Self.family(l, MillerIndex(1, 0, 0)))
        #expect(a.name == "{1 0 0}")
        #expect(a.multiplicity == 4)
        #expect(abs(a.d - 4.8412) <= 5e-5)
        let b = try #require(Self.family(l, MillerIndex(1, 1, 0)))
        #expect(b.name == "{1 1 0}")
        #expect(b.multiplicity == 8)
        let c = try #require(Self.family(l, MillerIndex(2, 0, 0)))
        #expect(c.name == "{2 0 0} + {2 0 -1}")
        #expect(c.multiplicity == 8)
        #expect(abs(c.d - 2.4206) <= 5e-5)
        #expect(l.count == 83)
        #expect(l.filter { $0.orbits.count > 1 }.count == 28)
        #expect(l.map(\.multiplicity).reduce(0, +) == 698)
    }
}
