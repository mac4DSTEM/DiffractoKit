//
//  PointSymmetryTests.swift
//  Planted truth for the groups (registration 5.2): the holohedry orders, closure, determinants and finite orders;
//  the Laue class of each group by the 3.5 rule; orbit sizes at general and special positions; the canonical forms
//  and names the apps use today (3.6). Group theory only; no reference code.
//

import XCTest
import DiffractoCore

final class PointSymmetryTests: XCTestCase {

    static let holohedries: [(UnitCell.MetricFamily, Int, LaueClass)] = [
        (.monoclinic(uniqueAxis: .b), 4, .twoOverM),
        (.monoclinic(uniqueAxis: .a), 4, .twoOverM),
        (.monoclinic(uniqueAxis: .c), 4, .twoOverM),
        (.cubic, 48, .mThreeBarM),
        (.hexagonal, 24, .sixOverMmm),
        (.tetragonal, 16, .fourOverMmm),
        (.rhombohedral, 12, .threeBarM),
        (.orthorhombic, 8, .mmm),
        (.triclinic, 2, .oneBar),
    ]

    /// Every holohedry has its registered order, is closed, every element has determinant ±1 and finite order, and
    /// its Laue class is the holohedry's. Mutant: a generator with a wrong entry (closure runs to 200 or the order moves).
    func testHolohedriesAreTheRegisteredGroups() {
        for (family, order, laue) in Self.holohedries {
            let group = PointSymmetry.holohedry(of: family)
            XCTAssertEqual(group.order, order, "\(family)")
            XCTAssertTrue(group.isClosed, "\(family) is not closed")
            XCTAssertTrue(group.containsInversion, "\(family)")
            for r in group.rotations {
                XCTAssertEqual(abs(r.determinant), 1, "\(family) \(r)")
                XCTAssertNotNil(r.order, "\(family) \(r) has no finite order")
            }
            XCTAssertEqual(group.laueClass, laue, "\(family)")
            XCTAssertEqual(group.properRotations.count, order / 2, "\(family)")
        }
    }

    static let c4z = IntegerMatrix([0, -1, 0], [1, 0, 0], [0, 0, 1])
    static let c6z = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
    static let c3z = c6z.multiplied(by: c6z)
    static let c2x = IntegerMatrix([1, 0, 0], [0, -1, 0], [0, 0, -1])
    static let c2y = IntegerMatrix([-1, 0, 0], [0, 1, 0], [0, 0, -1])
    static let c3diag = IntegerMatrix([0, 0, 1], [1, 0, 0], [0, 1, 0])
    static let c2aHex = IntegerMatrix([1, -1, 0], [0, -1, 0], [0, 0, -1])

    static let fourOverM = PointSymmetry(closureOf: [c4z, .inversion])
    static let threeBar = PointSymmetry(closureOf: [c3z, .inversion])
    static let threeBarM = PointSymmetry(closureOf: [c3z, c2aHex, .inversion])
    static let sixOverM = PointSymmetry(closureOf: [c6z, .inversion])
    static let mThreeBar = PointSymmetry(closureOf: [c2x, c2y, c3diag, .inversion])

    /// The Laue classes the holohedries do not reach: 4/m, -3, -3m (hexagonal axes), 6/m, m-3, each told apart by its
    /// order and its axis orders. Mutant: the 8-with-4-fold branch inverted (4/m and mmm swap).
    func testLaueClassesOfTheSmallerGroups() {
        XCTAssertEqual(Self.fourOverM.order, 8)
        XCTAssertEqual(Self.fourOverM.laueClass, .fourOverM)
        XCTAssertEqual(Self.threeBar.order, 6)
        XCTAssertEqual(Self.threeBar.laueClass, .threeBar)
        XCTAssertEqual(Self.threeBarM.order, 12)
        XCTAssertEqual(Self.threeBarM.laueClass, .threeBarM)
        XCTAssertEqual(Self.sixOverM.order, 12)
        XCTAssertEqual(Self.sixOverM.laueClass, .sixOverM)
        XCTAssertEqual(Self.mThreeBar.order, 24)
        XCTAssertEqual(Self.mThreeBar.laueClass, .mThreeBar)
        // A group without the inversion is no Laue class.
        XCTAssertNil(PointSymmetry(closureOf: [Self.c4z]).laueClass)
        // A set that is not closed is no Laue class either.
        XCTAssertNil(PointSymmetry(rotations: [.identity, .inversion, Self.c4z]).laueClass)
    }

    /// General-position orbit sizes (5.2). Mutant: the orbit as hkl R⁻¹ (the 6-fold's orbit changes size for no hkl but
    /// its members do: covered by the pins), the mirror dropped from 2/m (orbits halve).
    func testGeneralOrbitSizes() {
        let general = MillerIndex(1, 2, 3)
        let cases: [(PointSymmetry, Int, String)] = [
            (PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b)), 4, "2/m"),
            (PointSymmetry.holohedry(of: .triclinic), 2, "-1"),
            (PointSymmetry.holohedry(of: .orthorhombic), 8, "mmm"),
            (Self.fourOverM, 8, "4/m"),
            (PointSymmetry.holohedry(of: .tetragonal), 16, "4/mmm"),
            (Self.threeBar, 6, "-3"),
            (Self.threeBarM, 12, "-3m"),
            (Self.sixOverM, 12, "6/m"),
            (PointSymmetry.holohedry(of: .hexagonal), 24, "6/mmm"),
            (Self.mThreeBar, 24, "m-3"),
            (PointSymmetry.holohedry(of: .cubic), 48, "m-3m"),
        ]
        for (group, size, name) in cases {
            XCTAssertEqual(group.orbit(of: general).count, size, name)
        }
    }

    /// Special positions (5.2): the apps' multiplicities. Mutant: C₂ about a instead of b ({h0l} becomes 4).
    func testSpecialOrbitSizes() {
        let mono = PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))
        XCTAssertEqual(mono.orbit(of: MillerIndex(2, 0, 1)).count, 2)
        XCTAssertEqual(mono.orbit(of: MillerIndex(0, 2, 0)).count, 2)
        XCTAssertEqual(mono.orbit(of: MillerIndex(1, 1, 1)).count, 4)
        XCTAssertEqual(mono.orbit(of: MillerIndex(0, 1, 1)).count, 4)
        XCTAssertEqual(mono.orbit(of: MillerIndex(2, 0, 1)), [MillerIndex(2, 0, 1), MillerIndex(-2, 0, -1)])
        XCTAssertFalse(mono.orbit(of: MillerIndex(2, 0, 1)).contains(MillerIndex(2, 0, -1)), "(2 0 1) and (2 0 -1) are two orbits")

        let cubic = PointSymmetry.holohedry(of: .cubic)
        XCTAssertEqual(cubic.orbit(of: MillerIndex(1, 0, 0)).count, 6)
        XCTAssertEqual(cubic.orbit(of: MillerIndex(1, 1, 0)).count, 12)
        XCTAssertEqual(cubic.orbit(of: MillerIndex(1, 1, 1)).count, 8)
        XCTAssertEqual(cubic.orbit(of: MillerIndex(2, 1, 0)).count, 24)

        let hex = PointSymmetry.holohedry(of: .hexagonal)
        XCTAssertEqual(hex.orbit(of: MillerIndex(1, 0, 0)).count, 6)
        XCTAssertEqual(hex.orbit(of: MillerIndex(0, 0, 2)).count, 2)
        XCTAssertEqual(hex.orbit(of: MillerIndex(1, 0, 1)).count, 12)
        XCTAssertEqual(hex.orbit(of: MillerIndex(1, 1, 0)).count, 6)
        XCTAssertTrue(hex.orbit(of: MillerIndex(1, 1, 0)).contains(MillerIndex(-1, -1, 0)), "(1 1 0) and (-1 -1 0) are one orbit")

        // 4/m keeps {2 1 0} and {1 2 0} apart; 4/mmm joins them. hk0 is a special form (the mirror ⊥ c fixes it):
        // 4 under 4/m and 8 under 4/mmm; the general hkl are 8 and 16 (registration 5.2, amended 2026-10-14).
        XCTAssertEqual(Self.fourOverM.orbit(of: MillerIndex(2, 1, 0)).count, 4)
        XCTAssertEqual(Self.fourOverM.orbit(of: MillerIndex(1, 2, 0)).count, 4)
        XCTAssertFalse(Self.fourOverM.orbit(of: MillerIndex(2, 1, 0)).contains(MillerIndex(1, 2, 0)))
        XCTAssertEqual(Self.fourOverM.orbit(of: MillerIndex(2, 1, 3)).count, 8)
        XCTAssertEqual(PointSymmetry.holohedry(of: .tetragonal).orbit(of: MillerIndex(2, 1, 0)).count, 8)
        XCTAssertEqual(PointSymmetry.holohedry(of: .tetragonal).orbit(of: MillerIndex(2, 1, 3)).count, 16)
    }

    /// The one canonical rule reproduces the apps' names (3.6). Mutant: the key without the non-negative count
    /// (hexagonal {1 1 -2 0} becomes {2 -1 -1 0}); the sign of l not tied to h in 2/m.
    func testCanonicalFormsAndNamesMatchTheApps() {
        let cubic = PointSymmetry.holohedry(of: .cubic)
        XCTAssertEqual(cubic.canonical(MillerIndex(-3, 1, -1)), MillerIndex(3, 1, 1))
        XCTAssertEqual(cubic.orbitName(of: MillerIndex(1, -3, 1)), "{3 1 1}")
        XCTAssertEqual(cubic.orbitName(of: MillerIndex(0, 0, -2)), "{2 0 0}")

        let hex = PointSymmetry.holohedry(of: .hexagonal)
        XCTAssertEqual(hex.canonical(MillerIndex(-1, -1, 0)), MillerIndex(1, 1, 0))
        XCTAssertEqual(hex.orbitName(of: MillerIndex(-1, -1, 0)), "{1 1 -2 0}")
        XCTAssertEqual(hex.orbitName(of: MillerIndex(0, 1, -1)), "{1 0 -1 1}")
        XCTAssertEqual(hex.orbitName(of: MillerIndex(0, 0, -2)), "{0 0 0 2}")
        XCTAssertEqual(hex.orbitName(of: MillerIndex(-1, 1, 0)), "{1 0 -1 0}")
        XCTAssertEqual(hex.orbitName(of: MillerIndex(3, 0, 2)), "{3 0 -3 2}")
        XCTAssertEqual(hex.orbitName(of: MillerIndex(2, 1, 3)), "{2 1 -3 3}")

        let mono = PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))
        XCTAssertEqual(mono.canonical(MillerIndex(-2, 0, 1)), MillerIndex(2, 0, -1))
        XCTAssertEqual(mono.orbitName(of: MillerIndex(-2, 0, 1)), "{2 0 -1}")
        XCTAssertEqual(mono.orbitName(of: MillerIndex(2, 0, 1)), "{2 0 1}")
        XCTAssertEqual(mono.orbitName(of: MillerIndex(0, -2, 0)), "{0 2 0}")
        XCTAssertEqual(mono.orbitName(of: MillerIndex(0, 1, -1)), "{0 1 1}")
        XCTAssertEqual(mono.orbitName(of: MillerIndex(-1, -1, 1)), "{1 1 -1}")

        // 4/m: two orbits, two names.
        XCTAssertEqual(Self.fourOverM.orbitName(of: MillerIndex(-2, 1, 0)), "{1 2 0}")
        XCTAssertEqual(Self.fourOverM.orbitName(of: MillerIndex(2, 1, 0)), "{2 1 0}")
        // -3 in hexagonal axes prints four indices too.
        XCTAssertTrue(Self.threeBar.orbitName(of: MillerIndex(1, 0, 1)).hasPrefix("{"))
        XCTAssertEqual(Self.threeBar.orbitName(of: MillerIndex(1, 0, 1)).split(separator: " ").count, 4)
    }

    /// The canonical form is constant on an orbit and differs between orbits (bijective), on every class.
    func testCanonicalIsAnOrbitInvariant() {
        let groups = [PointSymmetry.holohedry(of: .cubic), PointSymmetry.holohedry(of: .hexagonal), Self.fourOverM,
                      Self.threeBar, Self.sixOverM, Self.mThreeBar, PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))]
        for group in groups {
            var seen: [MillerIndex: Set<MillerIndex>] = [:]
            for h in -2...2 { for k in -2...2 { for l in -2...2 where !(h == 0 && k == 0 && l == 0) {
                let hkl = MillerIndex(h, k, l)
                let orbit = group.orbit(of: hkl)
                let canonical = group.canonical(hkl)
                XCTAssertTrue(orbit.contains(canonical))
                for member in orbit { XCTAssertEqual(group.canonical(member), canonical) }
                if let other = seen[canonical] { XCTAssertEqual(other, orbit) } else { seen[canonical] = orbit }
            } } }
        }
    }

    /// The Laue group of C2/m's eight operations (the β″ CIF's loop) is 2/m, order 4.
    func testLaueGroupOfOperations() throws {
        let ops = try ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z",
                       "1/2+x,1/2-y,z"].map { try SymmetryOperation(xyz: $0) }
        let laue = PointSymmetry.laueGroup(of: ops)
        XCTAssertEqual(laue.order, 4)
        XCTAssertEqual(laue.laueClass, .twoOverM)
        XCTAssertEqual(laue, PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b)))
        // From a loop without the inversion (P2's two operations) the Laue group still has it.
        let p2 = try ["x,y,z", "-x,y,-z"].map { try SymmetryOperation(xyz: $0) }
        XCTAssertEqual(PointSymmetry.laueGroup(of: p2).laueClass, .twoOverM)
    }

    func testLaueClassTable() {
        XCTAssertEqual(LaueClass(pointGroupSymbol: "32"), .threeBarM)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "-42m"), .fourOverMmm)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "m-3"), .mThreeBar)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "6/m"), .sixOverM)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "2"), .twoOverM)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "m"), .twoOverM)
        XCTAssertEqual(LaueClass(pointGroupSymbol: "1"), .oneBar)
        XCTAssertNil(LaueClass(pointGroupSymbol: "5"))
        for laue in LaueClass.allCases {
            XCTAssertEqual(LaueClass(pointGroupSymbol: laue.symbol), laue)
            XCTAssertEqual(LaueClass(rawValue: laue.symbol), laue)
        }
        XCTAssertEqual(LaueClass.allCases.count, 11)
        XCTAssertEqual(Set(LaueClass.allCases.map(\.crystalSystem)).count, 7)
    }
}
