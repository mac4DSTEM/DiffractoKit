//
//  PointSymmetryTests.swift
//  Planted truth for the groups (registration 5.2): the holohedry orders, closure, determinants and finite orders;
//  the Laue class of each group by the 3.5 rule; orbit sizes at general and special positions; the canonical forms
//  and names the reference implementations use today (3.6). Group theory only; no reference code.
//

import DiffractoKitTestSupport
import Testing
import DiffractoCore

struct PointSymmetryTests {

    static let holohedries: [(UnitCell.MetricFamily, Int, LaueClass)] = [
        (.monoclinic(uniqueAxis: .b), 4, .twoOverM),
        (.monoclinic(uniqueAxis: .a), 4, .twoOverM),
        (.monoclinic(uniqueAxis: .c), 4, .twoOverM),
        (.cubic, 48, .mThreeBarM),
        (.hexagonal, 24, .sixOverMmm),
        (.tetragonal(uniqueAxis: .c), 16, .fourOverMmm),
        (.tetragonal(uniqueAxis: .a), 16, .fourOverMmm),
        (.tetragonal(uniqueAxis: .b), 16, .fourOverMmm),
        (.rhombohedral, 12, .threeBarM),
        (.orthorhombic, 8, .mmm),
        (.triclinic, 2, .oneBar),
    ]

    /// Every holohedry has its registered order, is closed, every element has determinant ±1 and finite order, and
    /// its Laue class is the holohedry's. Mutant: a generator with a wrong entry (closure runs to 200 or the order moves).
    @Test func testHolohedriesAreTheRegisteredGroups() {
        for (family, order, laue) in Self.holohedries {
            let group = PointSymmetry.holohedry(of: family)
            #expect(group.order == order, "\(family)")
            #expect(group.isClosed, "\(family) is not closed")
            #expect(group.containsInversion, "\(family)")
            for r in group.rotations {
                #expect(abs(r.determinant) == 1, "\(family) \(r)")
                #expect(r.order != nil, "\(family) \(r) has no finite order")
            }
            #expect(group.laueClass == laue, "\(family)")
            #expect(group.properRotations.count == order / 2, "\(family)")
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
    @Test func testLaueClassesOfTheSmallerGroups() {
        #expect(Self.fourOverM.order == 8)
        #expect(Self.fourOverM.laueClass == .fourOverM)
        #expect(Self.threeBar.order == 6)
        #expect(Self.threeBar.laueClass == .threeBar)
        #expect(Self.threeBarM.order == 12)
        #expect(Self.threeBarM.laueClass == .threeBarM)
        #expect(Self.sixOverM.order == 12)
        #expect(Self.sixOverM.laueClass == .sixOverM)
        #expect(Self.mThreeBar.order == 24)
        #expect(Self.mThreeBar.laueClass == .mThreeBar)
        // A group without the inversion is no Laue class.
        #expect(PointSymmetry(closureOf: [Self.c4z]).laueClass == nil)
        // A set that is not closed is no Laue class either.
        #expect(PointSymmetry(rotations: [.identity, .inversion, Self.c4z]).laueClass == nil)
    }

    /// General-position orbit sizes (5.2). Mutant: the orbit as hkl R⁻¹ (the 6-fold's orbit changes size for no hkl but
    /// its members do: covered by the pins), the mirror dropped from 2/m (orbits halve).
    @Test func testGeneralOrbitSizes() {
        let general = MillerIndex(1, 2, 3)
        let cases: [(PointSymmetry, Int, String)] = [
            (PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b)), 4, "2/m"),
            (PointSymmetry.holohedry(of: .triclinic), 2, "-1"),
            (PointSymmetry.holohedry(of: .orthorhombic), 8, "mmm"),
            (Self.fourOverM, 8, "4/m"),
            (PointSymmetry.holohedry(of: .tetragonal(uniqueAxis: .c)), 16, "4/mmm"),
            (Self.threeBar, 6, "-3"),
            (Self.threeBarM, 12, "-3m"),
            (Self.sixOverM, 12, "6/m"),
            (PointSymmetry.holohedry(of: .hexagonal), 24, "6/mmm"),
            (Self.mThreeBar, 24, "m-3"),
            (PointSymmetry.holohedry(of: .cubic), 48, "m-3m"),
        ]
        for (group, size, name) in cases {
            #expect(group.orbit(of: general).count == size, "\(name)")
        }
    }

    /// Special positions (5.2): the reference implementations' multiplicities. Mutant: C₂ about a instead of b ({h0l} becomes 4).
    @Test func testSpecialOrbitSizes() {
        let mono = PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))
        #expect(mono.orbit(of: MillerIndex(2, 0, 1)).count == 2)
        #expect(mono.orbit(of: MillerIndex(0, 2, 0)).count == 2)
        #expect(mono.orbit(of: MillerIndex(1, 1, 1)).count == 4)
        #expect(mono.orbit(of: MillerIndex(0, 1, 1)).count == 4)
        #expect(mono.orbit(of: MillerIndex(2, 0, 1)) == [MillerIndex(2, 0, 1), MillerIndex(-2, 0, -1)])
        #expect(!(mono.orbit(of: MillerIndex(2, 0, 1)).contains(MillerIndex(2, 0, -1))), "(2 0 1) and (2 0 -1) are two orbits")

        let cubic = PointSymmetry.holohedry(of: .cubic)
        #expect(cubic.orbit(of: MillerIndex(1, 0, 0)).count == 6)
        #expect(cubic.orbit(of: MillerIndex(1, 1, 0)).count == 12)
        #expect(cubic.orbit(of: MillerIndex(1, 1, 1)).count == 8)
        #expect(cubic.orbit(of: MillerIndex(2, 1, 0)).count == 24)

        let hex = PointSymmetry.holohedry(of: .hexagonal)
        #expect(hex.orbit(of: MillerIndex(1, 0, 0)).count == 6)
        #expect(hex.orbit(of: MillerIndex(0, 0, 2)).count == 2)
        #expect(hex.orbit(of: MillerIndex(1, 0, 1)).count == 12)
        #expect(hex.orbit(of: MillerIndex(1, 1, 0)).count == 6)
        #expect(hex.orbit(of: MillerIndex(1, 1, 0)).contains(MillerIndex(-1, -1, 0)), "(1 1 0) and (-1 -1 0) are one orbit")

        // 4/m keeps {2 1 0} and {1 2 0} apart; 4/mmm joins them. hk0 is a special form (the mirror ⊥ c fixes it):
        // 4 under 4/m and 8 under 4/mmm; the general hkl are 8 and 16 (registration 5.2, amended 2026-10-14).
        #expect(Self.fourOverM.orbit(of: MillerIndex(2, 1, 0)).count == 4)
        #expect(Self.fourOverM.orbit(of: MillerIndex(1, 2, 0)).count == 4)
        #expect(!(Self.fourOverM.orbit(of: MillerIndex(2, 1, 0)).contains(MillerIndex(1, 2, 0))))
        #expect(Self.fourOverM.orbit(of: MillerIndex(2, 1, 3)).count == 8)
        #expect(PointSymmetry.holohedry(of: .tetragonal(uniqueAxis: .c)).orbit(of: MillerIndex(2, 1, 0)).count == 8)
        #expect(PointSymmetry.holohedry(of: .tetragonal(uniqueAxis: .c)).orbit(of: MillerIndex(2, 1, 3)).count == 16)
    }

    /// The one canonical rule reproduces the reference implementations' names (3.6). Mutant: the key without the non-negative count
    /// (hexagonal {1 1 -2 0} becomes {2 -1 -1 0}); the sign of l not tied to h in 2/m.
    @Test func testCanonicalFormsAndNamesMatchTheApps() {
        let cubic = PointSymmetry.holohedry(of: .cubic)
        #expect(cubic.canonical(MillerIndex(-3, 1, -1)) == MillerIndex(3, 1, 1))
        #expect(cubic.orbitName(of: MillerIndex(1, -3, 1)) == "{3 1 1}")
        #expect(cubic.orbitName(of: MillerIndex(0, 0, -2)) == "{2 0 0}")

        let hex = PointSymmetry.holohedry(of: .hexagonal)
        #expect(hex.canonical(MillerIndex(-1, -1, 0)) == MillerIndex(1, 1, 0))
        #expect(hex.orbitName(of: MillerIndex(-1, -1, 0)) == "{1 1 -2 0}")
        #expect(hex.orbitName(of: MillerIndex(0, 1, -1)) == "{1 0 -1 1}")
        #expect(hex.orbitName(of: MillerIndex(0, 0, -2)) == "{0 0 0 2}")
        #expect(hex.orbitName(of: MillerIndex(-1, 1, 0)) == "{1 0 -1 0}")
        #expect(hex.orbitName(of: MillerIndex(3, 0, 2)) == "{3 0 -3 2}")
        #expect(hex.orbitName(of: MillerIndex(2, 1, 3)) == "{2 1 -3 3}")

        let mono = PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))
        #expect(mono.canonical(MillerIndex(-2, 0, 1)) == MillerIndex(2, 0, -1))
        #expect(mono.orbitName(of: MillerIndex(-2, 0, 1)) == "{2 0 -1}")
        #expect(mono.orbitName(of: MillerIndex(2, 0, 1)) == "{2 0 1}")
        #expect(mono.orbitName(of: MillerIndex(0, -2, 0)) == "{0 2 0}")
        #expect(mono.orbitName(of: MillerIndex(0, 1, -1)) == "{0 1 1}")
        #expect(mono.orbitName(of: MillerIndex(-1, -1, 1)) == "{1 1 -1}")

        // 4/m: two orbits, two names.
        #expect(Self.fourOverM.orbitName(of: MillerIndex(-2, 1, 0)) == "{1 2 0}")
        #expect(Self.fourOverM.orbitName(of: MillerIndex(2, 1, 0)) == "{2 1 0}")
        // -3 in hexagonal axes prints four indices too.
        #expect(Self.threeBar.orbitName(of: MillerIndex(1, 0, 1)).hasPrefix("{"))
        #expect(Self.threeBar.orbitName(of: MillerIndex(1, 0, 1)).split(separator: " ").count == 4)
    }

    /// The canonical form is constant on an orbit and differs between orbits (bijective), on every class.
    @Test func testCanonicalIsAnOrbitInvariant() {
        let groups = [PointSymmetry.holohedry(of: .cubic), PointSymmetry.holohedry(of: .hexagonal), Self.fourOverM,
                      Self.threeBar, Self.sixOverM, Self.mThreeBar, PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b))]
        for group in groups {
            var seen: [MillerIndex: Set<MillerIndex>] = [:]
            for h in -2...2 { for k in -2...2 { for l in -2...2 where !(h == 0 && k == 0 && l == 0) {
                let hkl = MillerIndex(h, k, l)
                let orbit = group.orbit(of: hkl)
                let canonical = group.canonical(hkl)
                #expect(orbit.contains(canonical))
                for member in orbit { #expect(group.canonical(member) == canonical) }
                if let other = seen[canonical] { #expect(other == orbit) } else { seen[canonical] = orbit }
            } } }
        }
    }

    /// The Laue group of C2/m's eight operations (the β″ CIF's loop) is 2/m, order 4.
    @Test func testLaueGroupOfOperations() throws {
        let ops = try ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z",
                       "1/2+x,1/2-y,z"].map { try SymmetryOperation(xyz: $0) }
        let laue = PointSymmetry.laueGroup(of: ops)
        #expect(laue.order == 4)
        #expect(laue.laueClass == .twoOverM)
        #expect(laue == PointSymmetry.holohedry(of: .monoclinic(uniqueAxis: .b)))
        // From a loop without the inversion (P2's two operations) the Laue group still has it.
        let p2 = try ["x,y,z", "-x,y,-z"].map { try SymmetryOperation(xyz: $0) }
        #expect(PointSymmetry.laueGroup(of: p2).laueClass == .twoOverM)
    }

    @Test func testLaueClassTable() {
        #expect(LaueClass(pointGroupSymbol: "32") == .threeBarM)
        #expect(LaueClass(pointGroupSymbol: "-42m") == .fourOverMmm)
        #expect(LaueClass(pointGroupSymbol: "m-3") == .mThreeBar)
        #expect(LaueClass(pointGroupSymbol: "6/m") == .sixOverM)
        #expect(LaueClass(pointGroupSymbol: "2") == .twoOverM)
        #expect(LaueClass(pointGroupSymbol: "m") == .twoOverM)
        #expect(LaueClass(pointGroupSymbol: "1") == .oneBar)
        #expect(LaueClass(pointGroupSymbol: "5") == nil)
        for laue in LaueClass.allCases {
            #expect(LaueClass(pointGroupSymbol: laue.symbol) == laue)
            #expect(LaueClass(rawValue: laue.symbol) == laue)
        }
        #expect(LaueClass.allCases.count == 11)
        #expect(Set(LaueClass.allCases.map(\.crystalSystem)).count == 7)
    }
}
