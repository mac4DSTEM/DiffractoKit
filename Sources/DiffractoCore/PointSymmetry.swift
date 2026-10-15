//
//  PointSymmetry.swift
//  Role: a finite group of integer rotation matrices in the fractional basis: the holohedry of a metric, the Laue
//        group of a structure, the orbit of a reflection and its canonical form and name. Registration
//        docs/registration-slice-1.md sections 3.4 (the holohedry generators), 3.5 (the Laue class of a group) and 3.6
//        (orbits: (h k l) R; the canonical form is the orbit member greatest under (number of non-negative indices,
//        lexicographic), the one rule that reproduces the reference implementations' cubic, hexagonal and monoclinic names).
//

import Foundation

public struct PointSymmetry: Sendable, Equatable {
    /// The group's elements, deduplicated, in a deterministic order (by their nine integers).
    public let rotations: [IntegerMatrix]

    /// The elements as given, deduplicated and ordered; closure is not enforced (`isClosed` says).
    public init(rotations: [IntegerMatrix]) {
        self.rotations = PointSymmetry.ordered(Set(rotations))
    }

    /// The closure of `generators` under multiplication (the identity included). Stops at 200 elements, which no
    /// crystallographic group reaches, so a generator of infinite order cannot run away.
    public init(closureOf generators: [IntegerMatrix]) {
        var elements: Set<IntegerMatrix> = [.identity]
        var frontier: [IntegerMatrix] = [.identity]
        while let current = frontier.popLast(), elements.count < 200 {
            for g in generators {
                let product = current.multiplied(by: g)
                if elements.insert(product).inserted { frontier.append(product) }
            }
        }
        self.rotations = PointSymmetry.ordered(elements)
    }

    private static func ordered(_ set: Set<IntegerMatrix>) -> [IntegerMatrix] {
        set.sorted { lhs, rhs in
            let l = lhs.rows.flatMap { [$0.x, $0.y, $0.z] }, r = rhs.rows.flatMap { [$0.x, $0.y, $0.z] }
            return l.lexicographicallyPrecedes(r)
        }
    }

    // MARK: Named groups

    /// The rotation parts of `operations`, closed under inversion: {±R}.
    public static func laueGroup(of operations: [SymmetryOperation]) -> PointSymmetry {
        laueGroup(ofRotations: operations.map(\.rotation))
    }

    public static func laueGroup(ofRotations rotations: [IntegerMatrix]) -> PointSymmetry {
        PointSymmetry(rotations: rotations + rotations.map(\.negated))
    }

    /// The holohedry (the full Laue group of the lattice) of a metric family (section 3.4): the conventional setting, with
    /// the 4-fold of a tetragonal metric about its odd axis and the unique axis of a monoclinic one where the metric puts it.
    public static func holohedry(of family: UnitCell.MetricFamily) -> PointSymmetry {
        let c2x = IntegerMatrix([1, 0, 0], [0, -1, 0], [0, 0, -1])
        let c2y = IntegerMatrix([-1, 0, 0], [0, 1, 0], [0, 0, -1])
        let c2z = IntegerMatrix([-1, 0, 0], [0, -1, 0], [0, 0, 1])
        let c4z = IntegerMatrix([0, -1, 0], [1, 0, 0], [0, 0, 1])
        let c3diag = IntegerMatrix([0, 0, 1], [1, 0, 0], [0, 1, 0])
        let c6z = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        // The 4-fold along a and along b (permuted c4z; the 2-fold perpendicular to it is c2y / c2z).
        let c4x = IntegerMatrix([1, 0, 0], [0, 0, -1], [0, 1, 0])
        let c4y = IntegerMatrix([0, 0, 1], [0, 1, 0], [-1, 0, 0])
        let c2aHex = IntegerMatrix([1, -1, 0], [0, -1, 0], [0, 0, -1])
        let c2rhombo = IntegerMatrix([0, -1, 0], [-1, 0, 0], [0, 0, -1])
        let generators: [IntegerMatrix]
        switch family {
        case .cubic: generators = [c4z, c3diag, .inversion]
        case .hexagonal: generators = [c6z, c2aHex, .inversion]
        case .tetragonal(let axis):
            switch axis {
            case .a: generators = [c4x, c2y, .inversion]
            case .b: generators = [c4y, c2z, .inversion]
            case .c: generators = [c4z, c2x, .inversion]
            }
        case .rhombohedral: generators = [c3diag, c2rhombo, .inversion]
        case .orthorhombic: generators = [c2x, c2y, .inversion]
        case .monoclinic(let axis):
            switch axis {
            case .a: generators = [c2x, .inversion]
            case .b: generators = [c2y, .inversion]
            case .c: generators = [c2z, .inversion]
            }
        case .triclinic: generators = [.inversion]
        }
        return PointSymmetry(closureOf: generators)
    }

    // MARK: Group properties

    public var order: Int { rotations.count }

    public func contains(_ r: IntegerMatrix) -> Bool { rotations.contains(r) }

    public var containsInversion: Bool { contains(.inversion) }

    /// Every product of two elements is an element.
    public var isClosed: Bool {
        let set = Set(rotations)
        for a in rotations {
            for b in rotations where !set.contains(a.multiplied(by: b)) { return false }
        }
        return true
    }

    /// The proper rotations (determinant +1).
    public var properRotations: [IntegerMatrix] { rotations.filter { $0.determinant == 1 } }

    /// The Laue class by the order and the axis orders present (section 3.5); nil when the group is not closed, lacks
    /// the inversion, or has an order no Laue group has.
    public var laueClass: LaueClass? {
        guard containsInversion, isClosed else { return nil }
        let orders = Set(properRotations.compactMap(\.order))
        let has4 = orders.contains(4), has6 = orders.contains(6)
        switch order {
        case 2: return .oneBar
        case 4: return .twoOverM
        case 8: return has4 ? .fourOverM : .mmm
        case 16: return .fourOverMmm
        case 6: return .threeBar
        case 12: return has6 ? .sixOverM : .threeBarM
        case 24: return has6 ? .sixOverMmm : .mThreeBar
        case 48: return .mThreeBarM
        default: return nil
        }
    }

    // MARK: Orbits (section 3.6)

    /// {(h k l) R : R in the group}.
    public func orbit(of hkl: MillerIndex) -> Set<MillerIndex> {
        Set(rotations.map { hkl.applying($0) })
    }

    /// The orbit member greatest under (number of non-negative indices, then lexicographic (h, k, l)).
    public func canonical(_ hkl: MillerIndex) -> MillerIndex {
        PointSymmetry.canonical(of: orbit(of: hkl))
    }

    public static func canonical(of orbit: Set<MillerIndex>) -> MillerIndex {
        orbit.max { lhs, rhs in
            (lhs.nonNegativeCount, lhs.h, lhs.k, lhs.l) < (rhs.nonNegativeCount, rhs.h, rhs.k, rhs.l)
        } ?? .zero
    }

    /// "{h k l}", or "{h k i l}" with i = −(h + k) when `millerBravais`.
    public static func orbitName(ofCanonical hkl: MillerIndex, millerBravais: Bool) -> String {
        millerBravais ? "{\(hkl.h) \(hkl.k) \(hkl.millerBravaisI) \(hkl.l)}" : "{\(hkl.h) \(hkl.k) \(hkl.l)}"
    }

    /// The name of `hkl`'s orbit, in Miller–Bravais form for the hexagonal-axes Laue classes.
    public func orbitName(of hkl: MillerIndex) -> String {
        PointSymmetry.orbitName(ofCanonical: canonical(hkl), millerBravais: laueClass?.usesMillerBravais ?? false)
    }
}
