//
//  MillerIndex.swift
//  Role: (h, k, l), the integer triple that names a reflection or a lattice plane, and the integer 3 × 3 matrices the
//        symmetry operations rotate it with. The covariant action is the one that matters here (registration section
//        3.6): a rotation R acting on fractional coordinates as p' = R p sends the reflection (h k l) to the row
//        vector (h k l) R, so |F| is preserved along that orbit and a reflection fixed by R is absent unless
//        (h k l) · t is an integer (section 3.8).
//

import Foundation

/// (h, k, l): a reflection, a plane, or an orbit's canonical form. Lexicographic order on (h, k, l).
public struct MillerIndex: Sendable, Hashable, Comparable, CustomStringConvertible {
    public let h: Int
    public let k: Int
    public let l: Int

    public init(_ h: Int, _ k: Int, _ l: Int) {
        self.h = h; self.k = k; self.l = l
    }

    public static let zero = MillerIndex(0, 0, 0)

    public var isZero: Bool { h == 0 && k == 0 && l == 0 }

    /// The Friedel mate (−h, −k, −l).
    public var negated: MillerIndex { MillerIndex(-h, -k, -l) }

    /// Lexicographic on (h, k, l).
    public static func < (a: MillerIndex, b: MillerIndex) -> Bool {
        (a.h, a.k, a.l) < (b.h, b.k, b.l)
    }

    /// The row vector (h k l) times `rotation`: where a symmetry operation with that rotation part sends this
    /// reflection (registration section 3.6).
    public func applying(_ rotation: IntegerMatrix) -> MillerIndex {
        MillerIndex(
            h * rotation.rows[0].x + k * rotation.rows[1].x + l * rotation.rows[2].x,
            h * rotation.rows[0].y + k * rotation.rows[1].y + l * rotation.rows[2].y,
            h * rotation.rows[0].z + k * rotation.rows[1].z + l * rotation.rows[2].z)
    }

    /// The fourth Miller–Bravais index i = −(h + k), for the hexagonal-axes Laue classes.
    public var millerBravaisI: Int { -(h + k) }

    /// h x + k y + l z.
    public func dot(_ v: SIMD3<Double>) -> Double {
        Double(h) * v.x + Double(k) * v.y + Double(l) * v.z
    }

    /// Of the orbit-canonical key (registration 3.6): how many indices are non-negative.
    public var nonNegativeCount: Int {
        (h >= 0 ? 1 : 0) + (k >= 0 ? 1 : 0) + (l >= 0 ? 1 : 0)
    }

    public var description: String { "(\(h) \(k) \(l))" }
}

/// A 3 × 3 integer matrix: the rotation part of a symmetry operation in the fractional basis, where every operation of
/// a space group in any setting is an integer matrix of determinant ±1 and finite order.
public struct IntegerMatrix: Sendable, Hashable, CustomStringConvertible {
    /// The three rows.
    public let rows: [SIMD3<Int>]

    public init(rows: [SIMD3<Int>]) {
        precondition(rows.count == 3, "an IntegerMatrix has three rows")
        self.rows = rows
    }

    public init(_ r0: SIMD3<Int>, _ r1: SIMD3<Int>, _ r2: SIMD3<Int>) {
        self.rows = [r0, r1, r2]
    }

    public static let identity = IntegerMatrix([1, 0, 0], [0, 1, 0], [0, 0, 1])
    public static let inversion = IntegerMatrix([-1, 0, 0], [0, -1, 0], [0, 0, -1])

    public subscript(row: Int, column: Int) -> Int {
        rows[row][column]
    }

    /// Matrix product (self · other), the composition that applies `other` first to a column vector.
    public func multiplied(by other: IntegerMatrix) -> IntegerMatrix {
        func row(_ r: SIMD3<Int>) -> SIMD3<Int> {
            SIMD3(r.x * other.rows[0].x + r.y * other.rows[1].x + r.z * other.rows[2].x,
                  r.x * other.rows[0].y + r.y * other.rows[1].y + r.z * other.rows[2].y,
                  r.x * other.rows[0].z + r.y * other.rows[1].z + r.z * other.rows[2].z)
        }
        return IntegerMatrix(rows: rows.map(row))
    }

    public var negated: IntegerMatrix { IntegerMatrix(rows: rows.map { SIMD3(-$0.x, -$0.y, -$0.z) }) }

    public var transposed: IntegerMatrix {
        IntegerMatrix(SIMD3(rows[0].x, rows[1].x, rows[2].x),
                      SIMD3(rows[0].y, rows[1].y, rows[2].y),
                      SIMD3(rows[0].z, rows[1].z, rows[2].z))
    }

    /// The determinant. Computed in wrapping arithmetic: it never traps, and it is the exact determinant whenever that fits
    /// in an `Int` (the value mod 2⁶⁴ is the value), as it does for every unimodular matrix however large its entries.
    public var determinant: Int {
        let m = rows
        let a: Int = m[1].y &* m[2].z &- m[1].z &* m[2].y
        let b: Int = m[1].x &* m[2].z &- m[1].z &* m[2].x
        let c: Int = m[1].x &* m[2].y &- m[1].y &* m[2].x
        return m[0].x &* a &- m[0].y &* b &+ m[0].z &* c
    }

    public var trace: Int { rows[0].x + rows[1].y + rows[2].z }

    public var isIdentity: Bool { self == .identity }

    /// The smallest n ≥ 1 with Rⁿ = I, or nil when none is found within 12 (a crystallographic rotation has order 1,
    /// 2, 3, 4 or 6; an improper one up to 6 as well, since (−I)² = I).
    public var order: Int? {
        var power = self
        for n in 1...12 {
            if power.isIdentity { return n }
            power = power.multiplied(by: self)
        }
        return nil
    }

    /// Applies the matrix to a column vector of doubles (p' = R p).
    public func apply(_ p: SIMD3<Double>) -> SIMD3<Double> {
        SIMD3(
            Double(rows[0].x) * p.x + Double(rows[0].y) * p.y + Double(rows[0].z) * p.z,
            Double(rows[1].x) * p.x + Double(rows[1].y) * p.y + Double(rows[1].z) * p.z,
            Double(rows[2].x) * p.x + Double(rows[2].y) * p.y + Double(rows[2].z) * p.z)
    }

    public var description: String {
        rows.map { "[\($0.x) \($0.y) \($0.z)]" }.joined(separator: " ")
    }
}
