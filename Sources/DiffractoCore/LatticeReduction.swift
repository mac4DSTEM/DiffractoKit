//
//  LatticeReduction.swift
//  Role: the Niggli-reduced cell of a lattice and the integer change of basis that reaches it (registration
//        docs/registration-slice-1b.md section 1). The algorithm is Křivý & Gruber (1976), Acta Cryst. A32, 297,
//        steps A1 to A8, with Grosse-Kunstleve, Sauter & Adams (2004), Acta Cryst. A60, 1-6: every comparison carries
//        an epsilon (x < y becomes x < y - ε, x = y becomes |x - y| ≤ ε), so a cell within ε of a boundary between
//        two reduced forms reduces to one of them and never cycles; and the metric is recomputed from the exact input
//        metric and the integer matrix at every step, so rounding does not accumulate.
//
//        ε = `relativeEpsilon` · V^(2/3) (Å², the unit of the metric entries); the default relative value is 1e-5,
//        the value spglib (`niggli_reduce(eps=1e-5)`) and pymatgen (`get_niggli_reduced_lattice(tol=1e-5)`) use.
//

import Foundation
import simd

/// A Niggli-reduced cell and the integer matrix that reaches it.
public struct ReducedCell: Sendable, Equatable {
    /// The reduced cell (a ≤ b ≤ c up to ε, the main conditions of the Niggli form).
    public let cell: UnitCell
    /// P with determinant +1: the reduced basis vectors are the COLUMNS of P in terms of the input basis, so
    /// x_input = P x_reduced and the reduced metric is Pᵀ G P.
    public let transformation: IntegerMatrix
    /// ε in Å² that the comparisons used.
    public let epsilon: Double
    /// Whether the algorithm stopped on its own conditions (false only on an invalid cell, a shear beyond 2⁴⁰ lattice
    /// units, or after 2000 steps; rounding in the six parameters can also keep a cell sheared by 1e5 or more lattice
    /// units from settling). `ConventionalCell.reductionConverged` says what the conventional route does then.
    public let converged: Bool
    /// The number of steps (each one swap, sign change or shear) the reduction took.
    public let steps: Int
}

extension UnitCell {
    /// The default relative ε of `niggliReduced` (spglib's and pymatgen's 1e-5).
    public static let niggliRelativeEpsilon = 1e-5

    /// The Niggli-reduced cell of this lattice (the cell is taken as a PRIMITIVE lattice: centring in a file's atoms
    /// is not a property of the six parameters) and the matrix P that reaches it.
    ///
    /// A `relativeEpsilon` that is NaN, negative or infinite is refused the way an invalid cell is: the cell comes back
    /// unreduced (identity transformation, epsilon 0) with `converged == false`, never as reduced. A finite epsilon so
    /// large that no comparison can fail also stops at once and reports what it reached (a cell reduced only to within
    /// that epsilon): keep it to the order of 1e-5.
    public func niggliReduced(relativeEpsilon: Double = UnitCell.niggliRelativeEpsilon) -> ReducedCell {
        guard isValid, relativeEpsilon.isFinite, relativeEpsilon >= 0 else {
            return ReducedCell(cell: self, transformation: .identity, epsilon: 0, converged: false, steps: 0)
        }
        let g = metricTensor
        let eps = relativeEpsilon * pow(volume, 2.0 / 3.0)
        var cols: [SIMD3<Int>] = [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, 0, 1)]
        func q(_ i: Int, _ j: Int) -> Double { UnitCell.bilinear(g, cols[i], cols[j]) }
        func sign(_ x: Double) -> Int { x > eps ? 1 : (x < -eps ? -1 : 0) }
        var steps = 0
        var converged = false
        while steps < 2000 {
            steps += 1
            let a = q(0, 0), b = q(1, 1), c = q(2, 2)
            let xi = 2 * q(1, 2), eta = 2 * q(0, 2), zeta = 2 * q(0, 1)
            // A1, A2: order the lengths (ties by |ξ| ≤ |η| ≤ |ζ|)
            if a > b + eps || (abs(a - b) <= eps && abs(xi) > abs(eta) + eps) { cols.swapAt(0, 1); continue }
            if b > c + eps || (abs(b - c) <= eps && abs(eta) > abs(zeta) + eps) { cols.swapAt(1, 2); continue }
            // A3, A4: all of ξ, η, ζ positive (type I) or all non-positive (type II). Scaling the columns by
            // (1, e₃, e₂) multiplies (ξ, η, ζ) by (e₁, e₂, e₃) with e₁e₂e₃ = 1.
            let s = [sign(xi), sign(eta), sign(zeta)]
            var e = [1, 1, 1]
            if s[0] * s[1] * s[2] == 1 {
                e = s
            } else {
                for i in 0..<3 { e[i] = s[i] == 0 ? 1 : -s[i] }
                if e[0] * e[1] * e[2] == -1, let z = s.firstIndex(of: 0) { e[z] = -1 }
            }
            if e != [1, 1, 1] { cols[1] = cols[1] &* e[2]; cols[2] = cols[2] &* e[1]; continue }
            // A5, A6, A7: a shear that makes |ξ| ≤ B, |η| ≤ A, |ζ| ≤ A (with the boundary cases of the paper). The shear is
            // by an integer MULTIPLE when the excess is large (`shearMultiple`: 1 up to 50 lengths, the paper's unit
            // shear, so a cell that needs a few dozen shears reduces exactly as before), so a cell sheared by thousands of
            // lattice units reduces in a few steps and not after thousands (Gate B on 4876b18, F1).
            if abs(xi) > b + eps || (abs(xi - b) <= eps && 2 * eta < zeta - eps) || (abs(xi + b) <= eps && zeta < -eps) {
                let n: Int = UnitCell.shearMultiple(abs(xi), b, eps)
                if n < 0 { break }
                guard let sheared = UnitCell.sheared(cols[2], by: xi >= 0 ? n : -n, along: cols[1]) else { break }
                cols[2] = sheared; continue
            }
            if abs(eta) > a + eps || (abs(eta - a) <= eps && 2 * xi < zeta - eps) || (abs(eta + a) <= eps && zeta < -eps) {
                let n: Int = UnitCell.shearMultiple(abs(eta), a, eps)
                if n < 0 { break }
                guard let sheared = UnitCell.sheared(cols[2], by: eta >= 0 ? n : -n, along: cols[0]) else { break }
                cols[2] = sheared; continue
            }
            if abs(zeta) > a + eps || (abs(zeta - a) <= eps && 2 * xi < eta - eps) || (abs(zeta + a) <= eps && eta < -eps) {
                let n: Int = UnitCell.shearMultiple(abs(zeta), a, eps)
                if n < 0 { break }
                guard let sheared = UnitCell.sheared(cols[1], by: zeta >= 0 ? n : -n, along: cols[0]) else { break }
                cols[1] = sheared; continue
            }
            // A8: c → a + b + c when that shortens it
            let sum = xi + eta + zeta + a + b
            if sum < -eps || (abs(sum) <= eps && 2 * (a + eta) + zeta > eps) {
                guard let summed = UnitCell.checkedSum(cols[0], cols[1], cols[2]) else { break }
                cols[2] = summed; continue
            }
            converged = true
            break
        }
        var p = IntegerMatrix(columns: cols)
        if p.determinant < 0 { p = p.negated }          // −I keeps the metric and makes the basis right-handed
        let reduced = UnitCell.cell(fromMetric: UnitCell.gram(of: p, g))
        return ReducedCell(cell: reduced, transformation: p, epsilon: eps, converged: converged, steps: steps)
    }

    // MARK: Helpers shared with the conventional cell

    /// How many times to shear in one step when the cross term `x` (ξ, η or ζ, in absolute value) exceeds the length
    /// `len` (B or A) that bounds it: 1 while |x| ≤ 50·len (the paper's unit shear: every cell of the 2827 pinned ones
    /// reduces by unit shears alone, bit for bit as before), else the number that leaves
    /// |x| − 2·len·n in [len, 3·len), so the unit shears after it finish the job. Negative (stop, not converged) when
    /// the multiple is not a number an Int basis can hold (a shear beyond 2⁴⁰, where the metric's rounding is larger
    /// than ε anyway).
    static func shearMultiple(_ x: Double, _ len: Double, _ eps: Double) -> Int {
        if !(x <= 50 * len + eps) {
            let n: Double = ((x - len) / (2 * len)).rounded(.down)
            if !(n < 1.0995e12) { return -1 }
            return max(1, Int(n))
        }
        return 1
    }

    /// `v - n * w` componentwise; nil if any component overflows `Int` (the reduction then stops unconverged).
    static func sheared(_ v: SIMD3<Int>, by n: Int, along w: SIMD3<Int>) -> SIMD3<Int>? {
        var out = SIMD3<Int>(0, 0, 0)
        for k in 0..<3 {
            let (product, o1) = n.multipliedReportingOverflow(by: w[k])
            if o1 { return nil }
            let (difference, o2) = v[k].subtractingReportingOverflow(product)
            if o2 { return nil }
            out[k] = difference
        }
        return out
    }

    /// `a + b + c` componentwise; nil on overflow.
    static func checkedSum(_ a: SIMD3<Int>, _ b: SIMD3<Int>, _ c: SIMD3<Int>) -> SIMD3<Int>? {
        var out = SIMD3<Int>(0, 0, 0)
        for k in 0..<3 {
            let (ab, o1) = a[k].addingReportingOverflow(b[k])
            if o1 { return nil }
            let (abc, o2) = ab.addingReportingOverflow(c[k])
            if o2 { return nil }
            out[k] = abc
        }
        return out
    }

    /// xᵀ G y.
    static func bilinear(_ g: [SIMD3<Double>], _ x: SIMD3<Int>, _ y: SIMD3<Int>) -> Double {
        let xd = SIMD3<Double>(Double(x.x), Double(x.y), Double(x.z))
        let yd = SIMD3<Double>(Double(y.x), Double(y.y), Double(y.z))
        return dot(xd, SIMD3(dot(g[0], yd), dot(g[1], yd), dot(g[2], yd)))
    }

    /// Pᵀ G P for the columns of P.
    static func gram(of p: IntegerMatrix, _ g: [SIMD3<Double>]) -> [SIMD3<Double>] {
        let c = (0..<3).map { p.column($0) }
        func row(_ i: Int) -> SIMD3<Double> {
            SIMD3(bilinear(g, c[i], c[0]), bilinear(g, c[i], c[1]), bilinear(g, c[i], c[2]))
        }
        return [row(0), row(1), row(2)]
    }

    /// The cell whose metric is `g` (lengths and angles read from it).
    static func cell(fromMetric g: [SIMD3<Double>]) -> UnitCell {
        let a = g[0].x.squareRoot(), b = g[1].y.squareRoot(), c = g[2].z.squareRoot()
        func angle(_ gij: Double, _ x: Double, _ y: Double) -> Double {
            acos(min(1, max(-1, gij / (x * y)))) * 180 / .pi
        }
        return UnitCell(a: a, b: b, c: c, alphaDeg: angle(g[1].z, b, c), betaDeg: angle(g[0].z, a, c),
                        gammaDeg: angle(g[0].y, a, b))
    }
}

// MARK: - Integer matrices as lattice transformations (internal)

extension IntegerMatrix {
    /// The matrix whose columns are `columns`.
    init(columns c: [SIMD3<Int>]) {
        self.init(rows: [SIMD3(c[0].x, c[1].x, c[2].x), SIMD3(c[0].y, c[1].y, c[2].y), SIMD3(c[0].z, c[1].z, c[2].z)])
    }

    func column(_ j: Int) -> SIMD3<Int> { SIMD3(rows[0][j], rows[1][j], rows[2][j]) }

    /// M v.
    func times(_ v: SIMD3<Int>) -> SIMD3<Int> {
        SIMD3(rows[0].x * v.x + rows[0].y * v.y + rows[0].z * v.z,
              rows[1].x * v.x + rows[1].y * v.y + rows[1].z * v.z,
              rows[2].x * v.x + rows[2].y * v.y + rows[2].z * v.z)
    }

    /// The classical adjoint: M⁻¹ = adjugate / det.
    var adjugate: IntegerMatrix {
        let m = rows
        return IntegerMatrix(
            SIMD3(m[1].y * m[2].z - m[1].z * m[2].y, m[0].z * m[2].y - m[0].y * m[2].z, m[0].y * m[1].z - m[0].z * m[1].y),
            SIMD3(m[1].z * m[2].x - m[1].x * m[2].z, m[0].x * m[2].z - m[0].z * m[2].x, m[0].z * m[1].x - m[0].x * m[1].z),
            SIMD3(m[1].x * m[2].y - m[1].y * m[2].x, m[0].y * m[2].x - m[0].x * m[2].y, m[0].x * m[1].y - m[0].y * m[1].x))
    }

    /// The product `self · other`, or nil if any entry overflows `Int`.
    func multipliedChecked(by other: IntegerMatrix) -> IntegerMatrix? {
        func entry(_ r: SIMD3<Int>, _ column: Int) -> Int? {
            var total = 0
            for k in 0..<3 {
                let (product, o1) = r[k].multipliedReportingOverflow(by: other.rows[k][column])
                if o1 { return nil }
                let (sum, o2) = total.addingReportingOverflow(product)
                if o2 { return nil }
                total = sum
            }
            return total
        }
        var rowsOut: [SIMD3<Int>] = []
        for r in rows {
            guard let x = entry(r, 0), let y = entry(r, 1), let z = entry(r, 2) else { return nil }
            rowsOut.append(SIMD3(x, y, z))
        }
        return IntegerMatrix(rows: rowsOut)
    }

    /// The classical adjoint, or nil if any minor overflows `Int`.
    var adjugateChecked: IntegerMatrix? {
        func minor(_ a: Int, _ b: Int, _ c: Int, _ d: Int) -> Int? {       // a b - c d
            let (p, o1) = a.multipliedReportingOverflow(by: b)
            let (q, o2) = c.multipliedReportingOverflow(by: d)
            if o1 || o2 { return nil }
            let (difference, o3) = p.subtractingReportingOverflow(q)
            return o3 ? nil : difference
        }
        let m = rows
        guard let a00 = minor(m[1].y, m[2].z, m[1].z, m[2].y), let a01 = minor(m[0].z, m[2].y, m[0].y, m[2].z),
              let a02 = minor(m[0].y, m[1].z, m[0].z, m[1].y),
              let a10 = minor(m[1].z, m[2].x, m[1].x, m[2].z), let a11 = minor(m[0].x, m[2].z, m[0].z, m[2].x),
              let a12 = minor(m[0].z, m[1].x, m[0].x, m[1].z),
              let a20 = minor(m[1].x, m[2].y, m[1].y, m[2].x), let a21 = minor(m[0].y, m[2].x, m[0].x, m[2].y),
              let a22 = minor(m[0].x, m[1].y, m[0].y, m[1].x)
        else { return nil }
        return IntegerMatrix(SIMD3(a00, a01, a02), SIMD3(a10, a11, a12), SIMD3(a20, a21, a22))
    }

    /// M⁻¹ v.
    func inverseApplied(to v: SIMD3<Double>) -> SIMD3<Double> {
        let adj = adjugate.apply(v)
        return adj / Double(determinant)
    }

    /// M⁻¹ N M-style conjugation helper: (self · other · self⁻¹) when every entry divides, else nil.
    func conjugating(_ other: IntegerMatrix) -> IntegerMatrix? {
        let det = determinant
        guard det != 0 else { return nil }
        guard let adjugate = adjugateChecked, let first = multipliedChecked(by: other),
              let product = first.multipliedChecked(by: adjugate) else { return nil }
        var out: [SIMD3<Int>] = []
        for r in product.rows {
            guard r.x % det == 0, r.y % det == 0, r.z % det == 0 else { return nil }
            out.append(SIMD3(r.x / det, r.y / det, r.z / det))
        }
        return IntegerMatrix(rows: out)
    }

    /// (self⁻¹ · other · self) when every entry divides, else nil.
    func conjugatedBy(_ other: IntegerMatrix) -> IntegerMatrix? {
        let det = determinant
        guard det != 0 else { return nil }
        guard let adjugate = adjugateChecked, let first = adjugate.multipliedChecked(by: other),
              let product = first.multipliedChecked(by: self) else { return nil }
        var out: [SIMD3<Int>] = []
        for r in product.rows {
            guard r.x % det == 0, r.y % det == 0, r.z % det == 0 else { return nil }
            out.append(SIMD3(r.x / det, r.y / det, r.z / det))
        }
        return IntegerMatrix(rows: out)
    }
}
