//
//  UnitCell.swift
//  Role: a, b, c, α, β, γ and everything the metric alone determines: the real and reciprocal lattices, the metric
//        tensors, the volume, d(hkl) and g(hkl), the crystal family of the metric at a stated tolerance. Ported from
//        mac4DSTEM 7448401b `Crystal.swift` (the py4DSTEM / VESTA convention: a along x, b in the xy plane; reciprocal
//        vectors without 2π, |g| in Å⁻¹, d = 1/|g| in Å) and generalised to the seven systems; registration
//        docs/registration-slice-1.md sections 3.1 and 3.2.
//

import Foundation
import simd

/// A lattice: the six parameters and the matrices they determine.
public struct UnitCell: Sendable, Equatable {
    /// Lengths in Å.
    public let a: Double
    public let b: Double
    public let c: Double
    /// Angles in degrees: α between b and c, β between a and c, γ between a and b.
    public let alphaDeg: Double
    public let betaDeg: Double
    public let gammaDeg: Double

    /// The real lattice, rows a₁, a₂, a₃ (Å): a₁ along x, a₂ in the xy plane.
    public let realLattice: [SIMD3<Double>]
    /// The reciprocal lattice, rows b₁, b₂, b₃ (Å⁻¹), aᵢ · bⱼ = δᵢⱼ (no 2π).
    public let reciprocalLattice: [SIMD3<Double>]
    /// G = A Aᵀ (Å²).
    public let metricTensor: [SIMD3<Double>]
    /// G* = G⁻¹ (Å⁻²): |g|² = (h k l) G* (h k l)ᵀ.
    public let reciprocalMetricTensor: [SIMD3<Double>]
    /// |det A| (Å³).
    public let volume: Double

    public init(a: Double, b: Double, c: Double, alphaDeg: Double, betaDeg: Double, gammaDeg: Double) {
        self.a = a; self.b = b; self.c = c
        self.alphaDeg = alphaDeg; self.betaDeg = betaDeg; self.gammaDeg = gammaDeg

        let alpha = alphaDeg * .pi / 180
        let beta = betaDeg * .pi / 180
        let gamma = gammaDeg * .pi / 180
        let f = cos(beta) * cos(gamma) - cos(alpha)
        let vol = a * b * c * (1
            + 2 * cos(alpha) * cos(beta) * cos(gamma)
            - cos(alpha) * cos(alpha)
            - cos(beta) * cos(beta)
            - cos(gamma) * cos(gamma)).squareRoot()

        let real = [
            SIMD3<Double>(a, 0, 0),
            SIMD3<Double>(b * cos(gamma), b * sin(gamma), 0),
            SIMD3<Double>(c * cos(beta), -c * f / sin(gamma), vol / (a * b * sin(gamma))),
        ]
        self.realLattice = real
        self.volume = abs(UnitCell.determinant(real))
        let metric = UnitCell.gram(real)
        self.metricTensor = metric
        let metricInverse = UnitCell.invert3x3(metric) ?? UnitCell.identity3
        self.reciprocalMetricTensor = metricInverse
        self.reciprocalLattice = UnitCell.matMul(metricInverse, real)
    }

    // MARK: Constructors per system (the free parameters only)

    public static func cubic(a: Double) -> UnitCell {
        UnitCell(a: a, b: a, c: a, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
    }
    public static func tetragonal(a: Double, c: Double) -> UnitCell {
        UnitCell(a: a, b: a, c: c, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
    }
    public static func orthorhombic(a: Double, b: Double, c: Double) -> UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
    }
    /// Hexagonal axes (γ = 120°), also the conventional cell of a trigonal structure.
    public static func hexagonal(a: Double, c: Double) -> UnitCell {
        UnitCell(a: a, b: a, c: c, alphaDeg: 90, betaDeg: 90, gammaDeg: 120)
    }
    /// Rhombohedral axes: a = b = c, α = β = γ.
    public static func rhombohedral(a: Double, alphaDeg: Double) -> UnitCell {
        UnitCell(a: a, b: a, c: a, alphaDeg: alphaDeg, betaDeg: alphaDeg, gammaDeg: alphaDeg)
    }
    /// Unique axis b: α = γ = 90°.
    public static func monoclinic(a: Double, b: Double, c: Double, betaDeg: Double) -> UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: 90, betaDeg: betaDeg, gammaDeg: 90)
    }
    public static func triclinic(a: Double, b: Double, c: Double,
                                 alphaDeg: Double, betaDeg: Double, gammaDeg: Double) -> UnitCell {
        UnitCell(a: a, b: b, c: c, alphaDeg: alphaDeg, betaDeg: betaDeg, gammaDeg: gammaDeg)
    }

    // MARK: Reciprocal space

    /// g = h b₁ + k b₂ + l b₃ (Å⁻¹, Cartesian).
    public func gVector(_ hkl: MillerIndex) -> SIMD3<Double> {
        Double(hkl.h) * reciprocalLattice[0] + Double(hkl.k) * reciprocalLattice[1]
            + Double(hkl.l) * reciprocalLattice[2]
    }

    /// |g| from the Cartesian vector (Å⁻¹).
    public func gLength(_ hkl: MillerIndex) -> Double {
        length(gVector(hkl))
    }

    /// d = 1/|g| from the reciprocal metric tensor (Å); the self-check against `gLength` is the tests'.
    public func d(_ hkl: MillerIndex) -> Double {
        1 / gSquaredFromMetric(hkl).squareRoot()
    }

    /// |g|² = (h k l) G* (h k l)ᵀ.
    public func gSquaredFromMetric(_ hkl: MillerIndex) -> Double {
        let v = SIMD3<Double>(Double(hkl.h), Double(hkl.k), Double(hkl.l))
        let gv = SIMD3<Double>(dot(reciprocalMetricTensor[0], v),
                               dot(reciprocalMetricTensor[1], v),
                               dot(reciprocalMetricTensor[2], v))
        return dot(v, gv)
    }

    /// The angle between the normals of two planes (degrees), from their g vectors.
    public func interplanarAngleDeg(_ p: MillerIndex, _ q: MillerIndex) -> Double {
        let gp = gVector(p), gq = gVector(q)
        let cosine = dot(gp, gq) / (length(gp) * length(gq))
        return acos(min(1, max(-1, cosine))) * 180 / .pi
    }

    // MARK: Validation

    public enum Issue: Sendable, Equatable, CustomStringConvertible {
        case invalidLengths
        case invalidAngles
        case degenerateVolume

        public var description: String {
            switch self {
            case .invalidLengths: return "Cell lengths must be finite and positive."
            case .invalidAngles: return "Cell angles must be finite and strictly between 0° and 180°."
            case .degenerateVolume: return "The lattice vectors do not form a finite positive-volume cell."
            }
        }
    }

    /// The apps' `CrystalModel.validationIssues` for the cell part, in that order.
    public var issues: [Issue] {
        var out: [Issue] = []
        if [a, b, c].contains(where: { !$0.isFinite || $0 <= 0 }) { out.append(.invalidLengths) }
        if [alphaDeg, betaDeg, gammaDeg].contains(where: { !$0.isFinite || $0 <= 0 || $0 >= 180 }) {
            out.append(.invalidAngles)
        }
        // A cell flatter than V = 1e-6 · a b c is coplanar up to rounding (α = β = γ = 120° gives V ≈ 3e-6 Å³ for 5 Å
        // edges, and structure factors of 1e6): refused (Gate B on slice 2, 2026-10-15).
        if !volume.isFinite || volume <= 0 || volume < 1e-6 * a * b * c { out.append(.degenerateVolume) }
        return out
    }

    public var isValid: Bool { issues.isEmpty }

    // MARK: The crystal family of the metric (registration 3.2)

    public enum Axis: Int, Sendable, Equatable, CaseIterable {
        case a = 0, b = 1, c = 2
    }

    /// What the six parameters alone say, first that applies in this order; a trigonal structure in hexagonal axes
    /// is `.hexagonal` here and `.trigonal` only by its symmetry.
    public enum MetricFamily: Sendable, Equatable {
        case cubic
        case hexagonal
        case tetragonal
        case rhombohedral
        case orthorhombic
        case monoclinic(uniqueAxis: Axis)
        case triclinic
    }

    /// Lengths equal within `lengthTolerance` (relative), angles within `angleToleranceDeg` (degrees). The defaults
    /// are macTEM's indexing tolerance (1e-9); a CIF reader that wants the apps' 1e-3 / 0.05° passes them.
    public func metricFamily(lengthTolerance: Double = 1e-9, angleToleranceDeg: Double = 1e-9) -> MetricFamily {
        func same(_ x: Double, _ y: Double) -> Bool { abs(x - y) <= lengthTolerance * max(abs(x), abs(y)) }
        func angle(_ x: Double, _ target: Double) -> Bool { abs(x - target) <= angleToleranceDeg }
        let ab = same(a, b), ac = same(a, c), bc = same(b, c)
        let ra = angle(alphaDeg, 90), rb = angle(betaDeg, 90), rc = angle(gammaDeg, 90)
        if ab && ac && ra && rb && rc { return .cubic }
        if ab && ra && rb && angle(gammaDeg, 120) { return .hexagonal }
        if ab && !ac && ra && rb && rc { return .tetragonal }
        if ab && ac && bc && same(alphaDeg, betaDeg) && same(alphaDeg, gammaDeg) && !ra { return .rhombohedral }
        if ra && rb && rc { return .orthorhombic }
        if rb && rc && !ra { return .monoclinic(uniqueAxis: .a) }
        if ra && rc && !rb { return .monoclinic(uniqueAxis: .b) }
        if ra && rb && !rc { return .monoclinic(uniqueAxis: .c) }
        return .triclinic
    }

    // MARK: - Small linear-algebra helpers (row-vector matrices), as the apps'

    static let identity3: [SIMD3<Double>] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    static func gram(_ m: [SIMD3<Double>]) -> [SIMD3<Double>] {
        [
            SIMD3(dot(m[0], m[0]), dot(m[0], m[1]), dot(m[0], m[2])),
            SIMD3(dot(m[1], m[0]), dot(m[1], m[1]), dot(m[1], m[2])),
            SIMD3(dot(m[2], m[0]), dot(m[2], m[1]), dot(m[2], m[2])),
        ]
    }

    static func matMul(_ aM: [SIMD3<Double>], _ bM: [SIMD3<Double>]) -> [SIMD3<Double>] {
        func row(_ r: SIMD3<Double>) -> SIMD3<Double> {
            r.x * bM[0] + r.y * bM[1] + r.z * bM[2]
        }
        return [row(aM[0]), row(aM[1]), row(aM[2])]
    }

    static func determinant(_ m: [SIMD3<Double>]) -> Double {
        m[0].x * (m[1].y * m[2].z - m[1].z * m[2].y)
      - m[0].y * (m[1].x * m[2].z - m[1].z * m[2].x)
      + m[0].z * (m[1].x * m[2].y - m[1].y * m[2].x)
    }

    static func invert3x3(_ m: [SIMD3<Double>]) -> [SIMD3<Double>]? {
        let det = determinant(m)
        guard det.isFinite, abs(det) > 1e-15 else { return nil }
        let inv = 1 / det
        return [
            SIMD3((m[1].y * m[2].z - m[1].z * m[2].y) * inv,
                  (m[0].z * m[2].y - m[0].y * m[2].z) * inv,
                  (m[0].y * m[1].z - m[0].z * m[1].y) * inv),
            SIMD3((m[1].z * m[2].x - m[1].x * m[2].z) * inv,
                  (m[0].x * m[2].z - m[0].z * m[2].x) * inv,
                  (m[0].z * m[1].x - m[0].x * m[1].z) * inv),
            SIMD3((m[1].x * m[2].y - m[1].y * m[2].x) * inv,
                  (m[0].y * m[2].x - m[0].x * m[2].y) * inv,
                  (m[0].x * m[1].y - m[0].y * m[1].x) * inv),
        ]
    }
}

/// Minimum-image separation along one fractional axis, in [0, 0.5] (the apps' `axisDelta`).
public func fractionalAxisDelta(_ x: Double, _ y: Double) -> Double {
    var d = abs(x - y).truncatingRemainder(dividingBy: 1)
    if d > 0.5 { d = 1 - d }
    return d
}

/// Wraps fractional coordinates into [0, 1), snapping values within 1e-9 of 0 or 1 to 0 (the apps' `wrap`).
public func wrapFractional(_ v: SIMD3<Double>) -> SIMD3<Double> {
    func wrapAxis(_ value: Double) -> Double {
        var m = value.truncatingRemainder(dividingBy: 1)
        if m < 0 { m += 1 }
        if abs(m) < 1e-9 || abs(m - 1) < 1e-9 { m = 0 }
        return m
    }
    return SIMD3(wrapAxis(v.x), wrapAxis(v.y), wrapAxis(v.z))
}
