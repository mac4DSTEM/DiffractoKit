//
//  Misorientation.swift
//  Role: orix 0.15.0's `Misorientation((Gl, Gr))` (`equivalent`, `reduce` = `map_into_symmetry_reduced_zone`) and `Orientation`
//        (`angle_with`, `in_euler_fundamental_region`) on the kit's rotation matrices and `CrystalOrientationSymmetry` groups
//        (docs/parity.md, ROADMAP task 35). Pins: Tests/DiffractoOrientationTests/Fixtures/orix-misorientation-pins.json.
//
//  CONVENTIONS. A `Misorientation` carries a rotation matrix R that maps vectors of the SOURCE crystal to vectors of the TARGET crystal
//  (v_target = R v_source). Its symmetrically equivalent forms are g_t R g_s, g_t a proper rotation of the target group, g_s of the source
//  group: orix's `Misorientation(M, symmetry=(Gl, Gr))` with Gl = source and Gr = target, whose `equivalent()` is Gr * M * Gl (the matrix
//  of a quaternion product is the product of the matrices; orix's `to_matrix` is the plain active rotation matrix). The kit's orientation
//  matrices (columns = lab axes in crystal coordinates) are orix's `Orientation` matrices, equivalents `S * M` (unchanged from
//  `CrystalOrientationSymmetry`).
//
//  DEVIATION (orix `Misorientation.reduce`): orix returns the equivalent inside its misorientation fundamental zone, which is documented as
//  "the smallest rotation angle", and is not for every pair of groups: in the pins, 1 of 41 draws of (m-3m, 6/mmm) and 3 of 41 of (m-3, -3)
//  come back with a larger angle than another equivalent (orix 2.18 against 0.92 rad on one draw). The kit returns the smallest angle
//  (brute force over the equivalents, ties by the canonical quaternion of `CrystalOrientationSymmetry.reduce`); the pins compare the kit with
//  orix where orix is minimal and with the minimum where it is not.
//

import Foundation
import simd
import DiffractoCore

/// A misorientation between two crystals with their proper point groups (the Laue class's rotation subgroup).
public struct Misorientation: Sendable {
    /// Maps source-crystal vectors to target-crystal vectors.
    public let matrix: simd_double3x3
    public let source: CrystalOrientationSymmetry
    public let target: CrystalOrientationSymmetry

    /// Throws `OrientationSymmetryError.invalidRotation` unless `matrix` is a finite proper rotation (tolerance 1e-8).
    public init(_ matrix: simd_double3x3, source: CrystalOrientationSymmetry, target: CrystalOrientationSymmetry) throws {
        guard Misorientation.isRotation(matrix) else { throw OrientationSymmetryError.invalidRotation }
        self.matrix = matrix; self.source = source; self.target = target
    }

    /// The rotation angle, radians, in [0, pi].
    public var angleRad: Double { Misorientation.angle(of: matrix) }

    /// orix `equivalent(grain_exchange)`: `g_t R g_s` for every target operator, then every source operator (orix's order, duplicates kept:
    /// |target| x |source| matrices). With `grainExchange` and equal groups the inverse and its equivalents follow.
    public func equivalent(grainExchange: Bool = false) -> [simd_double3x3] {
        var bases = [matrix]
        if grainExchange && source.sameGroup(as: target) { bases.append(matrix.transpose) }
        var out: [simd_double3x3] = []
        out.reserveCapacity(bases.count * target.operators.count * source.operators.count)
        for base in bases {
            for gt in target.operators { for gs in source.operators { out.append(gt * base * gs) } }
        }
        return out
    }

    /// orix `reduce()` / `map_into_symmetry_reduced_zone()`: the equivalent with the smallest rotation angle (see the DEVIATION in the file header).
    public func reduce() -> Misorientation {
        var best = matrix
        var bestAngle = Double.infinity
        var bestKey: [Double]?
        for candidate in equivalent() {
            let angle = Misorientation.angle(of: candidate)
            let key = Misorientation.canonicalKey(candidate)
            if angle < bestAngle - 1e-13 || (abs(angle - bestAngle) <= 1e-13 && bestKey.map { Misorientation.precedes(key, $0) } ?? true) {
                best = candidate; bestAngle = angle; bestKey = key
            }
        }
        return Misorientation(unchecked: best, source: source, target: target)
    }

    public func mapIntoSymmetryReducedZone() -> Misorientation { reduce() }

    private init(unchecked matrix: simd_double3x3, source: CrystalOrientationSymmetry, target: CrystalOrientationSymmetry) {
        self.matrix = matrix; self.source = source; self.target = target
    }

    // MARK: Helpers shared with `Orientation`

    static func isRotation(_ m: simd_double3x3) -> Bool {
        for c in 0 ..< 3 { for r in 0 ..< 3 where !m[c][r].isFinite { return false } }
        let gap = (m.transpose * m - matrix_identity_double3x3)
        var largest = 0.0
        for c in 0 ..< 3 { for r in 0 ..< 3 { largest = max(largest, abs(gap[c][r])) } }
        return largest <= 1e-8 && abs(simd_determinant(m) - 1) <= 1e-8
    }

    static func angle(of m: simd_double3x3) -> Double {
        let q = simd_quatd(m)
        return 2 * atan2(simd_length(q.imag), abs(q.real))
    }

    /// The quaternion `(-w, x, y, z)` with `w >= 0` (a half-turn: first nonzero imaginary component positive), the tie key of `reduce`.
    static func canonicalKey(_ m: simd_double3x3) -> [Double] {
        let q = simd_quatd(m)
        var v = [q.real, q.imag.x, q.imag.y, q.imag.z]
        let negate: Bool
        if abs(v[0]) < 1e-14 { negate = v.dropFirst().first(where: { abs($0) > 1e-14 }).map { $0 < 0 } == true } else { negate = v[0] < 0 }
        if negate { v = v.map { -$0 } }
        return [-v[0], v[1], v[2], v[3]]
    }

    static func precedes(_ a: [Double], _ b: [Double]) -> Bool {
        for (x, y) in zip(a, b) where abs(x - y) > 1e-13 { return x < y }
        return false
    }
}

extension CrystalOrientationSymmetry {
    /// True when both hold the same set of proper operators.
    func sameGroup(as other: CrystalOrientationSymmetry) -> Bool {
        guard operators.count == other.operators.count else { return false }
        return operators.allSatisfy { a in
            other.operators.contains { b in
                var gap = 0.0
                for c in 0 ..< 3 { for r in 0 ..< 3 { gap = max(gap, abs(a[c][r] - b[c][r])) } }
                return gap < 1e-9
            }
        }
    }
}

/// An orientation (the kit's matrix convention, columns = lab axes in crystal coordinates) with the crystal's symmetry attached.
public struct Orientation: Sendable {
    public let matrix: simd_double3x3
    public let symmetry: CrystalOrientationSymmetry

    public init(_ matrix: simd_double3x3, symmetry: CrystalOrientationSymmetry) throws {
        guard Misorientation.isRotation(matrix) else { throw OrientationSymmetryError.invalidRotation }
        self.matrix = matrix; self.symmetry = symmetry
    }

    /// orix `Orientation.reduce()`: the equivalent `S * M` with the smallest rotation angle (`CrystalOrientationSymmetry.reduce`).
    public func reduce() throws -> Orientation {
        let r = try symmetry.reduce(matrix)
        return Orientation(unchecked: r.matrix, symmetry: symmetry)
    }

    /// orix `Orientation.angle_with(other)`: the disorientation angle to `other` of the same crystal, radians. Throws `differentSymmetry` when
    /// the two orientations do not share a group.
    public func angle(with other: Orientation) throws -> Double {
        guard symmetry.sameGroup(as: other.symmetry) else { throw OrientationSymmetryError.differentSymmetry }
        return try symmetry.disorientation(from: matrix, to: other.matrix).angleRad
    }

    private init(unchecked matrix: simd_double3x3, symmetry: CrystalOrientationSymmetry) {
        self.matrix = matrix; self.symmetry = symmetry
    }

    // MARK: Euler fundamental region

    /// The proper subgroup's name in orix's convention for the kit's standard setting ("121" for 2/m with b unique, "312" for -31m, "32" for -3m1).
    var properGroupName: String {
        switch symmetry.laueClass {
        case .oneBar: return "1"
        case .twoOverM: return "121"
        case .mmm: return "222"
        case .fourOverM: return "4"
        case .fourOverMmm: return "422"
        case .threeBar: return "3"
        case .threeBarM: return symmetry.threeBarMSetting == .threeBarMOne ? "32" : "312"
        case .sixOverM: return "6"
        case .sixOverMmm: return "622"
        case .mThreeBar: return "23"
        case .mThreeBarM: return "432"
        }
    }

    /// orix `Symmetry.euler_fundamental_region`: the maximum (phi1, Phi, phi2) in degrees.
    var eulerRegionDegrees: SIMD3<Double> {
        switch properGroupName {
        case "121": return SIMD3(360, 90, 360)
        case "222": return SIMD3(360, 90, 180)
        case "4": return SIMD3(360, 180, 90)
        case "422": return SIMD3(360, 90, 90)
        case "3": return SIMD3(360, 180, 120)
        case "312", "32": return SIMD3(360, 90, 120)
        case "6": return SIMD3(360, 180, 60)
        case "622": return SIMD3(360, 90, 60)
        case "23": return SIMD3(360, 90, 180)
        case "432": return SIMD3(360, 90, 90)
        default: return SIMD3(360, 180, 360)
        }
    }

    /// orix `Symmetry._primary_axis_order` (nil in orix for a name it does not know; the 1-fold cases are 1 here).
    var primaryAxisOrder: Int {
        switch properGroupName {
        case "222", "23": return 2
        case "3", "32", "312": return 3
        case "4", "422", "432": return 4
        case "6", "622": return 6
        default: return 1
        }
    }

    /// orix `Symmetry._special_rotation` (from MTEX `Symmetry.rotation_special`): the proper operators other than rotations about the c axis,
    /// as the quaternions orix multiplies in front of the orientation.
    func specialRotations() -> [simd_quatd] {
        func about(_ axis: SIMD3<Double>, _ n: Int) -> [simd_quatd] {
            let v = simd_normalize(axis)
            return (0 ..< n).map { k in
                let theta = 2 * Double.pi * Double(k) / Double(n)
                return simd_quatd(real: cos(theta / 2), imag: sin(theta / 2) * v)
            }
        }
        // orix flattens `rot.outer(two_fold)` with the two-fold index running slowest: [r0 e0, r1 e0, r2 e0, r0 e1, ...], each r_i * e_j
        func outer(_ a: [simd_quatd], _ b: [simd_quatd]) -> [simd_quatd] { b.flatMap { y in a.map { x in x * y } } }
        switch properGroupName {
        case "121":
            return symmetry.operators.map { simd_quatd($0) }
        case "222", "422", "622", "32":
            return about(SIMD3(-1, 0, 0), 2)
        case "312":
            return about(SIMD3(-1, 1, 0), 2)
        case "23":
            return outer(about(SIMD3(-1, -1, -1), 3), about(SIMD3(-1, 0, 0), 2))
        case "432":
            return outer(about(SIMD3(-1, -1, -1), 3), about(SIMD3(-1, -1, 0), 2))
        default:
            return [simd_quatd(real: 1, imag: .zero)]
        }
    }

    /// orix `Orientation.in_euler_fundamental_region()`: the Bunge Euler angles (radians) of the first symmetry-equivalent orientation inside the
    /// fundamental region of the proper subgroup (MTEX `project2EulerFR`): the special rotations are tried in front of the orientation, phi2 is
    /// taken modulo 2 pi over the primary axis order, and the first triplet inside the region wins (the last candidate when none is: orix returns
    /// the first of `argmax` of an all-false row, which is index 0).
    public func inEulerFundamentalRegion() -> SIMD3<Double> {
        let q = simd_quatd(matrix)
        let region = eulerRegionDegrees * (Double.pi / 180)
        let order = Double(primaryAxisOrder)
        var triplets: [SIMD3<Double>] = []
        for s in specialRotations() {
            let e = Orientation.bungeEuler(simd_normalize(s * q))
            triplets.append(SIMD3(e.x, e.y, Orientation.positiveMod(e.z, 2 * Double.pi / order)))
        }
        // orix compares with <=, with no tolerance
        let inside = triplets.firstIndex { $0.x <= region.x && $0.y <= region.y && $0.z <= region.z } ?? 0
        return triplets[inside]
    }

    /// orix `qu2eu_single` (Rowenhorst 2015, Eq. A.14), radians in [0, 2 pi).
    static func bungeEuler(_ q: simd_quatd) -> SIMD3<Double> {
        let a = q.real, b = q.imag.x, c = q.imag.y, d = q.imag.z
        let qad = a * a + d * d, qbc = b * b + c * c
        let chi = (qad * qbc).squareRoot()
        var eu = SIMD3<Double>(0, 0, 0)
        let eps = 1e-9
        if chi < eps {
            if qbc < eps {
                eu.x = atan2(-2 * a * d, a * a - d * d)
            } else {
                eu.x = atan2(-2 * b * c, b * b - c * c)
                eu.y = Double.pi
            }
            return SIMD3(positiveMod(eu.x, 2 * .pi), positiveMod(eu.y, 2 * .pi), positiveMod(eu.z, 2 * .pi))
        }
        eu.x = atan2((b * d - a * c) / chi, (-a * b - c * d) / chi)
        eu.y = atan2(2 * chi, qad - qbc)
        eu.z = atan2((a * c + b * d) / chi, (c * d - a * b) / chi)
        for i in 0 ..< 3 where abs(eu[i]) < eps { eu[i] = 0 }
        return SIMD3(positiveMod(eu.x, 2 * .pi), positiveMod(eu.y, 2 * .pi), positiveMod(eu.z, 2 * .pi))
    }

    static func positiveMod(_ x: Double, _ m: Double) -> Double {
        let r = x - m * (x / m).rounded(.down)       // numpy.mod: the result has the divisor's sign
        return r
    }
}
