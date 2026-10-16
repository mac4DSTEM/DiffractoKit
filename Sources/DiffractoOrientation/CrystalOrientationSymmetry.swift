import Foundation
import simd
import DiffractoCore

/// Explicit colour conventions. Neither key performs py4DSTEM's direction-family reduction.
public enum IPFColorKey: Sendable, Equatable, CaseIterable {
    /// py4DSTEM f050d207 `plot_orientation_maps`, full mask (1), beam colour, rounded to Float.
    /// Supply its `family` pole, except for 2/m which uses the original pole (the NNLS branch).
    case py4DSTEM
    /// The original cubic, hexagonal and b-unique monoclinic keys. This convention folds the
    /// input pole itself; it is independent of py4DSTEM's colour convention.
    case native
}

/// Explicit refusals of invalid inputs or a colour convention with no defined key.
public enum OrientationSymmetryError: Error, Sendable, Equatable {
    /// Expected a finite, orthonormal matrix with determinant +1 (tolerance 1e-8).
    case invalidRotation
    /// Expected a finite nonzero pole with a representable length greater than 1e-15.
    case invalidDirection
    /// Sampling accepts counts in 0...4096; other counts are refused before allocation.
    case invalidSampleCount
    /// The native colour convention defines only m-3m, 6/mmm and b-unique 2/m.
    case nativeKeyUnavailable(LaueClass)
}

/// Proper orientation symmetry for every Laue class, without an unknown-class fallback.
///
/// Standard Cartesian settings: cubic/orthorhombic axes x,y,z; monoclinic unique b = y;
/// tetragonal c = z; trigonal/hexagonal a = x, c = z, b at 120° to a. The -3m table is
/// pymatgen's standard hexagonal setting. Transform an orientation into this frame before use.
/// Columns of an input orientation are lab axes expressed in crystal coordinates, as in
/// py4DSTEM. Equivalent rotations are `S * matrix`. Inversion is used only for pole directions.
public struct CrystalOrientationSymmetry: Sendable {
    /// The Laue class in the documented standard setting.
    public let laueClass: LaueClass
    /// Proper rotations only, deterministic order, verified against pymatgen's point-group tables.
    public let operators: [simd_double3x3]

    /// Constructs the complete proper subgroup of any of the eleven Laue classes.
    public init(laueClass: LaueClass) {
        self.laueClass = laueClass
        let x = IntegerMatrix([1,0,0],[0,-1,0],[0,0,-1])
        let y = IntegerMatrix([-1,0,0],[0,1,0],[0,0,-1])
        let four = IntegerMatrix([0,-1,0],[1,0,0],[0,0,1])
        let three = IntegerMatrix([0,-1,0],[1,-1,0],[0,0,1])
        let six = IntegerMatrix([1,-1,0],[1,0,0],[0,0,1])
        let diagonal = IntegerMatrix([0,0,1],[1,0,0],[0,1,0])
        let generators: [IntegerMatrix]
        switch laueClass {
        case .oneBar: generators = []
        case .twoOverM: generators = [y]
        case .mmm: generators = [x,y]
        case .fourOverM: generators = [four]
        case .fourOverMmm: generators = [four,x]
        case .threeBar: generators = [three]
        case .threeBarM:
            generators = [three,IntegerMatrix([0,-1,0],[-1,0,0],[0,0,-1])]
        case .sixOverM: generators = [six]
        case .sixOverMmm:
            generators = [six,IntegerMatrix([1,-1,0],[0,-1,0],[0,0,-1])]
        case .mThreeBar: generators = [x,diagonal]
        case .mThreeBarM: generators = [four,diagonal]
        }
        let fractional = PointSymmetry(closureOf:generators).rotations
        let hex = [.threeBar,.threeBarM,.sixOverM,.sixOverMmm].contains(laueClass)
        let frame = hex ? simd_double3x3(columns:(SIMD3(1,0,0),SIMD3(-0.5,sqrt(3)/2,0),SIMD3(0,0,1)))
                        : matrix_identity_double3x3
        self.operators = fractional.map { r in
            let rows = r.rows.map { SIMD3<Double>(Double($0.x),Double($0.y),Double($0.z)) }
            return frame * simd_double3x3(rows:rows) * frame.inverse
        }
    }

    /// Returns `S * matrix` with minimum rotation angle to identity (proper rotations only).
    /// Canonical quaternion components `(−w,x,y,z)` break ties, comparing differences over 1e-13;
    /// at a 180° rotation the first nonzero imaginary component fixes the quaternion sign.
    /// The angle is radians stored as Float, matching orientation-result diagnostics.
    /// This is a single rotation representative, not py4DSTEM's three independent pole families.
    public func reduce(_ matrix: simd_double3x3) throws ->
        (matrix: simd_double3x3, disorientationRad: Float) {
        guard Self.finite(matrix), Self.matrixGap(matrix.transpose * matrix,matrix_identity_double3x3) <= 1e-8,
              abs(simd_determinant(matrix)-1) <= 1e-8 else { throw OrientationSymmetryError.invalidRotation }
        var best = matrix
        var bestKey: [Double]?
        var angle = Double.pi
        for op in operators {
            let candidate = op * matrix
            let q = simd_quatd(candidate.transpose)
            var values = [q.real,q.imag.x,q.imag.y,q.imag.z]
            let negate: Bool
            if abs(values[0]) < 1e-14 {
                negate = values.dropFirst().first(where: { abs($0) > 1e-14 }).map { $0 < 0 } == true
            } else {
                negate = values[0] < 0
            }
            if negate {
                values = values.map { -$0 }
            }
            let key = [-values[0],values[1],values[2],values[3]]
            if bestKey == nil || Self.precedes(key,bestKey!) {
                bestKey = key
                best = candidate
                angle = 2 * acos(min(1,max(0,values[0])))
            }
        }
        return (best,Float(angle))
    }

    /// Folds a pole under the full Laue group (proper operators and inversion).
    /// Native sectors are preserved for m-3m (z ≥ x ≥ y ≥ 0), 6/mmm (z ≥ 0, azimuth 0...30°)
    /// and 2/m (y ≥ 0, azimuth atan2(z,x) in 0...180°). For all other classes the largest
    /// `(z,x,y)` orbit member is used, comparing differences over 1e-13.
    /// This deterministic sampling sector is not py4DSTEM's direction-family reduction.
    public func reduceDirection(_ direction: SIMD3<Double>) throws -> SIMD3<Double> {
        try Self.validate(direction)
        switch laueClass {
        case .mThreeBarM: return NativeCubicSymmetry.reduceDirection(direction)
        case .sixOverMmm: return NativeHexagonalSymmetry.reduceDirection(direction)
        case .twoOverM: return NativeMonoclinicSymmetry.reduceDirection(direction)
        default:
            let d = simd_normalize(direction)
            var best = d
            for op in operators {
                for sign in [-1.0,1.0] {
                    let v = sign * (op * d)
                    if Self.precedes([-v.z,-v.x,-v.y],[-best.z,-best.x,-best.y]) { best = v }
                }
            }
            return simd_normalize(best)
        }
    }

    /// Deterministic Fibonacci candidates folded into the direction sector, followed by
    /// spherical farthest-point selection. The three native-class samplers preserve their original
    /// seeds and arithmetic. Other classes start with folded z, x, y (deduplicated), using 16
    /// candidates per requested point (at least 192). Zero returns []; counts above 4096 or below
    /// zero throw before allocation. Sampling is a native policy, not py4DSTEM's angular grid.
    public func sampleFundamentalZone(count: Int) throws -> [SIMD3<Double>] {
        guard (0...4096).contains(count) else { throw OrientationSymmetryError.invalidSampleCount }
        guard count > 0 else { return [] }
        switch laueClass {
        case .mThreeBarM: return NativeCubicSymmetry.sampleFundamentalZone(count:count)
        case .sixOverMmm: return NativeHexagonalSymmetry.sampleFundamentalZone(count:count)
        case .twoOverM: return NativeMonoclinicSymmetry.sampleFundamentalZone(count:count)
        default: break
        }
        var selected: [SIMD3<Double>] = []
        for d in [SIMD3<Double>(0,0,1),SIMD3<Double>(1,0,0),SIMD3<Double>(0,1,0)] {
            let v = try reduceDirection(d)
            if !selected.contains(where:{simd_dot(v,$0)>1-1e-13}) { selected.append(v) }
        }
        if count <= selected.count { return Array(selected.prefix(count)) }
        let n = max(192,count*16)
        let golden = Double.pi * (3-sqrt(5))
        let candidates = try (0..<n).map { i -> SIMD3<Double> in
            let z = 1-2*(Double(i)+0.5)/Double(n)
            let r = sqrt(max(0,1-z*z)), a = golden*Double(i)
            return try reduceDirection([r*cos(a),r*sin(a),z])
        }
        var distances = candidates.map { v in selected.map { 1-simd_dot(v,$0) }.min()! }
        while selected.count < count {
            let best = distances.indices.max { distances[$0] < distances[$1] }!
            let v = candidates[best]
            selected.append(v)
            for i in candidates.indices { distances[i] = min(distances[i],1-simd_dot(candidates[i],v)) }
            distances[best] = -.infinity
        }
        return selected
    }

    /// Colours a pole with the explicitly selected convention, returning RGB Float components.
    ///
    /// `.py4DSTEM`: pass the already reduced `Orientation.family` pole in the pinned automatic
    /// sector. For 2/m, pass the original matrix pole instead (Python uses NNLS, not `family`).
    /// The full-mask beam key is reproduced; correlation masking and plotting belong to the caller.
    /// Inputs are not renormalized in this branch, matching Python's clipping and normalization.
    /// `.native`: pass the original pole; folding and the original three-class keys are applied.
    /// Other native classes throw `nativeKeyUnavailable`; no identity colour fallback exists.
    public func ipfColor(direction: SIMD3<Double>, key: IPFColorKey) throws -> SIMD3<Float> {
        try Self.validate(direction)
        switch key {
        case .native:
            switch laueClass {
            case .mThreeBarM: return NativeCubicSymmetry.ipfColor(direction:direction)
            case .sixOverMmm: return NativeHexagonalSymmetry.ipfColor(direction:direction)
            case .twoOverM: return NativeMonoclinicSymmetry.ipfColor(direction:direction)
            default: throw OrientationSymmetryError.nativeKeyUnavailable(laueClass)
            }
        case .py4DSTEM:
            return py4DSTEMColor(direction)
        }
    }

    private static func validate(_ d: SIMD3<Double>) throws {
        let length = simd_length(d)
        guard length.isFinite, length > 1e-15 else { throw OrientationSymmetryError.invalidDirection }
    }
    private static func finite(_ m: simd_double3x3) -> Bool {
        (0..<3).allSatisfy { c in (0..<3).allSatisfy { m[c][$0].isFinite } }
    }
    private static func matrixGap(_ a: simd_double3x3,_ b: simd_double3x3) -> Double {
        (0..<3).flatMap { c in (0..<3).map { abs(a[c][$0]-b[c][$0]) } }.max()!
    }
    private static func precedes(_ a: [Double],_ b: [Double]) -> Bool {
        for (x,y) in zip(a,b) where abs(x-y)>1e-13 { return x<y }
        return false
    }
}
