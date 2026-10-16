//
//  OrientationResult.swift
//  Role: the value types of an orientation result: `EulerAngles` (the py4DSTEM / SciPy Bunge convention, both directions),
//        `ZoneAxisFrame` (the detector frame of a zone axis and the orientation matrix built from it) and `OrientationResult`
//        (one position's result). Registration docs/registration-slice-4.md, sections 1 and 4. The map and its presentation
//        stay with the caller; symmetry primitives are in CrystalOrientationSymmetry.swift.
//
//  `ZoneAxisFrame` is the apps' application-specific frame type under a library name; the renaming is the only change.
//

import Foundation
import simd

/// py4DSTEM/orix-compatible Bunge Euler angles in radians. py4DSTEM stores an
/// orientation matrix whose columns are lab x/y/z in crystal coordinates,
/// then exports `matrix.T` through SciPy's extrinsic `zxz` decomposition.
public struct EulerAngles: Sendable, Equatable, Hashable {
    public var phi1: Float
    public var Phi: Float
    public var phi2: Float

    public init(phi1: Float, Phi: Float, phi2: Float) {
        self.phi1 = phi1
        self.Phi = Phi
        self.phi2 = phi2
    }

    public static let zero = EulerAngles(phi1: 0, Phi: 0, phi2: 0)

    /// Decompose in py4DSTEM's orientation-matrix convention: transpose
    /// crystal-coordinate lab axes to crystal→lab, then extrinsic zxz. At
    /// gimbal lock, match SciPy by folding the observable rotation into phi1
    /// and setting phi2 to zero.
    public init(py4DSTEMOrientationMatrix matrix: simd_double3x3) {
        let r = simd_transpose(matrix)
        let r22 = min(1.0, max(-1.0, r.columns.2.z))
        let middle = acos(r22)
        let sine = sin(middle)
        let first: Double
        let third: Double
        // The general branch needs the angle in the off-diagonal terms too: when R22 rounds to ±(1 − 2⁻⁵³) acos gives
        // Φ = 1.49e-8 (sin Φ above the threshold) while R20, R21 can be exactly 0, and atan2(0, 0) would lose the rotation
        // about the beam. Such a matrix is at the gimbal and is folded into φ1 below. Identical to the apps' original wherever it
        // returns the rotation.
        if abs(sine) > 1e-8, (r.columns.0.z * r.columns.0.z + r.columns.1.z * r.columns.1.z).squareRoot() > 1e-8 {
            first = atan2(r.columns.0.z, r.columns.1.z)   // atan2(R20,R21)
            third = atan2(r.columns.2.x, -r.columns.2.y) // atan2(R02,-R12)
        } else if r22 > 0 {
            first = atan2(r.columns.0.y, r.columns.0.x)   // Phi = 0: phi1+phi2
            third = 0
        } else {
            first = atan2(-r.columns.0.y, r.columns.0.x)  // Phi = pi: phi1-phi2
            third = 0
        }
        self.phi1 = Float(Self.signedAngle(first))
        self.Phi = Float(middle)
        self.phi2 = Float(Self.signedAngle(third))
    }

    /// Crystal→lab rotation matrix corresponding to SciPy `from_euler("zxz")`.
    public var crystalToLabMatrix: simd_double3x3 {
        let a = Double(phi1), b = Double(Phi), c = Double(phi2)
        let ca = cos(a), sa = sin(a), cb = cos(b), sb = sin(b)
        let cc = cos(c), sc = sin(c)
        return simd_double3x3(columns: (
            SIMD3(cc * ca - sc * cb * sa,
                  sc * ca + cc * cb * sa,
                  sb * sa),
            SIMD3(-cc * sa - sc * cb * ca,
                  -sc * sa + cc * cb * ca,
                  sb * ca),
            SIMD3(sc * sb, -cc * sb, cb)
        ))
    }

    public var py4DSTEMOrientationMatrix: simd_double3x3 {
        simd_transpose(crystalToLabMatrix)
    }

    private static func signedAngle(_ value: Double) -> Double {
        var angle = value.truncatingRemainder(dividingBy: 2 * .pi)
        if angle <= -.pi { angle += 2 * .pi }
        if angle > .pi { angle -= 2 * .pi }
        return abs(angle) < 1e-12 ? 0 : angle
    }

    /// (φ1, Φ, φ2) in degrees: each Float angle times `Float(180 / π)`.
    public var degrees: (Float, Float, Float) {
        let scale = Float(180.0 / Double.pi)
        return (phi1 * scale, Phi * scale, phi2 * scale)
    }
}

/// The detector frame of a zone axis and the py4DSTEM orientation matrix built from it.
/// Columns are detector x, detector y, and beam z expressed in crystal axes,
/// matching py4DSTEM's `Orientation.matrix` convention (the frame rotation of the
/// registration, section 4: e1 = normalise(n × r), e2 = n × e1, r = x̂ when |n_x| < 0.9 else ŷ).
public enum ZoneAxisFrame {
    public static func detectorBasis(zoneAxis n: SIMD3<Double>) -> simd_double3x3 {
        let reference: SIMD3<Double> = abs(n.x) < 0.9 ? [1, 0, 0] : [0, 1, 0]
        let e1 = simd_normalize(simd_cross(n, reference))
        let e2 = simd_cross(n, e1)
        return simd_double3x3(columns: (e1, e2, n))
    }

    /// py4DSTEM applies this in-plane matrix on the right of the zone basis.
    public static func matrix(detectorBasis basis: simd_double3x3,
                              inPlaneAngle angle: Float) -> simd_double3x3 {
        let c = cos(Double(angle)), s = sin(Double(angle))
        let inPlane = simd_double3x3(columns: (
            SIMD3(c, -s, 0), SIMD3(s, c, 0), SIMD3(0, 0, 1)
        ))
        return basis * inPlane
    }
}

/// One scan position's orientation result: the template, the Bunge angles, the in-plane angle and the scores.
public struct OrientationResult: Sendable, Equatable {
    /// Index into the template library (zone axis) for the best match.
    public var templateIndex: Int
    public var euler: EulerAngles
    /// Best in-plane rotation angle (radians), from azimuthal correlation.
    public var inPlaneAngle: Float = 0
    /// Best normalized cross-correlation score.
    public var score: Float
    /// Second-best score, used for reliability.
    public var secondScore: Float
    public var phaseID: Int
    /// Rotation angle of the symmetry-reduced representative, useful as a scalar
    /// orientation diagnostic.
    public var symmetryDisorientationRad: Float = 0
    /// The winner came from the matcher's conjugated (in-plane mirror) pass:
    /// the pattern is the x-mirror of template `templateIndex` rotated by
    /// `inPlaneAngle` (which then carries py4DSTEM's +π), and `euler` already
    /// carries the beam flip that makes it a proper rotation (py4DSTEM's
    /// `Orientation.mirror`; a review finding). A template spot at azimuth a
    /// then sits at π − (a + inPlaneAngle), not a + inPlaneAngle.
    public var mirrored: Bool = false

    /// EBSD-style reliability. Higher values indicate a clearer best match.
    public var reliability: Float {
        guard score > 0 else { return 0 }
        let value = 1 - (secondScore / score)
        return max(0, min(1, value))
    }

    public static let empty = OrientationResult(
        templateIndex: -1,
        euler: .zero,
        inPlaneAngle: 0,
        score: 0,
        secondScore: 0,
        phaseID: -1
    )

    // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
    public init(templateIndex: Int, euler: EulerAngles, inPlaneAngle: Float = 0, score: Float, secondScore: Float, phaseID: Int, symmetryDisorientationRad: Float = 0, mirrored: Bool = false) {
        self.templateIndex = templateIndex
        self.euler = euler
        self.inPlaneAngle = inPlaneAngle
        self.score = score
        self.secondScore = secondScore
        self.phaseID = phaseID
        self.symmetryDisorientationRad = symmetryDisorientationRad
        self.mirrored = mirrored
    }
}
