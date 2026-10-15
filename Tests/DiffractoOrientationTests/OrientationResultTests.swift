//
//  OrientationResultTests.swift
//  Role: the value types of an orientation result: the Euler angles in both directions, the zone-axis frame and its switch of
//        reference axis (registration slice 4, section 4). The comparison with SciPy on 1000 seeded rotations is gate G3 and
//        lives with the pins; these are the closed forms.
//

import Testing
import Foundation
import simd
import DiffractoOrientation

@Suite struct OrientationResultTests {

    private func maxDifference(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var m = 0.0
        for c in 0..<3 { for r in 0..<3 { m = max(m, abs(a[c][r] - b[c][r])) } }
        return m
    }

    /// Angles → matrix → angles, away from the gimbal: the Float angles come back within Float precision.
    /// Mutation: the transpose dropped (the angles come back negated and swapped).
    @Test func anglesSurviveTheRoundTripThroughTheOrientationMatrix() {
        let triples: [(Float, Float, Float)] = [
            (0.3, 0.7, -1.2), (2.9, 1.5, 0.4), (-2.5, 3.0, 2.2), (1.0, 0.01, 0.5), (-0.4, 2.0, -3.0),
        ]
        for (a, b, c) in triples {
            let e = EulerAngles(phi1: a, Phi: b, phi2: c)
            let back = EulerAngles(py4DSTEMOrientationMatrix: e.py4DSTEMOrientationMatrix)
            #expect(abs(back.phi1 - a) < 2e-6 && abs(back.Phi - b) < 2e-6 && abs(back.phi2 - c) < 2e-6,
                    "\(e) came back as \(back)")
            #expect(maxDifference(back.py4DSTEMOrientationMatrix, e.py4DSTEMOrientationMatrix) < 2e-6)
        }
        // The matrix is the transpose of the crystal → lab matrix.
        let e = EulerAngles(phi1: 0.3, Phi: 0.7, phi2: -1.2)
        #expect(maxDifference(e.py4DSTEMOrientationMatrix, simd_transpose(e.crystalToLabMatrix)) == 0)
        // A proper rotation.
        #expect(abs(simd_determinant(e.crystalToLabMatrix) - 1) < 1e-12)
    }

    /// SciPy's gimbal-lock rule: at Φ = 0 the observable rotation (φ1 + φ2) goes in φ1 and φ2 is 0; at Φ = π it is φ1 − φ2.
    /// The rotation itself is kept. Mutation: the gimbal threshold 1e-8 widened to 1e-4 (the second case, Φ = 5e-6, then
    /// loses φ2).
    @Test func theGimbalCasesFoldIntoTheFirstAngle() {
        let flat = EulerAngles(py4DSTEMOrientationMatrix: EulerAngles(phi1: 0.7, Phi: 0, phi2: 0.4).py4DSTEMOrientationMatrix)
        #expect(abs(flat.phi1 - 1.1) < 1e-6)
        #expect(flat.Phi == 0)
        #expect(flat.phi2 == 0)

        // Φ = π exactly (a Float cannot hold it: Float(π) is above π and is 9e-8 off the gimbal, so the matrix is built in
        // Double): R = Rz(φ2) Rx(π) Rz(φ1), the matrix its transpose.
        func rotZ(_ t: Double) -> simd_double3x3 {
            simd_double3x3(columns: (SIMD3(cos(t), sin(t), 0), SIMD3(-sin(t), cos(t), 0), SIMD3(0, 0, 1)))
        }
        func rotX(_ t: Double) -> simd_double3x3 {
            simd_double3x3(columns: (SIMD3(1, 0, 0), SIMD3(0, cos(t), sin(t)), SIMD3(0, -sin(t), cos(t))))
        }
        let flipped = simd_transpose(rotZ(0.4) * rotX(.pi) * rotZ(0.7))
        let upside = EulerAngles(py4DSTEMOrientationMatrix: flipped)
        #expect(abs(upside.Phi - Float.pi) < 1e-6)
        #expect(upside.phi2 == 0)
        #expect(maxDifference(upside.py4DSTEMOrientationMatrix, flipped) < 1e-6)

        // Just off the gimbal, above the 1e-8 threshold in sin Φ: the two angles are still separate.
        let nearMatrix = EulerAngles(phi1: 0.7, Phi: 5e-6, phi2: 0.4).py4DSTEMOrientationMatrix
        let nearBack = EulerAngles(py4DSTEMOrientationMatrix: nearMatrix)
        #expect(abs(nearBack.phi2 - 0.4) < 1e-3, "φ2 was folded away at Φ = 5e-6, above the gimbal threshold")
    }

    /// The angles are wrapped to (−π, π] and a value under 1e-12 is 0: the identity is all zeros.
    @Test func theIdentityIsAllZeroAndAnglesAreWrapped() {
        let identity = EulerAngles(py4DSTEMOrientationMatrix: matrix_identity_double3x3)
        #expect(identity == .zero)
        for e in [EulerAngles(phi1: 3.1, Phi: 1.0, phi2: -3.1), EulerAngles(phi1: -3.1, Phi: 2.0, phi2: 3.1)] {
            let back = EulerAngles(py4DSTEMOrientationMatrix: e.py4DSTEMOrientationMatrix)
            #expect(back.phi1 > -Float.pi && back.phi1 <= Float.pi)
            #expect(back.phi2 > -Float.pi && back.phi2 <= Float.pi)
        }
        let deg = EulerAngles(phi1: Float.pi, Phi: Float.pi / 2, phi2: -Float.pi / 2).degrees
        #expect(abs(deg.0 - 180) < 1e-4 && abs(deg.1 - 90) < 1e-4 && abs(deg.2 + 90) < 1e-4)
    }

    /// e1 = normalise(n × r), e2 = n × e1, the columns (e1, e2, n); r = x̂ when |n_x| < 0.9, else ŷ. Mutation 3: the switch
    /// at 0.5 instead of 0.9 (the n_x = 0.7 case then takes ŷ).
    @Test func theDetectorBasisFollowsTheDocumentedFrameAndSwitchesAtNineTenths() {
        // [001]
        let z = ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(0, 0, 1))
        #expect(z.columns.0 == SIMD3(0, 1, 0))
        #expect(z.columns.1 == SIMD3(-1, 0, 0))
        #expect(z.columns.2 == SIMD3(0, 0, 1))
        // [100]: |n_x| ≥ 0.9, reference ŷ.
        let x = ZoneAxisFrame.detectorBasis(zoneAxis: SIMD3(1, 0, 0))
        #expect(x.columns.0 == SIMD3(0, 0, 1))
        #expect(x.columns.1 == SIMD3(0, -1, 0))
        // n_x = 0.7 keeps x̂ as the reference: e1 ∝ n × x̂.
        let n7 = simd_normalize(SIMD3(0.7, 0.3, 0.6))
        let b7 = ZoneAxisFrame.detectorBasis(zoneAxis: n7)
        #expect(simd_length(b7.columns.0 - simd_normalize(simd_cross(n7, SIMD3(1, 0, 0)))) < 1e-15)
        // The switch is at 0.9: just below keeps x̂, at 0.9 takes ŷ.
        let below = SIMD3(0.8999999, 0.0, 0.43588990)
        let above = SIMD3(0.9, 0.0, 0.43588989)
        #expect(simd_length(ZoneAxisFrame.detectorBasis(zoneAxis: below).columns.0
                            - simd_normalize(simd_cross(below, SIMD3(1, 0, 0)))) < 1e-15)
        #expect(simd_length(ZoneAxisFrame.detectorBasis(zoneAxis: above).columns.0
                            - simd_normalize(simd_cross(above, SIMD3(0, 1, 0)))) < 1e-15)
        // Orthonormal, right-handed.
        for basis in [z, x, b7] {
            #expect(maxDifference(simd_transpose(basis) * basis, matrix_identity_double3x3) < 1e-14)
            #expect(abs(simd_determinant(basis) - 1) < 1e-14)
        }
    }

    /// matrix(detectorBasis: B, inPlaneAngle: a) = B · R_z(a) with R_z's columns (c, −s, 0), (s, c, 0), (0, 0, 1).
    @Test func theOrientationMatrixIsTheBasisTimesTheInPlaneRotation() {
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: simd_normalize(SIMD3(0.2, -0.5, 0.8)))
        #expect(maxDifference(ZoneAxisFrame.matrix(detectorBasis: basis, inPlaneAngle: 0), basis) == 0)
        let a: Float = 0.9
        let m = ZoneAxisFrame.matrix(detectorBasis: basis, inPlaneAngle: a)
        let c = cos(Double(a)), s = sin(Double(a))
        #expect(simd_length(m.columns.0 - (c * basis.columns.0 - s * basis.columns.1)) < 1e-15)
        #expect(simd_length(m.columns.1 - (s * basis.columns.0 + c * basis.columns.1)) < 1e-15)
        #expect(m.columns.2 == basis.columns.2)
    }

    /// The result value: its reliability is 1 − second / best, clamped, and 0 without a score.
    @Test func theResultReliabilityIsOneMinusTheRatioOfTheScores() {
        var r = OrientationResult.empty
        #expect(r.reliability == 0)
        r = OrientationResult(templateIndex: 3, euler: .zero, score: 0.8, secondScore: 0.2, phaseID: 1)
        #expect(abs(r.reliability - 0.75) < 1e-6)
        r.secondScore = 0.9
        #expect(r.reliability == 0)
        #expect(OrientationResult.empty.templateIndex == -1 && OrientationResult.empty.phaseID == -1)
    }
}
