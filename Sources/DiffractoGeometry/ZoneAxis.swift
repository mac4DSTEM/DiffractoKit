//
//  ZoneAxis.swift
//  Role: zone-axis primitives (registration docs/registration-slice-3.md §1, "ZoneAxis"): the 49 low-index beam
//        directions, a direction [u v w] as a Cartesian unit vector of a cell, the detector basis of a beam direction,
//        and the projection of a reflection list onto it with the excitation error. Seeds, read at mac4DSTEM 7d6fdbd1
//        by `git show`: `PhaseReferenceLibrary.lowIndexZoneAxes` and `cartesianZoneAxis`, `ACOMOrientation.detectorBasis`
//        (Core/Analysis/OrientationResult.swift), `OrientationPlan.project` (Core/Crystal/OrientationPlan.swift).
//        Units: Å⁻¹ without 2π, as the rest of the kit.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

/// A zone axis [u v w]: a REAL-space direction u a₁ + v a₂ + w a₃ of a cell.
public struct ZoneAxis: Sendable, Hashable, CustomStringConvertible {
    public let u: Int
    public let v: Int
    public let w: Int

    public init(_ u: Int, _ v: Int, _ w: Int) {
        self.u = u; self.v = v; self.w = w
    }

    public var description: String { "[\(u) \(v) \(w)]" }

    /// Low-index beam directions: |u|, |v|, |w| ≤ 2, coprime, one of each ± pair, in the order of mac4DSTEM 7d6fdbd1
    /// `PhaseReferenceLibrary.lowIndexZoneAxes` (u, v, w each from −2 to 2, the first of a pair met kept): 49 of them,
    /// by value. Deliberately not reduced by crystal symmetry.
    public static let lowIndex: [ZoneAxis] = [
        (-1, -1, -1), (-2, -2, -1), (-1, -1, 0), (-2, -2, 1), (-1, -1, 1), (-2, -1, -2), (-2, -1, -1), (-2, -1, 0),
        (-2, -1, 1), (-2, -1, 2), (-1, 0, -1), (-2, 0, -1), (-1, 0, 0), (-2, 0, 1), (-1, 0, 1), (-2, 1, -2),
        (-2, 1, -1), (-2, 1, 0), (-2, 1, 1), (-2, 1, 2), (-1, 1, -1), (-2, 2, -1), (-1, 1, 0), (-2, 2, 1),
        (-1, 1, 1), (-1, -2, -2), (-1, -2, -1), (-1, -2, 0), (-1, -2, 1), (-1, -2, 2), (-1, -1, -2), (-1, -1, 2),
        (-1, 0, -2), (-1, 0, 2), (-1, 1, -2), (-1, 1, 2), (-1, 2, -2), (-1, 2, -1), (-1, 2, 0), (-1, 2, 1),
        (-1, 2, 2), (0, -1, -1), (0, -2, -1), (0, -1, 0), (0, -2, 1), (0, -1, 1), (0, -1, -2), (0, -1, 2),
        (0, 0, -1),
    ].map { ZoneAxis($0.0, $0.1, $0.2) }

    /// The beam direction as a Cartesian unit vector: u a₁ + v a₂ + w a₃ over the cell's real lattice rows, normalised.
    /// A zone axis is a real-space direction, not the Cartesian triple (u, v, w): the two agree only for a cubic cell.
    /// nil when the direction has no length (0 0 0, or a length at most 1e-12 Å).
    public func cartesian(in cell: UnitCell) -> SIMD3<Double>? {
        let direction = Double(u) * cell.realLattice[0] + Double(v) * cell.realLattice[1] + Double(w) * cell.realLattice[2]
        let n = simd_length(direction)
        guard n.isFinite, n > 1e-12 else { return nil }
        return direction / n
    }

    /// The detector basis of a unit beam direction `n`: columns e₁, e₂, n, orthonormal and right-handed. The reference is
    /// x unless |n.x| ≥ 0.9, then y; e₁ = normalize(n × reference), e₂ = n × e₁ (mac4DSTEM `ACOMOrientation.detectorBasis`).
    public static func detectorBasis(_ n: SIMD3<Double>) -> simd_double3x3 {
        let reference: SIMD3<Double> = abs(n.x) < 0.9 ? [1, 0, 0] : [0, 1, 0]
        let e1 = simd_normalize(simd_cross(n, reference))
        let e2 = simd_cross(n, e1)
        return simd_double3x3(columns: (e1, e2, n))
    }

    /// One reflection seen down a beam direction.
    public struct Projected: Sendable, Equatable {
        public let reflection: KinematicReflection
        /// The excitation error s_g (Å⁻¹).
        public let sg: Double
        /// g · e₁ and g · e₂ (Å⁻¹): the position in the detector plane.
        public let x: Double
        public let y: Double
    }

    /// The reflections excited at the unit beam direction `n`, in the order given, with their positions in the detector
    /// plane of `detectorBasis(n)`.
    ///
    /// Without a wavelength the Ewald sphere is flat: s_g = g · n, odd under g → −g, exactly 0 on the zone law.
    /// With `wavelengthAngstroms` the sphere's curvature is included, s_g = g · n + (λ/2)|g|²: even in g, so g and −g
    /// are excited differently (mac4DSTEM `OrientationPlan.project`).
    ///
    /// DEVIATION (sign) from py4DSTEM f050d207, `Crystal.excitation_errors`: it writes s_g = (2 g_z − λ|g|²) / (2 − 2λ g_z)
    /// with k₀ = [0, 0, −1/λ], and its lab z is the NEGATIVE of this beam direction (the two lab frames differ by an
    /// x ↔ y swap plus a z flip), so in terms of g · n the curvature term has the opposite sign here. Both are the same
    /// physics; only the frame differs. The denominator (2 − 2λ g_z) is not applied: s_g is the first-order form.
    ///
    /// A reflection is kept when |s_g| ≤ `sgMax` (a non-finite s_g is never kept).
    public static func project(_ reflections: [KinematicReflection], onto n: SIMD3<Double>,
                               wavelengthAngstroms: Double? = nil, sgMax: Double) -> [Projected] {
        let basis = detectorBasis(n)
        let e1 = basis.columns.0
        let e2 = basis.columns.1
        let curvature = (wavelengthAngstroms ?? 0) / 2
        var out: [Projected] = []
        out.reserveCapacity(reflections.count)
        for r in reflections {
            let g = r.reflection.g
            let sg = simd_dot(g, n) + curvature * simd_length_squared(g)
            guard abs(sg) <= sgMax else { continue }
            out.append(Projected(reflection: r, sg: sg, x: simd_dot(g, e1), y: simd_dot(g, e2)))
        }
        return out
    }
}
