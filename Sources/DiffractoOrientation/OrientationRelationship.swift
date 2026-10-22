//
//  OrientationRelationship.swift
//  Role: the variants of an orientation relationship between a precipitate and its matrix (ROADMAP task 35, study step 6), for ACOM and
//        phase-map orientations: the matrix group enumerates the crystallographically equivalent placements of the precipitate, the
//        precipitate group identifies those that differ only by one of its own symmetry operations. beta'' (Mg5Si6, C2/m) in Al has
//        24 / 2 = 12 variants, four needle variants for each of the three <001> axes.
//
//  No Python reference computes this (orix has `Misorientation.equivalent`, which this builds on; the class count is group theory).
//  The planted case is tools/beta-variants/planted.py and Andersen et al. 1998 (Acta Mater. 46, 3283), `validation: "none"`.
//

import Foundation
import simd
import DiffractoCore

public struct OrientationVariantSet: Sendable {
    /// Maps precipitate Cartesian vectors to matrix Cartesian vectors (the lattice frame of `UnitCell`: a along x, b in the xy plane).
    public let matrix: simd_double3x3
    public let matrixSymmetry: CrystalOrientationSymmetry
    public let precipitateSymmetry: CrystalOrientationSymmetry
    /// One representative per variant, in the order of the matrix group's operators: `g_matrix * matrix` for the first operator of each class.
    public let variants: [simd_double3x3]

    /// Throws `OrientationSymmetryError.invalidRotation` unless `matrix` is a proper rotation.
    public init(matrix: simd_double3x3, matrixSymmetry: CrystalOrientationSymmetry, precipitateSymmetry: CrystalOrientationSymmetry) throws {
        guard Misorientation.isRotation(matrix) else { throw OrientationSymmetryError.invalidRotation }
        self.matrix = matrix; self.matrixSymmetry = matrixSymmetry; self.precipitateSymmetry = precipitateSymmetry
        self.variants = OrientationRelationship.variants(of: matrix, matrix: matrixSymmetry, precipitate: precipitateSymmetry)
    }

    /// The precipitate-to-matrix rotation of two measured orientations in the kit's matrix convention (columns = lab axes in crystal
    /// coordinates, `v_crystal = M v_lab`): `v_matrix = M_matrix M_precipitate^T v_precipitate`. The precipitate's crystal frame must be the
    /// `UnitCell` Cartesian frame, as for `matrix`.
    public static func observedRelation(matrixOrientation: simd_double3x3, precipitateOrientation: simd_double3x3) -> simd_double3x3 {
        matrixOrientation * precipitateOrientation.transpose
    }

    /// The variant nearest to an observed precipitate-to-matrix rotation `observed` by disorientation angle (the precipitate group on the right,
    /// `min over s of angle(R_k^T observed s)`), its angle, and the angle to the second nearest (the margin; a small margin is an ambiguous call).
    public func nearestVariant(of observed: simd_double3x3) -> (index: Int, disorientationRad: Double, nextRad: Double?)? {
        guard Misorientation.isRotation(observed), !variants.isEmpty else { return nil }
        var angles = variants.map { r -> Double in
            precipitateSymmetry.operators.map { s in Misorientation.angle(of: r.transpose * observed * s) }.min() ?? .pi
        }
        let best = angles.indices.min { angles[$0] < angles[$1] }!
        let first = angles[best]
        angles.remove(at: best)
        return (best, first, angles.min())
    }

    /// `nearestVariant(of:)` for a matrix orientation and a precipitate orientation of the kit's matrix convention.
    public func nearestVariant(matrixOrientation: simd_double3x3, precipitateOrientation: simd_double3x3) -> (index: Int, disorientationRad: Double, nextRad: Double?)? {
        nearestVariant(of: Self.observedRelation(matrixOrientation: matrixOrientation, precipitateOrientation: precipitateOrientation))
    }

    static func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
        var g = 0.0
        for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
        return g
    }
}

extension OrientationRelationship {
    /// The variants of `relation` (precipitate to matrix): the images `g * relation` over the matrix group's operators `g`, two of them
    /// being the same variant when they differ by a precipitate operator on the right (`R_j = R_i s`). Representatives are the first of each
    /// class in operator order. For a relation that is itself symmetric (an operator of the matrix group maps the relation to its own
    /// precipitate image, as beta'''s 180 degree rotation about b does) the count is |matrix group| / |such operators|.
    public static func variants(of relation: simd_double3x3, matrix: CrystalOrientationSymmetry,
                                precipitate: CrystalOrientationSymmetry) -> [simd_double3x3] {
        var reps: [simd_double3x3] = []
        for g in matrix.operators {
            let candidate = g * relation
            let known = reps.contains { r in precipitate.operators.contains { s in OrientationVariantSet.gap(r * s, candidate) < 1e-9 } }
            if !known { reps.append(candidate) }
        }
        return reps
    }

    /// The precipitate-to-matrix rotation from two statements of the orientation relationship, as a crystallographer writes it
    /// ("[010] of the precipitate is parallel to [001] of the matrix", a plane or a direction on each side): the first statement fixes an axis
    /// exactly, the second only the rotation about it (its component perpendicular to the first is used, so choose the exactly parallel pair
    /// first). Both crystals are in their `UnitCell` Cartesian frame.
    public static func rotation(_ first: OrientationRelationship, _ second: OrientationRelationship,
                                precipitate: CrystalStructure, matrix: CrystalStructure) throws -> simd_double3x3 {
        func triad(_ a: SIMD3<Double>, _ b: SIMD3<Double>) throws -> simd_double3x3 {
            let e1 = simd_normalize(a)
            let rest = b - simd_dot(b, e1) * e1
            guard simd_length(rest) > 1e-9 * simd_length(b), simd_length(a) > 0 else { throw OrientationSymmetryError.invalidDirection }
            let e2 = simd_normalize(rest)
            return simd_double3x3(columns: (e1, e2, simd_cross(e1, e2)))
        }
        return try triad(first.matrix.cartesian(in: matrix), second.matrix.cartesian(in: matrix))
            * triad(first.candidate.cartesian(in: precipitate), second.candidate.cartesian(in: precipitate)).transpose
    }
}
