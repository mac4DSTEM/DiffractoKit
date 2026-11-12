//  MatrixSpotDistortion.swift
//  The detector's distortion measured on a spot pattern (DiffractoKit#96).
//
//  A scan of one crystal has no ring to fit an ellipse to. With the matrix structure and its orientation per position known
//  (`MatrixOrientationField`), the matrix's predicted spots are the reference instead: the least-squares 2 × 2 matrix that
//  carries the predicted spots onto the measured vectors, over all positions, holds the detector's stretch, its axis and the
//  error of the scale.
//
//  Method: ordinary least squares of measured = A · predicted over the paired spots (the normal equations of a linear map
//  in the plane), with the pairs made twice: first against the predicted spots after the per-position in-plane fit
//  (`fittingInPlane`), then against the same spots carried by the first A, so that spots the stretch had moved out of reach
//  are paired too. The stretch is read from A's singular value decomposition (Golub and Van Loan, Matrix Computations,
//  4th ed., section 2.4). No upstream function: py4DSTEM fits an ellipse to a ring (`fit_ellipse_1D`,
//  `fit_ellipse_amorphous_ring`), which a single-crystal pattern does not have.

import Foundation
import simd
import DiffractoCore

/// The detector's distortion as fitted on the matrix's spots: measured vector = `affine` · predicted vector, in Å⁻¹, in the
/// detector's axes (x to the right, y down, as the peaks are given).
public struct MatrixSpotDistortion: Sendable, Equatable {
    /// The fitted linear map from predicted to measured vectors.
    public let affine: simd_double2x2
    /// Pairs of a predicted spot and a measured vector the fit used.
    public let pairCount: Int
    /// RMS distance of the measured vectors from the predicted spots carried by `affine`, Å⁻¹.
    public let residualInvAngstrom: Double
    /// The same without the stretch: from the predicted spots carried by the best rotation and scale alone, Å⁻¹.
    public let residualWithoutStretchInvAngstrom: Double

    /// The two singular values of `affine`, the larger first: how much a predicted length is stretched along the long and
    /// the short axis. Both 1 for a detector without distortion and with the right scale.
    public var principalScales: (long: Double, short: Double) {
        let s = Self.singularValues(affine)
        return (s.sigma1, s.sigma2)
    }

    /// Long axis over short axis, at least 1: 1.05 for a pattern stretched by 5 %.
    public var stretchRatio: Double { let s = principalScales; return s.short > 0 ? s.long / s.short : .infinity }

    /// The direction of the long axis in the detector's axes, radians in [0, π), from x towards y.
    public var stretchAxisRad: Double {
        let angle = Self.singularValues(affine).longAxisRad
        let folded = angle.truncatingRemainder(dividingBy: .pi)
        return folded < 0 ? folded + .pi : folded
    }

    /// The geometric mean of the two scales: the factor by which the Å⁻¹ per pixel in use is too small (above 1) or too
    /// large (below 1) for the matrix's lattice.
    public var meanScale: Double { let s = principalScales; return (s.long * s.short).squareRoot() }

    /// The map that takes the stretch out of a measured vector and leaves rotation and mean scale as they are:
    /// √(long · short) · S⁻¹, with S the symmetric stretch of `affine`'s polar decomposition on the measured side.
    public var removingStretch: simd_double2x2 {
        let s = Self.singularValues(affine)
        let c = cos(s.longAxisRad), sn = sin(s.longAxisRad)
        let u = simd_double2x2(columns: (SIMD2(c, sn), SIMD2(-sn, c)))
        let mean = (s.sigma1 * s.sigma2).squareRoot()
        let inverse = simd_double2x2(diagonal: SIMD2(mean / s.sigma1, mean / s.sigma2))
        return u * inverse * u.transpose
    }

    /// Singular values of a 2 × 2 matrix and the direction of the first left singular vector (the long axis on the
    /// measured side), closed form.
    static func singularValues(_ m: simd_double2x2) -> (sigma1: Double, sigma2: Double, longAxisRad: Double) {
        // M Mᵀ is symmetric positive semidefinite; its eigenvalues are the squared singular values and its first
        // eigenvector the long axis on the measured side.
        let a = m * m.transpose
        let p = a[0][0], q = a[1][0], r = a[1][1]          // [[p, q], [q, r]]
        let half = (p + r) / 2, diff = (p - r) / 2
        let root = (diff * diff + q * q).squareRoot()
        let l1 = max(half + root, 0), l2 = max(half - root, 0)
        let angle = root > 0 ? 0.5 * atan2(2 * q, p - r) : 0
        return (l1.squareRoot(), l2.squareRoot(), angle)
    }
}

extension PhaseVectorMatcher {

    /// The detector's distortion from the matrix's spots, over a scan.
    ///
    /// At every position with a matrix orientation the matrix's predicted spots are laid on the measured vectors by the
    /// in-plane fit of the matrix removal, and every usable measured vector is paired with the nearest predicted spot
    /// within `pairRadiusInvAngstrom`. One linear map is fitted to all pairs, the pairs are made again against the
    /// predicted spots carried by that map, and it is fitted once more.
    ///
    /// - Parameters:
    ///   - scan: the peaks of every position, in detector pixels.
    ///   - field: the matrix orientation of every position (`matrixOrientations`), from the same scan.
    ///   - reflections: the matrix crystal's reflections (the list its template library was built from).
    ///   - pairRadiusInvAngstrom: how far a measured vector may be from a predicted spot to be paired with it. It has to
    ///     hold the distortion at the longest vectors: 0.06 Å⁻¹ holds a stretch of 8 % at 0.75 Å⁻¹.
    ///   - minimumPairs: fewer pairs give no fit.
    /// - Returns: nil when the scan and the field do not have the same number of positions, with fewer than
    ///   `minimumPairs` pairs, or when the paired spots do not span two directions (one row of spots fixes no stretch).
    public static func matrixSpotDistortion(scan: [[SIMD2<Double>]], field: MatrixOrientationField,
                                            reflections: [ACOMReflection],
                                            originX: Float, originY: Float, invAngstromPerPixel: Double,
                                            pairRadiusInvAngstrom: Double = 0.06,
                                            excitationSlabInvAngstrom: Double = matrixRemovalSlabInvAngstrom,
                                            directBeamRadiusInvAngstrom: Double = 0.15,
                                            maximumVectorInvAngstrom: Double = 0,
                                            minimumPairs: Int = 50) -> MatrixSpotDistortion? {
        guard scan.count == field.orientations.count, pairRadiusInvAngstrom > 0 else { return nil }
        // Per position: the usable measured vectors and the predicted spots after the in-plane fit.
        var positions: [(vectors: [SIMD2<Double>], spots: [SIMD2<Double>])] = []
        for (i, peaks) in scan.enumerated() {
            guard let orientation = field.orientations[i] else { continue }
            let spots = matrixSpots(reflections: reflections, orientation: orientation, excitationSlabInvAngstrom: excitationSlabInvAngstrom)
            let vectors = peaks.compactMap {
                experimentalVector(peak: $0, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                   directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom, maximumVectorInvAngstrom: maximumVectorInvAngstrom)
            }
            guard spots.count >= 2, vectors.count >= 2 else { continue }
            let fitted = fittingInPlane(spots: spots, vectors: vectors, toleranceInvAngstrom: pairRadiusInvAngstrom / 3)
            positions.append((vectors, fitted.spots))
        }
        // Normal equations of measured = A · predicted over the pairs, the predicted spots carried by `carry` for the pairing.
        func fit(carry: simd_double2x2) -> simd_double2x2? {
            var us = simd_double2x2(0), ss = simd_double2x2(0)
            var pairs = 0
            for position in positions {
                for u in position.vectors {
                    var best: SIMD2<Double>?, bestDistance = Double.infinity
                    for s in position.spots { let d = simd_length(u - carry * s); if d < bestDistance { bestDistance = d; best = s } }
                    guard let s = best, bestDistance <= pairRadiusInvAngstrom else { continue }
                    us = us + simd_double2x2(columns: (u * s.x, u * s.y))      // u sᵀ
                    ss = ss + simd_double2x2(columns: (s * s.x, s * s.y))      // s sᵀ
                    pairs += 1
                }
            }
            guard pairs >= max(2, minimumPairs) else { return nil }
            // Two directions: s sᵀ summed must be well conditioned (the smaller eigenvalue at least 1 % of the larger).
            let eigen = MatrixSpotDistortion.singularValues(ss)
            guard eigen.sigma1 > 0, eigen.sigma2 >= 0.01 * eigen.sigma1 else { return nil }
            return us * ss.inverse
        }
        guard let first = fit(carry: matrix_identity_double2x2), let affine = fit(carry: first) else { return nil }
        // Residuals over the final pairs: with the affine, and with the best rotation and scale alone (a similarity:
        // measured ≈ k R predicted, least squares in closed form).
        var sum = 0.0, sumSimilarity = 0.0, pairs = 0
        var dot = 0.0, cross = 0.0, norm = 0.0
        var paired: [(u: SIMD2<Double>, s: SIMD2<Double>)] = []
        for position in positions {
            for u in position.vectors {
                var best: SIMD2<Double>?, bestDistance = Double.infinity
                for s in position.spots { let d = simd_length(u - affine * s); if d < bestDistance { bestDistance = d; best = s } }
                guard let s = best, bestDistance <= pairRadiusInvAngstrom else { continue }
                paired.append((u, s))
                sum += bestDistance * bestDistance
                dot += s.x * u.x + s.y * u.y
                cross += s.x * u.y - s.y * u.x
                norm += simd_length_squared(s)
                pairs += 1
            }
        }
        guard pairs >= max(2, minimumPairs), norm > 0 else { return nil }
        let a = dot / norm, b = cross / norm                       // k cos, k sin
        for pair in paired {
            let carried = SIMD2(a * pair.s.x - b * pair.s.y, b * pair.s.x + a * pair.s.y)
            sumSimilarity += simd_length_squared(pair.u - carried)
        }
        return MatrixSpotDistortion(affine: affine, pairCount: pairs,
                                    residualInvAngstrom: (sum / Double(pairs)).squareRoot(),
                                    residualWithoutStretchInvAngstrom: (sumSimilarity / Double(pairs)).squareRoot())
    }
}
