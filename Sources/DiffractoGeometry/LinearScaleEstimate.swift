//
//  LinearScaleEstimate.swift
//  Role: the linear reciprocal scale Δq (Å⁻¹ per pixel) from the innermost reflection of a known crystal, measured on the peak
//        radii of many scan positions; ported from mac4DSTEM `KnownCrystalQCalibration.estimate` with the app's types
//        (Bragg vectors, origin, probe radius) replaced by plain radii. It REPORTS and refuses nothing: the second-shell
//        ratio is a number beside the crystal's prediction, never a verdict (the app measured that no threshold on it separates
//        a healthy cube from a mis-assigned ring). No camera length is needed, unlike `RingCalibration.q`.
//
//  THE ESTIMATOR (all float64, exact copies of the app's arithmetic):
//    per position: radii above `minimumRadiusPixels`, sorted; the SAME-SHELL CLUSTER is the run from the smallest radius up to
//    (1 + band) times it, band = min(0.08, (r₂/r₁ − 1)/2) of the crystal's two innermost shells (0.08 without a second shell);
//    the position's r₁ is the cluster's mean (a minimum is biased low by the order statistic of k symmetry equivalents);
//    the second shell's r₂ is the mean of the same-shell cluster starting at the first radius beyond (1 + separation) r₁,
//    separation = (expected ratio − 1)/2. Over positions: observed r₁ = median of the r₁s, Δq = reference q / observed r₁,
//    MAD = median |r₁ − observed|, shell ratio = median r₂ / observed.
//

import Foundation

/// What the second shell measured, beside what the crystal predicts. Reports; judges nothing.
public enum LinearScaleShellCheck: Sendable, Equatable {
    /// `positions` is how many positions contributed an r₂.
    case measured(observedRatio: Double, expectedRatio: Double, positions: Int)
    /// Not a pass: the innermost detected peak was assumed to be the innermost allowed reflection.
    case notSelfChecked(String)

    /// |observed / expected − 1|, or nil when nothing was checked. No threshold is applied to it anywhere.
    public var mismatch: Double? {
        guard case .measured(let observed, let expected, _) = self, expected != 0 else { return nil }
        return abs(observed / expected - 1)
    }
}

public struct LinearScaleEstimate: Sendable, Equatable {
    /// Δq, Å⁻¹ per pixel: the reference q over the observed radius.
    public let invAngstromPerPixel: Double
    public let observedRadiusPixels: Double
    public let referenceRadiusInvAngstrom: Double
    public let medianAbsoluteDeviationPixels: Double
    /// Positions that contributed an r₁.
    public let sampleCount: Int
    public let secondShellRadiusPixels: Double?
    public let shellCheck: LinearScaleShellCheck
    /// Median number of same-shell equivalents averaged per position (1: the cluster brought nothing).
    public let sameShellPeaksPerPosition: Double

    /// MAD over the observed radius (infinite for a zero radius). No threshold is applied to it.
    public var shellConsistency: Double {
        observedRadiusPixels > 0 ? medianAbsoluteDeviationPixels / observedRadiusPixels : .infinity
    }
}

public enum LinearScale {

    /// The estimate from per-position peak radii (pixels from the origin; the beam itself is dropped by `minimumRadiusPixels`).
    /// `referenceRadiusInvAngstrom` is the innermost allowed reflection's 1/d, `secondShellRadiusInvAngstrom` the next shell's (nil:
    /// the shell check says `.notSelfChecked`). Nil unless the reference is finite and positive and some position has a radius
    /// above the minimum.
    public static func estimate(
        radiiPerPosition: [[Double]],
        referenceRadiusInvAngstrom: Double,
        secondShellRadiusInvAngstrom: Double?,
        minimumRadiusPixels: Double = 2
    ) -> LinearScaleEstimate? {
        guard referenceRadiusInvAngstrom.isFinite, referenceRadiusInvAngstrom > 0 else { return nil }

        let expectedRatio = secondShellRadiusInvAngstrom.flatMap { second -> Double? in
            let ratio = second / referenceRadiusInvAngstrom
            return ratio.isFinite && ratio > 1 ? ratio : nil
        }
        let separation = expectedRatio.map { ($0 - 1) / 2 }
        let sameShellBand = min(0.08, separation ?? 0.08)

        func shellMean(from radii: ArraySlice<Double>) -> (mean: Double, count: Int)? {
            guard let first = radii.first else { return nil }
            let cluster = radii.prefix { $0 <= first * (1 + sameShellBand) }
            return (cluster.reduce(0, +) / Double(cluster.count), cluster.count)
        }

        var firstRadii: [Double] = []
        var secondRadii: [Double] = []
        var equivalentsPerPosition: [Double] = []
        firstRadii.reserveCapacity(radiiPerPosition.count)
        for peaks in radiiPerPosition {
            let radii = peaks.filter { $0.isFinite && $0 > minimumRadiusPixels }.sorted()
            guard let first = shellMean(from: radii[...]) else { continue }
            firstRadii.append(first.mean)
            equivalentsPerPosition.append(Double(first.count))
            if let separation,
               let secondIndex = radii.firstIndex(where: { $0 > first.mean * (1 + separation) }),
               let second = shellMean(from: radii[secondIndex...]) {
                secondRadii.append(second.mean)
            }
        }
        guard !firstRadii.isEmpty else { return nil }
        let observed = median(firstRadii)
        guard observed.isFinite, observed > 0 else { return nil }
        let mad = median(firstRadii.map { abs($0 - observed) })

        let secondObserved = secondRadii.isEmpty ? nil : median(secondRadii)
        let shellCheck: LinearScaleShellCheck
        if let expectedRatio, let secondObserved, secondObserved.isFinite {
            shellCheck = .measured(observedRatio: secondObserved / observed, expectedRatio: expectedRatio,
                                   positions: secondRadii.count)
        } else if expectedRatio == nil {
            shellCheck = .notSelfChecked("this crystal has only one allowed shell within range")
        } else {
            shellCheck = .notSelfChecked("only one shell is detectable")
        }
        return LinearScaleEstimate(
            invAngstromPerPixel: referenceRadiusInvAngstrom / observed,
            observedRadiusPixels: observed,
            referenceRadiusInvAngstrom: referenceRadiusInvAngstrom,
            medianAbsoluteDeviationPixels: mad,
            sampleCount: firstRadii.count,
            secondShellRadiusPixels: secondObserved,
            shellCheck: shellCheck,
            sameShellPeaksPerPosition: median(equivalentsPerPosition))
    }

    private static func median(_ values: [Double]) -> Double {
        let sorted = values.sorted()
        let middle = sorted.count / 2
        return sorted.count.isMultiple(of: 2) ? (sorted[middle - 1] + sorted[middle]) / 2 : sorted[middle]
    }
}
