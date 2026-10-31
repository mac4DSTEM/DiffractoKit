//
//  PatternMeasurement.swift
//  Role: three single-pattern measurements the apps asked for (DiffractoKit 0.10.0, kit asks 1, 3 and 9 of macTEM's diffraction room): the centre of
//        the unscattered beam, an ellipse-corrected radial profile about a centre, and the subpixel centre of a spot. Coordinates are
//        (row, column) as everywhere in the kit. CPU Double; patterns are Float32 and read once.
//
//  The spot refiner and the beam finder are the apps' interims, ported and tested against planted truth (no Python reference computes a
//  single-pattern beam centre or a pattern spot refinement: py4DSTEM's `get_origin` is a per-scan reduction, CuboKit's); the radial profile is the
//  definition (py4DSTEM `radial_integral` after `EllipseTransform`'s correction), checked on a planted ellipse. All `validation: "none"`.
//

import Foundation

/// The centre of the unscattered beam found in one pattern.
public struct PatternCentreFit: Sendable, Equatable {
    public let row: Double
    public let column: Double
    /// The radius of the disc above half height, pixels: sqrt(area / pi).
    public let beamRadiusPixels: Double
    /// Pixels above half height.
    public let areaPixels: Int
    /// "none": no standard has validated the finder.
    public var validation: String { "none" }
}

public enum PatternCentre {
    /// The unscattered beam's centre: the intensity-weighted centroid of the connected region around the brightest (box-smoothed) pixel that lies above the
    /// midpoint between the pattern's median and that maximum. Nil when there is no beam to find: no finite pixel, a maximum less than
    /// `minimumContrastMADs` median absolute deviations (at least 1 grey level) above the median, or a region that reaches the pattern's edge.
    /// A pattern whose beam is blocked (a beam stop) is not centred on its brightest feature, and the finder does not detect it: it refuses, or
    /// returns the centroid of the brightest Bragg spot or of the stop's rim, pixels from the centre. Use `FriedelCentre.find` there (with the stop's mask
    /// when there is one), or a ring fit about a start centre, `RingCalibration.fitRing`.
    public static func find(pattern: DiffractionPattern, minimumContrastMADs: Double = 20) -> PatternCentreFit? {
        let w = pattern.width, h = pattern.height
        guard w >= 8, h >= 8, pattern.pixels.count == w * h, minimumContrastMADs.isFinite, minimumContrastMADs >= 0 else { return nil }
        var sample = pattern.pixels.filter { $0.isFinite }
        guard sample.count > 16 else { return nil }
        sample.sort()
        let median = Double(sample[sample.count / 2])
        let mad = Double(sample.map { abs($0 - Float(median)) }.sorted()[sample.count / 2])
        // 3 x 3 box mean, so one hot pixel is not a beam
        var smooth = [Double](repeating: -Double.infinity, count: w * h)
        var best = -Double.infinity, bestAt = -1
        for y in 1 ..< h - 1 {
            for x in 1 ..< w - 1 {
                var sum = 0.0
                var ok = true
                for dy in -1 ... 1 { for dx in -1 ... 1 {
                    let v = pattern.pixels[(y + dy) * w + x + dx]
                    if !v.isFinite { ok = false } else { sum += Double(v) }
                } }
                guard ok else { continue }
                smooth[y * w + x] = sum / 9
                if sum / 9 > best { best = sum / 9; bestAt = y * w + x }
            }
        }
        guard bestAt >= 0, best - median >= minimumContrastMADs * max(mad, 1) else { return nil }
        let threshold = median + 0.5 * (best - median)
        var inRegion = [Bool](repeating: false, count: w * h)
        var stack = [bestAt]
        inRegion[bestAt] = true
        var members: [Int] = []
        while let p = stack.popLast() {
            members.append(p)
            let y = p / w, x = p % w
            if x <= 1 || y <= 1 || x >= w - 2 || y >= h - 2 { return nil }          // reaches the edge (the smoothed image stops one pixel short of it)
            for q in [p - 1, p + 1, p - w, p + w] where !inRegion[q] && smooth[q] > threshold {
                inRegion[q] = true; stack.append(q)
            }
        }
        var sum = 0.0, sr = 0.0, sc = 0.0
        for p in members {
            let weight = Double(pattern.pixels[p]) - threshold
            guard weight > 0 else { continue }
            sum += weight; sr += weight * Double(p / w); sc += weight * Double(p % w)
        }
        guard sum > 0 else { return nil }
        return PatternCentreFit(row: sr / sum, column: sc / sum, beamRadiusPixels: (Double(members.count) / Double.pi).squareRoot(), areaPixels: members.count)
    }
}

/// A radial profile: the mean intensity at each radius about a centre.
public struct PatternRadialProfile: Sendable, Equatable {
    /// Bin centres, pixels (bin k covers [k, k + 1) times the bin width).
    public let radiusPixels: [Double]
    /// Mean of the finite pixels in the bin; NaN for an empty bin.
    public let mean: [Double]
    public let count: [Int]
    public let binWidthPixels: Double
    public var validation: String { "none" }

    /// The profile of the pixels about (`centerRow`, `centerColumn`) by distance after the correction of `ellipse` (the `EllipseTransform` of a fitted
    /// ring: the radius a circle of the same scale would have), up to `maximumRadius` (default: the nearest edge). Without an ellipse the distance is the
    /// plain one, and the bins and `maximumRadius` are in corrected radius (the default limit is the nearest edge in plain pixels). Returns nil for an empty
    /// pattern, a bin width that is not positive, a centre outside the pattern, or an `ellipse` that is not one (a, b finite and positive, theta finite):
    /// a profile is never silently uncorrected. `maximumRadius` is capped at the pattern's diagonal.
    public static func about(pattern: DiffractionPattern, centerRow: Double, centerColumn: Double,
                             ellipse: (a: Double, b: Double, theta: Double)? = nil, binWidth: Double = 1, maximumRadius: Double? = nil) -> PatternRadialProfile? {
        let w = pattern.width, h = pattern.height
        guard w > 0, h > 0, pattern.pixels.count == w * h, binWidth.isFinite, binWidth > 0, centerRow.isFinite, centerColumn.isFinite,
              centerRow >= 0, centerRow <= Double(h - 1), centerColumn >= 0, centerColumn <= Double(w - 1) else { return nil }
        var limit = maximumRadius ?? min(centerRow, centerColumn, Double(h - 1) - centerRow, Double(w - 1) - centerColumn)
        guard limit.isFinite, limit > 0 else { return nil }
        let transform = ellipse.flatMap { EllipseTransform(a: $0.a, b: $0.b, theta: $0.theta) }
        if ellipse != nil && transform == nil { return nil }
        // with a strongly elongated ellipse a corrected radius can exceed the plain diagonal; the cap keeps the allocation bounded
        limit = min(limit, 4 * Double(w + h))
        let bins = Int((limit / binWidth).rounded(.up))
        var sum = [Double](repeating: 0, count: bins), count = [Int](repeating: 0, count: bins)
        for y in 0 ..< h {
            for x in 0 ..< w {
                let v = pattern.pixels[y * w + x]
                guard v.isFinite else { continue }
                let dx = Double(x) - centerColumn, dy = Double(y) - centerRow
                let r: Double
                if let transform {
                    let c = transform.correctedOffset(dx: Float(dx), dy: Float(dy))
                    r = (Double(c.x) * Double(c.x) + Double(c.y) * Double(c.y)).squareRoot()
                } else { r = (dx * dx + dy * dy).squareRoot() }
                let k = Int(r / binWidth)
                guard k < bins else { continue }
                sum[k] += Double(v); count[k] += 1
            }
        }
        return PatternRadialProfile(radiusPixels: (0 ..< bins).map { (Double($0) + 0.5) * binWidth },
                                    mean: (0 ..< bins).map { count[$0] > 0 ? sum[$0] / Double(count[$0]) : .nan },
                                    count: count, binWidthPixels: binWidth)
    }
}

public enum PatternSpot {
    /// The subpixel centre of the spot nearest a click: the intensity-weighted centroid of the (2 half + 1) square window about it, the window's minimum
    /// removed as background, recentred and repeated three times. The window is clamped to the pattern. Nil when a click lies further than `half`
    /// outside the pattern, the clamped window is under two pixels wide, it holds no finite pixel or no signal above its minimum. The centroid of a
    /// symmetric spot is unbiased while the window holds the spot and no neighbour; a window that cuts the spot's tail (or is clamped by the pattern's
    /// edge) biases it toward the click or the pattern (`PatternMeasurementTests` pins sigma 2 with half 5 to 0.05 px; a smaller half than about 2.5 sigma is worse).
    public static func refine(pattern: DiffractionPattern, row: Double, column: Double, half: Int = 5) -> (row: Double, column: Double)? {
        let w = pattern.width, h = pattern.height
        guard w > 0, h > 0, pattern.pixels.count == w * h, half >= 1, half < 10_000, row.isFinite, column.isFinite,
              row >= -Double(half), row <= Double(h - 1 + half), column >= -Double(half), column <= Double(w - 1 + half) else { return nil }
        var r = row, c = column
        for _ in 0 ..< 3 {
            let r0 = max(Int(r.rounded()) - half, 0), r1 = min(Int(r.rounded()) + half, h - 1)
            let c0 = max(Int(c.rounded()) - half, 0), c1 = min(Int(c.rounded()) + half, w - 1)
            guard r0 < r1, c0 < c1 else { return nil }
            var low = Float.greatestFiniteMagnitude
            for y in r0 ... r1 { for x in c0 ... c1 where pattern.pixels[y * w + x].isFinite { low = min(low, pattern.pixels[y * w + x]) } }
            var sum = 0.0, sr = 0.0, sc = 0.0
            for y in r0 ... r1 {
                for x in c0 ... c1 {
                    let v = pattern.pixels[y * w + x]
                    guard v.isFinite else { continue }
                    let weight = Double(v - low)
                    sum += weight; sr += weight * Double(y); sc += weight * Double(x)
                }
            }
            guard sum > 0 else { return nil }
            r = sr / sum; c = sc / sum
        }
        return (r, c)
    }
}
