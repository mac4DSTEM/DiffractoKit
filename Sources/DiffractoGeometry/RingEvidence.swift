//
//  RingEvidence.swift
//  Role: what the data say about a ring at a fitted ellipse, read sector by sector and without a model in the azimuth:
//        whether each of the 36 sectors holds the ring, whether the ring is a peak, and whether the sectors put it where
//        the fit says. DiffractoKit 0.1.0 (docs/registration-0.0.4.md §1.3), the kit's answer to slice 3b's Gate B gaps F1
//        (noise in an angular gap counted as coverage) and F2 (a disk edge, a blob, a textured circle calibrated).
//
//  THE MEASUREMENT. Every finite, non-impulse annulus pixel (the annulus about the caller's start centre, as the fit read
//  it) gets its offset δ = ρR − R from the fitted ellipse (ρ the profile model's normalised elliptical radius, R = (a + b)/2)
//  and its azimuth about the fitted centre. Pixels with −(3 w_in + 2) ≤ δ ≤ 3 w_out + 2 are read, w the fitted ring sigmas
//  (a conic has none: the half-maximum sigma of its seed profile, else (outer − inner)/8).
//  - Per sector (36), a linear least-squares fit y = α g(δ) + τ g′(δ) + β + γ δ, g the fitted asymmetric-Gaussian ring of
//    unit peak: α is the sector's ring amplitude with its standard error, t = −τ/α its radial offset from the fit (first
//    order), β + γδ a local line. A sector needs 8 pixels.
//  - A = the median α over the sectors that have 8 pixels. A sector holds evidence when α ≥ 0.5 A and α ≥ 3 SE(α).
//  - The peak: the profile pooled over all read pixels in 0.5 px bins of δ falls from its mean over |δ| ≤ 0.5 by at least
//    0.25 A on each side (inner minimum over −(3 w_in + 2) ≤ δ ≤ −w_in, outer over w_out ≤ δ ≤ 3 w_out + 2; bins of at least
//    3 pixels).
//  - The geometry: over the sectors with evidence, t = c0 + c1 cos φ + s1 sin φ + c2 cos 2φ + s2 sin 2φ weighted by
//    1/(SE(t)² + 0.01²), its covariance scaled by max(1, χ²/(n − 5)); the centre shift √(c1² + s1²) and the ellipse shift
//    √(c2² + s2²) each fail beyond 0.25 px AND 4 standard errors.
//  - One ring (added after Gate B on the 0.1.0 candidate, F1: a neighbouring ring inside the annulus was fitted, or fitted between the
//    two): the profile pooled over the WHOLE annulus in 0.5 px bins of δ, smoothed over 3 bins, must not hold a second
//    peak more than 1 px from the fitted radius: a bin that rises above the lowest bin between it and the ring AND above the lowest bin
//    between it and the annulus's edge, by 0.25 A or more on both sides. A tail rising to the edge (a central beam) is
//    not a peak.
//  All numbers are stored without NaN (nil where undefined), so two evidences compare with `==`.
//

import Foundation

/// The ring evidence of one fit (`EllipseCalibration.ringEvidence`). `passes` is what a calibration requires.
public struct RingEvidence: Sendable, Equatable {
    /// The median sector ring amplitude A, in the pattern's units (0 when no sector could be read).
    public let amplitude: Double
    /// Per sector (36, by azimuth about the fitted centre from the row axis towards the column axis): the ring amplitude
    /// α and its standard error (nil with fewer than 8 pixels or a singular fit), and the radial offset t in pixels of
    /// the data's ring from the fitted ellipse (nil unless the sector holds evidence).
    public let sectorAmplitudes: [Double?]
    public let sectorAmplitudeErrors: [Double?]
    public let sectorShifts: [Double?]
    /// Sectors that hold the ring: α ≥ 0.5 A and α ≥ 3 SE(α).
    public let occupiedBins: Int
    /// How far the pooled profile falls from the ring's peak on each side, as a fraction of A (nil when a side has no bin).
    public let innerDrop: Double?
    public let outerDrop: Double?
    /// The consensus radial offset c0, the centre shift √(c1² + s1²) and the ellipse shift √(c2² + s2²), pixels, with the
    /// standard errors of the two shifts (nil with fewer than 6 sectors holding evidence or a singular fit).
    public let radiusShift: Double?
    public let centreShift: Double?
    public let centreShiftError: Double?
    public let ellipseShift: Double?
    public let ellipseShiftError: Double?
    /// The highest second peak of the pooled profile more than 1 px from the fitted radius inside the annulus, as a fraction of A (the
    /// smaller of its rises above the lowest bins on either side); nil when no such peak or no bins there.
    public let secondaryPeak: Double?

    public init(amplitude: Double, sectorAmplitudes: [Double?], sectorAmplitudeErrors: [Double?], sectorShifts: [Double?],
                occupiedBins: Int, innerDrop: Double?, outerDrop: Double?, radiusShift: Double?, centreShift: Double?,
                centreShiftError: Double?, ellipseShift: Double?, ellipseShiftError: Double?,
                secondaryPeak: Double? = nil) {
        self.secondaryPeak = secondaryPeak
        self.amplitude = amplitude
        self.sectorAmplitudes = sectorAmplitudes
        self.sectorAmplitudeErrors = sectorAmplitudeErrors
        self.sectorShifts = sectorShifts
        self.occupiedBins = occupiedBins
        self.innerDrop = innerDrop
        self.outerDrop = outerDrop
        self.radiusShift = radiusShift
        self.centreShift = centreShift
        self.centreShiftError = centreShiftError
        self.ellipseShift = ellipseShift
        self.ellipseShiftError = ellipseShiftError
    }

    // MARK: Registered bounds (registration-0.0.4 §1.3), now the fields of `Policy` (0.2.0)

    /// The bounds the ring evidence is measured and judged by (registration-0.0.4 §1.3): the suite's choices, with
    /// `Policy.default` the set every release before 0.2.0 used. An invalid policy (`isValid` false) gives empty evidence
    /// from `EllipseCalibration.ringEvidence` and `false` from every judgement below.
    public struct Policy: Sendable, Equatable {
        /// The fewest pixels a sector needs to be read. Default 8.
        public var minimumSectorPixels: Int
        /// A sector holds the ring when its amplitude is at least this fraction of the median sector amplitude. Default 0.5.
        public var minimumAmplitudeFraction: Double
        /// ... and at least this many standard errors. Default 3.
        public var minimumSignificance: Double
        /// The fewest sectors (of 36) that must hold the ring. Default 30, the same bound as the coverage rule.
        public var minimumOccupiedBins: Int
        /// The ring is a peak when the profile falls by at least this fraction of the amplitude on both sides. Default 0.25.
        public var minimumDrop: Double
        /// The centre and ellipse shifts fail beyond this many pixels ... Default 0.25.
        public var maximumShift: Double
        /// ... and this many standard errors. Default 4.
        public var shiftSignificance: Double
        /// A second peak of this height or more (a fraction of the amplitude) is another ring in the annulus. Default 0.25.
        public var maximumSecondaryPeak: Double

        public init(minimumSectorPixels: Int = RingEvidence.minimumSectorPixels,
                    minimumAmplitudeFraction: Double = RingEvidence.minimumAmplitudeFraction,
                    minimumSignificance: Double = RingEvidence.minimumSignificance,
                    minimumOccupiedBins: Int = RingEvidence.minimumOccupiedBins,
                    minimumDrop: Double = RingEvidence.minimumDrop,
                    maximumShift: Double = RingEvidence.maximumShift,
                    shiftSignificance: Double = RingEvidence.shiftSignificance,
                    maximumSecondaryPeak: Double = RingEvidence.maximumSecondaryPeak) {
            self.minimumSectorPixels = minimumSectorPixels
            self.minimumAmplitudeFraction = minimumAmplitudeFraction
            self.minimumSignificance = minimumSignificance
            self.minimumOccupiedBins = minimumOccupiedBins
            self.minimumDrop = minimumDrop
            self.maximumShift = maximumShift
            self.shiftSignificance = shiftSignificance
            self.maximumSecondaryPeak = maximumSecondaryPeak
        }

        /// The suite's values.
        public static let `default` = Policy()

        /// At least 5 pixels per sector (the local fit has four parameters and needs one more), 0 to 36 occupied
        /// sectors, non-negative finite fractions, shifts and significances.
        public var isValid: Bool {
            let counts: Bool = minimumSectorPixels >= 5 && minimumOccupiedBins >= 0
                && minimumOccupiedBins <= RingEvidence.sectorCount
            let fractions: Bool = minimumAmplitudeFraction.isFinite && minimumAmplitudeFraction >= 0
                && minimumDrop.isFinite && minimumDrop >= 0
                && maximumSecondaryPeak.isFinite && maximumSecondaryPeak >= 0
            let shifts: Bool = maximumShift.isFinite && maximumShift >= 0
                && shiftSignificance.isFinite && shiftSignificance >= 0
                && minimumSignificance.isFinite && minimumSignificance >= 0
            return counts && fractions && shifts
        }
    }

    // The suite's values, which `Policy.default` takes (and which stay public as the documented defaults).
    public static let sectorCount = 36
    public static let minimumSectorPixels = 8
    public static let minimumAmplitudeFraction = 0.5
    public static let minimumSignificance = 3.0
    /// The same bound as the coverage rule: 30 of 36 sectors.
    public static let minimumOccupiedBins = 30
    public static let minimumDrop = 0.25
    public static let maximumShift = 0.25
    public static let shiftSignificance = 4.0
    /// A second peak of this height or more (a fraction of A) is another ring in the annulus.
    public static let maximumSecondaryPeak = 0.25

    /// Both drops at least `minimumDrop`.
    public var isPeak: Bool { isPeak(.default) }

    /// Both drops at least the policy's `minimumDrop`.
    public func isPeak(_ policy: Policy) -> Bool {
        guard policy.isValid, let innerDrop, let outerDrop else { return false }
        return innerDrop >= policy.minimumDrop && outerDrop >= policy.minimumDrop
    }

    /// Neither shift beyond `maximumShift` and `shiftSignificance` standard errors; false when the shifts are undefined.
    public var geometryHolds: Bool { geometryHolds(.default) }

    public func geometryHolds(_ policy: Policy) -> Bool {
        guard policy.isValid, let centreShift, let centreShiftError, let ellipseShift, let ellipseShiftError else { return false }
        func exceeds(_ shift: Double, _ error: Double) -> Bool {
            shift > policy.maximumShift && shift > policy.shiftSignificance * error
        }
        return !exceeds(centreShift, centreShiftError) && !exceeds(ellipseShift, ellipseShiftError)
    }

    /// No second peak of `maximumSecondaryPeak` or more in the annulus.
    public var isSingleRing: Bool { isSingleRing(.default) }

    public func isSingleRing(_ policy: Policy) -> Bool {
        policy.isValid && (secondaryPeak ?? 0) < policy.maximumSecondaryPeak
    }

    /// What a calibration requires: 30 sectors holding the ring, a peak, one ring in the annulus, and the geometry.
    public var passes: Bool { passes(.default) }

    public func passes(_ policy: Policy) -> Bool {
        occupiedBins >= policy.minimumOccupiedBins && isPeak(policy) && isSingleRing(policy) && geometryHolds(policy)
    }
}

extension EllipseCalibration {

    /// The ring evidence of `fit` in `pattern`, over the annulus (innerRadius, outerRadius] about the start centre
    /// (row `centerQX`, column `centerQY`) the fit was given. A measurement: it never throws and never changes the fit.
    public static func ringEvidence(pattern: DiffractionPattern, fit: EllipseCalibrationFit,
                                    centerQX: Double, centerQY: Double,
                                    innerRadius: Double, outerRadius: Double,
                                    policy: RingEvidence.Policy = .default,
                                    impulseRangeFactor: Double = 2) -> RingEvidence {
        let sectors = RingEvidence.sectorCount
        let empty = RingEvidence(amplitude: 0, sectorAmplitudes: Array(repeating: nil, count: sectors),
                                 sectorAmplitudeErrors: Array(repeating: nil, count: sectors),
                                 sectorShifts: Array(repeating: nil, count: sectors), occupiedBins: 0,
                                 innerDrop: nil, outerDrop: nil, radiusShift: nil, centreShift: nil,
                                 centreShiftError: nil, ellipseShift: nil, ellipseShiftError: nil)
        let a = fit.a, b = fit.b
        guard policy.isValid, impulseRangeFactor.isFinite, impulseRangeFactor >= 0 else { return empty }
        guard patternSizeIsConsistent(pattern), pattern.width > 2, pattern.height > 2,
              a.isFinite, b.isFinite, a > 0, b > 0, fit.theta.isFinite, fit.centerQX.isFinite, fit.centerQY.isFinite,
              centerQX.isFinite, centerQY.isFinite, innerRadius >= 0, outerRadius > innerRadius,
              outerRadiusIsSane(outerRadius, pattern) else { return empty }

        let annulus = annulusPixels(pattern, centerQX: centerQX, centerQY: centerQY, innerRadius: innerRadius,
                                    outerRadius: outerRadius, includeOuter: true,
                                    impulseFactor: impulseRangeFactor).pixels
        guard !annulus.isEmpty else { return empty }

        // The ring's widths: the profile's sigmas, or for a conic the seed's half-maximum sigma.
        let widthIn: Double, widthOut: Double
        if let profile = fit.profile {
            widthIn = profile.innerSigma
            widthOut = profile.outerSigma
        } else {
            let w = halfMaximumRingSigma(values: annulus.map(\.value), xs: annulus.map { Double($0.qx) },
                                         ys: annulus.map { Double($0.qy) }, seed: fit, innerRadius: innerRadius,
                                         outerRadius: outerRadius) ?? (outerRadius - innerRadius) / 8
            widthIn = w
            widthOut = w
        }
        guard widthIn.isFinite, widthOut.isFinite, widthIn > 0, widthOut > 0 else { return empty }
        let lowerOffset = -(3 * widthIn + 2), upperOffset = 3 * widthOut + 2

        // The profile model's elliptical coordinate (EllipseCalibration.fitAmorphousRing's `predictions`).
        let sinTheta = sin(fit.theta), cosTheta = cos(fit.theta)
        let a2 = a * a, b2 = b * b
        let radius = (a + b) / 2
        let radius2 = radius * radius
        let coefficientA = (sinTheta * sinTheta / b2 + cosTheta * cosTheta / a2) * radius2
        let coefficientB = 2 * (b2 - a2) * sinTheta * cosTheta / (a2 * b2) * radius2
        let coefficientC = (cosTheta * cosTheta / b2 + sinTheta * sinTheta / a2) * radius2

        struct Read { let sector: Int; let offset: Double; let value: Double }
        var reads: [Read] = []
        var wholeBins: [Int: (sum: Double, count: Int)] = [:]
        for pixel in annulus {
            let dx = Double(pixel.qx) - fit.centerQX, dy = Double(pixel.qy) - fit.centerQY
            let rho2 = max(0, coefficientA * dx * dx + coefficientB * dx * dy + coefficientC * dy * dy)
            let offset = rho2.squareRoot() - radius
            let wholeBin = Int((offset / 0.5).rounded(.down))
            wholeBins[wholeBin, default: (0, 0)].sum += pixel.value
            wholeBins[wholeBin, default: (0, 0)].count += 1
            guard offset >= lowerOffset, offset <= upperOffset else { continue }
            var angle = atan2(dy, dx)
            if angle < 0 { angle += 2 * .pi }
            let sector = min(sectors - 1, Int(angle / (2 * .pi) * Double(sectors)))
            reads.append(Read(sector: sector, offset: offset, value: pixel.value))
        }

        func shape(_ offset: Double) -> (g: Double, slope: Double) {
            let w = offset < 0 ? widthIn : widthOut
            let g = exp(-offset * offset / (2 * w * w))
            return (g, -offset / (w * w) * g)
        }

        // Per sector: y = α g + τ g′ + β + γ δ.
        var bySector = [[Read]](repeating: [], count: sectors)
        for read in reads { bySector[read.sector].append(read) }
        var alphas = [Double?](repeating: nil, count: sectors)
        var alphaErrors = [Double?](repeating: nil, count: sectors)
        var shifts = [Double?](repeating: nil, count: sectors)
        var shiftErrors = [Double?](repeating: nil, count: sectors)
        for sector in 0..<sectors {
            let rows = bySector[sector]
            guard rows.count >= policy.minimumSectorPixels else { continue }
            let design = rows.map { row -> [Double] in
                let s = shape(row.offset)
                return [s.g, s.slope, 1, row.offset]
            }
            guard let fitted = leastSquares(design, rows.map(\.value)) else { continue }
            let alpha = fitted.coefficients[0], tau = fitted.coefficients[1]
            let varianceAlpha = fitted.covariance[0][0], varianceTau = fitted.covariance[1][1]
            let covariance = fitted.covariance[0][1]
            guard alpha.isFinite, varianceAlpha.isFinite, varianceAlpha >= 0 else { continue }
            alphas[sector] = alpha
            alphaErrors[sector] = varianceAlpha.squareRoot()
            if alpha != 0 {
                let t = -tau / alpha
                let varianceT = (varianceTau + t * t * varianceAlpha + 2 * t * covariance) / (alpha * alpha)
                if t.isFinite, varianceT.isFinite {
                    shifts[sector] = t
                    shiftErrors[sector] = max(0, varianceT).squareRoot()
                }
            }
        }

        let readable = alphas.compactMap { $0 }
        guard !readable.isEmpty, let amplitude = median(of: readable), amplitude > 0 else {
            return RingEvidence(amplitude: 0, sectorAmplitudes: alphas, sectorAmplitudeErrors: alphaErrors,
                                sectorShifts: Array(repeating: nil, count: sectors), occupiedBins: 0, innerDrop: nil,
                                outerDrop: nil, radiusShift: nil, centreShift: nil, centreShiftError: nil,
                                ellipseShift: nil, ellipseShiftError: nil)
        }
        var occupied = [Bool](repeating: false, count: sectors)
        for sector in 0..<sectors {
            guard let alpha = alphas[sector], let error = alphaErrors[sector] else { continue }
            occupied[sector] = alpha >= policy.minimumAmplitudeFraction * amplitude
                && alpha >= policy.minimumSignificance * error
        }
        let occupiedBins = occupied.filter { $0 }.count

        // The peak, on the pooled profile.
        let binWidth = 0.5
        var binSums: [Int: (sum: Double, count: Int)] = [:]
        var centreSum = 0.0, centreCount = 0
        for read in reads {
            let bin = Int((read.offset / binWidth).rounded(.down))
            binSums[bin, default: (0, 0)].sum += read.value
            binSums[bin, default: (0, 0)].count += 1
            if abs(read.offset) <= 0.5 { centreSum += read.value; centreCount += 1 }
        }
        var innerDrop: Double?, outerDrop: Double?
        if centreCount > 0 {
            let peak = centreSum / Double(centreCount)
            var innerMinimum: Double?, outerMinimum: Double?
            for (bin, entry) in binSums where entry.count >= 3 {
                let middle = (Double(bin) + 0.5) * binWidth
                let mean = entry.sum / Double(entry.count)
                if middle >= lowerOffset, middle <= -widthIn { innerMinimum = min(innerMinimum ?? mean, mean) }
                if middle >= widthOut, middle <= upperOffset { outerMinimum = min(outerMinimum ?? mean, mean) }
            }
            innerDrop = innerMinimum.map { (peak - $0) / amplitude }.flatMap { $0.isFinite ? $0 : nil }
            outerDrop = outerMinimum.map { (peak - $0) / amplitude }.flatMap { $0.isFinite ? $0 : nil }
        }

        // One ring: a second peak of the whole annulus's pooled profile beyond the ring's ±w.
        var secondaryPeak: Double?
        let keys = wholeBins.keys.filter { wholeBins[$0]!.count >= 3 }.sorted()
        if keys.count >= 3 {
            let means = keys.map { wholeBins[$0]!.sum / Double(wholeBins[$0]!.count) }
            let smooth = means.indices.map { i -> Double in
                let window = means[max(0, i - 1)...min(means.count - 1, i + 1)]
                return window.reduce(0, +) / Double(window.count)
            }
            let middles = keys.map { (Double($0) + 0.5) * 0.5 }
            func height(_ side: [Int]) -> Double? {   // indices ordered from the ring outwards to the annulus's edge
                var best: Double?
                for (position, i) in side.enumerated() {
                    let towardRing = side[..<position].map { smooth[$0] }.min()
                    let towardEdge = side[(position + 1)...].map { smooth[$0] }.min()
                    guard let towardRing, let towardEdge else { continue }
                    let rise = min(smooth[i] - towardRing, smooth[i] - towardEdge) / amplitude
                    if rise > 0, rise.isFinite { best = max(best ?? rise, rise) }
                }
                return best
            }
            // From 1 px off the fitted radius: a fit can swallow its neighbour into a wide sigma, and a flank that only
            // falls (or a tail that rises to the edge) never meets the two-sided rise, so starting close is safe.
            let outer = keys.indices.filter { middles[$0] > 1 }
            let inner = keys.indices.filter { middles[$0] < -1 }.reversed()
            secondaryPeak = [height(Array(outer)), height(Array(inner))].compactMap { $0 }.max()
        }

        // The geometry: harmonics 0, 1, 2 of the offsets of the sectors that hold the ring.
        var harmonicRows: [[Double]] = [], harmonicValues: [Double] = [], weights: [Double] = []
        var keptShifts = [Double?](repeating: nil, count: sectors)
        for sector in 0..<sectors where occupied[sector] {
            guard let t = shifts[sector], let error = shiftErrors[sector] else { continue }
            keptShifts[sector] = t
            let phi = (Double(sector) + 0.5) / Double(sectors) * 2 * .pi
            harmonicRows.append([1, cos(phi), sin(phi), cos(2 * phi), sin(2 * phi)])
            harmonicValues.append(t)
            weights.append(1 / (error * error + 0.01 * 0.01))
        }
        var radiusShift: Double?, centreShift: Double?, centreShiftError: Double?
        var ellipseShift: Double?, ellipseShiftError: Double?
        if harmonicRows.count >= 6,
           let harmonic = weightedLeastSquares(harmonicRows, harmonicValues, weights) {
            let c = harmonic.coefficients, v = harmonic.covariance
            func amplitudeAndError(_ i: Int, _ j: Int) -> (Double, Double)? {
                let value = (c[i] * c[i] + c[j] * c[j]).squareRoot()
                let variance = value > 0
                    ? (c[i] * c[i] * v[i][i] + c[j] * c[j] * v[j][j] + 2 * c[i] * c[j] * v[i][j]) / (value * value)
                    : max(v[i][i], v[j][j])
                guard value.isFinite, variance.isFinite else { return nil }
                return (value, max(0, variance).squareRoot())
            }
            if c[0].isFinite, let centre = amplitudeAndError(1, 2), let ellipse = amplitudeAndError(3, 4) {
                radiusShift = c[0]
                (centreShift, centreShiftError) = centre
                (ellipseShift, ellipseShiftError) = ellipse
            }
        }
        return RingEvidence(amplitude: amplitude, sectorAmplitudes: alphas, sectorAmplitudeErrors: alphaErrors,
                            sectorShifts: keptShifts, occupiedBins: occupiedBins, innerDrop: innerDrop,
                            outerDrop: outerDrop, radiusShift: radiusShift, centreShift: centreShift,
                            centreShiftError: centreShiftError, ellipseShift: ellipseShift,
                            ellipseShiftError: ellipseShiftError, secondaryPeak: secondaryPeak)
    }

    // MARK: Least squares for the evidence (small systems, normal equations)

    static func median(of values: [Double]) -> Double? {
        guard !values.isEmpty else { return nil }
        let sorted = values.sorted()
        let middle = sorted.count / 2
        return sorted.count.isMultiple(of: 2) ? (sorted[middle - 1] + sorted[middle]) / 2 : sorted[middle]
    }

    /// Ordinary least squares: the coefficients and their covariance s² (XᵀX)⁻¹, s² = RSS/(n − p). Nil when singular or
    /// when n ≤ p.
    static func leastSquares(_ design: [[Double]], _ values: [Double]) -> (coefficients: [Double], covariance: [[Double]])? {
        guard let p = design.first?.count, design.count > p else { return nil }
        guard let fitted = weightedLeastSquares(design, values, Array(repeating: 1, count: values.count),
                                                scaleByChiSquare: false) else { return nil }
        var rss = 0.0
        for (row, value) in zip(design, values) {
            var prediction = 0.0
            for k in 0..<p { prediction += row[k] * fitted.coefficients[k] }
            rss += (value - prediction) * (value - prediction)
        }
        let s2 = rss / Double(design.count - p)
        return (fitted.coefficients, fitted.covariance.map { $0.map { $0 * s2 } })
    }

    /// Weighted least squares: the coefficients and (AᵀWA)⁻¹, scaled by max(1, χ²/(n − p)) unless told not to.
    static func weightedLeastSquares(_ design: [[Double]], _ values: [Double], _ weights: [Double],
                                     scaleByChiSquare: Bool = true) -> (coefficients: [Double], covariance: [[Double]])? {
        guard let p = design.first?.count, design.count > p else { return nil }
        var normal = [[Double]](repeating: [Double](repeating: 0, count: p), count: p)
        var rhs = [Double](repeating: 0, count: p)
        for (index, row) in design.enumerated() {
            let w = weights[index]
            for i in 0..<p {
                rhs[i] += w * row[i] * values[index]
                for j in i..<p { normal[i][j] += w * row[i] * row[j] }
            }
        }
        for i in 0..<p { for j in 0..<i { normal[i][j] = normal[j][i] } }
        // Scale the columns to unit diagonal so `solve`'s absolute pivot floor is not met by small units.
        let column = (0..<p).map { normal[$0][$0] > 0 ? 1 / normal[$0][$0].squareRoot() : 1 }
        let scaled = (0..<p).map { i in (0..<p).map { j in normal[i][j] * column[i] * column[j] } }
        guard let solution = solve(scaled, rhs: (0..<p).map { rhs[$0] * column[$0] }) else { return nil }
        let coefficients = (0..<p).map { solution[$0] * column[$0] }
        var inverse = [[Double]](repeating: [Double](repeating: 0, count: p), count: p)
        for k in 0..<p {
            var unit = [Double](repeating: 0, count: p)
            unit[k] = 1
            guard let x = solve(scaled, rhs: unit) else { return nil }
            for i in 0..<p { inverse[i][k] = x[i] * column[i] * column[k] }
        }
        guard coefficients.allSatisfy(\.isFinite) else { return nil }
        guard scaleByChiSquare else { return (coefficients, inverse) }
        var chiSquare = 0.0
        for (index, row) in design.enumerated() {
            var prediction = 0.0
            for k in 0..<p { prediction += row[k] * coefficients[k] }
            chiSquare += weights[index] * (values[index] - prediction) * (values[index] - prediction)
        }
        let factor = max(1, chiSquare / Double(design.count - p))
        return (coefficients, inverse.map { $0.map { $0 * factor } })
    }
}
