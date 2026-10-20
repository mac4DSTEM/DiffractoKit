//
//  ACOMMatcher.swift
//  Role: ACOM orientation matching on the CPU (the numbers of the apps' `OrientationMatcher`): a pattern's detected Bragg peaks
//        become a polar image that is correlated with every template of an `ACOMTemplateLibrary`, which gives the best beam
//        direction, its in-plane rotation, the score and the runner-up.
//
//          corr(Δφ) = IFFT_φ( Σ_ring conj(FFT(exp_ring)) · FFT(template_ring) )
//
//  The peak of corr over Δφ is the in-plane angle, its height the score. A second pass on the unconjugated product correlates
//  the azimuth-reversed image, covering the mirror half of the proper-rotation zone that the sampled zone does not reach
//  (py4DSTEM `inversion_symmetry`, `Orientation.mirror`). The scan loop, selection, progress and cancellation are in ACOMScan.swift; the
//  Metal option (`ACOMBackend.metal`) is in ACOMMetal.swift; the CPU path in this file is the reference.
//

import Foundation
import simd
import DiffractoCore
import ImagoFFT

/// Matches measured peaks against an ``ACOMTemplateLibrary``. A value: it holds the library and an immutable FFT plan, and
/// every call owns its scratch, so one matcher can serve many threads.
public struct ACOMMatcher: Sendable {
    public let library: ACOMTemplateLibrary
    private let fft: FFT1D
    private let symmetry: CrystalOrientationSymmetry

    /// - Returns: nil when the library's azimuthal length has no FFT plan (it cannot happen for a library `generate` made).
    public init?(library: ACOMTemplateLibrary) {
        guard let fft = FFT1D(n: library.geometry.nAzimuthal) else { return nil }
        self.library = library
        self.fft = fft
        self.symmetry = CrystalOrientationSymmetry(laueClass: library.laueClass)
    }

    /// The least number of peaks a pattern needs.
    public static let minimumPeakCount = 2

    /// Matches one pattern.
    ///
    /// - Parameters:
    ///   - peaks: detected peaks in detector pixels.
    ///   - originX: the pattern centre, pixels (the unscattered beam).
    ///   - originY: likewise.
    ///   - invAngstromPerPixel: the Q pixel size; the peak radii are converted with it to the library's Å⁻¹.
    /// - Returns: `OrientationResult.empty` for fewer than two peaks or an empty library. Otherwise `score` is the best
    ///   normalised correlation (clamped at 0), `secondScore` the best score among templates farther than the library's
    ///   distinct angle from the winner (0 when there is none), `euler` the symmetry-reduced orientation, `inPlaneAngle` the
    ///   rotation before reduction and `mirrored` whether the conjugated pass won.
    public func match(peaks: [ACOMPeak], originX: Float, originY: Float, invAngstromPerPixel: Double = 1)
        -> OrientationResult {
        guard let experimental = experimentalFFT(peaks: peaks, originX: originX, originY: originY,
                                                 invAngstromPerPixel: invAngstromPerPixel) else { return .empty }
        let all = correlate(experimentalRe: experimental.re, experimentalIm: experimental.im)
        let selected = Self.selectOrientation(zoneAxes: library.zoneAxes, scores: all.scores, bins: all.bins,
                                              distinctOrientationRadians: library.distinctOrientationRadians)
        guard selected.template >= 0 else { return .empty }
        return makeResult(bestTemplate: selected.template, bestBin: selected.bin, bestScore: selected.score,
                          secondScore: selected.secondScore, mirrored: all.mirrored[selected.template])
    }

    /// Every template's best score for one pattern, in library order: the larger of the forward and the mirrored pass.
    /// Empty when `match` would return `.empty` for lack of peaks.
    public func templateScores(peaks: [ACOMPeak], originX: Float, originY: Float, invAngstromPerPixel: Double = 1) -> [Float] {
        guard let experimental = experimentalFFT(peaks: peaks, originX: originX, originY: originY,
                                                 invAngstromPerPixel: invAngstromPerPixel) else { return [] }
        return correlate(experimentalRe: experimental.re, experimentalIm: experimental.im).scores
    }

    /// The measured pattern's polar image after the deposition, blur, per-ring mean subtraction and unit normalisation of the
    /// templates (the same exponent, kernel and grid): the image whose ring FFTs are correlated. nil for fewer than two peaks.
    public func experimentalPolarImage(peaks: [ACOMPeak], originX: Float, originY: Float,
                                       invAngstromPerPixel: Double = 1) -> [Float]? {
        guard peaks.count >= Self.minimumPeakCount else { return nil }
        let power = library.settings.intensityPower
        var spots: [ACOMTemplateLibrary.Spot] = []
        spots.reserveCapacity(peaks.count)
        for p in peaks {
            let dx = Double(p.x - originX), dy = Double(p.y - originY)
            let r = (dx * dx + dy * dy).squareRoot() * invAngstromPerPixel
            let intensity = Double(max(p.intensity, 0))
            spots.append((r: r, azimuth: atan2(dy, dx), weight: power == 1 ? intensity : pow(intensity, power)))
        }
        var polar = ACOMTemplateLibrary.buildPolar(spots: spots, geometry: library.geometry, azimuthalBlurBins: 1.5,
                                                   radialKernel: library.settings.radialKernel)
        ACOMTemplateLibrary.normalizeUnit(&polar)
        return polar
    }

    /// The forward ring FFTs of the measured pattern's polar image (ring major), what the correlation multiplies the templates by.
    public func experimentalFFT(peaks: [ACOMPeak], originX: Float, originY: Float, invAngstromPerPixel: Double = 1)
        -> (re: [Float], im: [Float])? {
        guard let polar = experimentalPolarImage(peaks: peaks, originX: originX, originY: originY,
                                                 invAngstromPerPixel: invAngstromPerPixel), library.count > 0
        else { return nil }
        return ACOMTemplateLibrary.ringFFTs(polar, geometry: library.geometry, fft: fft)
    }

    /// Both passes for every template: per template the larger score (the mirrored pass only when strictly larger), its
    /// azimuthal bin and which pass won.
    private func correlate(experimentalRe: [Float], experimentalIm: [Float])
        -> (scores: [Float], bins: [UInt32], mirrored: [Bool]) {
        let na = library.geometry.nAzimuthal, nr = library.geometry.nRadial
        let invScale = 1 / Float(na)    // the FFT's inverse is unnormalised
        let count = library.count
        var scores = [Float](repeating: 0, count: count)
        var bins = [UInt32](repeating: 0, count: count)
        var mirrored = [Bool](repeating: false, count: count)
        var sumRe = [Float](repeating: 0, count: na)
        var sumIm = [Float](repeating: 0, count: na)
        for pass in 0..<2 {
            // pass 0: conj(E)·T, the forward correlation; pass 1: E·T, the correlation of the azimuth-reversed image.
            let imaginarySign: Float = pass == 0 ? -1 : 1
            for t in 0..<count {
                let offset = t * nr * na
                for a in 0..<na { sumRe[a] = 0; sumIm[a] = 0 }
                for r in 0..<nr {
                    let base = r * na
                    for a in 0..<na {
                        let er = experimentalRe[base + a], eIm = imaginarySign * experimentalIm[base + a]
                        let tr = library.templateFFTRe[offset + base + a], ti = library.templateFFTIm[offset + base + a]
                        sumRe[a] += er * tr - eIm * ti
                        sumIm[a] += er * ti + eIm * tr
                    }
                }
                fft.transform(re: &sumRe, im: &sumIm, forward: false)
                var best: Float = -.greatestFiniteMagnitude
                var bestBin = 0
                for a in 0..<na where sumRe[a] > best { best = sumRe[a]; bestBin = a }
                best *= invScale
                if pass == 0 || best > scores[t] {
                    scores[t] = best
                    bins[t] = UInt32(bestBin)
                    mirrored[t] = pass == 1
                }
            }
        }
        return (scores, bins, mirrored)
    }

    /// The winner and the best *distinct* runner-up. `1 − second/best` is a confidence only if "second" is a different
    /// orientation: on a dense library the runner-up over all templates is the neighbouring beam direction, whose image is
    /// near-identical to the winner's. Templates within the distinct angle of the winner are therefore skipped. (py4DSTEM's
    /// `match_single_pattern` zeroes the correlation near each taken match, `min_angle_between_matches_deg`, and repeats; this
    /// is the same rule applied once, for one runner-up.)
    package static func selectOrientation(zoneAxes: [SIMD3<Double>], scores: [Float], bins: [UInt32],
                                          distinctOrientationRadians: Double)
        -> (template: Int, bin: Int, score: Float, secondScore: Float) {
        var best: Float = -.greatestFiniteMagnitude
        var bestTemplate = -1
        for t in 0..<scores.count where scores[t] > best { best = scores[t]; bestTemplate = t }
        guard bestTemplate >= 0 else { return (-1, 0, 0, 0) }
        let winner = zoneAxes[bestTemplate]
        let cosLimit = cos(distinctOrientationRadians)
        var second: Float = -.greatestFiniteMagnitude
        for t in 0..<scores.count where t != bestTemplate {
            if simd_dot(zoneAxes[t], winner) > cosLimit { continue }
            if scores[t] > second { second = scores[t] }
        }
        if second == -.greatestFiniteMagnitude { second = 0 }
        return (bestTemplate, Int(bins[bestTemplate]), best, second)
    }

    func makeResult(bestTemplate: Int, bestBin: Int, bestScore: Float, secondScore: Float, mirrored: Bool)
        -> OrientationResult {
        let na = library.geometry.nAzimuthal
        // The correlation peaks at −θ; negate to report the pattern rotation. A mirrored match carries py4DSTEM's extra half
        // turn (crystal_ACOM.py: `+ np.pi`) and the beam flip of `orientationMatrix(...mirrored:)`.
        let bin = (na - bestBin) % na
        var angle = Float(Double(bin) / Double(na) * 2 * Double.pi)
        if mirrored { angle = (angle + .pi).truncatingRemainder(dividingBy: 2 * .pi) }
        let matrix = library.orientationMatrix(templateIndex: bestTemplate, inPlaneAngle: angle, mirrored: mirrored)
        // The matrix is a proper rotation by construction; the fallback keeps `match` total.
        let reduced = (try? symmetry.reduce(matrix)) ?? (matrix: matrix, disorientationRad: Float(0))
        return OrientationResult(templateIndex: bestTemplate,
                                 euler: EulerAngles(py4DSTEMOrientationMatrix: reduced.matrix),
                                 inPlaneAngle: angle, score: max(bestScore, 0), secondScore: max(secondScore, 0), phaseID: 0,
                                 symmetryDisorientationRad: reduced.disorientationRad, mirrored: mirrored)
    }
}
