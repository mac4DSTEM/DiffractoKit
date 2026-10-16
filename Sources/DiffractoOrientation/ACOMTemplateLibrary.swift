//
//  ACOMTemplateLibrary.swift
//  Role: the ACOM template library of the simplified polar-correlation route (the numbers of the apps' `OrientationPlan`): for a
//        crystal's kinematic reflections, sample the beam directions of the fundamental zone of its Laue class, project the
//        excited reflections of each to a polar (radius x azimuth) image, and keep that image's per-ring azimuthal FFT.
//        The matcher is `ACOMMatcher` (ACOMMatcher.swift).
//
//  This is NOT py4DSTEM's `crystal_ACOM`: there is no SO(3) correlogram, no shell grouping and no sub-pixel refinement. The
//  principle is the same (a polar template per beam direction, the in-plane rotation recovered by a 1D FFT along the azimuth),
//  and docs/parity.md (P20 in its ACOM row) says which quantities are compared with py4DSTEM and which are not. The route is
//  held to the apps' output; `ACOMExcitationModel.py4DSTEM` is this kit's option that follows py4DSTEM's excitation error.
//
//  What a library holds is Float, as the apps' plan: the images and their FFTs are Float32, the deposition weights Double.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import ImagoFFT

/// A reflection as the template library needs it: a reciprocal-lattice vector and the intensity it scatters with.
public struct ACOMReflection: Sendable, Equatable {
    /// The reciprocal-lattice vector in Å⁻¹ (no 2π), in the crystal's Cartesian frame.
    public var g: SIMD3<Double>
    /// The kinematic intensity |F|² (Å⁻⁴), before the excitation envelope.
    public var intensity: Double

    public init(g: SIMD3<Double>, intensity: Double) {
        self.g = g
        self.intensity = intensity
    }
}

/// How a reflection's excitation error s is computed for a beam direction `n` (a unit vector in the crystal's Cartesian frame).
public enum ACOMExcitationModel: Sendable, Equatable {
    /// s = g·n: a flat Ewald sphere. g and -g are then excited alike, every template is exactly π-periodic in azimuth and the
    /// in-plane angle is only determined modulo π.
    case flat
    /// s = g·n + λ|g|²/2: the apps' curved sphere (their frame's z axis is the beam direction itself).
    case appCurvature(wavelengthAngstrom: Double)
    /// s = (2 g·n − λ|g|²) / (2 − 2 λ g·n): py4DSTEM f050d207 `Crystal.excitation_errors` for `n` the third column of the
    /// orientation matrix. Use this model to match patterns that py4DSTEM generated or py4DSTEM's conventions describe.
    case py4DSTEM(wavelengthAngstrom: Double)
}

/// The settings of a template library (the apps' `OrientationPlan.generate` defaults).
public struct ACOMTemplateSettings: Sendable, Equatable {
    /// Radial bins of the polar image: bin `i` is at radius `i * kMax / nRadial` (Å⁻¹). Default 32.
    public var nRadial: Int
    /// Azimuthal bins, a power of two (the FFT length). Default 128, a step of 2.8°.
    public var nAzimuthal: Int
    /// Beam directions sampled in the fundamental zone (at most 4096). Default 600.
    public var zoneAxisCount: Int
    /// The width w of the excitation envelope exp(−(s/w)²), Å⁻¹ (py4DSTEM's σ√2). Default 0.03.
    public var sgWidth: Double
    /// Reflections with |s| above this (Å⁻¹) are not excited. Default 0.1.
    public var sgMax: Double
    /// σ of the circular Gaussian blur along the azimuth, in bins. Default 1.5; 0 for none.
    public var azimuthalBlurBins: Double
    /// The exponent applied to the excited intensity before deposition: py4DSTEM's `power_intensity`, 0.25. The matcher applies
    /// the same exponent to the measured intensities.
    public var intensityPower: Double
    /// Each template spot is also multiplied by its radius to this power (py4DSTEM's `power_radial`, which is 1 there). The
    /// apps measured their own matcher to be better without it and use 0, so it is 0 here.
    public var radialPower: Double
    /// The radial Gaussian σ (Å⁻¹) a spot is spread over: py4DSTEM's `corr_kernel_size`, 0.08. 0 deposits at the nearest bin.
    public var radialKernel: Double
    /// Templates closer than this (degrees) to the winner are not alternatives for the runner-up that the reliability is
    /// measured against. Default 10.
    public var distinctOrientationDegrees: Double
    /// The excitation error model. Default `.flat`, as the apps' default.
    public var excitation: ACOMExcitationModel

    public init(nRadial: Int = 32, nAzimuthal: Int = 128, zoneAxisCount: Int = 600, sgWidth: Double = 0.03,
                sgMax: Double = 0.1, azimuthalBlurBins: Double = 1.5, intensityPower: Double = 0.25,
                radialPower: Double = 0, radialKernel: Double = 0.08, distinctOrientationDegrees: Double = 10,
                excitation: ACOMExcitationModel = .flat) {
        self.nRadial = nRadial
        self.nAzimuthal = nAzimuthal
        self.zoneAxisCount = zoneAxisCount
        self.sgWidth = sgWidth
        self.sgMax = sgMax
        self.azimuthalBlurBins = azimuthalBlurBins
        self.intensityPower = intensityPower
        self.radialPower = radialPower
        self.radialKernel = radialKernel
        self.distinctOrientationDegrees = distinctOrientationDegrees
        self.excitation = excitation
    }
}

/// Why a library or a match was refused.
public enum ACOMError: Error, Sendable, Equatable, CustomStringConvertible {
    /// A setting is outside its range; the text names it.
    case invalidSettings(String)
    /// The reflection list is empty (or holds only g = 0) or holds a non-finite value.
    case noUsableReflections
    /// The structure's reflection list was refused (an element outside the scattering table, an invalid cell, ...).
    case kinematics(String)
    /// The structure has no Laue class (its symmetry could not be determined), so the fundamental zone is undefined.
    case noLaueClass
    /// The zone sampling was refused (the count is outside 0...4096).
    case zoneSampling(OrientationSymmetryError)

    public var description: String {
        switch self {
        case .invalidSettings(let s): return "Invalid template settings: \(s)."
        case .noUsableReflections: return "There is no finite reflection to build templates from."
        case .kinematics(let s): return "The structure's reflections were refused: \(s)"
        case .noLaueClass: return "The structure has no Laue class, so there is no fundamental zone to sample."
        case .zoneSampling(let e): return "The beam directions could not be sampled: \(e)."
        }
    }
}

/// One projected reflection of a template: radius (Å⁻¹) and azimuth (radians) in the template's polar frame, and the weight
/// deposited. Retained so a matched template can be drawn back onto the pattern.
public struct ACOMTemplateSpot: Sendable, Equatable {
    public let radius: Float
    public let azimuth: Float
    public let weight: Float

    public init(radius: Float, azimuth: Float, weight: Float) {
        self.radius = radius
        self.azimuth = azimuth
        self.weight = weight
    }
}

/// A detected diffraction peak: detector position (pixels) and intensity.
public struct ACOMPeak: Sendable, Equatable {
    public var x: Float
    public var y: Float
    public var intensity: Float

    public init(x: Float, y: Float, intensity: Float) {
        self.x = x
        self.y = y
        self.intensity = intensity
    }
}

/// The ACOM template library of the apps' polar-correlation route, with every template's azimuthal ring FFT precomputed.
///
/// For each sampled beam direction `n` of the fundamental zone of the crystal's Laue class (``CrystalOrientationSymmetry``):
/// the reflections with |s| ≤ `sgMax` are projected on the plane perpendicular to `n` in the detector frame of
/// ``ZoneAxisFrame/detectorBasis(zoneAxis:)``, weighted by (|F|² exp(−(s/w)²))^`intensityPower`, deposited on a polar grid
/// (a Gaussian of σ = `radialKernel` across the radial bins, then a circular Gaussian blur along the azimuth), each ring's mean
/// subtracted and the image scaled to unit norm. Because the in-plane rotation is recovered by the azimuthal FFT at match time,
/// only the two degrees of freedom of the beam direction are sampled.
///
/// Where it differs from py4DSTEM's `orientation_plan` (all measured or read from its source, see docs/parity.md): py4DSTEM
/// keeps the image's mean (its mean-subtraction lines are commented out) and groups reflections into shells of equal |g| where
/// this route uses a uniform radial grid; it samples an angular grid of the zone where this route samples a Fibonacci sphere
/// folded into the sector; `radialPower` is 1 there and 0 here (the apps measured the omission to be better for them).
public struct ACOMTemplateLibrary: Sendable {
    /// The polar grid: `nRadial` rings of `nAzimuthal` bins, `radialScale` Å⁻¹ per ring.
    public struct Geometry: Sendable, Equatable {
        public let nRadial: Int
        public let nAzimuthal: Int
        public let radialScale: Double

        public init(nRadial: Int, nAzimuthal: Int, radialScale: Double) {
            self.nRadial = nRadial
            self.nAzimuthal = nAzimuthal
            self.radialScale = radialScale
        }
    }

    public let geometry: Geometry
    public let laueClass: LaueClass
    public let settings: ACOMTemplateSettings
    /// The beam direction of each template, unit vectors in the crystal's Cartesian frame.
    public let zoneAxes: [SIMD3<Double>]
    /// The detector frame of each template: columns detector x, detector y and the beam, in crystal coordinates.
    public let detectorBases: [simd_double3x3]
    /// Per template: the spots its image was built from.
    public let templateSpots: [[ACOMTemplateSpot]]
    /// Per template: the real polar image, `nRadial * nAzimuthal` Floats, ring major, ring-mean-free and of unit norm.
    public let templates: [[Float]]
    /// Every template's per-ring azimuthal FFT, `[template][ring][azimuth]` contiguous (forward, unnormalised).
    public let templateFFTRe: [Float]
    public let templateFFTIm: [Float]
    /// The angle (radians) under which two beam directions count as one orientation for the runner-up.
    public let distinctOrientationRadians: Double

    public var count: Int { zoneAxes.count }

    // MARK: Generation

    /// Builds the library from reflections.
    ///
    /// - Parameters:
    ///   - reflections: the reflections to project (`g` ≠ 0; a zero vector has no azimuth and is skipped).
    ///   - kMax: the radius (Å⁻¹) of the outermost ring; reflections beyond it are not deposited.
    ///   - laueClass: the Laue class whose fundamental zone the beam directions sample.
    ///   - settings: see ``ACOMTemplateSettings``.
    /// - Throws: ``ACOMError``.
    public static func generate(reflections: [ACOMReflection], kMax: Double, laueClass: LaueClass,
                                settings: ACOMTemplateSettings = ACOMTemplateSettings()) throws(ACOMError)
        -> ACOMTemplateLibrary {
        try validate(settings, kMax: kMax)
        guard reflections.allSatisfy({ $0.g.x.isFinite && $0.g.y.isFinite && $0.g.z.isFinite && $0.intensity.isFinite })
        else { throw .noUsableReflections }
        let usable = reflections.filter { simd_length($0.g) > 0 }
        guard !usable.isEmpty else { throw .noUsableReflections }
        guard let fft = FFT1D(n: settings.nAzimuthal) else {
            throw .invalidSettings("nAzimuthal \(settings.nAzimuthal) is not a power of two")
        }
        let axes: [SIMD3<Double>]
        do {
            axes = try CrystalOrientationSymmetry(laueClass: laueClass).sampleFundamentalZone(count: settings.zoneAxisCount)
        } catch let error as OrientationSymmetryError { throw .zoneSampling(error) }
        catch { throw .invalidSettings("\(error)") }
        let geo = Geometry(nRadial: settings.nRadial, nAzimuthal: settings.nAzimuthal,
                           radialScale: kMax / Double(settings.nRadial))
        var bases: [simd_double3x3] = []
        var allSpots: [[ACOMTemplateSpot]] = []
        var templates: [[Float]] = []
        var fftRe: [Float] = [], fftIm: [Float] = []
        let ringCount = geo.nRadial * geo.nAzimuthal
        fftRe.reserveCapacity(axes.count * ringCount)
        fftIm.reserveCapacity(axes.count * ringCount)
        for n in axes {
            let spots = project(reflections: usable, zoneAxis: n, settings: settings)
            var polar = buildPolar(spots: spots, geometry: geo, azimuthalBlurBins: settings.azimuthalBlurBins,
                                   radialKernel: settings.radialKernel)
            normalizeUnit(&polar)
            let (re, im) = ringFFTs(polar, geometry: geo, fft: fft)
            bases.append(ZoneAxisFrame.detectorBasis(zoneAxis: n))
            allSpots.append(spots.map { ACOMTemplateSpot(radius: Float($0.r), azimuth: Float($0.azimuth), weight: Float($0.weight)) })
            templates.append(polar)
            fftRe.append(contentsOf: re)
            fftIm.append(contentsOf: im)
        }
        return ACOMTemplateLibrary(geometry: geo, laueClass: laueClass, settings: settings, zoneAxes: axes,
                                   detectorBases: bases, templateSpots: allSpots, templates: templates,
                                   templateFFTRe: fftRe, templateFFTIm: fftIm,
                                   distinctOrientationRadians: settings.distinctOrientationDegrees * .pi / 180)
    }

    /// Builds the library of a crystal: the kit's kinematic reflections up to `kMax` (after the |F| cut) and the crystal's own
    /// Laue class.
    ///
    /// - Throws: ``ACOMError`` (`.kinematics` for the kit's refusal of the structure, `.noLaueClass` when its symmetry is
    ///   undetermined).
    public static func generate(crystal: CrystalStructure, kMax: Double,
                                settings: ACOMTemplateSettings = ACOMTemplateSettings()) throws(ACOMError)
        -> ACOMTemplateLibrary {
        guard let laue = crystal.laueClass else { throw .noLaueClass }
        let list: KinematicReflectionList
        do { list = try crystal.kinematicReflections(kMax: kMax, retainingAll: false) }
        catch { throw .kinematics(error.description) }
        let reflections = list.kept.map { ACOMReflection(g: $0.reflection.g, intensity: $0.structureFactor.intensity) }
        return try generate(reflections: reflections, kMax: kMax, laueClass: laue, settings: settings)
    }

    private static func validate(_ s: ACOMTemplateSettings, kMax: Double) throws(ACOMError) {
        func bad(_ text: String) -> ACOMError { .invalidSettings(text) }
        guard kMax.isFinite, kMax > 0 else { throw bad("kMax must be a positive finite number") }
        guard s.nRadial >= 1, s.nRadial <= 4096 else { throw bad("nRadial must be in 1...4096") }
        guard s.nAzimuthal >= 4, s.nAzimuthal <= 4096, s.nAzimuthal & (s.nAzimuthal - 1) == 0 else {
            throw bad("nAzimuthal must be a power of two in 4...4096")
        }
        guard (0...4096).contains(s.zoneAxisCount) else { throw bad("zoneAxisCount must be in 0...4096") }
        guard s.sgWidth.isFinite, s.sgWidth > 0, s.sgMax.isFinite, s.sgMax >= 0 else {
            throw bad("sgWidth must be positive and sgMax non-negative")
        }
        guard s.azimuthalBlurBins.isFinite, s.azimuthalBlurBins >= 0, s.radialKernel.isFinite, s.radialKernel >= 0,
              s.intensityPower.isFinite, s.intensityPower > 0, s.radialPower.isFinite,
              s.distinctOrientationDegrees.isFinite, s.distinctOrientationDegrees >= 0 else {
            throw bad("blur, kernel, powers and the distinct angle must be finite, the powers' exponent positive")
        }
        switch s.excitation {
        case .flat: break
        case .appCurvature(let l), .py4DSTEM(let l):
            guard l.isFinite, l > 0 else { throw bad("the wavelength must be positive and finite") }
        }
    }

    // MARK: Projection

    /// One projected spot (Double): radius, azimuth in (−π, π] and weight.
    package typealias Spot = (r: Double, azimuth: Double, weight: Double)

    /// Projects reflections for the beam direction `n` (unit vector): the excited ones (|s| ≤ `sgMax`) map to (radius, azimuth)
    /// in the plane perpendicular to the beam, in the detector frame of ``ZoneAxisFrame/detectorBasis(zoneAxis:)``, with the
    /// weight (|F|² exp(−(s/w)²))^`intensityPower` (times radius^`radialPower`).
    ///
    /// The excitation model decides s. With `.flat`, s = g·n is odd under g → −g while the envelope is even, so every template
    /// is exactly π-periodic in azimuth and the in-plane angle is only defined modulo π; either curved model breaks the tie.
    /// The two curved models differ in the sign of the curvature term: the apps' frame has the beam along +z, py4DSTEM's wave
    /// vector along −z (see ``ACOMExcitationModel``).
    public static func project(reflections: [ACOMReflection], zoneAxis n: SIMD3<Double>,
                               settings: ACOMTemplateSettings) -> [(r: Double, azimuth: Double, weight: Double)] {
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: n)
        let e1 = basis.columns.0, e2 = basis.columns.1
        var spots: [Spot] = []
        spots.reserveCapacity(reflections.count)
        for refl in reflections {
            let gz = simd_dot(refl.g, n)
            let g2 = simd_length_squared(refl.g)
            let sg: Double
            switch settings.excitation {
            case .flat: sg = gz
            case .appCurvature(let l): sg = gz + (l / 2) * g2
            case .py4DSTEM(let l): sg = (2 * gz - l * g2) / (2 - 2 * l * gz)
            }
            if abs(sg) > settings.sgMax { continue }
            let excited = refl.intensity * exp(-(sg / settings.sgWidth) * (sg / settings.sgWidth))
            var w = settings.intensityPower == 1 ? excited : pow(excited, settings.intensityPower)
            let x = simd_dot(refl.g, e1), y = simd_dot(refl.g, e2)
            let radius = (x * x + y * y).squareRoot()
            if settings.radialPower != 0 { w *= pow(radius, settings.radialPower) }
            spots.append((r: radius, azimuth: atan2(y, x), weight: w))
        }
        return spots
    }

    // MARK: Polar image (shared with the matcher)

    /// Deposits spots into a `[nRadial × nAzimuthal]` image: each spot goes to its nearest azimuthal bin and is spread over the
    /// radial bins within ±3σ (σ = `radialKernel` / `radialScale`, at least one bin) by a Gaussian of its unrounded radius (nearest
    /// bin when `radialKernel` is 0); then a circular Gaussian blur along the azimuth (σ = `azimuthalBlurBins`, radius
    /// ⌈3σ⌋ bins); then each ring's mean is subtracted (py4DSTEM keeps the mean; this route removes the DC term of every ring).
    package static func buildPolar(spots: [Spot], geometry geo: Geometry, azimuthalBlurBins: Double,
                                   radialKernel: Double) -> [Float] {
        var img = [Float](repeating: 0, count: geo.nRadial * geo.nAzimuthal)
        let twoPi = 2 * Double.pi
        let sigmaBins = radialKernel / geo.radialScale
        let reach = sigmaBins > 0 ? max(1, Int((3 * sigmaBins).rounded())) : 0
        for s in spots {
            let exact = s.r / geo.radialScale
            let rb = Int(exact.rounded())
            if rb + reach < 0 || rb - reach >= geo.nRadial { continue }
            var a = s.azimuth.truncatingRemainder(dividingBy: twoPi)
            if a < 0 { a += twoPi }
            let ab = Int((a / twoPi * Double(geo.nAzimuthal)).rounded()) % geo.nAzimuthal
            if reach == 0 {
                if rb >= 0 && rb < geo.nRadial { img[rb * geo.nAzimuthal + ab] += Float(s.weight) }
                continue
            }
            for bin in (rb - reach)...(rb + reach) where bin >= 0 && bin < geo.nRadial {
                let d = (Double(bin) - exact) / sigmaBins
                img[bin * geo.nAzimuthal + ab] += Float(s.weight * exp(-0.5 * d * d))
            }
        }
        if azimuthalBlurBins > 0 { blurAzimuthal(&img, geometry: geo, sigmaBins: azimuthalBlurBins) }
        for r in 0..<geo.nRadial {
            let base = r * geo.nAzimuthal
            var mean: Float = 0
            for a in 0..<geo.nAzimuthal { mean += img[base + a] }
            mean /= Float(geo.nAzimuthal)
            for a in 0..<geo.nAzimuthal { img[base + a] -= mean }
        }
        return img
    }

    /// The forward FFT along the azimuth of every ring, ring major.
    package static func ringFFTs(_ polar: [Float], geometry geo: Geometry, fft: FFT1D) -> (re: [Float], im: [Float]) {
        var outRe = [Float](repeating: 0, count: geo.nRadial * geo.nAzimuthal)
        var outIm = [Float](repeating: 0, count: geo.nRadial * geo.nAzimuthal)
        var re = [Float](repeating: 0, count: geo.nAzimuthal)
        var im = [Float](repeating: 0, count: geo.nAzimuthal)
        for r in 0..<geo.nRadial {
            let base = r * geo.nAzimuthal
            for a in 0..<geo.nAzimuthal { re[a] = polar[base + a]; im[a] = 0 }
            fft.transform(re: &re, im: &im, forward: true)
            for a in 0..<geo.nAzimuthal { outRe[base + a] = re[a]; outIm[base + a] = im[a] }
        }
        return (outRe, outIm)
    }

    private static func blurAzimuthal(_ img: inout [Float], geometry geo: Geometry, sigmaBins: Double) {
        let radius = max(1, Int((3 * sigmaBins).rounded()))
        var taps = [Float](repeating: 0, count: 2 * radius + 1)
        var s: Float = 0
        for i in -radius...radius {
            let v = exp(-Float(i * i) / Float(2 * sigmaBins * sigmaBins))
            taps[i + radius] = v
            s += v
        }
        for i in 0..<taps.count { taps[i] /= s }
        let na = geo.nAzimuthal
        var ring = [Float](repeating: 0, count: na)
        for r in 0..<geo.nRadial {
            let base = r * na
            for a in 0..<na {
                var acc: Float = 0
                for t in -radius...radius {
                    let aa = ((a + t) % na + na) % na
                    acc += img[base + aa] * taps[t + radius]
                }
                ring[a] = acc
            }
            for a in 0..<na { img[base + a] = ring[a] }
        }
    }

    package static func normalizeUnit(_ img: inout [Float]) {
        var norm: Float = 0
        for v in img { norm += v * v }
        norm = norm.squareRoot()
        if norm > 0 { for i in 0..<img.count { img[i] /= norm } }
    }

    /// The orientation matrix (columns: lab x, y, z in crystal coordinates, py4DSTEM's convention) of a template at an in-plane
    /// angle: the template's detector basis times the in-plane rotation. For a mirrored match the detector-y and beam columns
    /// are negated (py4DSTEM's `Orientation.mirror`, a half turn about x), still a proper rotation. This is the matrix before
    /// the symmetry reduction that `OrientationResult.euler` carries.
    public func orientationMatrix(templateIndex: Int, inPlaneAngle: Float, mirrored: Bool) -> simd_double3x3 {
        var m = ZoneAxisFrame.matrix(detectorBasis: detectorBases[templateIndex], inPlaneAngle: inPlaneAngle)
        if mirrored {
            m.columns.1 = -m.columns.1
            m.columns.2 = -m.columns.2
        }
        return m
    }
}
