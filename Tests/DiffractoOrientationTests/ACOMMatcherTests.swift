//
//  ACOMMatcherTests.swift
//  Role: the ACOM template library and matcher (docs/parity.md P25). Pins from tools/pins/acom-matcher/pins.py:
//        (1) the apps' route, a NumPy TRANSCRIPTION of the app source (OrientationPlan, OrientationMatcher at mac4DSTEM c809679c):
//            beam directions, template images and spots, polar images, ring FFTs, every template's score, the winner, the
//            runner-up, the in-plane angle and the orientation matrix, for m-3m, 6/mmm and 2/m;
//        (2) py4DSTEM f050d207, run: the spots of a zone axis (flat at lambda -> 0 and curved at 200 kV) against
//            `ACOMTemplateLibrary.project`, and planted orientations against py4DSTEM's own match_single_pattern.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

// MARK: Fixture

private struct Pins: Decodable {
    struct Provenance: Decodable { let py4dstem: String; let py4dstem_lock: String; let numpy: String; let scipy: String }
    struct Template: Decodable {
        let index: Int
        let image: String
        let spots: String
        let ringFFTRe: String?
        let ringFFTIm: String?
    }
    struct Case: Decodable {
        let peaks: String
        let origin: [Double]
        let invAngstromPerPixel: Double
        let trueMatrix: [Double]
        let template: Int
        let bin: Int
        let mirrored: Bool
        let score: Double
        let second: Double
        let angle: Double
        let matrix: [Double]
        let scores: String
        let errorToTruthDeg: Double
        let marginTemplate: Double
        let marginInPlane: Double
        let euler: [Double]?
        let polar: String?
        let ringFFTRe: String?
        let ringFFTIm: String?
    }
    struct App: Decodable {
        let zoneAxisCount: Int
        let axes: String
        let templates: [Template]
        let cases: [Case]
    }
    struct Spot: Decodable {
        let mode: String
        let wavelength: Double
        let axis: [Double]
        let sgMax: Double
        let qx: String
        let qy: String
        let intensity: String
        let hkl: String
    }
    struct Transcription: Decodable {
        let matrix: [Double]
        let score: Double
        let template: Int
        let errorToTruthDeg: Double
        let errorToPy4dstemDeg: Double
    }
    struct Planted: Decodable {
        let trueMatrix: [Double]
        let qx: String
        let qy: String
        let intensity: String
        let py4dstemMatrix: [Double]
        let py4dstemCorr: Double
        let py4dstemErrorToTruthDeg: Double
        let transcriptionPy4dstemModel: Transcription?
    }
    struct Crystal: Decodable {
        let name: String
        let laue: String
        let kmax: Double
        let wavelength: Double
        let g: String
        let intensity: String
        let count: Int
        let appRoute: App?
        let py4dstemSpots: [Spot]
        let planted: [Planted]
    }
    let provenance: Provenance
    let sigma: Double
    let sgMax: Double
    let crystals: [String: Crystal]
}

private func fixture() throws -> Pins {
    let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/acom-matcher-pins.json")
    return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
}

private func decode<T>(_ s: String, as: T.Type) -> [T] {
    let data = Data(base64Encoded: s)!
    return data.withUnsafeBytes { raw in
        (0..<(data.count / MemoryLayout<T>.size)).map { raw.loadUnaligned(fromByteOffset: $0 * MemoryLayout<T>.size, as: T.self) }
    }
}
private func f8(_ s: String) -> [Double] { decode(s, as: Double.self) }
private func f4(_ s: String) -> [Float] { decode(s, as: Float.self) }

private func vectors(_ flat: [Double]) -> [SIMD3<Double>] {
    stride(from: 0, to: flat.count, by: 3).map { SIMD3(flat[$0], flat[$0 + 1], flat[$0 + 2]) }
}
private func matrix(_ columnMajor: [Double]) -> simd_double3x3 {
    simd_double3x3(columns: (SIMD3(columnMajor[0], columnMajor[1], columnMajor[2]),
                             SIMD3(columnMajor[3], columnMajor[4], columnMajor[5]),
                             SIMD3(columnMajor[6], columnMajor[7], columnMajor[8])))
}
private func reflections(_ c: Pins.Crystal) -> [ACOMReflection] {
    zip(vectors(f8(c.g)), f8(c.intensity)).map { ACOMReflection(g: $0, intensity: $1) }
}
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    (0..<3).flatMap { c in (0..<3).map { abs(a[c][$0] - b[c][$0]) } }.max()!
}
private func maxAbs(_ a: [Float], _ b: [Float]) -> Double {
    #expect(a.count == b.count)
    return zip(a, b).map { abs(Double($0) - Double($1)) }.max() ?? 0
}
/// The smallest rotation angle (degrees) between two orientations over the crystal's proper operations.
private func misorientationDegrees(_ a: simd_double3x3, _ b: simd_double3x3, laue: LaueClass) -> Double {
    CrystalOrientationSymmetry(laueClass: laue).operators.map { op in
        let c = (simd_trace(b.transpose * op * a) - 1) / 2
        return acos(min(1, max(-1, c))) * 180 / .pi
    }.min()!
}
private func peaks(_ c: Pins.Case) -> [ACOMPeak] {
    let flat = f4(c.peaks)
    return stride(from: 0, to: flat.count, by: 3).map { ACOMPeak(x: flat[$0], y: flat[$0 + 1], intensity: flat[$0 + 2]) }
}

private let appKeys = ["al", "mg", "mono"]

private func appLibrary(_ key: String, _ crystal: Pins.Crystal) throws -> ACOMTemplateLibrary {
    let app = try #require(crystal.appRoute)
    let settings = ACOMTemplateSettings(zoneAxisCount: app.zoneAxisCount,
                                        excitation: .appCurvature(wavelengthAngstrom: crystal.wavelength))
    return try ACOMTemplateLibrary.generate(reflections: reflections(crystal), kMax: crystal.kmax,
                                            laueClass: try #require(LaueClass(rawValue: crystal.laue)), settings: settings)
}

// MARK: The apps' route (transcription of the app source)

struct ACOMAppRouteTests {
    @Test func provenanceNamesTheOraclesAndTheFixtureIsSmall() throws {
        let p = try fixture()
        #expect(p.provenance.py4dstem == "0.14.19")
        #expect(p.provenance.py4dstem_lock.hasPrefix("f050d207"))
        let size = try FileManager.default.attributesOfItem(
            atPath: URL(fileURLWithPath: #filePath).deletingLastPathComponent()
                .appendingPathComponent("Fixtures/acom-matcher-pins.json").path)[.size] as! Int
        #expect(size < 1_000_000)
    }

    /// The beam directions are the apps' fundamental-zone samples; the templates are their polar images and spots.
    @Test(arguments: appKeys) func zoneAxesTemplatesAndSpotsEqualTheTranscription(key: String) throws {
        let c = try #require(try fixture().crystals[key])
        let lib = try appLibrary(key, c)
        let app = try #require(c.appRoute)
        let axes = vectors(f8(app.axes))
        #expect(lib.zoneAxes.count == axes.count)
        let axisGap = zip(lib.zoneAxes, axes).map { simd_length($0 - $1) }.max()!
        #expect(axisGap < 2e-14)
        var imageGap = 0.0, spotGap = 0.0
        for t in app.templates {
            imageGap = max(imageGap, maxAbs(lib.templates[t.index], f4(t.image)))
            let spots = f8(t.spots)             // n rows of (r, azimuth, weight)
            #expect(lib.templateSpots[t.index].count == spots.count / 3)
            for (i, s) in lib.templateSpots[t.index].enumerated() {
                spotGap = max(spotGap, abs(Double(s.radius) - spots[3 * i]), abs(Double(s.azimuth) - spots[3 * i + 1]),
                              abs(Double(s.weight) - spots[3 * i + 2]) / abs(spots[3 * i + 2]))
            }
            if let re = t.ringFFTRe, let im = t.ringFFTIm {
                let n = lib.geometry.nRadial * lib.geometry.nAzimuthal
                let off = t.index * n
                let fre = Array(lib.templateFFTRe[off..<(off + n)]), fim = Array(lib.templateFFTIm[off..<(off + n)])
                let g = max(maxAbs(fre, f4(re)), maxAbs(fim, f4(im)))
                print("ACOM \(key) template FFT gap \(g)")
                #expect(g < 1e-6)
            }
        }
        print("ACOM \(key) axes gap \(axisGap), template image gap \(imageGap), spot gap \(spotGap)")
        #expect(imageGap < 1e-6)
        #expect(spotGap < 1e-6)
    }

    /// Per planted pattern: the experimental polar image and ring FFT (one pattern), every template's score, the winner, the
    /// runner-up, the in-plane angle, the orientation matrix and (m-3m) the reduced Euler angles.
    @Test(arguments: appKeys) func matchesEqualTheTranscription(key: String) throws {
        let c = try #require(try fixture().crystals[key])
        let lib = try appLibrary(key, c)
        let matcher = try #require(ACOMMatcher(library: lib))
        var worstScores = 0.0, worstScore = 0.0, worstMatrix = 0.0, worstPolar = 0.0
        for (k, cs) in try #require(c.appRoute).cases.enumerated() {
            let pk = peaks(cs)
            let ox = Float(cs.origin[0]), oy = Float(cs.origin[1])
            let r = matcher.match(peaks: pk, originX: ox, originY: oy, invAngstromPerPixel: cs.invAngstromPerPixel)
            // The pinned choices are not near ties (marginTemplate, marginInPlane > 1e-4 in the generator).
            #expect(cs.marginTemplate > 1e-4 && cs.marginInPlane > 1e-4)
            #expect(r.templateIndex == cs.template, "\(key) case \(k)")
            #expect(r.mirrored == cs.mirrored, "\(key) case \(k)")
            #expect(abs(Double(r.inPlaneAngle) - cs.angle) < 2e-6, "\(key) case \(k): bin \(cs.bin)")
            worstScore = max(worstScore, abs(Double(r.score) - cs.score), abs(Double(r.secondScore) - cs.second))
            let scores = matcher.templateScores(peaks: pk, originX: ox, originY: oy, invAngstromPerPixel: cs.invAngstromPerPixel)
            worstScores = max(worstScores, maxAbs(scores, f4(cs.scores)))
            let m = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
            worstMatrix = max(worstMatrix, gap(m, matrix(cs.matrix)))
            // reduction: the result's Euler angles are the symmetry-reduced matrix's
            let ops = CrystalOrientationSymmetry(laueClass: try #require(LaueClass(rawValue: c.laue)))
            let reduced = try ops.reduce(m)
            #expect(r.euler == EulerAngles(py4DSTEMOrientationMatrix: reduced.matrix))
            #expect(r.symmetryDisorientationRad == reduced.disorientationRad)
            if let e = cs.euler {
                let want = [e[0], e[1], e[2]]
                let got = [Double(r.euler.phi1), Double(r.euler.Phi), Double(r.euler.phi2)]
                for (a, b) in zip(got, want) {
                    var d = abs(a - b).truncatingRemainder(dividingBy: 2 * .pi)
                    d = min(d, 2 * .pi - d)
                    #expect(d < 2e-6, "\(key) case \(k) Euler")
                }
            }
            #expect(r.score >= 0 && r.secondScore >= 0 && r.reliability >= 0 && r.reliability <= 1)
            if let polar = cs.polar {
                let image = try #require(matcher.experimentalPolarImage(peaks: pk, originX: ox, originY: oy,
                                                                        invAngstromPerPixel: cs.invAngstromPerPixel))
                worstPolar = max(worstPolar, maxAbs(image, f4(polar)))
                if let re = cs.ringFFTRe, let im = cs.ringFFTIm {
                    let f = try #require(matcher.experimentalFFT(peaks: pk, originX: ox, originY: oy,
                                                                 invAngstromPerPixel: cs.invAngstromPerPixel))
                    let g = max(maxAbs(f.re, f4(re)), maxAbs(f.im, f4(im)))
                    print("ACOM \(key) experimental ring FFT gap \(g)")
                    #expect(g < 1e-6)
                }
            }
        }
        print("ACOM \(key) worst: scores \(worstScores), winner/second \(worstScore), matrix \(worstMatrix), polar \(worstPolar)")
        #expect(worstScores < 1e-6)
        #expect(worstScore < 1e-6)
        #expect(worstMatrix < 2e-6)
        #expect(worstPolar < 1e-6)
    }
}

// MARK: py4DSTEM

struct ACOMPy4DSTEMTests {
    /// py4DSTEM `generate_diffraction_pattern` for orientation_matrix = the kit's detector frame of a zone axis: the spots at
    /// lambda -> 0 (`.flat`) and at 200 kV (`.py4DSTEM(wavelength:)`, py4DSTEM's own excitation error), elementwise.
    @Test(arguments: ["al", "mg", "mono", "rutile"]) func spotsEqualPy4DSTEMsGeneratedPattern(key: String) throws {
        let c = try #require(try fixture().crystals[key])
        let p = try fixture()
        let refl = reflections(c)
        var worst = 0.0, compared = 0
        for s in c.py4dstemSpots {
            let model: ACOMExcitationModel = s.mode == "flat" ? .flat : .py4DSTEM(wavelengthAngstrom: s.wavelength)
            let settings = ACOMTemplateSettings(sgWidth: p.sigma * 2.0.squareRoot(), sgMax: s.sgMax, intensityPower: 1,
                                                excitation: model)
            let spots = ACOMTemplateLibrary.project(reflections: refl, zoneAxis: SIMD3(s.axis[0], s.axis[1], s.axis[2]),
                                                    settings: settings)
            let qx = f8(s.qx), qy = f8(s.qy), inten = f8(s.intensity)
            #expect(spots.count == qx.count, "\(key) \(s.mode) \(s.axis): spot count")
            guard spots.count == qx.count else { continue }
            for (i, sp) in spots.enumerated() {
                let r = hypot(qx[i], qy[i])
                var da = abs(sp.azimuth - atan2(qy[i], qx[i]))
                da = min(da, 2 * .pi - da)
                worst = max(worst, abs(sp.r - r) / r, da, abs(sp.weight - inten[i]) / inten[i])
                compared += 1
            }
        }
        print("ACOM \(key) spots vs py4DSTEM: \(compared) spots, worst relative gap \(worst)")
        #expect(compared > 50)
        #expect(worst < 1e-6)   // py4DSTEM's intensities are float32 (struct_factors_int); everything else is double
    }

    /// Planted orientations: py4DSTEM's own peaks, its own match, and this kit's route (with py4DSTEM's excitation model).
    @Test(arguments: ["al", "mg", "mono", "rutile"]) func plantedOrientationsAgreeWithPy4DSTEMsMatch(key: String) throws {
        let p = try fixture()
        let c = try #require(p.crystals[key])
        let laue = try #require(LaueClass(rawValue: c.laue))
        let count = c.appRoute?.zoneAxisCount ?? 300
        let settings = ACOMTemplateSettings(zoneAxisCount: count, sgWidth: p.sigma * 2.0.squareRoot(),
                                            excitation: .py4DSTEM(wavelengthAngstrom: c.wavelength))
        let lib = try ACOMTemplateLibrary.generate(reflections: reflections(c), kMax: c.kmax, laueClass: laue, settings: settings)
        let matcher = try #require(ACOMMatcher(library: lib))
        var rows: [String] = []
        for (k, pl) in c.planted.enumerated() {
            let qx = f8(pl.qx), qy = f8(pl.qy), inten = f8(pl.intensity)
            let pk = zip(zip(qx, qy), inten).map { ACOMPeak(x: Float($0.0), y: Float($0.1), intensity: Float($1)) }
            let r = matcher.match(peaks: pk, originX: 0, originY: 0)
            let m = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
            let truth = matrix(pl.trueMatrix), py = matrix(pl.py4dstemMatrix)
            let toTruth = misorientationDegrees(m, truth, laue: laue)
            let toPy = misorientationDegrees(m, py, laue: laue)
            let pyToTruth = misorientationDegrees(py, truth, laue: laue)
            rows.append("\(key) \(k): kit->truth \(String(format: "%.2f", toTruth)), kit->py4DSTEM \(String(format: "%.2f", toPy)), "
                        + "py4DSTEM->truth \(String(format: "%.2f", pyToTruth)) (score \(String(format: "%.3f", r.score)), "
                        + "py4DSTEM corr \(String(format: "%.3f", pl.py4dstemCorr)))")
            #expect(abs(pyToTruth - pl.py4dstemErrorToTruthDeg) < 1e-6)
            if let t = pl.transcriptionPy4dstemModel {
                // the kit equals its NumPy transcription with the same excitation model, whichever orientation that finds
                #expect(r.templateIndex == t.template, "\(key) \(k)")
                #expect(gap(m, matrix(t.matrix)) < 1e-6, "\(key) \(k)")
                #expect(abs(toTruth - t.errorToTruthDeg) < 1e-4)
                #expect(abs(toPy - t.errorToPy4dstemDeg) < 1e-4)
                // The measured behaviour, stated: where the route finds the planted orientation (every pattern but one of 18
                // in the four crystals; the exception is a monoclinic pattern whose mirrored match to another sampled beam
                // direction scores 0.94 against 0.90 for the planted one) it is within 5 degrees of the truth and of py4DSTEM.
                if t.errorToTruthDeg < 5 { #expect(toTruth < 5 && toPy < 5, "\(key) \(k)") }
            } else {
                // rutile (4/mmm): no transcription of the kit's generic zone sampler exists; the bounds are the measured gaps
                // (largest 1.9 and 2.2 degrees) with margin
                #expect(toTruth < 3 && toPy < 4, "\(key) \(k)")
            }
        }
        rows.forEach { print("ACOM planted " + $0) }
    }
}

// MARK: The public contract

struct ACOMContractTests {
    private func smallLibrary(_ settings: ACOMTemplateSettings = ACOMTemplateSettings(zoneAxisCount: 40)) throws
        -> (ACOMTemplateLibrary, [ACOMReflection]) {
        let c = try #require(try fixture().crystals["al"])
        let r = reflections(c)
        return (try ACOMTemplateLibrary.generate(reflections: r, kMax: c.kmax, laueClass: .mThreeBarM, settings: settings), r)
    }

    @Test func settingsAreRefusedBeforeAnyWork() throws {
        let c = try #require(try fixture().crystals["al"])
        let r = reflections(c)
        func refused(_ s: ACOMTemplateSettings, kMax: Double = 1.6, _ reflections: [ACOMReflection]? = nil) -> Bool {
            do {
                _ = try ACOMTemplateLibrary.generate(reflections: reflections ?? r, kMax: kMax, laueClass: .mThreeBarM, settings: s)
                return false
            } catch { return true }
        }
        #expect(refused(ACOMTemplateSettings(nAzimuthal: 100)))
        #expect(refused(ACOMTemplateSettings(nAzimuthal: 2)))
        #expect(refused(ACOMTemplateSettings(nRadial: 0)))
        #expect(refused(ACOMTemplateSettings(zoneAxisCount: 4097)))
        #expect(refused(ACOMTemplateSettings(zoneAxisCount: -1)))
        #expect(refused(ACOMTemplateSettings(sgWidth: 0)))
        #expect(refused(ACOMTemplateSettings(excitation: .py4DSTEM(wavelengthAngstrom: 0))))
        #expect(refused(ACOMTemplateSettings(), kMax: .nan))
        #expect(refused(ACOMTemplateSettings(), []))
        #expect(refused(ACOMTemplateSettings(), [ACOMReflection(g: [.nan, 0, 0], intensity: 1)]))
        #expect(refused(ACOMTemplateSettings(), [ACOMReflection(g: [0, 0, 0], intensity: 1)]))
        #expect(throws: ACOMError.noUsableReflections) {
            try ACOMTemplateLibrary.generate(reflections: [], kMax: 1, laueClass: .mThreeBarM)
        }
    }

    @Test func tooFewPeaksGiveTheEmptyResult() throws {
        let (lib, _) = try smallLibrary()
        let matcher = try #require(ACOMMatcher(library: lib))
        #expect(matcher.match(peaks: [], originX: 0, originY: 0) == .empty)
        #expect(matcher.match(peaks: [ACOMPeak(x: 1, y: 0, intensity: 1)], originX: 0, originY: 0) == .empty)
        #expect(matcher.templateScores(peaks: [], originX: 0, originY: 0).isEmpty)
        #expect(matcher.experimentalPolarImage(peaks: [], originX: 0, originY: 0) == nil)
    }

    /// With a flat Ewald sphere every template is exactly pi-periodic in azimuth (the apps' own note): the image is
    /// unchanged by a half turn, so the in-plane angle is only defined modulo pi.
    @Test func flatTemplatesArePiPeriodicAndCurvedOnesAreNot() throws {
        let (flat, _) = try smallLibrary()
        let (curved, _) = try smallLibrary(ACOMTemplateSettings(
            zoneAxisCount: 40, excitation: .appCurvature(wavelengthAngstrom: 0.0251)))
        let na = flat.geometry.nAzimuthal, nr = flat.geometry.nRadial
        func asymmetry(_ lib: ACOMTemplateLibrary) -> Double {
            var worst = 0.0
            for t in lib.templates {
                for r in 0..<nr {
                    for a in 0..<na { worst = max(worst, abs(Double(t[r * na + a]) - Double(t[r * na + (a + na / 2) % na]))) }
                }
            }
            return worst
        }
        #expect(asymmetry(flat) < 1e-5)
        #expect(asymmetry(curved) > 1e-3)
    }

    /// A template's own spots, rotated in plane, are matched to that template (or an equivalent neighbour) with a score near 1,
    /// including a calibrated pixel size and an off-centre origin.
    @Test func aPatternBuiltFromATemplatesSpotsMatchesIt() throws {
        let (lib, _) = try smallLibrary(ACOMTemplateSettings(zoneAxisCount: 60,
                                                             excitation: .appCurvature(wavelengthAngstrom: 0.0251)))
        let matcher = try #require(ACOMMatcher(library: lib))
        let t = 17, phi = 0.9
        let ipp = 0.004
        let pk = lib.templateSpots[t].map { s -> ACOMPeak in
            let a = Double(s.azimuth) + phi, r = Double(s.radius) / ipp
            return ACOMPeak(x: Float(100 + r * cos(a)), y: Float(-40 + r * sin(a)), intensity: pow(s.weight, 4))
        }
        let r = matcher.match(peaks: pk, originX: 100, originY: -40, invAngstromPerPixel: ipp)
        #expect(r.templateIndex == t)
        #expect(r.score > 0.93)
        #expect(!r.mirrored)
        let want = lib.orientationMatrix(templateIndex: t, inPlaneAngle: Float(phi), mirrored: false)
        let got = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
        #expect(misorientationDegrees(got, want, laue: .mThreeBarM) < 3.5)
    }

    @Test func libraryAndMatcherAreSendable() throws {
        let (lib, _) = try smallLibrary()
        let matcher = try #require(ACOMMatcher(library: lib))
        let s: any Sendable = matcher
        _ = s
        #expect(matcher.library.count == 40)
    }
}

// MARK: Scan matching

struct ACOMScanTests {
    private func setup() throws -> (ACOMMatcher, ACOMTemplateLibrary) {
        let c = try #require(try fixture().crystals["al"])
        let lib = try ACOMTemplateLibrary.generate(
            reflections: reflections(c), kMax: c.kmax, laueClass: .mThreeBarM,
            settings: ACOMTemplateSettings(zoneAxisCount: 40, excitation: .appCurvature(wavelengthAngstrom: 0.0251)))
        return (try #require(ACOMMatcher(library: lib)), lib)
    }

    /// A scan whose pattern at (x, y) is template (3 * x + y) % 40 rotated by a position-dependent angle; every fifth is empty.
    private func scan(_ lib: ACOMTemplateLibrary, width: Int, height: Int) -> [[ACOMPeak]] {
        (0..<(width * height)).map { i in
            if i % 5 == 4 { return [] }
            let phi = 0.1 * Double(i % 7), ipp = 0.004
            return lib.templateSpots[(3 * i + 1) % lib.count].map { s in
                let a = Double(s.azimuth) + phi, r = Double(s.radius) / ipp
                return ACOMPeak(x: Float(64 + r * cos(a)), y: Float(64 + r * sin(a)), intensity: pow(s.weight, 4))
            }
        }
    }

    @Test func selectionsListTheirPositions() {
        #expect(ACOMScanSelection.full.sourceIndices(width: 3, height: 2) == [0, 1, 2, 3, 4, 5])
        #expect(ACOMScanSelection.full.sourceIndices(width: 0, height: 2).isEmpty)
        // 5 x 3 into at most 2 blocks a side: step 3, so columns 0 and 3 of row 0.
        #expect(ACOMScanSelection.preview(maxDimension: 2).sourceIndices(width: 5, height: 3) == [0, 3])
        #expect(ACOMScanSelection.preview(maxDimension: 2).previewStep(width: 5, height: 3) == 3)
        #expect(ACOMScanSelection.preview(maxDimension: 0).previewStep(width: 5, height: 3) == 5)
        #expect(ACOMScanSelection.full.previewStep(width: 5, height: 3) == 1)
        #expect(ACOMScanSelection.square(centerX: 0, centerY: 0, radius: 1).sourceIndices(width: 4, height: 4) == [0, 1, 4, 5])
        #expect(ACOMScanSelection.square(centerX: 9, centerY: 9, radius: 1).sourceIndices(width: 4, height: 4).isEmpty)
        #expect(ACOMScanSelection.square(centerX: 1, centerY: 1, radius: -3).sourceIndices(width: 4, height: 4) == [5])
        #expect(ACOMScanSelection.square(centerX: 2, centerY: 1, radius: 1).positionCount(width: 4, height: 4) == 9)
    }

    @Test func fullScanEqualsPatternByPatternForAnyWorkerCount() throws {
        let (matcher, lib) = try setup()
        let (w, h) = (6, 5)
        let peaks = scan(lib, width: w, height: h)
        let want = peaks.map { matcher.match(peaks: $0, originX: 64, originY: 64, invAngstromPerPixel: 0.004) }
        #expect(want.contains { $0.templateIndex >= 0 } && want.contains { $0.templateIndex < 0 })
        for workers in [1, 3, nil] as [Int?] {
            let got = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64,
                                                     invAngstromPerPixel: 0.004, maxWorkers: workers))
            #expect(got.results == want)
            #expect(got.width == w && got.height == h && got.templateCount == lib.count)
            #expect(got[2, 1] == want[1 * w + 2])
        }
    }

    @Test func progressEndsAtOneAndNeverExceedsIt() throws {
        let (matcher, lib) = try setup()
        let peaks = scan(lib, width: 8, height: 8)
        let seen = LockedBox<[Double]>([])
        let done = matcher.matchAll(peaks: peaks, width: 8, height: 8, originX: 64, originY: 64, invAngstromPerPixel: 0.004,
                                    maxWorkers: 4, progress: { f in seen.update { $0.append(f) } })
        #expect(done != nil)
        let values = seen.value
        #expect(!values.isEmpty)
        #expect(values.allSatisfy { $0 > 0 && $0 <= 1 })
        #expect(values.max() == 1)
        #expect(values.last == 1)
        #expect(zip(values, values.dropFirst()).allSatisfy { $0 <= $1 })
    }

    @Test func hugeSquareAndPreviewArgumentsClipInsteadOfTrapping() {
        #expect(ACOMScanSelection.square(centerX: 2, centerY: 2, radius: Int.max).positionCount(width: 4, height: 3) == 12)
        #expect(ACOMScanSelection.square(centerX: Int.min, centerY: 0, radius: 1).positionCount(width: 4, height: 3) == 0)
        #expect(ACOMScanSelection.preview(maxDimension: Int.max).sourceIndices(width: 4, height: 3).count == 12)
    }

    @Test func squareLeavesTheOutsideEmpty() throws {
        let (matcher, lib) = try setup()
        let (w, h) = (6, 5)
        let peaks = scan(lib, width: w, height: h)
        let full = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64, invAngstromPerPixel: 0.004))
        let sq = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64, invAngstromPerPixel: 0.004,
                                               selection: .square(centerX: 2, centerY: 2, radius: 1)))
        for y in 0..<h {
            for x in 0..<w {
                let inside = abs(x - 2) <= 1 && abs(y - 2) <= 1
                #expect(sq[x, y] == (inside ? full[x, y] : .empty))
            }
        }
    }

    @Test func previewCopiesEachBlocksCornerOverTheBlock() throws {
        let (matcher, lib) = try setup()
        let (w, h) = (6, 4)
        let peaks = scan(lib, width: w, height: h)
        let full = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64, invAngstromPerPixel: 0.004))
        let sel = ACOMScanSelection.preview(maxDimension: 3)    // step 2
        let pv = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64, invAngstromPerPixel: 0.004,
                                               selection: sel))
        for y in 0..<h { for x in 0..<w { #expect(pv[x, y] == full[(x / 2) * 2, (y / 2) * 2]) } }
    }

    @Test func cancellationReturnsNilAndStopsStartingPositions() throws {
        let (matcher, lib) = try setup()
        let peaks = scan(lib, width: 10, height: 10)
        let upFront = matcher.matchAll(peaks: peaks, width: 10, height: 10, originX: 64, originY: 64, isCancelled: { true })
        #expect(upFront == nil)
        // Cancel after about 20 polls: the run is nil, and far fewer than 100 positions were started (one worker, so the polls
        // are the positions).
        let polls = LockedBox(0)
        let r = matcher.matchAll(peaks: peaks, width: 10, height: 10, originX: 64, originY: 64, invAngstromPerPixel: 0.004, maxWorkers: 1,
                                 isCancelled: { polls.update { $0 += 1 }; return polls.value > 20 })
        #expect(r == nil)
        #expect(polls.value <= 23)
    }

    @Test func badInputIsRefused() throws {
        let (matcher, lib) = try setup()
        let peaks = scan(lib, width: 3, height: 3)
        #expect(matcher.matchAll(peaks: peaks, width: 3, height: 4, originX: 64, originY: 64) == nil)
        #expect(matcher.matchAll(peaks: [], width: 0, height: 0, originX: 64, originY: 64) == nil)
        #expect(matcher.matchAll(peaks: peaks, width: 3, height: 3, originX: 64, originY: 64,
                                 selection: .square(centerX: 9, centerY: 9, radius: 1)) == nil)
    }
}

// MARK: Metal backend

/// `ACOMBackend.metal` against the CPU reference. Written tolerances: scores and runner-up scores (normalised correlations in
/// [0, 1], Float32 sums over 32 rings and a 128-term inverse transform) within 1e-6 absolute (the CPU's own FFT and the GPU's
/// exact DFT round differently; measured worst gaps are printed); template, in-plane angle, mirror flag, Euler angles and
/// disorientation identical (they follow from the template and the bin). Cases where two bins or templates tie within the score
/// tolerance could differ; on the pinned fixtures the pinned margins (> 1e-4) exclude them. Skipped where no Metal device exists.
struct ACOMMetalTests {
    static let scoreTolerance: Float = 1e-6

    private func setup() throws -> (ACOMMatcher, ACOMTemplateLibrary) {
        let c = try #require(try fixture().crystals["al"])
        let lib = try ACOMTemplateLibrary.generate(
            reflections: reflections(c), kMax: c.kmax, laueClass: .mThreeBarM,
            settings: ACOMTemplateSettings(zoneAxisCount: 40, excitation: .appCurvature(wavelengthAngstrom: 0.0251)))
        return (try #require(ACOMMatcher(library: lib)), lib)
    }

    private func scan(_ lib: ACOMTemplateLibrary, width: Int, height: Int) -> [[ACOMPeak]] {
        (0..<(width * height)).map { i in
            if i % 5 == 4 { return [] }
            let phi = 0.1 * Double(i % 7), ipp = 0.004
            // every third pattern is mirrored (azimuth reversed) so the second pass wins sometimes
            return lib.templateSpots[(3 * i + 1) % lib.count].map { s in
                let a = (i % 3 == 0 ? -Double(s.azimuth) : Double(s.azimuth)) + phi, r = Double(s.radius) / ipp
                return ACOMPeak(x: Float(64 + r * cos(a)), y: Float(64 + r * sin(a)), intensity: pow(s.weight, 4))
            }
        }
    }

    /// Same result apart from the scores, which agree within the tolerance; returns the gap. With `tiesAllowed`, a position where
    /// the GPU took another route (template, bin or pass) to a top score and runner-up that equal the CPU's within the
    /// tolerance is a tie within Float32 rounding (typically a pattern that is its own mirror image: both passes score alike): it
    /// is counted in `ties` instead of failing.
    @discardableResult
    private func same(_ gpu: OrientationResult, _ cpu: OrientationResult, _ what: String,
                      tiesAllowed ties: LockedBox<Int>? = nil) -> Float {
        let gap = max(abs(gpu.score - cpu.score), abs(gpu.secondScore - cpu.secondScore))
        #expect(gap <= Self.scoreTolerance, "\(what) score gap \(gap)")
        let sameRoute = gpu.templateIndex == cpu.templateIndex && gpu.inPlaneAngle == cpu.inPlaneAngle && gpu.mirrored == cpu.mirrored
        if !sameRoute, let ties { ties.update { $0 += 1 }; return gap }
        #expect(gpu.templateIndex == cpu.templateIndex, "\(what) template")
        #expect(gpu.inPlaneAngle == cpu.inPlaneAngle, "\(what) in-plane angle")
        #expect(gpu.mirrored == cpu.mirrored, "\(what) mirrored")
        #expect(gpu.euler == cpu.euler && gpu.symmetryDisorientationRad == cpu.symmetryDisorientationRad && gpu.phaseID == cpu.phaseID,
                "\(what) euler")
        return gap
    }

    @Test func metalEqualsTheCPUOnTheAppRoutePins() throws {
        var worst: Float = 0, cases = 0
        for key in appKeys {
            let c = try #require(try fixture().crystals[key])
            let lib = try appLibrary(key, c)
            let matcher = try #require(ACOMMatcher(library: lib))
            guard matcher.isMetalAvailable else { return }
            let pinned = try #require(c.appRoute).cases
            for (k, cs) in pinned.enumerated() {
                let one = [peaks(cs)], ox = Float(cs.origin[0]), oy = Float(cs.origin[1])
                let cpu = try #require(matcher.matchAll(peaks: one, width: 1, height: 1, originX: ox, originY: oy,
                                                         invAngstromPerPixel: cs.invAngstromPerPixel, backend: .cpu)).results[0]
                let gpu = try #require(matcher.matchAll(peaks: one, width: 1, height: 1, originX: ox, originY: oy,
                                                         invAngstromPerPixel: cs.invAngstromPerPixel, backend: .metal)).results[0]
                #expect(cs.marginTemplate > 1e-4 && cs.marginInPlane > 1e-4)
                #expect(gpu.templateIndex == cs.template && gpu.mirrored == cs.mirrored, "\(key) \(k) versus the pin")
                #expect(abs(Double(gpu.inPlaneAngle) - cs.angle) < 2e-6, "\(key) \(k) angle versus the pin")
                worst = max(worst, same(gpu, cpu, "\(key) case \(k)"))   // no ties allowed on the pins
                cases += 1
            }
        }
        print("ACOM Metal versus CPU on the app-route pins: \(cases) patterns, worst score gap \(worst)")
    }

    @Test func metalScanEqualsCPUAcrossBatchesSelectionsAndEmptyPatterns() throws {
        let (matcher, lib) = try setup()
        guard matcher.isMetalAvailable else { return }
        let (w, h) = (13, 12)   // 156 positions: three GPU batches of at most 64
        let peaks = scan(lib, width: w, height: h)
        var worst: Float = 0
        let ties = LockedBox(0)
        for selection in [ACOMScanSelection.full, .preview(maxDimension: 5), .square(centerX: 6, centerY: 5, radius: 2)] {
            let cpu = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64,
                                                     invAngstromPerPixel: 0.004, selection: selection, backend: .cpu))
            let gpu = try #require(matcher.matchAll(peaks: peaks, width: w, height: h, originX: 64, originY: 64,
                                                     invAngstromPerPixel: 0.004, selection: selection, maxWorkers: 3, backend: .metal))
            #expect(gpu.width == w && gpu.height == h && gpu.templateCount == lib.count)
            for i in cpu.results.indices { worst = max(worst, same(gpu.results[i], cpu.results[i], "\(selection) position \(i)", tiesAllowed: ties)) }
            if selection == .full {
                #expect(cpu.results.contains { $0.templateIndex < 0 } && cpu.results.contains { $0.mirrored } && cpu.results.contains { !$0.mirrored && $0.templateIndex >= 0 })
            }
        }
        print("ACOM Metal versus CPU scan (156 positions, three selections): worst score gap \(worst), ties within rounding \(ties.value)")
        #expect(ties.value <= 3, "ties are the exception, not the rule: \(ties.value)")
    }

    @Test func metalProgressCancellationAndRefusals() throws {
        let (matcher, lib) = try setup()
        guard matcher.isMetalAvailable else { return }
        let peaks = scan(lib, width: 9, height: 9)
        let seen = LockedBox<[Double]>([])
        let r = matcher.matchAll(peaks: peaks, width: 9, height: 9, originX: 64, originY: 64, invAngstromPerPixel: 0.004,
                                 backend: .metal, progress: { f in seen.update { $0.append(f) } })
        #expect(r != nil && seen.value.last == 1 && zip(seen.value, seen.value.dropFirst()).allSatisfy { $0 <= $1 } && seen.value.allSatisfy { $0 <= 1 })
        #expect(matcher.matchAll(peaks: peaks, width: 9, height: 9, originX: 64, originY: 64, backend: .metal,
                                 isCancelled: { true }) == nil)
        let calls = LockedBox(0)
        #expect(matcher.matchAll(peaks: peaks, width: 9, height: 9, originX: 64, originY: 64, backend: .metal,
                                 isCancelled: { calls.update { $0 += 1 }; return calls.value > 3 }) == nil)
        #expect(matcher.matchAll(peaks: peaks, width: 9, height: 8, originX: 64, originY: 64, backend: .metal) == nil)
        #expect(matcher.matchAll(peaks: peaks, width: 9, height: 9, originX: 64, originY: 64,
                                 selection: .square(centerX: 99, centerY: 99, radius: 1), backend: .metal) == nil)
    }

    @Test func cpuStaysTheDefaultAndTheBackendCasesAreStable() throws {
        #expect(ACOMBackend.allCases == [.cpu, .metal])
        let (matcher, lib) = try setup()
        let peaks = scan(lib, width: 4, height: 3)
        let implicit = try #require(matcher.matchAll(peaks: peaks, width: 4, height: 3, originX: 64, originY: 64, invAngstromPerPixel: 0.004))
        let explicit = try #require(matcher.matchAll(peaks: peaks, width: 4, height: 3, originX: 64, originY: 64, invAngstromPerPixel: 0.004, backend: .cpu))
        #expect(implicit.results == explicit.results)
    }
}

private final class LockedBox<T: Sendable>: @unchecked Sendable {
    private let lock = NSLock()
    private var stored: T
    init(_ value: T) { stored = value }
    var value: T { lock.lock(); defer { lock.unlock() }; return stored }
    func update(_ body: (inout T) -> Void) { lock.lock(); defer { lock.unlock() }; body(&stored) }
}

// MARK: Off-grid self-recovery (diagnosis 2026-10-21)

/// The apps' open item "off-grid self-recovery fails 37-43 of 200 at 1.4 degrees": 200 seeded patterns from beam directions
/// perturbed off the library grid (up to 2.9 degrees from a template) and a random in-plane angle, built from the kit's own
/// excitation model, matched, compared with the planted orientation over the Laue operations. DIAGNOSIS: the failures are
/// the in-plane bin, not the beam-direction grid. 128 bins are 2.81 degrees apart, so the angle alone is off by up to
/// 1.41; going to 512 bins takes the failures from 35 to 9, quadrupling the zone axes (600 to 2400) does not move the
/// in-plane term, and fitting a parabola through the correlation peak (`refinesInPlaneAngle`) takes 35 to 14 at 128 bins.
/// The rest (7 of 200 above 10 degrees, 36 or 60 degrees, mirrored or not) are coincidence twins of the truth, not a defect
/// of the matcher (diagnosis for issue 8, `theLargeFailuresAreCoincidenceTwinsOfTheTruth`): the winner lies within 1.6
/// degrees of the truth turned by exactly 60 degrees about a <111> (sigma 3, 5 cases, beams near <112>) or 36.87 degrees about
/// a <100> (sigma 5, 2 cases, beams near <013>). That axis lies in the zone plane, so both orientations excite the same
/// zero-order spots; the flat-sphere part of the excitation error is even in g, and only the curvature term (lambda |g|^2 / 2,
/// about 0.003 per inverse angstrom at 200 kV, then taken to the power 0.25) tells them apart. The two exact patterns overlap
/// 0.95 to 1.0000 in the matcher's own polar image, more than the grid and bin quantisation loses (the best orientation within
/// 3 degrees of the truth scores 0.0007 to 0.021 lower). More zone axes (2400) or more bins (512) still leave 3 to 9 cases
/// above 10 degrees (the three left at 512 bins and 2400 axes are twins too); ten times the wavelength leaves none. py4DSTEM f050d207 (its own patterns, model and match_single_pattern)
/// picks the twin in 7 of the 14 truth and twin patterns.
struct ACOMOffGridRecoveryTests {
    private func failures(bins: Int, zoneAxes: Int, refine: Bool) throws -> (over1p4: Int, over10: Int, median: Double) {
        let c = try #require(try fixture().crystals["al"])
        let refl = reflections(c)
        let lib = try ACOMTemplateLibrary.generate(
            reflections: refl, kMax: c.kmax, laueClass: .mThreeBarM,
            settings: ACOMTemplateSettings(nAzimuthal: bins, zoneAxisCount: zoneAxes,
                                           excitation: .appCurvature(wavelengthAngstrom: c.wavelength)))
        let matcher = try #require(ACOMMatcher(library: lib, refinesInPlaneAngle: refine))
        var seed: UInt64 = 12345
        func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
        var errors: [Double] = []
        for _ in 0..<200 {
            let n = simd_normalize(lib.zoneAxes[Int(rnd() * Double(lib.count))]
                                   + SIMD3(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5) * 0.05)
            let phi = rnd() * 2 * .pi
            let ipp = 0.004
            let peaks = ACOMTemplateLibrary.project(reflections: refl, zoneAxis: n, settings: lib.settings).map { s in
                ACOMPeak(x: Float(s.r * cos(s.azimuth + phi) / ipp), y: Float(s.r * sin(s.azimuth + phi) / ipp),
                         intensity: Float(pow(s.weight, 1 / lib.settings.intensityPower)))
            }
            let r = matcher.match(peaks: peaks, originX: 0, originY: 0, invAngstromPerPixel: ipp)
            let truth = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: n), inPlaneAngle: Float(phi))
            let got = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
            errors.append(misorientationDegrees(got, truth, laue: .mThreeBarM))
        }
        errors.sort()
        return (errors.filter { $0 > 1.4 }.count, errors.filter { $0 > 10 }.count, errors[100])
    }

    @Test func theInPlaneBinIsTheFloorAndTheParabolaMovesIt() throws {
        let coarse = try failures(bins: 128, zoneAxes: 600, refine: false)
        let refined = try failures(bins: 128, zoneAxes: 600, refine: true)
        let fine = try failures(bins: 512, zoneAxes: 600, refine: false)
        let denseAxes = try failures(bins: 128, zoneAxes: 2400, refine: false)
        print("ACOM off-grid >1.4 deg of 200: 128 bins \(coarse), refined \(refined), 512 bins \(fine), 2400 axes \(denseAxes)")
        #expect(coarse.over1p4 >= 25)              // the app's 37-43 of 200, reproduced
        #expect(refined.over1p4 <= 20)
        #expect(refined.median < coarse.median)
        #expect(fine.over1p4 <= 12)
        #expect(denseAxes.over1p4 >= 20)           // more beam directions do not remove the in-plane floor
        #expect(coarse.over10 <= 10)
    }

    /// Every failure above 10 degrees of the 200 off-grid plants is a coincidence twin of the truth (60 degrees about <111>
    /// or 36.87 degrees about <100>) whose exact pattern the truth pattern overlaps by at least 0.95 in the matcher's polar
    /// image; with ten times the wavelength (the curvature term that alone separates the twins, ten times larger) none is left.
    @Test func theLargeFailuresAreCoincidenceTwinsOfTheTruth() throws {
        let c = try #require(try fixture().crystals["al"])
        let refl = reflections(c)
        func twinRotation(_ axis: SIMD3<Double>, _ degrees: Double) -> simd_double3x3 {
            simd_double3x3(simd_quatd(angle: degrees * .pi / 180, axis: simd_normalize(axis)))
        }
        let twins: [simd_double3x3] = [SIMD3<Double>(1, 1, 1), [-1, 1, 1], [1, -1, 1], [1, 1, -1]].flatMap {
            [twinRotation($0, 60), twinRotation($0, -60)]
        } + [SIMD3<Double>(1, 0, 0), [0, 1, 0], [0, 0, 1]].flatMap { [twinRotation($0, 36.87), twinRotation($0, -36.87)] }
        func run(wavelength: Double) throws -> (failures: Int, twinned: Int, worstGap: Double, leastOverlap: Double) {
            let settings = ACOMTemplateSettings(nAzimuthal: 128, zoneAxisCount: 600,
                                                excitation: .appCurvature(wavelengthAngstrom: wavelength))
            let lib = try ACOMTemplateLibrary.generate(reflections: refl, kMax: c.kmax, laueClass: .mThreeBarM, settings: settings)
            let matcher = try #require(ACOMMatcher(library: lib))
            let ipp = 0.004
            // the pattern of an orientation (columns lab x, lab y, beam in crystal coordinates), as `failures` builds it
            func pattern(_ o: simd_double3x3) -> [ACOMPeak] {
                let n = o.columns.2
                return refl.compactMap { r in
                    let sg = simd_dot(r.g, n) + wavelength / 2 * simd_length_squared(r.g)
                    guard abs(sg) <= settings.sgMax else { return nil }
                    let w = r.intensity * exp(-(sg / settings.sgWidth) * (sg / settings.sgWidth))
                    return ACOMPeak(x: Float(simd_dot(r.g, o.columns.0) / ipp), y: Float(simd_dot(r.g, o.columns.1) / ipp),
                                    intensity: Float(w))
                }
            }
            func polar(_ o: simd_double3x3) throws -> [Float] {
                try #require(matcher.experimentalPolarImage(peaks: pattern(o), originX: 0, originY: 0, invAngstromPerPixel: ipp))
            }
            var seed: UInt64 = 12345
            func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
            var failures = 0, twinned = 0, worstGap = 0.0, leastOverlap = 1.0
            for _ in 0..<200 {
                let n = simd_normalize(lib.zoneAxes[Int(rnd() * Double(lib.count))]
                                       + SIMD3(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5) * 0.05)
                let phi = rnd() * 2 * .pi
                let truth = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: n), inPlaneAngle: Float(phi))
                let r = matcher.match(peaks: pattern(truth), originX: 0, originY: 0, invAngstromPerPixel: ipp)
                let got = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
                guard misorientationDegrees(got, truth, laue: .mThreeBarM) > 10 else { continue }
                failures += 1
                let twin = twins.map { $0 * truth }.min {
                    misorientationDegrees(got, $0, laue: .mThreeBarM) < misorientationDegrees(got, $1, laue: .mThreeBarM)
                }!
                let gap = misorientationDegrees(got, twin, laue: .mThreeBarM)
                let overlap = zip(try polar(truth), try polar(twin)).reduce(0.0) { $0 + Double($1.0) * Double($1.1) }
                worstGap = max(worstGap, gap)
                leastOverlap = min(leastOverlap, overlap)
                if gap < 2.5 && overlap > 0.95 { twinned += 1 }
            }
            return (failures, twinned, worstGap, leastOverlap)
        }
        let at200kV = try run(wavelength: c.wavelength)
        print("ACOM off-grid failures above 10 deg: \(at200kV.failures), coincidence twins \(at200kV.twinned), winner to exact twin at most \(at200kV.worstGap) deg, truth-twin overlap at least \(at200kV.leastOverlap)")
        #expect(at200kV.failures >= 1 && at200kV.failures <= 10)
        #expect(at200kV.twinned == at200kV.failures)
        let curved = try run(wavelength: c.wavelength * 10)
        #expect(curved.failures == 0)
    }

    @Test func refinementMovesOnlyTheInPlaneAngleByAtMostHalfABin() throws {
        let c = try #require(try fixture().crystals["al"])
        let lib = try ACOMTemplateLibrary.generate(
            reflections: reflections(c), kMax: c.kmax, laueClass: .mThreeBarM,
            settings: ACOMTemplateSettings(zoneAxisCount: 60, excitation: .appCurvature(wavelengthAngstrom: c.wavelength)))
        let plain = try #require(ACOMMatcher(library: lib)), fine = try #require(ACOMMatcher(library: lib, refinesInPlaneAngle: true))
        #expect(!plain.refinesInPlaneAngle && fine.refinesInPlaneAngle)
        let ipp = 0.004
        let peaks = lib.templateSpots[17].map { s -> ACOMPeak in
            let a = Double(s.azimuth) + 0.9, r = Double(s.radius) / ipp
            return ACOMPeak(x: Float(r * cos(a)), y: Float(r * sin(a)), intensity: pow(s.weight, 4))
        }
        let a = plain.match(peaks: peaks, originX: 0, originY: 0, invAngstromPerPixel: ipp)
        let b = fine.match(peaks: peaks, originX: 0, originY: 0, invAngstromPerPixel: ipp)
        #expect(a.templateIndex == b.templateIndex && a.mirrored == b.mirrored && a.score == b.score)
        let halfBin = Float.pi / Float(lib.geometry.nAzimuthal)
        var d = abs(a.inPlaneAngle - b.inPlaneAngle); d = min(d, 2 * .pi - d)
        #expect(d <= halfBin * 1.0001)
        #expect(abs(Double(b.inPlaneAngle) - 0.9) <= abs(Double(a.inPlaneAngle) - 0.9) + 1e-6)
    }
}

/// Issue 8, the sparse corner: is there a corner of the m-3m zone ([001]-[101]-[111]) where the farthest-point sampling leaves
/// a hole larger than elsewhere, so that a pattern between templates lands far away? DIAGNOSIS: no. The cause that was
/// suspected (a sparse corner) would be refuted by holes near the three vertices and along the edges no larger than inside,
/// and by patterns at the largest holes recovering to within the hole plus the in-plane bin. Both hold: at 600 axes the
/// largest hole anywhere is 1.06 degrees (an interior point); within 5 degrees of [001], [101] and [111] it is 0.92, 0.97 and
/// 1.03; at 2 400 axes 0.53. Planted at the 30 largest holes (8 in-plane angles each), 230 of 240 patterns land within the
/// hole plus one in-plane bin and the other 10 are coincidence twins of the truth whose axis lies in the zone plane (8 sigma 3,
/// 60 degrees about <111>; 2 sigma 7, 38.21 degrees about <111>, beam 1 degree from <213>), the ambiguity of
/// `theLargeFailuresAreCoincidenceTwinsOfTheTruth`. py4DSTEM f050d207 with its own orientation_plan at 2 degree steps has the
/// same largest hole (1.07 degrees, 1.03 to 1.05 near the vertices); its match_single_pattern on its own patterns at its 12
/// largest holes (4 in-plane angles each) is off by the hole, at most 2.2 degrees, and on the sigma 7 beam above it picks the
/// twin on all 6 truth patterns (20 of 48 truth and twin patterns at four beams near <213>). Both libraries behave alike: a
/// pattern between templates is recovered to the covering distance, unless a coincidence twin wins. "Zone axis misses up to
/// 12.8 degrees beyond the bank sampling" is not reproduced. Results stay `validation: "none"`.
struct ACOMSparseCornerTests {
    @Test func noCornerOfTheZoneIsSparserAndAPatternInTheLargestHoleIsRecoveredToTheHole() throws {
        let c = try #require(try fixture().crystals["al"])
        let refl = reflections(c)
        let bins = 128
        let lib = try ACOMTemplateLibrary.generate(
            reflections: refl, kMax: c.kmax, laueClass: .mThreeBarM,
            settings: ACOMTemplateSettings(nAzimuthal: bins, zoneAxisCount: 600,
                                           excitation: .appCurvature(wavelengthAngstrom: c.wavelength)))
        let matcher = try #require(ACOMMatcher(library: lib))
        func degrees(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double { acos(min(1, simd_dot(a, b))) * 180 / .pi }
        // dense points of the zone (z >= x >= y >= 0): a folded Fibonacci sphere plus the three edges
        let n = 100_000
        let golden = Double.pi * (3 - 5.0.squareRoot())
        var points: [SIMD3<Double>] = (0..<n).map { i in
            let z = 1 - 2 * (Double(i) + 0.5) / Double(n)
            let r = max(0, 1 - z * z).squareRoot(), a = golden * Double(i)
            let v = [abs(r * cos(a)), abs(r * sin(a)), abs(z)].sorted()
            return simd_normalize(SIMD3(v[1], v[0], v[2]))
        }
        for k in 0...400 {
            let t = Double(k) / 400
            points += [simd_normalize(SIMD3(t, 0, 1)), simd_normalize(SIMD3(t, t, 1)), simd_normalize(SIMD3(1, t, 1))]
        }
        // both a point and the templates lie in the zone, so the nearest symmetry image of a template is the template itself
        let hole = points.map { p in acos(min(1, lib.zoneAxes.map { simd_dot($0, p) }.max()!)) * 180 / .pi }
        let largest = hole.max()!
        let vertices: [SIMD3<Double>] = [[0, 0, 1], simd_normalize([1, 0, 1]), simd_normalize([1, 1, 1])]
        let atVertices = vertices.map { v in points.indices.filter { degrees(points[$0], v) < 5 }.map { hole[$0] }.max()! }
        #expect(largest < 1.2)
        #expect(atVertices.allSatisfy { $0 <= largest })
        // the 30 largest holes, at least 1.5 degrees apart
        var centres: [Int] = []
        for i in hole.indices.sorted(by: { hole[$0] > hole[$1] }) where centres.count < 30 {
            if centres.allSatisfy({ degrees(points[$0], points[i]) > 1.5 }) { centres.append(i) }
        }
        func rotation(_ axis: SIMD3<Double>, _ deg: Double) -> simd_double3x3 {
            simd_double3x3(simd_quatd(angle: deg * .pi / 180, axis: simd_normalize(axis)))
        }
        // coincidence rotations (sigma 3, 5, 7, 9, 11) about every variant of their axis: when the axis lies in the zone plane
        // the rotation keeps the zero-order spots, and only the curvature term tells it from the truth
        var twins: [(name: String, axis: SIMD3<Double>, rotation: simd_double3x3)] = []
        let families: [(String, [SIMD3<Double>], Double)] = [
            ("sigma 3", [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], 60),
            ("sigma 5", [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 36.87),
            ("sigma 7", [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], 38.21),
            ("sigma 9", [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]], 38.94),
            ("sigma 11", [[1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]], 50.48),
        ]
        for (name, family, angle) in families {
            for u in family { twins += [(name, u, rotation(u, angle)), (name, u, rotation(u, -angle))] }
        }
        var seed: UInt64 = 777
        func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
        let ipp = 0.004, bin = 360 / Double(bins)
        var within = 0, twinned: [String] = [], other: [String] = []
        for i in centres {
            let p = points[i]
            let spots = ACOMTemplateLibrary.project(reflections: refl, zoneAxis: p, settings: lib.settings)
            for _ in 0..<8 {
                let phi = rnd() * 2 * .pi
                let peaks = spots.map { s in
                    ACOMPeak(x: Float(s.r * cos(s.azimuth + phi) / ipp), y: Float(s.r * sin(s.azimuth + phi) / ipp),
                             intensity: Float(pow(s.weight, 1 / lib.settings.intensityPower)))
                }
                let r = matcher.match(peaks: peaks, originX: 0, originY: 0, invAngstromPerPixel: ipp)
                let truth = ZoneAxisFrame.matrix(detectorBasis: ZoneAxisFrame.detectorBasis(zoneAxis: p), inPlaneAngle: Float(phi))
                let got = lib.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
                let error = misorientationDegrees(got, truth, laue: .mThreeBarM)
                if error <= hole[i] + bin {
                    within += 1
                } else if let t = twins.first(where: { t in
                    abs(degrees(simd_normalize(t.axis), p) - 90) < 3
                        && misorientationDegrees(got, t.rotation * truth, laue: .mThreeBarM) < hole[i] + bin
                }) {
                    twinned.append(String(format: "%.1f deg %@", error, t.name))
                } else {
                    other.append(String(format: "%.2f deg at hole %.2f", error, hole[i]))
                }
            }
        }
        print("ACOM sparse corner, 600 axes: largest hole \(largest) deg, near [001] [101] [111] \(atVertices); at the 30 largest holes x 8 angles: within hole + one bin \(within), coincidence twins with the axis in the zone plane \(twinned.count) \(twinned), other \(other)")
        #expect(other.isEmpty)
        #expect(within >= 200)
    }
}

struct OrientationMapMirrorPlaneTests {
    /// The flag plane and the pattern's own in-plane angle: a mirrored result at 0.30 turns (0.80 raw) reads 0.30, a forward
    /// one at 0.10 reads 0.10, a mirrored one raw 0.20 wraps to 0.70. RED UNDER: the planes absent (compile error).
    @Test func mirroredPlaneAndPatternAngle() {
        func result(_ angleTurns: Float, mirrored: Bool) -> OrientationResult {
            OrientationResult(templateIndex: 0, euler: .init(phi1: 0, Phi: 0, phi2: 0), inPlaneAngle: angleTurns * 2 * .pi,
                              score: 1, secondScore: 0, phaseID: 0, mirrored: mirrored)
        }
        let rs = [result(0.80, mirrored: true), result(0.10, mirrored: false), result(0.20, mirrored: true)]
        #expect(OrientationMapAnalysis.mirroredPlane(rs) == [1, 0, 1])
        let a = (DeprecatedPatternAngle.self as PatternAnglePlane.Type).plane(rs)
        #expect(abs(a[0] - 0.30) < 1e-6 && abs(a[1] - 0.10) < 1e-6 && abs(a[2] - 0.70) < 1e-6)
    }
}

/// `patternInPlaneAnglePlane` is deprecated but still public, so its test stays; the call goes through a protocol so the test
/// build does not warn. Delete with the function.
private protocol PatternAnglePlane { static func plane(_ results: [OrientationResult]) -> [Float] }
private enum DeprecatedPatternAngle: PatternAnglePlane {
    @available(*, deprecated)
    static func plane(_ results: [OrientationResult]) -> [Float] { OrientationMapAnalysis.patternInPlaneAnglePlane(results) }
}
