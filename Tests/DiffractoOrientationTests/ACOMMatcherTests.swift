//
//  ACOMMatcherTests.swift
//  Role: the ACOM template library and matcher (docs/parity.md P20). Pins from tools/acom-matcher/pins.py:
//        (1) the apps' route, a NumPy TRANSCRIPTION of the app source (OrientationPlan, OrientationMatcher at mac4DSTEM c809679c):
//            beam directions, template images and spots, polar images, ring FFTs, every template's score, the winner, the
//            runner-up, the in-plane angle and the orientation matrix, for m-3m, 6/mmm and 2/m;
//        (2) py4DSTEM f050d207, run: the spots of a zone axis (flat at lambda -> 0 and curved at 200 kV) against
//            `ACOMTemplateLibrary.project`, and planted orientations against py4DSTEM's own match_single_pattern.
//

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
