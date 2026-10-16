//
//  ACOMChecks.swift
//  Role: the ACOM matcher's public API (DiffractoOrientation), used from outside the package: a library from a crystal (gold) and
//        from hand-written reflections, a pattern built from a template's own spots matched back, the refusals. `runACOMChecks`
//        is called from main.swift.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

func runACOMChecks(_ report: inout Report) {
    report.section("ACOM matcher") { r in
        let gold = CrystalStructure(
            cell: UnitCell(a: 4.0782, b: 4.0782, c: 4.0782, alphaDeg: 90, betaDeg: 90, gammaDeg: 90),
            sites: [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]].map {
                AtomSite(z: 79, fractional: SIMD3<Double>($0[0], $0[1], $0[2]))
            }, operations: [.identity])
        let settings = ACOMTemplateSettings(zoneAxisCount: 30, excitation: .appCurvature(wavelengthAngstrom: 0.025079))
        // gold's file has the identity only, so its detected Laue class is -1 (all sampled beam directions are distinct)
        let lib = try ACOMTemplateLibrary.generate(crystal: gold, kMax: 1.5, settings: settings)
        r.equal("the library has the requested number of beam directions", lib.count, 30)
        r.equal("each template is nRadial x nAzimuthal", lib.templates[0].count, 32 * 128)
        r.equal("template FFTs are contiguous per template", lib.templateFFTRe.count, 30 * 32 * 128)
        r.equal("a template has unit norm", abs(lib.templates[3].reduce(0) { $0 + Double($1) * Double($1) } - 1) < 1e-4, true)
        let matcher = try unwrap(ACOMMatcher(library: lib))
        let t = 5
        let peaks = lib.templateSpots[t].map { s -> ACOMPeak in
            let a = Double(s.azimuth) + 0.6, rad = Double(s.radius) / 0.01
            return ACOMPeak(x: Float(200 + rad * cos(a)), y: Float(190 + rad * sin(a)), intensity: pow(s.weight, 4))
        }
        let m = matcher.match(peaks: peaks, originX: 200, originY: 190, invAngstromPerPixel: 0.01)
        r.equal("a template's own spots, rotated and calibrated, match that template", m.templateIndex, t)
        r.check("the score is near 1", m.score > 0.9, "score \(m.score)")
        r.check("the reliability is in [0, 1]", m.reliability >= 0 && m.reliability <= 1)
        r.equal("fewer than two peaks give the empty result", matcher.match(peaks: [peaks[0]], originX: 0, originY: 0) == .empty, true)
        let spots = ACOMTemplateLibrary.project(
            reflections: [ACOMReflection(g: [1, 0, 0], intensity: 1)], zoneAxis: [0, 0, 1], settings: ACOMTemplateSettings())
        r.equal("a reflection in the zone is a spot at its radius", spots.count == 1 && abs(spots[0].r - 1) < 1e-12, true)
        r.equal("an empty reflection list is refused", { () -> Bool in
            do { _ = try ACOMTemplateLibrary.generate(reflections: [], kMax: 1, laueClass: .mThreeBarM); return false }
            catch { return (error as? ACOMError) == .noUsableReflections }
        }(), true)
        let mat = lib.orientationMatrix(templateIndex: t, inPlaneAngle: 0, mirrored: false)
        r.near("the orientation matrix is a rotation", simd_determinant(mat), 1, 1e-12)
        let s: any Sendable = matcher
        _ = s
    }
}
