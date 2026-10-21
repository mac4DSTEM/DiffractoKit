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

        // The scan layer: 4 x 3 positions from template spots, a cancelled run, a square selection, a preview.
        let scanPeaks: [[ACOMPeak]] = (0..<12).map { i in
            i % 4 == 3 ? [] : lib.templateSpots[(i * 2 + 1) % lib.count].map { s in
                let a = Double(s.azimuth) + 0.1 * Double(i), rad = Double(s.radius) / 0.01
                return ACOMPeak(x: Float(200 + rad * cos(a)), y: Float(190 + rad * sin(a)), intensity: pow(s.weight, 4))
            }
        }
        let scan = matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, invAngstromPerPixel: 0.01,
                                    selection: .full, maxWorkers: 2, isCancelled: { false }, progress: { _ in })
        r.equal("a scan has one result per position", scan?.results.count, 12)
        r.equal("a scan position equals the single-pattern match",
                scan?[1, 1] == matcher.match(peaks: scanPeaks[5], originX: 200, originY: 190, invAngstromPerPixel: 0.01), true)
        r.equal("a position without peaks stays empty", scan?[3, 0] == .empty, true)
        r.equal("a cancelled scan is nil",
                matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, isCancelled: { true }) == nil, true)
        r.equal("the square selection lists its positions",
                ACOMScanSelection.square(centerX: 1, centerY: 1, radius: 1).sourceIndices(width: 4, height: 3), [0, 1, 2, 4, 5, 6, 8, 9, 10])
        let preview = matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, invAngstromPerPixel: 0.01,
                                       selection: .preview(maxDimension: 2))
        r.equal("a preview copies the block corner over the block", preview?[1, 1] == preview?[0, 0], true)
        r.equal("the preview step fits the longer side", ACOMScanSelection.preview(maxDimension: 2).previewStep(width: 4, height: 3), 2)
        let cpuScan = matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, invAngstromPerPixel: 0.01, backend: .cpu)
        r.equal("the CPU backend is the default", cpuScan?.results == scan?.results, true)
        if matcher.isMetalAvailable {
            let gpu = matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, invAngstromPerPixel: 0.01,
                                       maxWorkers: 2, backend: .metal)
            r.equal("the Metal scan has one result per position", gpu?.results.count, 12)
            r.equal("the Metal scan finds the same templates as the CPU",
                    zip(gpu?.results ?? [], scan?.results ?? []).filter { $0.templateIndex != $1.templateIndex }.count <= 1, true)
        } else {
            r.check("Metal is unavailable: the Metal scan is nil", matcher.matchAll(peaks: scanPeaks, width: 4, height: 3, originX: 200, originY: 190, backend: .metal) == nil)
        }
        // 0.7.0: a peaks view, one-peak polar image, result from scores, generate cancellation.
        let viewedScan = matcher.matchAll(width: 4, height: 3, originX: 200, originY: 190, invAngstromPerPixel: 0.01,
                                          selection: .full, maxWorkers: 2, peaksAt: { scanPeaks[$0].lazy.map { $0 } })
        r.equal("a peaks view gives the array form's results", viewedScan?.results == scan?.results, true)
        r.equal("one peak draws the polar image",
                matcher.experimentalPolarImage(peaks: [ACOMPeak(x: 230, y: 190, intensity: 1)], originX: 200, originY: 190,
                                               invAngstromPerPixel: 0.01) != nil, true)
        r.equal("malformed scores give the empty result", matcher.result(scores: [1], bins: [0], mirrored: [false]) == .empty, true)
        r.equal("a result from the pass scores is the empty result for an empty pass",
                matcher.result(forwardScores: [], forwardBins: [], mirrorScores: [], mirrorBins: []) == .empty, true)
        do {
            _ = try ACOMTemplateLibrary.generate(reflections: [ACOMReflection(g: [1, 0, 0], intensity: 1)], kMax: 1.5,
                                                 laueClass: .mThreeBarM, isCancelled: { true })
            r.check("generate honours cancellation", false)
        } catch { r.equal("generate honours cancellation", (error as? ACOMError) == .cancelled, true) }
        let empty = ACOMScanResults(width: 2, height: 2, templateCount: 0)
        r.equal("an empty map expands to itself", empty.expandedPreview(selection: .full) == empty, true)
    }
}
