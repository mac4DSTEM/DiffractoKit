//
//  ACOMMatcherCallerAPITests.swift
//  Role: the four ACOM asks of the apps' kit switch (0.7.0): `generate` cancellation and progress, `result(scores:...)` from
//        scores a caller computed, `experimentalPolarImage` with one peak, and `matchAll` taking a peaks view. Each is checked
//        against the existing CPU route (`match`, `matchAll(peaks:)`), which is the reference; nothing here changes a number.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

private func cubicReflections() -> [ACOMReflection] {
    var out: [ACOMReflection] = []
    for h in -3...3 { for k in -3...3 { for l in -3...3 where (h, k, l) != (0, 0, 0) {
        out.append(ACOMReflection(g: SIMD3(Double(h), Double(k), Double(l)) / 4, intensity: 1 / Double(h * h + k * k + l * l)))
    } } }
    return out
}

private let settings = ACOMTemplateSettings(nRadial: 32, nAzimuthal: 64, zoneAxisCount: 24,
                                            excitation: .appCurvature(wavelengthAngstrom: 0.025))

private func library() throws -> ACOMTemplateLibrary {
    try ACOMTemplateLibrary.generate(reflections: cubicReflections(), kMax: 1.5, laueClass: .mThreeBarM, settings: settings)
}

private func pattern(_ lib: ACOMTemplateLibrary, template: Int, rotate: Double = 0.3) -> [ACOMPeak] {
    lib.templateSpots[template].map {
        let a = Double($0.azimuth) + rotate
        return ACOMPeak(x: Float(Double($0.radius) * cos(a)), y: Float(Double($0.radius) * sin(a)), intensity: $0.weight)
    }
}

struct ACOMMatcherCallerAPITests {
    @Test func generateReportsProgressAndHonoursCancellation() throws {
        final class Box: @unchecked Sendable { var fractions: [Double] = [] }
        let box = Box()
        let lib = try ACOMTemplateLibrary.generate(reflections: cubicReflections(), kMax: 1.5, laueClass: .mThreeBarM,
                                                   settings: settings, progress: { box.fractions.append($0) })
        #expect(box.fractions.count == lib.count)
        #expect(box.fractions.last == 1)
        #expect(box.fractions == box.fractions.sorted())
        // Cancel after the first direction: no library, the named error.
        final class Counter: @unchecked Sendable { var n = 0 }
        let counter = Counter()
        #expect(throws: ACOMError.cancelled) {
            try ACOMTemplateLibrary.generate(reflections: cubicReflections(), kMax: 1.5, laueClass: .mThreeBarM, settings: settings,
                                             isCancelled: { counter.n += 1; return counter.n > 2 })
        }
        #expect(throws: ACOMError.cancelled) {
            try ACOMTemplateLibrary.generate(reflections: cubicReflections(), kMax: 1.5, laueClass: .mThreeBarM, settings: settings,
                                             isCancelled: { true })
        }
    }

    @Test func resultFromScoresEqualsMatch() throws {
        let lib = try library()
        let matcher = try #require(ACOMMatcher(library: lib))
        let peaks = pattern(lib, template: 7)
        let direct = matcher.match(peaks: peaks, originX: 0, originY: 0)
        #expect(direct.templateIndex == 7 || direct.score > 0.9)
        // The scores a caller would hold: the same per-template scores, with the bin of the winner.
        let scores = matcher.templateScores(peaks: peaks, originX: 0, originY: 0)
        #expect(scores.count == lib.count)
        let winner = direct.templateIndex
        let na = lib.geometry.nAzimuthal
        // makeResult's bin is the correlation peak at -theta: recover it from the reported angle.
        var bins = [UInt32](repeating: 0, count: lib.count)
        bins[winner] = UInt32((na - Int((Double(direct.inPlaneAngle - (direct.mirrored ? .pi : 0)) / (2 * .pi) * Double(na)).rounded())
                               % na + na) % na)
        var mirrored = [Bool](repeating: false, count: lib.count)
        mirrored[winner] = direct.mirrored
        let rebuilt = matcher.result(scores: scores, bins: bins, mirrored: mirrored)
        #expect(rebuilt.templateIndex == direct.templateIndex)
        #expect(rebuilt.score == direct.score)
        #expect(rebuilt.secondScore == direct.secondScore)
        #expect(rebuilt.mirrored == direct.mirrored)
        #expect(abs(rebuilt.inPlaneAngle - direct.inPlaneAngle) < 1e-5)
        #expect(rebuilt.euler == direct.euler)
    }

    @Test func resultTwoPassFormTakesTheMirrorOnlyWhenStrictlyLarger() throws {
        let lib = try library()
        let matcher = try #require(ACOMMatcher(library: lib))
        let n = lib.count
        var fwd = [Float](repeating: 0.1, count: n), mir = fwd
        let fb = [UInt32](repeating: 3, count: n), mb = [UInt32](repeating: 5, count: n)
        fwd[2] = 0.8; mir[2] = 0.8          // tie: forward
        let tie = matcher.result(forwardScores: fwd, forwardBins: fb, mirrorScores: mir, mirrorBins: mb)
        #expect(tie.templateIndex == 2 && tie.mirrored == false)
        mir[2] = 0.8001                      // strictly larger: mirrored
        let win = matcher.result(forwardScores: fwd, forwardBins: fb, mirrorScores: mir, mirrorBins: mb)
        #expect(win.templateIndex == 2 && win.mirrored == true)
        // Malformed input is the empty result, never a trap.
        #expect(matcher.result(scores: [], bins: [], mirrored: []) == .empty)
        #expect(matcher.result(scores: [1], bins: [0], mirrored: [false]) == .empty)
        #expect(matcher.result(scores: [Float](repeating: 1, count: n), bins: [UInt32](repeating: 9999, count: n),
                               mirrored: [Bool](repeating: false, count: n)) == .empty)
    }

    @Test func oneBraggPeakDrawsThePolarImageButDoesNotMatch() throws {
        let lib = try library()
        let matcher = try #require(ACOMMatcher(library: lib))
        let one = [ACOMPeak(x: 20, y: 5, intensity: 1)]
        let image = try #require(matcher.experimentalPolarImage(peaks: one, originX: 0, originY: 0, invAngstromPerPixel: 0.02))
        #expect(image.count == lib.geometry.nRadial * lib.geometry.nAzimuthal)
        #expect(image.contains { $0 != 0 })
        // Two peaks give the same image code (the first of them unchanged).
        #expect(matcher.experimentalPolarImage(peaks: [], originX: 0, originY: 0) == nil)
        #expect(matcher.match(peaks: one, originX: 0, originY: 0) == .empty)
        #expect(matcher.experimentalFFT(peaks: one, originX: 0, originY: 0) == nil)
        #expect(matcher.templateScores(peaks: one, originX: 0, originY: 0).isEmpty)
    }

    @Test func matchAllTakesAPeaksViewAndEqualsTheArrayForm() throws {
        let lib = try library()
        let matcher = try #require(ACOMMatcher(library: lib))
        let width = 4, height = 3
        let grid = (0..<(width * height)).map { pattern(lib, template: ($0 * 5) % lib.count, rotate: 0.1 * Double($0)) }
        let reference = try #require(matcher.matchAll(peaks: grid, width: width, height: height, originX: 0, originY: 0,
                                                      maxWorkers: 2))
        // A view over another element type, no [[ACOMPeak]] built.
        struct Foreign: Sendable { var px: Float, py: Float, w: Float }
        let foreign = grid.map { $0.map { Foreign(px: $0.x, py: $0.y, w: $0.intensity) } }
        let viewedRun = matcher.matchAll(width: width, height: height, originX: 0, originY: 0, maxWorkers: 2,
                                         peaksAt: { i in foreign[i].lazy.map { ACOMPeak(x: $0.px, y: $0.py, intensity: $0.w) } })
        let viewed = try #require(viewedRun)
        #expect(viewed == reference)
        // Only the selected positions are asked for.
        final class Seen: @unchecked Sendable { let lock = NSLock(); var indices = Set<Int>() }
        let seen = Seen()
        _ = matcher.matchAll(width: width, height: height, originX: 0, originY: 0, selection: .square(centerX: 1, centerY: 1, radius: 0),
                             peaksAt: { i in seen.lock.lock(); seen.indices.insert(i); seen.lock.unlock(); return grid[i] })
        #expect(seen.indices == [1 * width + 1])
    }
}
