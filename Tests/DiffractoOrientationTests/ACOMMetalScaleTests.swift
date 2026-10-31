//
//  ACOMMetalScaleTests.swift
//  Role: the Metal option of `ACOMMatcher.matchAll` (P31) at other azimuthal lengths and scan sizes than `ACOMMetalTests`,
//        and its timing.
//        (1) `ACOMMetalAzimuthalLengthTests` (runs by default): the Metal scores equal the CPU reference within 1e-6, and the
//            template, in-plane angle and mirror flag are identical, at 4 to 1 024 azimuthal bins (the inverse transform is a
//            radix-2 FFT on the GPU) and at a scan whose size is not a multiple of the 8 positions one spectrum thread serves.
//        (2) `ACOMMetalBenchmarkTests` (skipped by default): ms per pattern at a 128 x 128 scan (16 384 patterns) and at 600
//            and 2 400 beam directions. Run on a quiet machine with
//              DIFFRACTOKIT_ACOM_BENCHMARK=1 swift test -c release --filter ACOMMetalBenchmarkTests
//            The CPU path is timed on the first 1 024 patterns (each position is independent, so its per-pattern cost does
//            not depend on the scan size). `DIFFRACTOKIT_ACOM_BENCHMARK_WORKERS` sets the CPU workers of both paths (default 1).
//

import Foundation
import Testing
import DiffractoCore
import DiffractoOrientation

private let benchmarkEnabled = ProcessInfo.processInfo.environment["DIFFRACTOKIT_ACOM_BENCHMARK"] != nil

/// The fcc crystal of tools/acom-benchmark: a = 4.05 angstrom, |g| < 1.6, intensity 1 / |g|^2.
private func fccReflections() -> [ACOMReflection] {
    var reflections: [ACOMReflection] = []
    for h in -4...4 {
        for k in -4...4 {
            for l in -4...4 where (h, k, l) != (0, 0, 0) && (h + k) % 2 == 0 && (k + l) % 2 == 0 {
                let g = SIMD3(Double(h), Double(k), Double(l)) / 4.05
                let norm = (g * g).sum().squareRoot()
                if norm < 1.6 { reflections.append(ACOMReflection(g: g, intensity: 1 / (norm * norm))) }
            }
        }
    }
    return reflections
}

/// `count` patterns from the library's own templates (every seventh), each turned by a seeded random in-plane angle.
private func scanPeaks(_ library: ACOMTemplateLibrary, count: Int, ipp: Double) -> [[ACOMPeak]] {
    var seed: UInt64 = 99
    func rnd() -> Double { seed = seed &* 6364136223846793005 &+ 1442695040888963407; return Double(seed >> 11) / Double(1 << 53) }
    return (0..<count).map { i in
        let phi = rnd() * 2 * .pi
        return library.templateSpots[(i * 7) % library.count].map { s in
            let az = Double(s.azimuth) + phi, r = Double(s.radius) / ipp
            return ACOMPeak(x: Float(128 + r * cos(az)), y: Float(128 + r * sin(az)), intensity: Float(s.weight))
        }
    }
}

struct ACOMMetalAzimuthalLengthTests {
    @Test func metalEqualsTheCPUAtOtherAzimuthalLengthsAndAPartialPositionBlock() throws {
        var passTies = 0, compared = 0
        defer {
            print("ACOM Metal versus CPU at 4 to 1 024 bins: \(compared) positions, pass ties \(passTies)")
            #expect(passTies <= 6, "pass ties are the exception: \(passTies)")
        }
        for bins in [4, 8, 16, 64, 256, 512, 1024] {
            let library = try ACOMTemplateLibrary.generate(reflections: fccReflections(), kMax: 1.6, laueClass: .mThreeBarM,
                                                           settings: ACOMTemplateSettings(nAzimuthal: bins, zoneAxisCount: 40))
            let matcher = try #require(ACOMMatcher(library: library))
            guard matcher.isMetalAvailable else { return }
            let peaks = scanPeaks(library, count: 13, ipp: 0.004)   // 13 = one full block of 8 and a block of 5
            let cpu = try #require(matcher.matchAll(peaks: peaks, width: 13, height: 1, originX: 128, originY: 128,
                                                    invAngstromPerPixel: 0.004, maxWorkers: 1))
            let metal = try #require(matcher.matchAll(peaks: peaks, width: 13, height: 1, originX: 128, originY: 128,
                                                      invAngstromPerPixel: 0.004, maxWorkers: 1, backend: .metal))
            for (a, b) in zip(cpu.results, metal.results) {
                compared += 1
                #expect(abs(a.score - b.score) <= 1e-6, "\(bins) bins: score \(a.score) against \(b.score)")
                #expect(abs(a.secondScore - b.secondScore) <= 1e-6, "\(bins) bins: runner-up \(a.secondScore) against \(b.secondScore)")
                #expect(a.templateIndex == b.templateIndex, "\(bins) bins: template \(a.templateIndex) against \(b.templateIndex)")
                // a wrong FFT output order would keep the best score but move its bin. A mirror-symmetric pattern scores the
                // same in both passes (P31), so a pass flip at an equal score is a tie; within one pass the bin must agree.
                if a.mirrored == b.mirrored {
                    #expect(a.inPlaneAngle == b.inPlaneAngle, "\(bins) bins: angle \(a.inPlaneAngle) against \(b.inPlaneAngle)")
                } else {
                    passTies += 1
                }
            }
        }
    }
}

struct ACOMMetalBenchmarkTests {
    private func seconds(_ body: () -> Void) -> Double {
        let t0 = DispatchTime.now().uptimeNanoseconds
        body()
        return Double(DispatchTime.now().uptimeNanoseconds - t0) / 1e9
    }

    @Test(.enabled(if: benchmarkEnabled)) func metalAgainstCPUAt16384PatternsAnd600Or2400Axes() throws {
        let workers = Int(ProcessInfo.processInfo.environment["DIFFRACTOKIT_ACOM_BENCHMARK_WORKERS"] ?? "") ?? 1
        let reflections = fccReflections()
        let side = 128, ipp = 0.004, cpuCount = 1024
        for axes in [600, 2400] {
            let library = try ACOMTemplateLibrary.generate(reflections: reflections, kMax: 1.6, laueClass: .mThreeBarM,
                                                           settings: ACOMTemplateSettings(zoneAxisCount: axes))
            let matcher = try #require(ACOMMatcher(library: library))
            guard matcher.isMetalAvailable else { print("ACOM benchmark: no Metal device"); return }
            let peaks = scanPeaks(library, count: side * side, ipp: ipp)
            let head = Array(peaks.prefix(cpuCount))
            var cpu: ACOMScanResults?, metalHead: ACOMScanResults?, metal: ACOMScanResults?
            let tCPU = seconds {
                cpu = matcher.matchAll(peaks: head, width: 32, height: 32, originX: 128, originY: 128, invAngstromPerPixel: ipp,
                                       maxWorkers: workers)
            }
            _ = matcher.matchAll(peaks: Array(peaks.prefix(64)), width: 8, height: 8, originX: 128, originY: 128,
                                 invAngstromPerPixel: ipp, maxWorkers: workers, backend: .metal)   // warm-up: shader compile
            let tMetalHead = seconds {
                metalHead = matcher.matchAll(peaks: head, width: 32, height: 32, originX: 128, originY: 128,
                                             invAngstromPerPixel: ipp, maxWorkers: workers, backend: .metal)
            }
            let tMetal = seconds {
                metal = matcher.matchAll(peaks: peaks, width: side, height: side, originX: 128, originY: 128,
                                         invAngstromPerPixel: ipp, maxWorkers: workers, backend: .metal)
            }
            let a = try #require(cpu).results, b = try #require(metalHead).results
            #expect(try #require(metal).results.count == side * side)
            let worst = zip(a, b).map { abs(Double($0.score) - Double($1.score)) }.max() ?? 0
            print(String(format: "ACOM benchmark %d axes, %d CPU workers: CPU %.3f ms/pattern (1 024), Metal %.3f ms/pattern (1 024), Metal %.3f ms/pattern (16 384); Metal-CPU worst score gap %.2e",
                         axes, workers, tCPU / Double(cpuCount) * 1000, tMetalHead / Double(cpuCount) * 1000,
                         tMetal / Double(side * side) * 1000, worst))
            #expect(worst <= 1e-6)
        }
    }
}
