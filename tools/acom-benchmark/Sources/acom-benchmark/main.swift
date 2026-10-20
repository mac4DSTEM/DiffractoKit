// Times ACOMMatcher.matchAll on a synthetic fcc crystal (a = 4.05 Å, the 111/200/220/311/222/400/331/420 families), 600 beam
// directions (the default library), 32 x 128 polar bins, patterns built from a template's own spots with a random in-plane
// angle. Speed only: the numbers it returns are pinned by the tests, not here.
import Foundation
import DiffractoCore
import DiffractoOrientation

func now() -> Double { Double(DispatchTime.now().uptimeNanoseconds) / 1e9 }

let args = CommandLine.arguments
let width = args.count > 1 ? Int(args[1]) ?? 64 : 64
let height = args.count > 2 ? Int(args[2]) ?? 64 : 64

let a = 4.05
var reflections: [ACOMReflection] = []
for h in -4...4 { for k in -4...4 { for l in -4...4 {
    guard (h, k, l) != (0, 0, 0), (h + k) % 2 == 0, (k + l) % 2 == 0 else { continue }    // fcc allowed
    let g = SIMD3(Double(h), Double(k), Double(l)) / a
    let norm = (g * g).sum().squareRoot()
    guard norm < 1.6 else { continue }
    reflections.append(ACOMReflection(g: g, intensity: 1 / (norm * norm)))
} } }

var t0 = now()
let library = try ACOMTemplateLibrary.generate(reflections: reflections, kMax: 1.6, laueClass: .mThreeBarM,
                                               settings: ACOMTemplateSettings())
let genTime = now() - t0
guard let matcher = ACOMMatcher(library: library) else { fatalError("no FFT plan") }
print("library: \(library.count) templates, \(reflections.count) reflections, generated in \(String(format: "%.2f", genTime)) s")

var rng = SystemRandomNumberGenerator()
let ipp = 0.004
let peaks: [[ACOMPeak]] = (0..<(width * height)).map { i in
    let phi = Double.random(in: 0..<(2 * .pi), using: &rng)
    return library.templateSpots[(i * 7) % library.count].map { s in
        let az = Double(s.azimuth) + phi, r = Double(s.radius) / ipp
        return ACOMPeak(x: Float(128 + r * cos(az)), y: Float(128 + r * sin(az)), intensity: Float(s.weight))
    }
}
let meanPeaks = Double(peaks.reduce(0) { $0 + $1.count }) / Double(peaks.count)
print("scan: \(width) x \(height) = \(width * height) patterns, \(String(format: "%.1f", meanPeaks)) peaks each, cores: \(ProcessInfo.processInfo.activeProcessorCount)")

@MainActor func run(_ label: String, workers: Int?, backend: ACOMBackend = .cpu, cancel: (@Sendable () -> Bool)? = nil) {
    t0 = now()
    let r = matcher.matchAll(peaks: peaks, width: width, height: height, originX: 128, originY: 128,
                             invAngstromPerPixel: ipp, maxWorkers: workers, backend: backend, isCancelled: cancel)
    let dt = now() - t0
    precondition(r != nil)
    let per = dt / Double(width * height) * 1000
    print(String(format: "%-28@ %7.2f s   %6.2f ms/pattern wall   %7.1f patterns/s   256x256 scan: %.1f min",
                 label as NSString, dt, per, Double(width * height) / dt, 65536 * per / 1000 / 60))
}
run("1 worker", workers: 1)
run("2 workers", workers: 2)
run("4 workers", workers: 4)
run("8 workers", workers: 8)
run("8 workers + cancel poll", workers: 8, cancel: { false })
if matcher.isMetalAvailable {
    run("Metal, 8 prep workers", workers: 8, backend: .metal)
    run("Metal, 4 prep workers", workers: 4, backend: .metal)
} else {
    print("Metal: not available on this machine")
}
