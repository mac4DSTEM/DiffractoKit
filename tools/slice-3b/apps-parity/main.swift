// DiffractoKit slice 3b, oracle (registration docs/registration-slice-3.md section 3, last item): a scratch tool that links the
// apps' own EllipseCalibration.swift and DiffractionPattern.swift (mac4DSTEM 7d6fdbd1, fetched by run.sh with `git show`,
// compiled beside this file with swiftc -O) and runs `fit1D`, `fitAmorphousRing` and `fitBestAvailable` on planted rings
// rendered HERE by the SplitMix64 + Box-Muller renderer of mac4DSTEM's EllipseNoisyRingTests (the same code as
// Tests/DiffractoGeometryTests/EllipseNoisyRingTests.swift, copied by hand; the dump carries the pixels' FNV-1a hash so the
// test proves it renders the same Float pixels). Every output field is dumped as Double / Float bit patterns.
//
//   apps-parity <out.json> [<sweep-dir>]       (a sweep dir gets the clean sweep rings' raw Float32 files, for sweep.py)
//
// Written from the registration and the seed, not from the kit. No kit code is read or linked.

import Foundation

// MARK: The plant (identical in the kit's test file)

struct Plant: Codable {
    var kind: String
    var name: String
    var size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, noise: Double, seed: Int
    var tail: Double
    var b: Double { 2 * radius / (2 + eps) }
    var a: Double { b * (1 + eps) }
    var centerQX: Double { Double(size) / 2 + 0.37 + 2.5 * sin(1.3 * Double(seed) + 0.2) }
    var centerQY: Double { Double(size) / 2 - 0.41 + 2.5 * cos(1.9 * Double(seed) + 0.7) }

    func render() -> DiffractionPattern {
        let theta = thetaDegrees * .pi / 180
        let sin2 = sin(theta) * sin(theta), cos2 = cos(theta) * cos(theta)
        let A = sin2 / (b * b) + cos2 / (a * a)
        let C = cos2 / (b * b) + sin2 / (a * a)
        let B = 2 * (b * b - a * a) * sin(theta) * cos(theta) / (a * a * b * b)
        let mean = (a + b) / 2
        let gapStart = Double((47 * seed) % 360)
        var state = UInt64(seed) &+ 1
        func uniform() -> Double {
            state &+= 0x9E37_79B9_7F4A_7C15
            var z = state
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return (Double((z ^ (z >> 31)) >> 11) + 0.5) / 9_007_199_254_740_992.0
        }
        var pixels = [Float](repeating: 0, count: size * size)
        var spare = 0.0
        for qx in 0..<size {
            for qy in 0..<size {
                let dx = Double(qx) - centerQX, dy = Double(qy) - centerQY
                let rho = max(0, A * dx * dx + B * dx * dy + C * dy * dy).squareRoot()
                let distance = (rho - 1) * mean
                var value = exp(-(distance * distance) / (2 * 3.0 * 3.0))
                if gap > 0 {
                    var angle = atan2(dy, dx) * 180 / .pi
                    angle = (angle - gapStart).truncatingRemainder(dividingBy: 360)
                    if angle < 0 { angle += 360 }
                    if angle < gap { value = 0 }
                }
                if tail > 0 {
                    let r = max(hypot(dx, dy), 1)
                    value += tail * (radius / r) * (radius / r)
                }
                if noise > 0 {
                    if (qx * size + qy) % 2 == 0 {
                        let u1 = uniform(), u2 = uniform()
                        let r = (-2 * log(u1)).squareRoot()
                        value += noise * r * cos(2 * .pi * u2)
                        spare = noise * r * sin(2 * .pi * u2)
                    } else {
                        value += spare
                    }
                }
                pixels[qx * size + qy] = Float(value)
            }
        }
        return DiffractionPattern(qy: size, qx: size, pixels: pixels)
    }
}

// MARK: Dump format (decoded by Tests/DiffractoGeometryTests/EllipseAppsParityTests.swift)

func hex(_ x: Double) -> String {
    let s = String(x.bitPattern, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}

func fnv(_ pixels: [Float]) -> String {
    var h: UInt64 = 0xcbf2_9ce4_8422_2325
    for p in pixels {
        var v = p.bitPattern
        for _ in 0..<4 {
            h = (h ^ UInt64(v & 0xff)) &* 0x0000_0100_0000_01b3
            v >>= 8
        }
    }
    let s = String(h, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}

struct FitOut: Codable {
    var centerQX: String, centerQY: String, a: String, b: String, theta: String
    var normalizedResidual: String, conicResidual: String
    var sampleCount: Int, occupiedAngularBins: Int
    var model: String
    /// centralIntensity, ringIntensity, centralSigma, innerSigma, outerSigma, background
    var profile: [String]?
    var profileFallbackReason: String?
    var sparseCoverage: Bool
}

struct Result: Codable {
    var error: String?
    var fit: FitOut?
}

struct Case: Codable {
    var plant: Plant
    var pixelHash: String
    var startQX: Double, startQY: Double, innerRadius: Double, outerRadius: Double
    /// fit1D (default residual 0.2), fit1DWide (0.4), fit1DSparse (0.4, acceptSparseCoverage), profile (from fit1DWide),
    /// best (fitBestAvailable defaults), bestSparse (fitBestAvailable, acceptSparseCoverage).
    var results: [String: Result]
}

func describe(_ error: Error) -> String {
    switch error as? EllipseCalibration.FitError {
    case .invalidInput(let m)?: return "invalidInput(\(m))"
    case .insufficientSignal?: return "insufficientSignal"
    case .insufficientAngularCoverage(let n)?: return "insufficientAngularCoverage(\(n))"
    case .moreThanOneRing(let lo, let hi)?: return "moreThanOneRing(\(hex(lo)),\(hex(hi)))"
    case .didNotConverge?: return "didNotConverge"
    case .invalidEllipse?: return "invalidEllipse"
    case .excessiveResidual(let r)?: return "excessiveResidual(\(hex(r)))"
    case .ringWidthAtFloor?: return "ringWidthAtFloor"
    case nil: return "other(\(error))"
    }
}

func out(_ fit: EllipseCalibrationFit) -> FitOut {
    FitOut(centerQX: hex(fit.centerQX), centerQY: hex(fit.centerQY), a: hex(fit.a), b: hex(fit.b), theta: hex(fit.theta),
           normalizedResidual: hex(fit.normalizedResidual), conicResidual: hex(fit.conicResidual),
           sampleCount: fit.sampleCount, occupiedAngularBins: fit.occupiedAngularBins, model: fit.model.rawValue,
           profile: fit.profile.map { [hex($0.centralIntensity), hex($0.ringIntensity), hex($0.centralSigma),
                                      hex($0.innerSigma), hex($0.outerSigma), hex($0.background)] },
           profileFallbackReason: fit.profileFallbackReason, sparseCoverage: fit.sparseCoverage)
}

func run(_ body: () throws -> EllipseCalibrationFit) -> Result {
    do { return Result(error: nil, fit: out(try body())) } catch { return Result(error: describe(error), fit: nil) }
}

// MARK: The cases

var plants: [Plant] = []

// The three noisy rings of EllipseNoisyRingTests (seed: mac4DSTEM 7d6fdbd1).
plants.append(Plant(kind: "noisy", name: "dead ring n192-R70-e10-t0.0-g60-s100-k3", size: 192, radius: 70, eps: 0.01,
                    thetaDegrees: 0, gap: 60, noise: 0.10, seed: 3, tail: 0))
plants.append(Plant(kind: "noisy", name: "drifting ring n192-R70-e10-t112.5-g0-s100-k2", size: 192, radius: 70, eps: 0.01,
                    thetaDegrees: 112.5, gap: 0, noise: 0.10, seed: 2, tail: 0))
plants.append(Plant(kind: "noisy", name: "tail ring n192-R65-e10-t112.5-g0-s100-k5-T80", size: 192, radius: 65, eps: 0.01,
                    thetaDegrees: 112.5, gap: 0, noise: 0.10, seed: 5, tail: 0.8))

// The four clean rings of EllipseStopRuleTests (96 x 96, R 30, s 3; index = seed).
for (index, eps, thetaDegrees) in [(1, 0.0, 0.0), (4, 0.01, 20.0), (37, 0.03, 112.5), (43, 0.03, 160.0)] {
    plants.append(Plant(kind: "stoprule", name: "stoprule ring \(index), eps \(eps), theta \(thetaDegrees)", size: 96,
                        radius: 30, eps: eps, thetaDegrees: thetaDegrees, gap: 0, noise: 0, seed: index, tail: 0))
}

// 20 clean sweep rings, the shape of mac4DSTEM tools/ellipse-calibration-test's clean sweep at a size the kit's Debug tests
// can fit: square 112, 128 and 144 px, R 32-52 px, eps 0 / 1 / 3 %, unbroken or a 60 degree arc missing, theta from the sweep's
// list, s 3 px, a deterministic seed each (the seed places the centre and the gap).
let sweepThetas = [20.0, 45.0, 67.5, 90.0, 112.5, 135.0, 160.0]
var sweepIndex = 0
for (size, radius) in [(112, 32.0), (128, 40.0), (144, 52.0)] {
    for eps in [0.0, 0.01, 0.03] {
        for gap in [0.0, 60.0] {
            if sweepIndex >= 20 { continue }
            let theta = eps > 0 ? sweepThetas[sweepIndex % sweepThetas.count] : 0
            plants.append(Plant(kind: "sweep", name: "sweep \(sweepIndex) n\(size) R\(Int(radius)) e\(eps) t\(theta) g\(Int(gap))",
                                size: size, radius: radius, eps: eps, thetaDegrees: theta, gap: gap, noise: 0,
                                seed: 100 + sweepIndex, tail: 0))
            sweepIndex += 1
        }
    }
}
// two more with eps 3 % and a different theta so the count is 20
while sweepIndex < 20 {
    let size = 128, radius = 44.0, theta = sweepThetas[sweepIndex % sweepThetas.count]
    plants.append(Plant(kind: "sweep", name: "sweep \(sweepIndex) n\(size) R\(Int(radius)) e0.03 t\(theta) g0", size: size,
                        radius: radius, eps: 0.03, thetaDegrees: theta, gap: 0, noise: 0, seed: 100 + sweepIndex, tail: 0))
    sweepIndex += 1
}

// Sparse and refused rings: a 120 and a 150 degree arc missing (24 and 19 bins: refused by default, fitted on the anyway path
// with the one-ring check), and strongly noisy 128 px rings, whose outcomes (a profile, the conic fallback, a refusal) are
// whatever the seed's arithmetic gives; the dump records each.
plants.append(Plant(kind: "gap", name: "gap 120 n128 R40", size: 128, radius: 40, eps: 0.02, thetaDegrees: 30, gap: 120,
                    noise: 0, seed: 201, tail: 0))
plants.append(Plant(kind: "gap", name: "gap 150 n128 R40", size: 128, radius: 40, eps: 0.02, thetaDegrees: 80, gap: 150,
                    noise: 0, seed: 202, tail: 0))
plants.append(Plant(kind: "gap", name: "gap 240 n128 R40", size: 128, radius: 40, eps: 0.02, thetaDegrees: 10, gap: 240,
                    noise: 0, seed: 204, tail: 0))
plants.append(Plant(kind: "gap", name: "gap 90 n128 R44 noisy", size: 128, radius: 44, eps: 0.03, thetaDegrees: 140, gap: 90,
                    noise: 0.1, seed: 203, tail: 0))
for seed in 1...4 {
    plants.append(Plant(kind: "noisy", name: "noisy n128 R40 noise 0.3 seed \(seed)", size: 128, radius: 40, eps: 0.02,
                        thetaDegrees: Double(25 * seed), gap: 0, noise: 0.3, seed: seed, tail: 0))
}

let args = CommandLine.arguments
let sweepDir = args.count > 2 ? args[2] : nil
if let d = sweepDir { try FileManager.default.createDirectory(atPath: d, withIntermediateDirectories: true) }

var cases: [Case] = []
for (n, plant) in plants.enumerated() {
    let pattern = plant.render()
    if plant.kind == "sweep", let d = sweepDir {
        let data = pattern.pixels.withUnsafeBufferPointer { Data(buffer: $0) }
        try data.write(to: URL(fileURLWithPath: d).appendingPathComponent("sweep-\(n).f32"))
    }
    let sx = plant.centerQX + 1.3, sy = plant.centerQY - 0.9
    let inner = plant.b - 20, outer = plant.a + 20
    var results: [String: Result] = [:]
    results["fit1D"] = run {
        try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner, outerRadius: outer)
    }
    var wide: EllipseCalibrationFit?
    results["fit1DWide"] = run {
        let f = try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                             outerRadius: outer, maximumResidual: 0.4)
        wide = f
        return f
    }
    results["fit1DSparse"] = run {
        try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner, outerRadius: outer,
                                     maximumResidual: 0.4, acceptSparseCoverage: true)
    }
    if let wide {
        results["profile"] = run {
            try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: wide, centerQX: sx, centerQY: sy,
                                                    innerRadius: inner, outerRadius: outer)
        }
    }
    results["best"] = run {
        try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                                outerRadius: outer)
    }
    results["bestSparse"] = run {
        try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                                outerRadius: outer, acceptSparseCoverage: true)
    }
    cases.append(Case(plant: plant, pixelHash: fnv(pattern.pixels), startQX: sx, startQY: sy, innerRadius: inner,
                      outerRadius: outer, results: results))
    var summary = "apps-parity: \(plant.name):"
    for key in results.keys.sorted() {
        let r = results[key]!
        summary += " \(key) " + (r.error ?? (r.fit!.model == "Asymmetric Gaussian ring" ? "profile" : "conic")) + ";"
    }
    FileHandle.standardError.write(Data((summary + "\n").utf8))
}

let encoder = JSONEncoder()
encoder.outputFormatting = [.sortedKeys, .prettyPrinted]
let path = args.count > 1 ? args[1] : "slice-3b-apps.json"
try encoder.encode(cases).write(to: URL(fileURLWithPath: path))
