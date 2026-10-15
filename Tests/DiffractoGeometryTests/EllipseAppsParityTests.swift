//
//  EllipseAppsParityTests.swift
//  DiffractoKit slice 3b (docs/registration-slice-3.md sections 2 and 3, last item): the kit's `EllipseCalibration` against the
//  apps' own compiled code, BIT FOR BIT. tools/slice-3b/apps-parity/run.sh fetches mac4DSTEM 7d6fdbd1's
//  `EllipseCalibration.swift` and `DiffractionPattern.swift` by `git show`, compiles them with swiftc -O beside a main.swift that
//  renders the planted rings with the SplitMix64 + Box-Muller renderer of EllipseNoisyRingTests, runs `fit1D`, `fitAmorphousRing`
//  and `fitBestAvailable` on each (six calls per ring: fit1D at the default 0.2 residual cap, at 0.4, at 0.4 with the anyway
//  path, the profile fit from the 0.4 conic, `fitBestAvailable` by default and with the anyway path) and dumps every output
//  field, Doubles as bit patterns, to Fixtures/slice-3b-apps.json.
//
//  This test (1) renders the same plants here from the parameters in the dump and requires the Float pixels to hash to the dump's
//  FNV-1a hash (so the renderer in this file IS the dump's), (2) runs the kit on them, (3) requires every Double equal in
//  bitPattern, every integer, flag, model and fallback reason equal, and every thrown error equal (its payload's bits included).
//  The ring sets: the three noisy plants of the apps' EllipseNoisyRingTests and four noisy 128 px rings (noise 0.3); the four
//  EllipseStopRuleTests rings; 20 clean sweep rings (square 112-144 px, R 32-52, eps 0 / 1 / 3 %, unbroken or a 60 degree arc
//  missing); four sparse rings (120, 150, 240 degree arcs missing and a noisy 90 degree one).
//
//  A difference is a failure of the port (operation order, Float against Double, a constant), never a tolerance to widen: the
//  test reports the difference in ulps and the first field.
//
//  RED UNDER: any constant, boundary or operation of EllipseCalibration changed (the break-once run of tools/slice-3b/mutants.py
//  records which of these tests, and which of the arithmetic ones, turn red for each mutant).
//

import Testing
import Foundation
import DiffractoGeometry

private struct Dump: Decodable {
    struct Plant: Decodable {
        let kind: String, name: String
        let size: Int, radius: Double, eps: Double, thetaDegrees: Double, gap: Double, noise: Double, seed: Int, tail: Double
    }
    struct Fit: Decodable {
        let centerQX: String, centerQY: String, a: String, b: String, theta: String
        let normalizedResidual: String, conicResidual: String
        let sampleCount: Int, occupiedAngularBins: Int
        let model: String
        let profile: [String]?
        let profileFallbackReason: String?
        let sparseCoverage: Bool
    }
    struct Result: Decodable { let error: String?, fit: Fit? }
    struct Case: Decodable {
        let plant: Plant
        let pixelHash: String
        let startQX: Double, startQY: Double, innerRadius: Double, outerRadius: Double
        let results: [String: Result]
    }

    /// The apps' code compiled -O (the Release configuration, `tools/test.sh`) or -Onone (a debug `swift test`): the two differ in
    /// the last bits of some profile fits, so each build of the kit is held bit for bit to the apps built the same way.
    #if DEBUG
    static let file = "slice-3b-apps-onone.json"
    #else
    static let file = "slice-3b-apps.json"
    #endif
    static let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/\(file)")
    static let cases: [Case] = {
        do { return try JSONDecoder().decode([Case].self, from: Data(contentsOf: url)) } catch { fatalError("\(file): \(error)") }
    }()
}

private func hex(_ x: Double) -> String {
    let s = String(x.bitPattern, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}

private func fnv(_ pixels: [Float]) -> String {
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

private func describe(_ error: Error) -> String {
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

struct EllipseAppsParityTests {

    private static func plant(_ p: Dump.Plant) -> EllipsePlant {
        EllipsePlant(size: p.size, radius: p.radius, eps: p.eps, thetaDegrees: p.thetaDegrees, gap: p.gap, noise: p.noise,
                     seed: p.seed, tail: p.tail)
    }

    /// The difference of two bit patterns' values, in ulps of the first, for the report.
    private static func ulps(_ bits: String, _ value: Double) -> String {
        guard let raw = UInt64(bits, radix: 16) else { return "unparsable" }
        let apps = Double(bitPattern: raw)
        guard apps.isFinite, value.isFinite, apps != 0 else { return "apps \(apps), kit \(value)" }
        return "apps \(apps), kit \(value), \((value - apps) / apps.ulp) ulps"
    }

    /// Every field equal to the dump's; returns the mismatches (empty when bit identical).
    private static func compare(_ name: String, _ fit: EllipseCalibrationFit, _ expected: Dump.Fit) -> [String] {
        var bad: [String] = []
        func check(_ field: String, _ value: Double, _ bits: String) {
            if hex(value) != bits { bad.append("\(name) \(field): \(ulps(bits, value))") }
        }
        check("centerQX", fit.centerQX, expected.centerQX)
        check("centerQY", fit.centerQY, expected.centerQY)
        check("a", fit.a, expected.a)
        check("b", fit.b, expected.b)
        check("theta", fit.theta, expected.theta)
        check("normalizedResidual", fit.normalizedResidual, expected.normalizedResidual)
        check("conicResidual", fit.conicResidual, expected.conicResidual)
        if fit.sampleCount != expected.sampleCount { bad.append("\(name) sampleCount \(fit.sampleCount) against \(expected.sampleCount)") }
        if fit.occupiedAngularBins != expected.occupiedAngularBins {
            bad.append("\(name) occupiedAngularBins \(fit.occupiedAngularBins) against \(expected.occupiedAngularBins)")
        }
        if fit.model.rawValue != expected.model { bad.append("\(name) model \(fit.model.rawValue) against \(expected.model)") }
        if fit.sparseCoverage != expected.sparseCoverage { bad.append("\(name) sparseCoverage") }
        if fit.profileFallbackReason != expected.profileFallbackReason {
            bad.append("\(name) profileFallbackReason \(fit.profileFallbackReason ?? "nil") against \(expected.profileFallbackReason ?? "nil")")
        }
        switch (fit.profile, expected.profile) {
        case (nil, nil): break
        case (let p?, let e?):
            let values = [p.centralIntensity, p.ringIntensity, p.centralSigma, p.innerSigma, p.outerSigma, p.background]
            let names = ["centralIntensity", "ringIntensity", "centralSigma", "innerSigma", "outerSigma", "background"]
            for i in 0..<6 { check("profile.\(names[i])", values[i], e[i]) }
        default: bad.append("\(name) profile presence")
        }
        return bad
    }

    /// Runs the six calls on one case and returns the mismatches.
    private static func run(_ c: Dump.Case) -> [String] {
        let plant = Self.plant(c.plant)
        let pattern = plant.render()
        var bad: [String] = []
        if fnv(pattern.pixels) != c.pixelHash { return ["\(c.plant.name): the rendered pixels differ from the dump's (hash)"] }
        let sx = c.startQX, sy = c.startQY, inner = c.innerRadius, outer = c.outerRadius
        func judge(_ key: String, _ body: () throws -> EllipseCalibrationFit) {
            guard let expected = c.results[key] else { bad.append("\(c.plant.name): no \(key) in the dump"); return }
            let label = "\(c.plant.name) \(key)"
            do {
                let fit = try body()
                if let e = expected.error { bad.append("\(label): kit fitted, apps threw \(e)"); return }
                guard let f = expected.fit else { bad.append("\(label): dump has neither fit nor error"); return }
                bad += compare(label, fit, f)
            } catch {
                guard let e = expected.error else { bad.append("\(label): kit threw \(describe(error)), apps fitted"); return }
                if describe(error) != e { bad.append("\(label): kit threw \(describe(error)), apps threw \(e)") }
            }
        }
        func judgeAgainst(_ expected: Dump.Result?, _ label: String, _ body: () throws -> EllipseCalibrationFit) {
            guard let expected else { bad.append("\(label): nothing expected"); return }
            do {
                let fit = try body()
                if let e = expected.error { bad.append("\(label): kit fitted, apps threw \(e)"); return }
                guard let f = expected.fit else { bad.append("\(label): dump has neither fit nor error"); return }
                bad += compare(label, fit, f)
            } catch {
                guard let e = expected.error else { bad.append("\(label): kit threw \(describe(error)), apps fitted"); return }
                if describe(error) != e { bad.append("\(label): kit threw \(describe(error)), apps threw \(e)") }
            }
        }
        judge("fit1D") {
            try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner, outerRadius: outer)
        }
        var wide: EllipseCalibrationFit?
        judge("fit1DWide") {
            let f = try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                                 outerRadius: outer, maximumResidual: 0.4)
            wide = f
            return f
        }
        judge("fit1DSparse") {
            try EllipseCalibration.fit1D(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner, outerRadius: outer,
                                         maximumResidual: 0.4, acceptSparseCoverage: true)
        }
        if let wide {
            judge("profile") {
                try EllipseCalibration.fitAmorphousRing(pattern: pattern, initial: wide, centerQX: sx, centerQY: sy,
                                                        innerRadius: inner, outerRadius: outer)
            }
        } else if c.results["profile"] != nil {
            bad.append("\(c.plant.name): the apps ran the profile fit, the kit's wide conic threw")
        }
        judge("best") {
            try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner, outerRadius: outer)
        }
        judge("bestSparse") {
            try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                                    outerRadius: outer, acceptSparseCoverage: true)
        }
        // DiffractoKit 0.0.4 (registration-0.0.4 §4): progress and a never-cancelling closure change no bit of "best".
        if c.results["best"] != nil {
            judgeAgainst(c.results["best"], "\(c.plant.name) bestWithClosures") {
                try EllipseCalibration.fitBestAvailable(pattern: pattern, centerQX: sx, centerQY: sy, innerRadius: inner,
                                                        outerRadius: outer, progress: { _ in }, isCancelled: { false })
            }
        }
        return bad
    }

    private static func indices(_ kind: String) -> [Int] { Dump.cases.indices.filter { Dump.cases[$0].plant.kind == kind } }

    /// The dump holds the planned sets, every call in each case, both outcomes (fits and throws) and both models.
    /// RED UNDER: a regenerated dump with a ring or a call missing.
    @Test func testTheDumpHoldsTheRegisteredSetsAndBothOutcomes() {
        #expect(Self.indices("noisy").count == 7)
        #expect(Self.indices("stoprule").count == 4)
        #expect(Self.indices("sweep").count == 20)
        #expect(Self.indices("gap").count == 4)
        let keys = Set(["fit1D", "fit1DWide", "fit1DSparse", "best", "bestSparse"])
        var errors = Set<String>(), models = Set<String>(), sparse = 0
        for c in Dump.cases {
            #expect(keys.isSubset(of: Set(c.results.keys)), "\(c.plant.name)")
            for r in c.results.values {
                if let e = r.error { errors.insert(String(e.prefix { $0 != "(" })) }
                if let f = r.fit { models.insert(f.model); if f.sparseCoverage { sparse += 1 } }
            }
        }
        #expect(errors.isSuperset(of: ["insufficientAngularCoverage", "excessiveResidual", "didNotConverge"]), "\(errors)")
        #expect(models == ["Intensity-weighted conic", "Asymmetric Gaussian ring"])
        #expect(sparse > 0)
    }

    /// The three noisy plants of the apps' EllipseNoisyRingTests and four strongly noisy 128 px rings (one argument each, so
    /// the rings run side by side).
    @Test(arguments: EllipseAppsParityTests.indices("noisy"))
    func testNoisyPlantsAreBitIdenticalWithTheApps(_ index: Int) {
        let bad = Self.run(Dump.cases[index])
        #expect(bad.isEmpty, "\(bad.joined(separator: "\n"))")
    }

    /// The four rings of the apps' EllipseStopRuleTests.
    @Test(arguments: EllipseAppsParityTests.indices("stoprule"))
    func testStopRuleRingsAreBitIdenticalWithTheApps(_ index: Int) {
        let bad = Self.run(Dump.cases[index])
        #expect(bad.isEmpty, "\(bad.joined(separator: "\n"))")
    }

    /// The 20 clean sweep rings.
    @Test(arguments: EllipseAppsParityTests.indices("sweep"))
    func testCleanSweepRingsAreBitIdenticalWithTheApps(_ index: Int) {
        let bad = Self.run(Dump.cases[index])
        #expect(bad.isEmpty, "\(bad.joined(separator: "\n"))")
    }

    /// The sparse rings: refused by coverage, fitted on the anyway path with the one-ring check, `didNotConverge`.
    @Test(arguments: EllipseAppsParityTests.indices("gap"))
    func testSparseRingsAreBitIdenticalWithTheApps(_ index: Int) {
        let bad = Self.run(Dump.cases[index])
        #expect(bad.isEmpty, "\(bad.joined(separator: "\n"))")
    }
}
