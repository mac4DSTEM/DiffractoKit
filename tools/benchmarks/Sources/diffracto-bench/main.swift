// Times the DiffractoKit functions that have a Python counterpart (docs/benchmarks.md). Inputs are the files gen_inputs.py wrote
// under inputs/; the Python side (bench_pymatgen_orix.py, bench_py4dstem.py) reads the same files. One JSON line per
// measurement: {"side":"kit","name","size","n","backend","median_ms","min_ms","max_ms","note"}. Timing: one warm-up run, then
// five timed runs; median and the spread (min, max) are reported. Large sizes run only with DIFFRACTO_BENCH_LARGE=1 (the nightly).
import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry
import DiffractoOrientation

let inputs = URL(fileURLWithPath: CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "inputs")
let large = ProcessInfo.processInfo.environment["DIFFRACTO_BENCH_LARGE"] == "1"
let only = ProcessInfo.processInfo.environment["DIFFRACTO_BENCH_ONLY"]

func now() -> Double { Double(DispatchTime.now().uptimeNanoseconds) / 1e6 }

nonisolated(unsafe) var sink = 0.0

func measure(_ name: String, size: String, n: Int, backend: String = "cpu", note: String = "", _ body: () -> Double) {
    if let only, !name.hasPrefix(only) { return }
    if size == "large" && !large {
        print(#"{"side":"kit","name":"\#(name)","size":"large","n":\#(n),"backend":"\#(backend)","median_ms":null,"note":"nightly"}"#)
        return
    }
    sink += body()       // warm-up
    var t: [Double] = []
    for _ in 0..<5 { let t0 = now(); sink += body(); t.append(now() - t0) }
    t.sort()
    print(#"{"side":"kit","name":"\#(name)","size":"\#(size)","n":\#(n),"backend":"\#(backend)","median_ms":\#(t[2]),"min_ms":\#(t[0]),"max_ms":\#(t[4]),"note":"\#(note)"}"#)
}

func text(_ f: String) -> String { try! String(contentsOf: inputs.appendingPathComponent(f), encoding: .utf8) }

// ---- Niggli reduction of N random cells (pymatgen / spglib: get_niggli_reduced_lattice / niggli_reduce) -------------------------
let cells = try! JSONSerialization.jsonObject(with: Data(contentsOf: inputs.appendingPathComponent("cells.json"))) as! [[Double]]
for (size, n) in [("small", 100), ("typical", 1000), ("large", 10000)] {
    measure("niggli_reduce", size: size, n: n) {
        var s = 0.0
        for c in cells.prefix(n) {
            s += UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5]).niggliReduced().cell.a
        }
        return s
    }
}

// ---- Reflection list inside k_max (pymatgen: get_points_in_sphere on the reciprocal lattice) --------------------------------
let al1 = try! CIFImport.read(text("al_super1.cif"), fileBaseName: "al").structure
for (size, kMax) in [("small", 1.0), ("typical", 2.0), ("large", 4.0)] {
    measure("reflections", size: size, n: Int(kMax * 10), note: "k_max \(kMax)") {
        Double(try! al1.reflections(kMax: kMax).reflections.count)
    }
}

// ---- Kinematic reflections: F, |F|^2 cut (py4DSTEM: Crystal.calculate_structure_factors) ------------------------------------
let c2m = try! CIFImport.read(text("c2m.cif"), fileBaseName: "c2m").structure
for (size, kMax) in [("small", 1.0), ("typical", 2.0), ("large", 3.0)] {
    measure("kinematic_reflections", size: size, n: Int(kMax * 10), note: "C2/m cell, k_max \(kMax)") {
        Double(try! c2m.kinematicReflections(kMax: kMax, retainingAll: false).kept.count)
    }
}

// ---- Symmetry detection from the atoms (pymatgen: SpacegroupAnalyzer.get_symmetry_operations) --------------------------------
for (size, n) in [("small", 1), ("typical", 2), ("large", 3)] {
    let s = try! CIFImport.read(text("al_super\(n).cif"), fileBaseName: "al").structure
    measure("symmetry_detection", size: size, n: s.sites.count, note: "fcc Al supercell \(n)^3 as P1 (\(s.sites.count) atoms)") {
        Double(CrystalStructure.detectingSymmetry(cell: s.cell, sites: s.sites).operations.count)
    }
}

// ---- CIF import (pymatgen: CifParser.from_str(...).parse_structures) ----------------------------------------------------
for (size, f) in [("small", "c2m.cif"), ("typical", "al_super2.cif"), ("large", "al_super3.cif")] {
    let t = text(f)
    measure("cif_import", size: size, n: t.count, note: f) {
        Double(try! CIFImport.read(t, fileBaseName: "x").structure.sites.count)
    }
}

// ---- Orientation symmetry reduction of N rotations, m-3m (orix: Orientation.map_into_symmetry_reduced_zone) ---------------------
let rotData = try! Data(contentsOf: inputs.appendingPathComponent("rotations.f64"))
let rot: [Double] = rotData.withUnsafeBytes { Array($0.bindMemory(to: Double.self)) }
let sym = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
for (size, n) in [("small", 1000), ("typical", 10000), ("large", 100000)] {
    measure("orientation_reduce", size: size, n: n, note: "m-3m, 24 proper operators") {
        var s = 0.0
        for i in 0..<n {
            let r = Array(rot[(i * 9)..<(i * 9 + 9)])
            // rows of the row-major file are the matrix rows; simd columns are the transposed rows
            let m = simd_double3x3(columns: (SIMD3(r[0], r[3], r[6]), SIMD3(r[1], r[4], r[7]), SIMD3(r[2], r[5], r[8])))
            s += Double(try! sym.reduce(m).disorientationRad)
        }
        return s
    }
}

// ---- Ellipse fit of one clean ring (py4DSTEM: fit_ellipse_1D) -------------------------------------------------------------------
for (size, n) in [("small", 128), ("typical", 256), ("large", 512)] {
    let d = try! Data(contentsOf: inputs.appendingPathComponent("ring_\(n).f32"))
    let px: [Float] = d.withUnsafeBytes { Array($0.bindMemory(to: Float.self)) }
    let pattern = DiffractionPattern(width: n, height: n, pixels: px)
    let c = Double(n - 1) / 2 + 0.3
    measure("ellipse_fit_1d", size: size, n: n, note: "\(n) x \(n) pattern") {
        (try? EllipseCalibration.fit1D(pattern: pattern, centerQX: c, centerQY: c, innerRadius: 0.15 * Double(n),
                                       outerRadius: 0.45 * Double(n)).a) ?? -1
    }
}

// ---- ACOM: template library, and one pattern / a scan (py4DSTEM: orientation_plan, match_single_pattern) --------------------
do {
    let crystal = CrystalPresets.aluminum
    let refl = try! crystal.kinematicReflections(kMax: 1.6, retainingAll: false).kept.map {
        ACOMReflection(g: $0.reflection.g, intensity: $0.structureFactor.intensity)
    }
    for (size, count) in [("small", 100), ("typical", 600), ("large", 2400)] {
        measure("acom_library", size: size, n: count, note: "\(refl.count) reflections, 32 x 128 polar bins") {
            Double(try! ACOMTemplateLibrary.generate(reflections: refl, kMax: 1.6, laueClass: .mThreeBarM,
                settings: ACOMTemplateSettings(zoneAxisCount: count)).count)
        }
    }
    let lib = try! ACOMTemplateLibrary.generate(reflections: refl, kMax: 1.6, laueClass: .mThreeBarM,
        settings: ACOMTemplateSettings(zoneAxisCount: 600, excitation: .appCurvature(wavelengthAngstrom: 0.0251)))
    let matcher = ACOMMatcher(library: lib)!
    let ipp = 0.004
    func scan(_ n: Int) -> [[ACOMPeak]] {
        var seed: UInt64 = 99
        return (0..<n).map { i in
            seed = seed &* 6364136223846793005 &+ 1442695040888963407
            let phi = Double(seed >> 11) / Double(1 << 53) * 2 * .pi
            return lib.templateSpots[(i * 7) % lib.count].map { s in
                let az = Double(s.azimuth) + phi, r = Double(s.radius) / ipp
                return ACOMPeak(x: Float(128 + r * cos(az)), y: Float(128 + r * sin(az)), intensity: Float(pow(s.weight, 4)))
            }
        }
    }
    for (size, w) in [("small", 8), ("typical", 32), ("large", 128)] {
        let peaks = scan(w * w)
        measure("acom_match_scan", size: size, n: w * w, backend: "cpu 1 worker", note: "600 templates, \(w) x \(w) scan") {
            Double(matcher.matchAll(peaks: peaks, width: w, height: w, originX: 128, originY: 128, invAngstromPerPixel: ipp,
                                    maxWorkers: 1, backend: .cpu)?.results.count ?? 0)
        }
        measure("acom_match_scan", size: size, n: w * w, backend: "cpu all cores", note: "600 templates, \(w) x \(w) scan") {
            Double(matcher.matchAll(peaks: peaks, width: w, height: w, originX: 128, originY: 128, invAngstromPerPixel: ipp,
                                    maxWorkers: nil, backend: .cpu)?.results.count ?? 0)
        }
        if matcher.isMetalAvailable {
            measure("acom_match_scan", size: size, n: w * w, backend: "metal", note: "600 templates, \(w) x \(w) scan") {
                Double(matcher.matchAll(peaks: peaks, width: w, height: w, originX: 128, originY: 128, invAngstromPerPixel: ipp,
                                        maxWorkers: nil, backend: .metal)?.results.count ?? 0)
            }
        }
    }
}
