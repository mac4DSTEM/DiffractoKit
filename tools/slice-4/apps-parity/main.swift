// DiffractoKit slice 4, oracle (registration docs/registration-slice-4.md sections 5, 6 and gate G1): a scratch tool that links the
// apps' own PhaseReferenceLibrary.swift and PhaseVectorMatching.swift (mac4DSTEM 7d6fdbd1, fetched by run.sh with `git show`,
// compiled beside this file) with Crystal.swift, ScatteringFactors.swift, AnalysisCancellationToken.swift, the BraggPeak /
// BraggVectors heads cut by line range from DiskDetection.swift, and EulerAngles / ACOMOrientation cut by line range from
// OrientationResult.swift. It builds the libraries and the planted scans below, runs `PhaseVectorMatcher.map` and the other
// entry points, and dumps every number as Double / Float bit patterns (b / hex strings) beside the shortest round-trip decimal.
//
//   apps-parity <out.json>
//
// Written from the registration and the apps' tests (mac4DSTEMTests/PhaseVectorMatchingTests.swift at 7d6fdbd1 for the θ′ cell and
// the Al / β″ phase definitions), not from the kit. No kit code is read or linked. The planted scans come from a SplitMix64
// generator written here; the dump carries every input vector (pixel Float bit patterns, origin, scale) so a replay needs no
// generator.
//
// Dump layout (all hex strings are lower-case bit patterns, 16 digits for Double, 8 for Float):
//   D(x) = {"v": decimal|null, "b": bits}   F(x) the same for Float.
//   reflections[]   : the app's Crystal.reflections(kMax: 1.6) in its order, per phase
//   libraries[]     : every PhaseOrientationReference (phaseIndex, zoneAxis, rotation, every vector as [h,k,l,qx,qy,len,rel] bits)
//   scans[]         : the planted input, the matcher's experimental vectors, matrix fit, planted truth
//   runs[]          : one PhaseVectorMatcher.map run (scan x library x settings x matrix argument), every PhaseVectorResult
//   fieldEffects[]  : per PhaseVectorSettings field, how many positions each scan moved
//   zoneAxisFits[], claims[], resolution[], inPlaneSteps[], nearest[], euler[], detectorBases[]

import Foundation
import simd

// MARK: JSON helpers

func hex(_ x: Double) -> String {
    let s = String(x.bitPattern, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}
func hex(_ x: UInt64) -> String {
    let s = String(x, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}
func hex(_ x: Float) -> String {
    let s = String(x.bitPattern, radix: 16)
    return String(repeating: "0", count: 8 - s.count) + s
}
func dec(_ x: Double) -> String { x.isFinite ? "\(x)" : "null" }
func dec(_ x: Float) -> String { x.isFinite ? "\(x)" : "null" }
func D(_ x: Double) -> String { "{\"v\":\(dec(x)),\"b\":\"\(hex(x))\"}" }
func F(_ x: Float) -> String { "{\"v\":\(dec(x)),\"b\":\"\(hex(x))\"}" }
func arr(_ items: [String]) -> String { "[" + items.joined(separator: ",") + "]" }
func obj(_ pairs: [(String, String)]) -> String {
    "{" + pairs.map { "\"\($0.0)\":\($0.1)" }.joined(separator: ",") + "}"
}
func str(_ s: String) -> String {
    var out = "\""
    for ch in s.unicodeScalars {
        switch ch {
        case "\"": out += "\\\""
        case "\\": out += "\\\\"
        case "\n": out += "\\n"
        default: out += String(ch)
        }
    }
    return out + "\""
}
func i3(_ v: SIMD3<Int>) -> String { "[\(v.x),\(v.y),\(v.z)]" }
func mat(_ m: simd_double3x3) -> String {
    let c = [m.columns.0, m.columns.1, m.columns.2]
    let flat = c.flatMap { [$0.x, $0.y, $0.z] }     // column-major
    return obj([("order", str("column-major")), ("b", arr(flat.map { "\"\(hex($0))\"" })), ("v", arr(flat.map { dec($0) }))])
}

// The inputs this tool generates must not depend on the optimisation level. At -O the compiler fuses a cos(x) and a sin(x) of the
// same argument into one __sincos_stret call, whose last bit can differ from the separate calls -Onone makes (measured: the
// app's own PhaseReferenceLibrary.rotate differs that way at one in-plane step, see the README; a lone sin(x) can also be lowered
// to it). These two wrappers are compiled without optimisation and never inlined, so the inputs are the same in both builds and
// a difference between the dumps is the apps' code alone.
@_optimize(none) @inline(never) func cosI(_ x: Double) -> Double { Foundation.cos(x) }
@_optimize(none) @inline(never) func sinI(_ x: Double) -> Double { Foundation.sin(x) }

// MARK: Seeded generator (SplitMix64; written here, the dump carries every value it produced)

struct Rng {
    var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9E37_79B9_7F4A_7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
        z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
        return z ^ (z >> 31)
    }
    mutating func uniform() -> Double { Double(next() >> 11) * (1.0 / 9_007_199_254_740_992.0) }
    mutating func range(_ a: Double, _ b: Double) -> Double { a + (b - a) * uniform() }
    mutating func int(_ n: Int) -> Int { Int(next() % UInt64(n)) }
    mutating func gauss() -> Double {
        let u1 = max(uniform(), 1e-300), u2 = uniform()
        return (-2 * log(u1)).squareRoot() * cosI(2 * Double.pi * u2)
    }
}

// MARK: Crystals

func crystalJSON(_ c: Crystal) -> String {
    obj([("a", D(c.a)), ("b", D(c.b)), ("c", D(c.c)),
         ("alphaDeg", D(c.alphaDeg)), ("betaDeg", D(c.betaDeg)), ("gammaDeg", D(c.gammaDeg)),
         ("volume", D(c.volume)),
         ("sites", arr(c.sites.map { "[\($0.z),\(D($0.fractional.x)),\(D($0.fractional.y)),\(D($0.fractional.z)),\(D($0.occupancy))]" })),
         ("latReal", arr(c.latReal.map { arr([D($0.x), D($0.y), D($0.z)]) })),
         ("latInv", arr(c.latInv.map { arr([D($0.x), D($0.y), D($0.z)]) }))])
}

// θ′ (Al2Cu), the six-site cell of PhaseVectorMatchingTests.testOrientationRelationship... at 7d6fdbd1 (lines 874-881)
let thetaPrime = Crystal(a: 4.04, b: 4.04, c: 5.80, sites: [
    AtomSite(z: 13, fractional: [0.000000, 0.000000, 0.000000]),
    AtomSite(z: 13, fractional: [0.500000, 0.500000, 0.500000]),
    AtomSite(z: 13, fractional: [0.000000, 0.000000, 0.500000]),
    AtomSite(z: 13, fractional: [0.500000, 0.500000, 0.000000]),
    AtomSite(z: 29, fractional: [0.000000, 0.500000, 0.250000]),
    AtomSite(z: 29, fractional: [0.500000, 0.000000, 0.750000]),
])
// an oblique cell with two species, for the hexagonal-free general case (not in the apps' tests)
let triclinic = Crystal(a: 5.1, b: 6.3, c: 7.9, alphaDeg: 82, betaDeg: 96, gammaDeg: 104, sites: [
    AtomSite(z: 13, fractional: [0, 0, 0]),
    AtomSite(z: 29, fractional: [0.31, 0.42, 0.17]),
])

let relOR = OrientationRelationship(candidate: .plane(SIMD3(0, 0, 2)), matrix: .plane(SIMD3(2, 0, 0)))
let z001 = SIMD3<Int>(0, 0, 1), z011 = SIMD3<Int>(0, 1, 1), z111 = SIMD3<Int>(1, 1, 1)
let z100 = SIMD3<Int>(1, 0, 0), z010 = SIMD3<Int>(0, 1, 0)

func mainDefinitions(withOR: Bool) -> [PhaseDefinition] {
    let common = [z001, z011, z111]
    return [
        PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: common),
        PhaseDefinition(id: "au", displayName: "Au", crystal: .gold, role: .candidate, zoneAxes: common),
        PhaseDefinition(id: "theta-prime", displayName: "θ′-Al2Cu", crystal: thetaPrime, role: .candidate,
                        zoneAxes: common + [z100],
                        orientationRelationships: withOR ? [relOR] : []),
        PhaseDefinition(id: "beta-double-prime", displayName: "β″-Mg5Si6", crystal: .betaDoublePrime, role: .candidate,
                        zoneAxes: common + [z010]),
    ]
}

func alAllDefinitions() -> [PhaseDefinition] {
    [PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: []),
     PhaseDefinition(id: "au", displayName: "Au", crystal: .gold, role: .candidate, zoneAxes: [z001])]
}

func obliqueDefinitions() -> [PhaseDefinition] {
    [PhaseDefinition(id: "al", displayName: "Al", crystal: .aluminum, role: .matrix, zoneAxes: [z001]),
     PhaseDefinition(id: "mg", displayName: "Mg", crystal: .magnesium, role: .candidate, zoneAxes: [z001, z100]),
     PhaseDefinition(id: "tri", displayName: "triclinic", crystal: triclinic, role: .candidate,
                     zoneAxes: [z001, z011, SIMD3(1, 1, 2)])]
}

func latticeVectorJSON(_ v: LatticeVector) -> String {
    switch v {
    case .plane(let h): return "[\"plane\",\(i3(h))]"
    case .direction(let u): return "[\"direction\",\(i3(u))]"
    }
}

func phaseJSON(_ p: PhaseDefinition) -> String {
    obj([("id", str(p.id)), ("displayName", str(p.displayName)),
         ("role", str(p.role.rawValue)),
         ("zoneAxes", arr(p.zoneAxes.map(i3))),
         ("orientationRelationships", arr(p.orientationRelationships.map {
             "{\"candidate\":\(latticeVectorJSON($0.candidate)),\"matrix\":\(latticeVectorJSON($0.matrix))}" })),
         ("excitationSlabInvAngstrom", p.excitationSlabInvAngstrom.map(D) ?? "null"),
         ("crystal", crystalJSON(p.crystal))])
}

func referenceSettingsJSON(_ s: PhaseReferenceSettings) -> String {
    obj([("kMaxInvAngstrom", D(s.kMaxInvAngstrom)), ("excitationSlabInvAngstrom", D(s.excitationSlabInvAngstrom)),
         ("excitationWidthInvAngstrom", D(s.excitationWidthInvAngstrom)),
         ("minimumIntensityFraction", D(s.minimumIntensityFraction)),
         ("maximumVectorsPerEntry", "\(s.maximumVectorsPerEntry)"), ("inPlaneStepDeg", D(s.inPlaneStepDeg)),
         ("maximumEntries", "\(s.maximumEntries)")])
}

func settingsJSON(_ s: PhaseVectorSettings) -> String {
    obj([("directBeamRadiusInvAngstrom", D(s.directBeamRadiusInvAngstrom)),
         ("maximumVectorInvAngstrom", D(s.maximumVectorInvAngstrom)),
         ("matrixToleranceInvAngstrom", D(s.matrixToleranceInvAngstrom)),
         ("pairRadiusInvAngstrom", D(s.pairRadiusInvAngstrom)),
         ("minimumVectors", "\(s.minimumVectors)"),
         ("matrixFallbackExplainedFraction", D(s.matrixFallbackExplainedFraction)),
         ("minimumMatchedVectors", "\(s.minimumMatchedVectors)"),
         ("friedelPairMinimumMatchedVectors", "\(s.friedelPairMinimumMatchedVectors)"),
         ("chanceMatchMultiple", D(s.chanceMatchMultiple)),
         ("notIndexedAboveInvAngstrom", D(s.notIndexedAboveInvAngstrom)),
         ("minimumPhaseContrastInvAngstrom", D(s.minimumPhaseContrastInvAngstrom)),
         ("orientationRelationshipToleranceDeg", D(s.orientationRelationshipToleranceDeg)),
         ("completenessAwareCrossPhaseRanking", s.completenessAwareCrossPhaseRanking ? "true" : "false"),
         ("classificationRule", str(s.classificationRule.rawValue)),
         ("directMatrixMaximumVectors", "\(s.directMatrixMaximumVectors)"),
         ("residualCutoffInvAngstrom", D(s.residualCutoffInvAngstrom)),
         ("knownVariantsMinimumSpecificReflections", "\(s.knownVariantsMinimumSpecificReflections)")])
}

// MARK: Libraries

/// FNV-1a 64 over the little-endian bytes of these 64-bit words: phaseIndex, u, v, w, the rotation's Double bits, the vector count,
/// then per vector h, k, l (as Int64 bit patterns), q.x, q.y, length and relativeIntensity (Double bits). Every entry of every library
/// carries it ("hash"), so an entry whose vectors are not dumped in full is still checked bit for bit.
func entryHash(_ e: PhaseOrientationReference) -> UInt64 {
    var h: UInt64 = 0xcbf2_9ce4_8422_2325
    func feed(_ w: UInt64) {
        var x = w
        for _ in 0..<8 { h = (h ^ (x & 0xff)) &* 0x0000_0100_0000_01b3; x >>= 8 }
    }
    func feed(_ i: Int) { feed(UInt64(bitPattern: Int64(i))) }
    feed(e.phaseIndex); feed(e.zoneAxis.x); feed(e.zoneAxis.y); feed(e.zoneAxis.z)
    feed(e.inPlaneRotationRad.bitPattern); feed(e.vectors.count)
    for v in e.vectors {
        feed(v.h); feed(v.k); feed(v.l)
        feed(v.q.x.bitPattern); feed(v.q.y.bitPattern); feed(v.length.bitPattern); feed(v.relativeIntensity.bitPattern)
    }
    return h
}

/// Rotation-0 entries also carry decimals ("vd"); the rotated ones are hex only (decimals are derived by `rotate`).
/// `full` false drops the vector lists and keeps the hash (and the count) only.
func entryJSON(_ e: PhaseOrientationReference, chance: Double, rotationIndex: Int, full: Bool) -> String {
    var pairs: [(String, String)] = [
        ("p", "\(e.phaseIndex)"), ("u", i3(e.zoneAxis)), ("k", "\(rotationIndex)"), ("r", "\"\(hex(e.inPlaneRotationRad))\""),
        ("n", "\(e.vectors.count)"), ("hash", "\"\(hex(entryHash(e)))\""), ("chance", "\"\(hex(chance))\""),
    ]
    if full {
        pairs.append(("v", arr(e.vectors.map {
            "[\($0.h),\($0.k),\($0.l),\"\(hex($0.q.x))\",\"\(hex($0.q.y))\",\"\(hex($0.length))\",\"\(hex($0.relativeIntensity))\"]"
        })))
        if e.inPlaneRotationRad == 0 {
            pairs.append(("vd", arr(e.vectors.map { "[\(dec($0.q.x)),\(dec($0.q.y)),\(dec($0.length)),\(dec($0.relativeIntensity))]" })))
        }
    }
    return obj(pairs)
}

/// `fullRotations`: nil dumps every entry's vectors; otherwise only the entries whose in-plane step index k is in the set (all
/// entries keep their hash). k counts the entries of one (phase, zone axis) from 0 in library order.
func libraryJSON(_ name: String, _ lib: PhaseReferenceLibrary, fullRotations: Set<Int>? = nil) -> String {
    var counter: [String: Int] = [:]
    let rows: [String] = lib.entries.map { e in
        let key = "\(e.phaseIndex)/\(e.zoneAxis.x),\(e.zoneAxis.y),\(e.zoneAxis.z)"
        let k = counter[key, default: 0]
        counter[key] = k + 1
        return entryJSON(e, chance: e.chanceMatchFraction(pairRadius: 0.02, accessibleRadius: 0.9),
                         rotationIndex: k, full: fullRotations.map { $0.contains(k) } ?? true)
    }
    return obj([
        ("name", str(name)), ("referenceSettings", referenceSettingsJSON(lib.settings)),
        ("phases", arr(lib.phases.map(phaseJSON))),
        ("matrixPhaseIndex", "\(lib.matrixPhaseIndex)"), ("entryCount", "\(lib.entries.count)"),
        ("matrixEntryIndices", arr(lib.matrixEntryIndices.map { "\($0)" })),
        ("candidateEntryCount", "\(lib.candidateEntryIndices.count)"),
        ("chanceBasis", obj([("pairRadius", D(0.02)), ("accessibleRadius", D(0.9))])),
        ("fullRotationIndices", fullRotations.map { arr($0.sorted().map { "\($0)" }) } ?? "\"all\""),
        ("entries", "[\n" + rows.joined(separator: ",\n") + "\n]"),
    ])
}

func firstEntry(_ lib: PhaseReferenceLibrary, phase: Int, axis: SIMD3<Int>) -> Int {
    lib.entries.firstIndex { $0.phaseIndex == phase && $0.zoneAxis == axis }!
}

func strongVectors(_ e: PhaseOrientationReference, limit: Int, minLength: Double, maxLength: Double) -> [SIMD2<Double>] {
    let kept = e.vectors.enumerated().filter { $0.element.length >= minLength && $0.element.length <= maxLength }
    // strongest first, ties on the index (total order, so the plant never depends on a sort's stability)
    let ordered = kept.sorted { $0.element.relativeIntensity != $1.element.relativeIntensity
        ? $0.element.relativeIntensity > $1.element.relativeIntensity : $0.offset < $1.offset }
    return ordered.prefix(limit).map { $0.element.q }
}

// MARK: Scans

struct Scan {
    let name: String
    let width: Int, height: Int
    let originX: Float, originY: Float
    let scale: Double
    var peaks: [[BraggPeak]]
    var truth: [String]
    var bragg: BraggVectors { BraggVectors(scanWidth: width, scanHeight: height, peaks: peaks) }
}

func peak(_ q: SIMD2<Double>, originX: Float, originY: Float, scale: Double, jitter: Double, rng: inout Rng) -> BraggPeak {
    let jx = jitter > 0 ? rng.range(-jitter, jitter) : 0
    let jy = jitter > 0 ? rng.range(-jitter, jitter) : 0
    return BraggPeak(x: Float(Double(originX) + (q.x + jx) / scale),
                     y: Float(Double(originY) + (q.y + jy) / scale), intensity: 1)
}

func scanJSON(_ s: Scan) -> String {
    let defaults = PhaseVectorSettings()
    func vecs(_ radius: Double, _ maxV: Double) -> String {
        arr(s.peaks.map { p in
            arr(PhaseVectorMatcher.experimentalVectors(
                peaks: p, originX: s.originX, originY: s.originY, invAngstromPerPixel: s.scale,
                directBeamRadiusInvAngstrom: radius, maximumVectorInvAngstrom: maxV)
                .map { "[\"\(hex($0.x))\",\"\(hex($0.y))\"]" })
        })
    }
    return obj([("name", str(s.name)), ("scanWidth", "\(s.width)"), ("scanHeight", "\(s.height)"),
                ("originX", F(s.originX)), ("originY", F(s.originY)), ("invAngstromPerPixel", D(s.scale)),
                ("truth", arr(s.truth.map(str))),
                ("peaks", "[\n" + s.peaks.map { p in
                    arr(p.map { "[\(dec($0.x)),\(dec($0.y)),\"\(hex($0.x))\",\"\(hex($0.y))\"]" }) }.joined(separator: ",\n") + "\n]"),
                ("vectorsDefault", vecs(defaults.directBeamRadiusInvAngstrom, defaults.maximumVectorInvAngstrom)),
                ("vectorsRadius0Max0p9", vecs(0, 0.9))])
}

// MARK: Runs

struct RunOut {
    let json: String
    let signature: [String]
    let matrixEntryIndex: Int
    let returnedNil: Bool
    let verdictCounts: [Int]
    let map: PhaseMap?
}

func resultJSON(_ r: PhaseVectorResult) -> String {
    "[\(r.verdict.rawValue),\(r.phaseIndex),\(r.entryIndex),\"\(hex(r.score))\",\(dec(r.score)),\(r.matchedCount),\(r.survivingCount),\(r.removedCount),\(r.runnerUpPhaseIndex),\"\(hex(r.runnerUpScore))\",\(dec(r.runnerUpScore))]"
}
func resultSignature(_ r: PhaseVectorResult) -> String {
    "\(r.verdict.rawValue),\(r.phaseIndex),\(r.entryIndex),\(hex(r.score)),\(r.matchedCount),\(r.survivingCount),\(r.removedCount),\(r.runnerUpPhaseIndex),\(hex(r.runnerUpScore))"
}

func execute(id: String, label: String, scan: Scan, libraryName: String, library: PhaseReferenceLibrary,
             settings: PhaseVectorSettings, matrixArg: Int?, field: String? = nil) -> RunOut {
    let map = PhaseVectorMatcher.map(
        bragg: scan.bragg, library: library, settings: settings,
        originX: scan.originX, originY: scan.originY, invAngstromPerPixel: scan.scale, matrixEntryIndex: matrixArg)
    var counts = [0, 0, 0, 0]
    var sig: [String] = []
    var rows = "null"
    if let map {
        sig = map.results.map(resultSignature)
        for r in map.results { counts[Int(r.verdict.rawValue)] += 1 }
        rows = "[\n" + map.results.map(resultJSON).joined(separator: ",\n") + "\n]"
    }
    let json = obj([
        ("id", str(id)), ("label", str(label)), ("scan", str(scan.name)), ("library", str(libraryName)),
        ("field", field.map(str) ?? "null"), ("settings", settingsJSON(settings)),
        ("matrixArgument", matrixArg.map { "\($0)" } ?? "null"),
        ("returnedNil", map == nil ? "true" : "false"),
        ("matrixEntryIndex", map.map { "\($0.matrixEntryIndex)" } ?? "null"),
        ("matrixPhaseIndex", map.map { "\($0.matrixPhaseIndex)" } ?? "null"),
        ("phaseNames", map.map { arr($0.phaseNames.map(str)) } ?? "null"),
        ("verdictCounts", arr(counts.map { "\($0)" })),
        ("phaseCounts", map.map { arr($0.phaseCounts.map { "\($0)" }) } ?? "null"),
        ("results", rows),
    ])
    return RunOut(json: json, signature: sig, matrixEntryIndex: map?.matrixEntryIndex ?? -2,
                  returnedNil: map == nil, verdictCounts: counts, map: map)
}

// MARK: Euler

func rotZ(_ a: Double) -> simd_double3x3 {
    let c = cosI(a), s = sinI(a)
    return simd_double3x3(columns: (SIMD3(c, s, 0), SIMD3(-s, c, 0), SIMD3(0, 0, 1)))
}
func rotX(_ a: Double) -> simd_double3x3 {
    let c = cosI(a), s = sinI(a)
    return simd_double3x3(columns: (SIMD3(1, 0, 0), SIMD3(0, c, s), SIMD3(0, -s, c)))
}
func quaternionMatrix(_ w: Double, _ x: Double, _ y: Double, _ z: Double) -> simd_double3x3 {
    let n = (w * w + x * x + y * y + z * z).squareRoot()
    let (a, b, c, d) = (w / n, x / n, y / n, z / n)
    return simd_double3x3(rows: [
        SIMD3(1 - 2 * (c * c + d * d), 2 * (b * c - a * d), 2 * (b * d + a * c)),
        SIMD3(2 * (b * c + a * d), 1 - 2 * (b * b + d * d), 2 * (c * d - a * b)),
        SIMD3(2 * (b * d - a * c), 2 * (c * d + a * b), 1 - 2 * (b * b + c * c)),
    ])
}

func eulerCase(_ label: String, _ m: simd_double3x3) -> String {
    let e = EulerAngles(py4DSTEMOrientationMatrix: m)
    let again = EulerAngles(py4DSTEMOrientationMatrix: e.py4DSTEMOrientationMatrix)
    let deg = e.degrees
    return obj([("label", str(label)), ("matrix", mat(m)),
                ("phi1", F(e.phi1)), ("Phi", F(e.Phi)), ("phi2", F(e.phi2)),
                ("crystalToLab", mat(e.crystalToLabMatrix)), ("py4DSTEMMatrix", mat(e.py4DSTEMOrientationMatrix)),
                ("degrees", arr([F(deg.0), F(deg.1), F(deg.2)])),
                ("roundTrip", arr([F(again.phi1), F(again.Phi), F(again.phi2)]))])
}

// MARK: Main

func runAll(outPath: String) throws {
    let started = Date()
    let defaults = PhaseVectorSettings()
    var sections: [(String, String)] = []
    sections.append(("tool", obj([("source", str("mac4DSTEM 7d6fdbd1")),
                                  ("note", str("apps-parity-4: the apps' PhaseReferenceLibrary / PhaseVectorMatching / Crystal compiled by swiftc; see tools/slice-4/apps-parity/README.md")),
                                  ("defaultPhaseVectorSettings", settingsJSON(defaults)),
                                  ("defaultPhaseReferenceSettings", referenceSettingsJSON(PhaseReferenceSettings()))])))

    // ---- 0. the reflection lists the library reads (the app's order, kMax 1.6, default tolerance 1e-4)
    let crystalList: [(String, Crystal)] = [("Al", .aluminum), ("Au", .gold), ("theta-prime", thetaPrime),
                                            ("beta-double-prime", .betaDoublePrime), ("Mg", .magnesium), ("triclinic", triclinic)]
    var reflectionBlocks: [String] = []
    var tieNotes: [String] = []
    for (name, c) in crystalList {
        let refl = c.reflections(kMax: PhaseReferenceSettings().kMaxInvAngstrom)
        var exactTies = 0, tiesInLoopOrder = 0
        for i in 1..<max(1, refl.count) where refl[i].gLength.bitPattern == refl[i - 1].gLength.bitPattern {
            exactTies += 1
            let a = refl[i - 1], b = refl[i]
            if (a.h, a.k, a.l) < (b.h, b.k, b.l) { tiesInLoopOrder += 1 }
        }
        tieNotes.append("\(name): \(refl.count) reflections, \(exactTies) adjacent exact-|g| ties, \(tiesInLoopOrder) of them in (h,k,l) ascending order")
        reflectionBlocks.append(obj([
            ("phase", str(name)), ("kMax", D(1.6)), ("tolerance", D(1e-4)), ("count", "\(refl.count)"),
            ("crystal", crystalJSON(c)),
            ("hkl", arr(refl.map { "[\($0.h),\($0.k),\($0.l)]" })),
            ("gLength", arr(refl.map { "\"\(hex($0.gLength))\"" })),
            ("g", arr(refl.map { "[\"\(hex($0.g.x))\",\"\(hex($0.g.y))\",\"\(hex($0.g.z))\"]" })),
            ("re", arr(refl.map { "\"\(hex($0.structureFactorRe))\"" })),
            ("im", arr(refl.map { "\"\(hex($0.structureFactorIm))\"" })),
            ("intensity", arr(refl.map { "\"\(hex($0.intensity))\"" })),
            ("exactTies", "\(exactTies)"), ("tiesInHklAscendingOrder", "\(tiesInLoopOrder)"),
        ]))
    }
    sections.append(("reflections", "[\n" + reflectionBlocks.joined(separator: ",\n") + "\n]"))
    for t in tieNotes { print("apps-parity-4: reflections  \(t)") }

    // ---- 1. libraries
    let main = try PhaseReferenceLibrary.build(phases: mainDefinitions(withOR: false))
    let mainOR = try PhaseReferenceLibrary.build(phases: mainDefinitions(withOR: true))
    var sameEntries = main.entries.count == mainOR.entries.count
    if sameEntries {
        for (a, b) in zip(main.entries, mainOR.entries) {
            if a.phaseIndex != b.phaseIndex || a.zoneAxis != b.zoneAxis || a.inPlaneRotationRad.bitPattern != b.inPlaneRotationRad.bitPattern
                || a.vectors.count != b.vectors.count { sameEntries = false; break }
            for (u, v) in zip(a.vectors, b.vectors)
            where u.q.x.bitPattern != v.q.x.bitPattern || u.q.y.bitPattern != v.q.y.bitPattern
                || u.relativeIntensity.bitPattern != v.relativeIntensity.bitPattern || u.h != v.h || u.k != v.k || u.l != v.l {
                sameEntries = false
            }
        }
    }
    print("apps-parity-4: main library \(main.entries.count) entries; with the orientation relationship the entries are \(sameEntries ? "identical" : "DIFFERENT")")

    var alAllSettings = PhaseReferenceSettings()
    alAllSettings.maximumEntries = 20000
    let alAll = try PhaseReferenceLibrary.build(phases: alAllDefinitions(), settings: alAllSettings)
    let oblique = try PhaseReferenceLibrary.build(phases: obliqueDefinitions())
    print("apps-parity-4: al-all-axes library \(alAll.entries.count) entries, oblique \(oblique.entries.count)")

    var libraryBlocks: [String] = []
    libraryBlocks.append(libraryJSON("main", main))
    libraryBlocks.append(obj([("name", str("main-or")), ("sameEntriesAsMain", sameEntries ? "true" : "false"),
                              ("phases", arr(mainOR.phases.map(phaseJSON)))]))
    libraryBlocks.append(libraryJSON("al-all-axes", alAll, fullRotations: [0, 1, 17, 90, 179]))
    libraryBlocks.append(libraryJSON("oblique", oblique, fullRotations: [0, 1, 17, 90, 179]))
    // the 62 low-index axes, in the app's order
    let axes62 = PhaseReferenceLibrary.lowIndexZoneAxes
    sections.append(("lowIndexZoneAxes", arr(axes62.map(i3))))
    sections.append(("libraries", "[\n" + libraryBlocks.joined(separator: ",\n") + "\n]"))
    // in-plane steps for a set of step values
    let stepValues: [Double] = [2, 45, 7, 360, 359.9, 0.7, 0, -3, 361, .nan, .infinity, 90, 1, 0.5, 100]
    sections.append(("inPlaneSteps", arr(stepValues.map { step in
        var s = PhaseReferenceSettings(); s.inPlaneStepDeg = step
        let r = PhaseReferenceLibrary.inPlaneSteps(s)
        return obj([("step", D(step)), ("count", "\(r.count)"), ("radians", arr(r.map { "\"\(hex($0))\"" }))])
    })))

    // ---- 2. the planted scans
    let matrixIdx = firstEntry(main, phase: 0, axis: z001) + 3                   // Al [001] at 6 degrees
    let alAxis111Idx = firstEntry(main, phase: 0, axis: z111)
    let thetaFirst100 = firstEntry(main, phase: 2, axis: z100)
    // the θ′ [100] entry in the middle of the orientation-relationship window around the matrix entry
    let consistent = (0..<180).map { thetaFirst100 + $0 }.filter {
        PhaseVectorMatcher.orientationConsistent(
            candidate: main.entries[$0], candidateCrystal: thetaPrime, matrix: main.entries[matrixIdx], matrixCrystal: .aluminum,
            relationships: [relOR], toleranceDeg: 10)
    }
    precondition(!consistent.isEmpty, "no theta-prime [100] entry satisfies the orientation relationship")
    let thetaOR = consistent[consistent.count / 2]
    let thetaFirst011 = firstEntry(main, phase: 2, axis: z011)
    let betaFirst010 = firstEntry(main, phase: 3, axis: z010)
    let betaFirst001 = firstEntry(main, phase: 3, axis: z001)
    let auFirst011 = firstEntry(main, phase: 1, axis: z011)
    let auFirst001 = firstEntry(main, phase: 1, axis: z001)
    let alFirst001 = firstEntry(main, phase: 0, axis: z001)
    let thetaFirst001 = firstEntry(main, phase: 2, axis: z001)

    func describe(_ idx: Int) -> String {
        let e = main.entries[idx]
        return "\(main.phases[e.phaseIndex].id) \(e.zoneAxis.x)\(e.zoneAxis.y)\(e.zoneAxis.z) entry \(idx) @ \(e.inPlaneRotationRad * 180 / .pi) deg"
    }

    // (a) 16 x 16 Al [001] matrix scan with θ′ and β″ planted
    let (oxA, oyA, scaleA, jitter) = (Float(512.37), Float(511.81), 0.0031, 0.004)
    var rngA = Rng(seed: 20_261_015)
    let matrixVecs = strongVectors(main.entries[matrixIdx], limit: 8, minLength: 0.2, maxLength: 1.2)
    var scanA = Scan(name: "a-matrix-with-precipitates", width: 16, height: 16, originX: oxA, originY: oyA, scale: scaleA, peaks: [], truth: [])
    for y in 0..<16 {
        for x in 0..<16 {
            var qs = matrixVecs
            var truth = "matrix: \(describe(matrixIdx))"
            var extraCandidate: Int? = nil
            if x >= 8 && y < 8 { extraCandidate = thetaOR }
            else if x >= 8 && y >= 8 { extraCandidate = ((x * 3 + y * 5) % 7 == 0) ? auFirst011 + 12 : betaFirst010 + 25 }
            else if (2...4).contains(x) && (10...12).contains(y) { extraCandidate = thetaFirst011 + 40 }
            else if (x == 6 && y >= 8) { extraCandidate = betaFirst001 + 40 }
            if let c = extraCandidate {
                qs += strongVectors(main.entries[c], limit: 10, minLength: 0.2, maxLength: 1.2)
                truth += " + \(describe(c))"
            }
            var plantedPeaks: [BraggPeak] = []
            if x == 0 && y == 0 { qs = []; truth = "empty position" }
            else if x == 1 && y == 0 { qs = [SIMD2(0.02, -0.03)]; truth = "direct beam only (inside the beam radius)" }
            else {
                if (x * 7 + y * 3) % 13 == 0 {
                    for _ in 0..<2 { let r = rngA.range(0.2, 1.2), a = rngA.range(0, 2 * Double.pi); qs.append(SIMD2(r * cosI(a), r * sinI(a))) }
                    truth += " + 2 random vectors"
                }
            }
            for q in qs { plantedPeaks.append(peak(q, originX: oxA, originY: oyA, scale: scaleA, jitter: jitter, rng: &rngA)) }
            scanA.peaks.append(plantedPeaks); scanA.truth.append(truth)
        }
    }

    // (b) random vectors at every position
    var rngB = Rng(seed: 20_261_016)
    var scanB = Scan(name: "b-random-vectors", width: 16, height: 16, originX: oxA, originY: oyA, scale: scaleA, peaks: [], truth: [])
    for _ in 0..<256 {
        let n = 3 + rngB.int(10)
        var ps: [BraggPeak] = []
        for _ in 0..<n {
            let r = rngB.range(0.2, 1.3), a = rngB.range(0, 2 * Double.pi)
            ps.append(peak(SIMD2(r * cosI(a), r * sinI(a)), originX: oxA, originY: oyA, scale: scaleA, jitter: 0, rng: &rngB))
        }
        scanB.peaks.append(ps); scanB.truth.append("\(n) random vectors")
    }

    // (c) one Friedel pair per position (a random candidate entry's reference vector and its negative, jittered)
    var rngC = Rng(seed: 20_261_017)
    var scanC = Scan(name: "c-friedel-pairs", width: 16, height: 16, originX: oxA, originY: oyA, scale: scaleA, peaks: [], truth: [])
    let candidates = main.candidateEntryIndices
    for i in 0..<256 {
        var cIdx = candidates[rngC.int(candidates.count)]
        var pool = main.entries[cIdx].vectors.filter { $0.length >= 0.25 && $0.length <= 1.0 }
        while pool.isEmpty { cIdx = candidates[rngC.int(candidates.count)]; pool = main.entries[cIdx].vectors.filter { $0.length >= 0.25 && $0.length <= 1.0 } }
        let q = pool[rngC.int(pool.count)].q
        var ps = [peak(q, originX: oxA, originY: oyA, scale: scaleA, jitter: jitter, rng: &rngC),
                  peak(-q, originX: oxA, originY: oyA, scale: scaleA, jitter: jitter, rng: &rngC)]
        var truth = "Friedel pair of \(describe(cIdx))"
        if i % 4 == 1 {
            for m in matrixVecs.prefix(2) { ps.append(peak(m, originX: oxA, originY: oyA, scale: scaleA, jitter: jitter, rng: &rngC)) }
            truth += " + 2 matrix vectors"
        }
        scanC.peaks.append(ps); scanC.truth.append(truth)
    }

    // (d) special positions, integer origin (128, 128), 0.01 per pixel
    let (oxD, oyD, scaleD) = (Float(128), Float(128), 0.01)
    var rngD = Rng(seed: 20_261_018)
    var scanD = Scan(name: "d-special-positions", width: 4, height: 3, originX: oxD, originY: oyD, scale: scaleD, peaks: [], truth: [])
    func exact(_ qs: [SIMD2<Double>]) -> [BraggPeak] { qs.map { peak($0, originX: oxD, originY: oyD, scale: scaleD, jitter: 0, rng: &rngD) } }
    let m0 = firstEntry(main, phase: 0, axis: z001)           // Al [001] at 0 degrees
    let matrixExact = strongVectors(main.entries[m0], limit: 8, minLength: 0.2, maxLength: 1.2)
    let thetaExact = strongVectors(main.entries[thetaFirst100], limit: 10, minLength: 0.2, maxLength: 1.2)
    let betaExact = strongVectors(main.entries[betaFirst010], limit: 10, minLength: 0.2, maxLength: 1.2)
    let friedelQ = main.entries[betaFirst010].vectors.first { $0.length >= 0.3 }!.q
    scanD.peaks = [
        [],
        [BraggPeak(x: 128.3, y: 127.9, intensity: 1)],
        [BraggPeak(x: 143, y: 128, intensity: 1)],
        [BraggPeak(x: 143, y: 128, intensity: 1), BraggPeak(x: 128, y: 143, intensity: 1), BraggPeak(x: 128, y: 113, intensity: 1)],
        [BraggPeak(x: .nan, y: 130, intensity: 1), BraggPeak(x: .infinity, y: 130, intensity: 1), BraggPeak(x: 170, y: 150, intensity: 1), BraggPeak(x: 90, y: 100, intensity: 1)],
        [BraggPeak(x: 200, y: 140, intensity: 1)],
        exact(matrixExact),
        exact(matrixExact + thetaExact),
        exact(betaExact),
        exact([friedelQ, -friedelQ]),
        exact([friedelQ, friedelQ, matrixExact[0], matrixExact[0]]),
        exact([SIMD2(5.0, 0.2), SIMD2(-4.0, 3.0), matrixExact[1], betaExact[2], betaExact[3]]),
    ]
    scanD.truth = ["empty", "direct beam only", "one vector at |q| = 0.15 (the beam-radius boundary)",
                   "three vectors: two at the 0.15 boundary, one at 0.15 on the other side",
                   "NaN and infinite positions among two usable vectors",
                   "one usable vector", "matrix Al [001] 0 deg exactly", "matrix + theta-prime [100] 0 deg exactly",
                   "beta double prime [010] 0 deg exactly", "one exact Friedel pair of beta''",
                   "duplicated beta'' and matrix vectors", "far vectors (|q| 5) among valid ones"]
    // (e) malformed: count != width * height; (f) zero width
    let scanE = Scan(name: "e-malformed-count", width: 3, height: 3, originX: oxD, originY: oyD, scale: scaleD,
                     peaks: Array(scanD.peaks.prefix(8)), truth: ["8 positions for a 3 x 3 scan"])
    let scanF = Scan(name: "f-zero-width", width: 0, height: 4, originX: oxD, originY: oyD, scale: scaleD, peaks: [], truth: ["width 0"])
    let allScans = [scanA, scanB, scanC, scanD, scanE, scanF]

    // matrix fits
    var fitBlocks: [String] = []
    for s in allScans {
        let fit = PhaseVectorMatcher.fitMatrixOrientation(
            bragg: s.bragg, library: main, settings: defaults, originX: s.originX, originY: s.originY, invAngstromPerPixel: s.scale)
        fitBlocks.append(obj([("scan", str(s.name)),
                              ("fit", fit.map { obj([("entryIndex", "\($0.entryIndex)"), ("matched", "\($0.matched)"), ("meanDistance", D($0.meanDistance))]) } ?? "null")]))
    }
    sections.append(("scans", "[\n" + allScans.map(scanJSON).joined(separator: ",\n") + "\n]"))
    sections.append(("matrixFits", arr(fitBlocks)))
    sections.append(("planted", obj([("matrixEntryIndexScanA", "\(matrixIdx)"), ("thetaPrimeORIndex", "\(thetaOR)"),
                                     ("thetaPrimeORWindow", arr(consistent.map { "\($0)" })),
                                     ("jitterInvAngstrom", D(jitter)),
                                     ("entryIndexAl001Rot0", "\(m0)"), ("entryIndexAl111Rot0", "\(alAxis111Idx)"),
                                     ("entryIndexAu001Rot0", "\(auFirst001)"), ("entryIndexTheta001Rot0", "\(thetaFirst001)"),
                                     ("entryIndexAl001FirstOfAxis", "\(alFirst001)")])))

    // ---- runs
    var runs: [RunOut] = []
    var runIndex: [String: Int] = [:]
    func add(_ id: String, _ label: String, _ scan: Scan, _ libName: String, _ lib: PhaseReferenceLibrary,
             _ settings: PhaseVectorSettings, matrix: Int? = nil, field: String? = nil) {
        runIndex[id] = runs.count
        runs.append(execute(id: id, label: label, scan: scan, libraryName: libName, library: lib, settings: settings, matrixArg: matrix, field: field))
    }
    var known = defaults; known.classificationRule = .knownVariants

    for s in [scanA, scanB, scanC] {
        add("\(s.name)/main/default", "default settings", s, "main", main, defaults)
        add("\(s.name)/main/knownVariants", "classificationRule .knownVariants", s, "main", main, known)
    }
    for s in [scanD] {
        add("\(s.name)/main/default/fit", "default, matrix fitted", s, "main", main, defaults)
        add("\(s.name)/main/default/explicit", "default, matrixEntryIndex = Al [001] 0 deg", s, "main", main, defaults, matrix: m0)
        add("\(s.name)/main/knownVariants/fit", ".knownVariants, matrix fitted", s, "main", main, known)
        add("\(s.name)/main/knownVariants/explicit", ".knownVariants, matrixEntryIndex = Al [001] 0 deg", s, "main", main, known, matrix: m0)
    }
    add("\(scanA.name)/main/default/explicit-fitted", "default, explicit matrix entry = the planted 6 deg entry", scanA, "main", main, defaults, matrix: matrixIdx)
    add("\(scanA.name)/main/default/explicit-wrong", "default, explicit matrix entry = Al [111] 0 deg (wrong axis)", scanA, "main", main, defaults, matrix: alAxis111Idx)
    add("\(scanB.name)/main/default/explicit", "default, explicit matrix entry = Al [001] 0 deg", scanB, "main", main, defaults, matrix: m0)
    add("\(scanC.name)/main/default/explicit", "default, explicit matrix entry = Al [001] 0 deg", scanC, "main", main, defaults, matrix: m0)
    add("\(scanE.name)/main/default", "malformed scan", scanE, "main", main, defaults)
    add("\(scanF.name)/main/default", "zero-width scan", scanF, "main", main, defaults)
    add("\(scanA.name)/main-or/default", "default settings, theta-prime carries the orientation relationship", scanA, "main-or", mainOR, defaults)
    add("\(scanB.name)/main-or/default", "default settings, theta-prime carries the orientation relationship", scanB, "main-or", mainOR, defaults)
    add("\(scanC.name)/main-or/default", "default settings, theta-prime carries the orientation relationship", scanC, "main-or", mainOR, defaults)
    add("\(scanA.name)/main-or/tol30", "orientation relationship tolerance 30 deg", scanA, "main-or", mainOR,
        { var s = defaults; s.orientationRelationshipToleranceDeg = 30; return s }(), field: "orientationRelationshipToleranceDeg=30")

    var claimsSettings = defaults
    claimsSettings.directBeamRadiusInvAngstrom = 0.4; claimsSettings.pairRadiusInvAngstrom = 0.008; claimsSettings.matrixToleranceInvAngstrom = 0.003
    add("\(scanA.name)/main/claims-custom", "custom radii, used by the claims block", scanA, "main", main, claimsSettings)

    // one run per field: (field name, value text, apply, base kind)
    struct Variant { let name: String; let value: String; let apply: (inout PhaseVectorSettings) -> Void; let base: String }
    let variants: [Variant] = [
        Variant(name: "directBeamRadiusInvAngstrom", value: "0.4", apply: { $0.directBeamRadiusInvAngstrom = 0.4 }, base: "default"),
        Variant(name: "maximumVectorInvAngstrom", value: "0.9", apply: { $0.maximumVectorInvAngstrom = 0.9 }, base: "default"),
        Variant(name: "matrixToleranceInvAngstrom", value: "0.003", apply: { $0.matrixToleranceInvAngstrom = 0.003 }, base: "default"),
        Variant(name: "pairRadiusInvAngstrom", value: "0.008", apply: { $0.pairRadiusInvAngstrom = 0.008 }, base: "default"),
        Variant(name: "minimumVectors", value: "12", apply: { $0.minimumVectors = 12 }, base: "default"),
        Variant(name: "matrixFallbackExplainedFraction", value: "0.3", apply: { $0.matrixFallbackExplainedFraction = 0.3 }, base: "default"),
        Variant(name: "minimumMatchedVectors", value: "7", apply: { $0.minimumMatchedVectors = 7 }, base: "default"),
        Variant(name: "friedelPairMinimumMatchedVectors", value: "3", apply: { $0.friedelPairMinimumMatchedVectors = 3 }, base: "default"),
        Variant(name: "chanceMatchMultiple", value: "40", apply: { $0.chanceMatchMultiple = 40 }, base: "default"),
        Variant(name: "notIndexedAboveInvAngstrom", value: "0.003", apply: { $0.notIndexedAboveInvAngstrom = 0.003 }, base: "default"),
        Variant(name: "minimumPhaseContrastInvAngstrom", value: "0.02", apply: { $0.minimumPhaseContrastInvAngstrom = 0.02 }, base: "default"),
        Variant(name: "orientationRelationshipToleranceDeg", value: "1", apply: { $0.orientationRelationshipToleranceDeg = 1 }, base: "default-or"),
        Variant(name: "completenessAwareCrossPhaseRanking", value: "true", apply: { $0.completenessAwareCrossPhaseRanking = true }, base: "default"),
        Variant(name: "classificationRule", value: "knownVariants", apply: { $0.classificationRule = .knownVariants }, base: "default"),
        Variant(name: "directMatrixMaximumVectors", value: "4", apply: { $0.directMatrixMaximumVectors = 4 }, base: "knownVariants"),
        Variant(name: "residualCutoffInvAngstrom", value: "0.01", apply: { $0.residualCutoffInvAngstrom = 0.01 }, base: "knownVariants"),
        Variant(name: "knownVariantsMinimumSpecificReflections", value: "6", apply: { $0.knownVariantsMinimumSpecificReflections = 6 }, base: "knownVariants"),
    ]
    precondition(variants.count == 17)
    var effects: [String] = []
    var unchangedNotes: [String] = []
    for v in variants {
        var perScan: [(String, Int, Int, Bool)] = []                   // scan, changed positions, positions, matrix entry moved
        for s in [scanA, scanB, scanC] {
            var baseSettings = defaults
            var baseID = "\(s.name)/main/default"
            var libName = "main", lib = main
            switch v.base {
            case "knownVariants": baseSettings = known; baseID = "\(s.name)/main/knownVariants"
            case "default-or": libName = "main-or"; lib = mainOR; baseID = "\(s.name)/main-or/default"
            default: break
            }
            var settings = baseSettings
            v.apply(&settings)
            let id = "\(s.name)/\(libName)/field/\(v.name)"
            add(id, "\(v.name) = \(v.value) (base: \(v.base))", s, libName, lib, settings, field: v.name)
            let a = runs[runIndex[baseID]!], b = runs[runIndex[id]!]
            var changed = 0
            for (x, y) in zip(a.signature, b.signature) where x != y { changed += 1 }
            perScan.append((s.name, changed, a.signature.count, a.matrixEntryIndex != b.matrixEntryIndex))
        }
        effects.append(obj([("field", str(v.name)), ("value", str(v.value)), ("base", str(v.base)),
                            ("scans", arr(perScan.map { obj([("scan", str($0.0)), ("changedPositions", "\($0.1)"), ("positions", "\($0.2)"), ("matrixEntryMoved", $0.3 ? "true" : "false")]) }))]))
        let line = perScan.map { "\($0.0.prefix(1)):\($0.1)\($0.3 ? "+matrix" : "")" }.joined(separator: " ")
        print("apps-parity-4: field \(v.name)=\(v.value)  changed positions  \(line)")
        if perScan.allSatisfy({ $0.1 == 0 && !$0.3 }) { unchangedNotes.append(v.name) }
    }
    sections.append(("fieldEffects", "[\n" + effects.joined(separator: ",\n") + "\n]"))
    print("apps-parity-4: fields that changed nothing on scans a, b and c: \(unchangedNotes.isEmpty ? "none" : unchangedNotes.joined(separator: ", "))")
    for (id, i) in runIndex.sorted(by: { $0.value < $1.value }).prefix(12) {
        print("apps-parity-4: run \(id): verdictCounts (matrix, indexed, notIndexed, noData) \(runs[i].verdictCounts), matrix entry \(runs[i].matrixEntryIndex)\(runs[i].returnedNil ? " (nil)" : "")")
    }
    sections.append(("runs", "[\n" + runs.map(\.json).joined(separator: ",\n") + "\n]"))

    // ---- 3. fitZoneAxis, claims, resolution
    func fitJSON(_ f: PhaseVectorMatcher.ZoneAxisFit) -> String {
        obj([("zoneAxis", i3(f.zoneAxis)), ("inPlaneRotationRad", D(f.inPlaneRotationRad)), ("matchedVectors", "\(f.matchedVectors)"),
             ("totalVectors", "\(f.totalVectors)"), ("meanDistance", D(f.meanDistance)),
             ("chanceMatchedVectors", D(f.chanceMatchedVectors)), ("sweepMedianFraction", D(f.sweepMedianFraction))])
    }
    var fitZoneBlocks: [String] = []
    let fitA = PhaseVectorMatcher.fitZoneAxis(
        bragg: scanA.bragg, crystal: .aluminum, referenceSettings: PhaseReferenceSettings(), settings: defaults,
        originX: scanA.originX, originY: scanA.originY, invAngstromPerPixel: scanA.scale)
    fitZoneBlocks.append(obj([("scan", str(scanA.name)), ("arguments", str("Al, default PhaseReferenceSettings, default PhaseVectorSettings, 62 axes, step 5, sampleLimit 400")),
                              ("count", "\(fitA.count)"), ("fits", arr(fitA.map(fitJSON)))]))
    print("apps-parity-4: fitZoneAxis scan a top 3: " + fitA.prefix(3).map { "[\($0.zoneAxis.x)\($0.zoneAxis.y)\($0.zoneAxis.z)] \($0.matchedVectors)/\($0.totalVectors)" }.joined(separator: "; "))
    let fitA2 = PhaseVectorMatcher.fitZoneAxis(
        bragg: scanA.bragg, crystal: .aluminum, referenceSettings: PhaseReferenceSettings(), settings: defaults,
        originX: scanA.originX, originY: scanA.originY, invAngstromPerPixel: scanA.scale,
        candidateAxes: [z001, z011, z111, SIMD3(1, 1, 2)], inPlaneStepDeg: 10, sampleLimit: 100)
    fitZoneBlocks.append(obj([("scan", str(scanA.name)), ("arguments", str("Al, 4 axes, step 10, sampleLimit 100")),
                              ("count", "\(fitA2.count)"), ("fits", arr(fitA2.map(fitJSON)))]))
    let fitB = PhaseVectorMatcher.fitZoneAxis(
        bragg: scanB.bragg, crystal: .betaDoublePrime, referenceSettings: PhaseReferenceSettings(), settings: defaults,
        originX: scanB.originX, originY: scanB.originY, invAngstromPerPixel: scanB.scale,
        candidateAxes: Array(PhaseReferenceLibrary.lowIndexZoneAxes.prefix(20)), inPlaneStepDeg: 15, sampleLimit: 50)
    fitZoneBlocks.append(obj([("scan", str(scanB.name)), ("arguments", str("beta'' monoclinic, first 20 axes, step 15, sampleLimit 50")),
                              ("count", "\(fitB.count)"), ("fits", arr(fitB.map(fitJSON)))]))
    sections.append(("zoneAxisFits", arr(fitZoneBlocks)))

    // claims: every position of scans a, c, d and b from the default runs
    func claimCode(_ c: PhaseDiskClaim) -> Int {
        switch c {
        case .notMatched: return -1
        case .matrix: return -2
        case .unexplained: return -3
        case .phase(let p): return p
        }
    }
    var claimBlocks: [String] = []
    for s in [scanA, scanB, scanC, scanD] {
        let runID = s.name == scanD.name ? "\(s.name)/main/default/fit" : "\(s.name)/main/default"
        let run = runs[runIndex[runID]!]
        guard let map = run.map else { continue }
        let matrixEntry: PhaseOrientationReference? = map.matrixEntryIndex >= 0 ? main.entries[map.matrixEntryIndex] : nil
        var perPosition: [String] = []
        for (i, p) in s.peaks.enumerated() {
            let claims = PhaseVectorMatcher.claims(
                peaks: p, originX: s.originX, originY: s.originY, invAngstromPerPixel: s.scale, settings: defaults,
                matrixEntry: matrixEntry, result: map.results[i], library: main)
            perPosition.append(arr(claims.map { "\(claimCode($0))" }))
        }
        claimBlocks.append(obj([("scan", str(s.name)), ("run", str(runID)),
                                ("encoding", str("-1 notMatched, -2 matrix, -3 unexplained, >= 0 phase index")),
                                ("claims", "[\n" + perPosition.joined(separator: ",\n") + "\n]")]))
    }
    // and one run with a non-default radius so the claim logic sees different settings
    do {
        let s = scanA
        let runID = "\(s.name)/main/claims-custom"
        let run = runs[runIndex[runID]!]
        if let map = run.map {
            let matrixEntry: PhaseOrientationReference? = map.matrixEntryIndex >= 0 ? main.entries[map.matrixEntryIndex] : nil
            let perPosition = s.peaks.enumerated().map { (i, p) in
                arr(PhaseVectorMatcher.claims(peaks: p, originX: s.originX, originY: s.originY, invAngstromPerPixel: s.scale, settings: claimsSettings,
                                              matrixEntry: matrixEntry, result: map.results[i], library: main).map { "\(claimCode($0))" })
            }
            claimBlocks.append(obj([("scan", str(s.name)), ("run", str(runID)), ("settings", settingsJSON(claimsSettings)),
                                    ("claims", "[\n" + perPosition.joined(separator: ",\n") + "\n]")]))
        }
    }
    sections.append(("claims", "[\n" + claimBlocks.joined(separator: ",\n") + "\n]"))

    sections.append(("resolution", arr([0.0457, 0.01, 0.002, 0, -1, Double.nan].map { scale in
        let r = PhaseVectorResolution(settings: defaults, invAngstromPerPixel: scale)
        let scaled = r.scaledToDetector(defaults)
        return obj([("invAngstromPerPixel", D(scale)), ("pairRadiusPixels", D(r.pairRadiusPixels)),
                    ("matrixRemovalPixels", D(r.matrixRemovalPixels)), ("notIndexedAbovePixels", D(r.notIndexedAbovePixels)),
                    ("advice", r.advice.map(str) ?? "null"), ("scaledToDetector", settingsJSON(scaled))])
    })))

    // ---- nearest / score on library entries
    var rngN = Rng(seed: 20_261_019)
    var nearestBlocks: [String] = []
    let stride = max(1, main.entries.count / 30)
    let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
    for e in Swift.stride(from: 0, to: main.entries.count, by: stride) {
        let entry = main.entries[e]
        guard !entry.vectors.isEmpty else { continue }
        var us: [SIMD2<Double>] = []
        for k in 0..<12 {
            if k % 2 == 0 {
                let ref = entry.vectors[rngN.int(entry.vectors.count)].q
                us.append(ref + SIMD2(rngN.range(-0.03, 0.03), rngN.range(-0.03, 0.03)))
            } else {
                let r = rngN.range(0.2, 1.3), a = rngN.range(0, 2 * Double.pi)
                us.append(SIMD2(r * cosI(a), r * sinI(a)))
            }
        }
        func nn(_ radius: Double) -> String {
            arr(us.map { u in
                if let h = PhaseVectorMatcher.nearest(u, in: entry.vectors, radius: radius) { return "[\(h.index),\(D(h.distance))]" }
                return "null"
            })
        }
        let sc = PhaseVectorMatcher.score(vectors: us, against: entry, pairRadius: 0.02, scratch: scratch)
        nearestBlocks.append(obj([
            ("entryIndex", "\(e)"), ("us", arr(us.map { "[\"\(hex($0.x))\",\"\(hex($0.y))\"]" })),
            ("nearest0p02", nn(0.02)), ("nearest0p05", nn(0.05)),
            ("score0p02", sc.map { obj([("score", D($0.score)), ("matched", "\($0.matched)"), ("uniqueReferences", "\($0.uniqueReferences)")]) } ?? "null"),
            ("friedel0p02", PhaseVectorMatcher.containsFriedelPair(us, radius: 0.02) ? "true" : "false"),
        ]))
    }
    sections.append(("nearest", "[\n" + nearestBlocks.joined(separator: ",\n") + "\n]"))

    // ---- 4. Euler
    var rngE = Rng(seed: 20_261_020)
    var eulerBlocks: [String] = []
    for i in 0..<200 {
        let m = quaternionMatrix(rngE.gauss(), rngE.gauss(), rngE.gauss(), rngE.gauss())
        eulerBlocks.append(eulerCase("random \(i)", m))
    }
    let gimbalPhi: [Double] = [0, 1e-9, -1e-9, 5e-9, 1e-8, 2e-8, 1e-7, .pi, .pi - 1e-9, .pi - 1e-8, .pi - 2e-8, .pi - 1e-7, 1.0, .pi / 2]
    let anglePairs: [(Double, Double)] = [(0.3, 1.1), (-2.5, 0.7), (3.0, -3.0), (0, 0), (1.5, -0.4)]
    for (pi, phi) in gimbalPhi.enumerated() {
        for (ai, (p1, p2)) in anglePairs.enumerated() {
            let r = rotZ(p2) * rotX(phi) * rotZ(p1)             // crystal to lab, extrinsic zxz
            eulerBlocks.append(eulerCase("gimbal Phi[\(pi)] angles[\(ai)] Phi=\(phi) phi1=\(p1) phi2=\(p2)", simd_transpose(r)))
        }
    }
    sections.append(("euler", "[\n" + eulerBlocks.joined(separator: ",\n") + "\n]"))

    // detector bases for the 62 axes of Al and of the triclinic cell, and the in-plane matrix
    var basisBlocks: [String] = []
    for (name, crystal) in [("Al", Crystal.aluminum), ("triclinic", triclinic), ("beta-double-prime", Crystal.betaDoublePrime)] {
        var items: [String] = []
        for axis in axes62 {
            guard let n = PhaseReferenceLibrary.cartesianZoneAxis(axis, crystal: crystal) else { continue }
            let basis = ACOMOrientation.detectorBasis(zoneAxis: n)
            let angles: [Float] = [0, 0.5, -2.2]
            items.append(obj([("axis", i3(axis)), ("n", arr([D(n.x), D(n.y), D(n.z)])), ("basis", mat(basis)),
                              ("matrices", arr(angles.map { obj([("angle", F($0)), ("matrix", mat(ACOMOrientation.matrix(detectorBasis: basis, inPlaneAngle: $0)))]) }))]))
        }
        basisBlocks.append(obj([("crystal", str(name)), ("axes", "[\n" + items.joined(separator: ",\n") + "\n]")]))
    }
    sections.append(("detectorBases", "[\n" + basisBlocks.joined(separator: ",\n") + "\n]"))

    let out = "{\n" + sections.map { "\"\($0.0)\":\($0.1)" }.joined(separator: ",\n") + "\n}\n"
    try out.write(toFile: outPath, atomically: true, encoding: .utf8)
    print("apps-parity-4: wrote \(outPath), \(out.utf8.count) bytes, \(runs.count) runs, \(Int(Date().timeIntervalSince(started))) s (not in the dump: it must be deterministic)")
}

let outPath = CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "slice-4-apps.json"
try runAll(outPath: outPath)
