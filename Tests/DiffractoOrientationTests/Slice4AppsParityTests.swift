//
//  Slice4AppsParityTests.swift
//  Role: gate G1 of slice 4 (registration docs/registration-slice-4.md, section 7, and the amendment of 2026-10-15 16:08):
//        the kit reproduces, bit for bit, the dump of the apps' own compiled code (tools/slice-4/apps-parity, written from
//        mac4DSTEM 7d6fdbd1) in Fixtures/slice-4-apps.json. The dump is the -O build and these tests run at -c release, so
//        there is no exception: every Double is compared by `bitPattern`, every Float by `bitPattern`, every integer and
//        string by equality. A difference is a finding of the port, never a tolerance to widen.
//
//        Blocks of the dump and the tests that read them:
//          reflections    the reflection list a library is made from (order, hkl, g, |g|, Re F, Im F, |F|^2), 6 phases
//          libraries      `main`, `main-or` (against main's entries), `al-all-axes`, `oblique`: every entry's phase,
//                         axis, rotation, chance fraction, FNV-1a hash (the harness's, reproduced) and the vectors dumped
//          scans          the experimental vectors of every position (two radius / reach pairs)
//          runs           the 72 `map` runs (every PhaseVectorResult field, the fitted matrix entry, the counts); the
//                         malformed scans give nil; the default runs again at maximumWorkers 1 and 8
//          matrixFits, zoneAxisFits, claims, resolution, inPlaneSteps, nearest/score/Friedel, euler, detectorBases,
//          lowIndexZoneAxes
//
//        Inputs are rebuilt from the dump alone (cells and sites, Float pixel bits, origin, scale, the settings of each
//        run); no generator is replayed. Public API only.
//
//        Fail-first hook (test-only): the environment variable S4G1_FLIP=<block>[:n] flips the n-th expected value
//        (default 5) of the named block in memory, so the block must turn red. Unset, the tests are the plain G1 tests.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

// MARK: - Reading the dump

private typealias JSON = [String: Any]

private func u64(_ x: Any?) -> UInt64 {
    if let s = x as? String { return UInt64(s, radix: 16)! }
    if let j = x as? JSON { return u64(j["b"]) }
    preconditionFailure("not a bit pattern: \(String(describing: x))")
}
private func u32(_ x: Any?) -> UInt32 {
    if let s = x as? String { return UInt32(s, radix: 16)! }
    if let j = x as? JSON { return u32(j["b"]) }
    preconditionFailure("not a bit pattern: \(String(describing: x))")
}
private func dbl(_ x: Any?) -> Double { Double(bitPattern: u64(x)) }
private func flt(_ x: Any?) -> Float { Float(bitPattern: u32(x)) }
private func jint(_ x: Any?) -> Int { (x as! NSNumber).intValue }
private func bool(_ x: Any?) -> Bool { (x as! NSNumber).boolValue }
private func str(_ x: Any?) -> String { x as! String }
private func arr(_ x: Any?) -> [Any] { x as! [Any] }
private func objs(_ x: Any?) -> [JSON] { x as! [JSON] }
private func ints3(_ x: Any?) -> SIMD3<Int> { let a = arr(x); return SIMD3(jint(a[0]), jint(a[1]), jint(a[2])) }
private func isNull(_ x: Any?) -> Bool { x == nil || x is NSNull }
private func hex(_ x: UInt64) -> String {
    let s = String(x, radix: 16)
    return String(repeating: "0", count: 16 - s.count) + s
}
private func hex(_ x: UInt32) -> String {
    let s = String(x, radix: 16)
    return String(repeating: "0", count: 8 - s.count) + s
}
/// A column-major 3 x 3 matrix as the harness wrote it (`order`, nine Double bit patterns).
private func matrixBits(_ x: Any?) -> [UInt64] { arr((x as! JSON)["b"]).map { u64($0) } }
private func matrix3(_ x: Any?) -> simd_double3x3 {
    let m = matrixBits(x).map { Double(bitPattern: $0) }
    return simd_double3x3(columns: (SIMD3(m[0], m[1], m[2]), SIMD3(m[3], m[4], m[5]), SIMD3(m[6], m[7], m[8])))
}
private func flat(_ m: simd_double3x3) -> [Double] {
    [m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.1.x, m.columns.1.y, m.columns.1.z,
     m.columns.2.x, m.columns.2.y, m.columns.2.z]
}

// MARK: - The tally: every comparison counted, the first differences named

private struct Tally {
    let block: String
    var compared = 0
    var failures = 0
    var messages: [String] = []
    private let flipAt: Int?

    init(_ block: String) {
        self.block = block
        var flip: Int?
        if let spec = ProcessInfo.processInfo.environment["S4G1_FLIP"] {
            let parts = spec.split(separator: ":", maxSplits: 1).map(String.init)
            if parts.first == block { flip = parts.count > 1 ? Int(parts[1]) : 5 }
        }
        flipAt = flip
    }

    private var flipNow: Bool { compared == flipAt }

    private mutating func fail(_ message: String) {
        failures += 1
        if messages.count < 6 { messages.append(message) }
    }

    mutating func bits(_ got: UInt64, _ want: UInt64, _ path: @autoclosure () -> String) {
        compared += 1
        let w = flipNow ? want ^ 1 : want
        if got != w { fail("\(path()): kit \(hex(got)) dump \(hex(w))") }
    }
    mutating func bits(_ got: UInt32, _ want: UInt32, _ path: @autoclosure () -> String) {
        compared += 1
        let w = flipNow ? want ^ 1 : want
        if got != w { fail("\(path()): kit \(hex(got)) dump \(hex(w))") }
    }
    /// A Double against the dump's `{"v", "b"}` object or hex string.
    mutating func d(_ got: Double, _ want: Any?, _ path: @autoclosure () -> String) {
        bits(got.bitPattern, u64(want), path())
    }
    mutating func f(_ got: Float, _ want: Any?, _ path: @autoclosure () -> String) {
        bits(got.bitPattern, u32(want), path())
    }
    mutating func eq<T: Equatable>(_ got: T, _ want: T, _ path: @autoclosure () -> String) {
        compared += 1
        if got != want || flipNow { fail("\(path()): kit \(got) dump \(want)" + (flipNow ? " (expected value flipped by S4G1_FLIP)" : "")) }
    }
    mutating func int(_ got: Int, _ want: Any?, _ path: @autoclosure () -> String) {
        eq(got, jint(want), path())
    }
    mutating func note(_ message: String) { fail(message) }

    var report: Comment {
        Comment(rawValue: "\(block): \(failures) of \(compared) values differ. " + messages.joined(separator: " | "))
    }
    func printed(_ label: String) {
        print("S4G1 \(block) \(label): compared \(compared) values, \(failures) differences")
    }
}

// MARK: - Structures, settings and scans rebuilt from the dump

private func makeStructure(_ c: JSON) -> CrystalStructure {
    let cell = UnitCell(a: dbl(c["a"]), b: dbl(c["b"]), c: dbl(c["c"]),
                        alphaDeg: dbl(c["alphaDeg"]), betaDeg: dbl(c["betaDeg"]), gammaDeg: dbl(c["gammaDeg"]))
    let sites = arr(c["sites"]).map { row -> AtomSite in
        let s = arr(row)
        return AtomSite(z: jint(s[0]), fractional: SIMD3(dbl(s[1]), dbl(s[2]), dbl(s[3])), occupancy: dbl(s[4]))
    }
    return CrystalStructure(cell: cell, sites: sites, operations: [.identity])
}

private func latticeVector(_ x: Any?) -> LatticeVector {
    let a = arr(x)
    let v = ints3(a[1])
    return str(a[0]) == "plane" ? .plane(v) : .direction(v)
}

private func definitions(_ phases: [JSON]) -> [PhaseDefinition] {
    phases.map { p in
        PhaseDefinition(
            id: str(p["id"]), displayName: str(p["displayName"]), crystal: makeStructure(p["crystal"] as! JSON),
            role: PhaseDefinition.Role(rawValue: str(p["role"]))!,
            zoneAxes: arr(p["zoneAxes"]).map(ints3),
            orientationRelationships: objs(p["orientationRelationships"]).map {
                OrientationRelationship(candidate: latticeVector($0["candidate"]), matrix: latticeVector($0["matrix"]))
            },
            excitationSlabInvAngstrom: isNull(p["excitationSlabInvAngstrom"]) ? nil : dbl(p["excitationSlabInvAngstrom"]))
    }
}

private func referenceSettings(_ j: JSON) -> PhaseReferenceSettings {
    var s = PhaseReferenceSettings()
    s.kMaxInvAngstrom = dbl(j["kMaxInvAngstrom"])
    s.excitationSlabInvAngstrom = dbl(j["excitationSlabInvAngstrom"])
    s.excitationWidthInvAngstrom = dbl(j["excitationWidthInvAngstrom"])
    s.minimumIntensityFraction = dbl(j["minimumIntensityFraction"])
    s.maximumVectorsPerEntry = jint(j["maximumVectorsPerEntry"])
    s.inPlaneStepDeg = dbl(j["inPlaneStepDeg"])
    s.maximumEntries = jint(j["maximumEntries"])
    return s
}

private func vectorSettings(_ j: JSON) -> PhaseVectorSettings {
    var s = PhaseVectorSettings()
    s.directBeamRadiusInvAngstrom = dbl(j["directBeamRadiusInvAngstrom"])
    s.maximumVectorInvAngstrom = dbl(j["maximumVectorInvAngstrom"])
    s.matrixToleranceInvAngstrom = dbl(j["matrixToleranceInvAngstrom"])
    s.pairRadiusInvAngstrom = dbl(j["pairRadiusInvAngstrom"])
    s.minimumVectors = jint(j["minimumVectors"])
    s.matrixFallbackExplainedFraction = dbl(j["matrixFallbackExplainedFraction"])
    s.minimumMatchedVectors = jint(j["minimumMatchedVectors"])
    s.friedelPairMinimumMatchedVectors = jint(j["friedelPairMinimumMatchedVectors"])
    s.chanceMatchMultiple = dbl(j["chanceMatchMultiple"])
    s.notIndexedAboveInvAngstrom = dbl(j["notIndexedAboveInvAngstrom"])
    s.minimumPhaseContrastInvAngstrom = dbl(j["minimumPhaseContrastInvAngstrom"])
    s.orientationRelationshipToleranceDeg = dbl(j["orientationRelationshipToleranceDeg"])
    s.completenessAwareCrossPhaseRanking = bool(j["completenessAwareCrossPhaseRanking"])
    s.classificationRule = PhaseVectorSettings.ClassificationRule(rawValue: str(j["classificationRule"]))!
    s.directMatrixMaximumVectors = jint(j["directMatrixMaximumVectors"])
    s.residualCutoffInvAngstrom = dbl(j["residualCutoffInvAngstrom"])
    s.knownVariantsMinimumSpecificReflections = jint(j["knownVariantsMinimumSpecificReflections"])
    return s
}

/// One dumped scan: the detector-pixel positions as the apps held them (Float bits), widened to Double as the registration's
/// substitution 2 says an application does, plus the origin and the scale.
private struct DumpedScan: Sendable {
    let name: String
    let width: Int
    let height: Int
    let originX: Float
    let originY: Float
    let scale: Double
    let peaks: [[SIMD2<Double>]]
    var coordinates: PeakCoordinates {
        .detectorPixels(originX: originX, originY: originY, invAngstromPerPixel: scale)
    }
    var fixture: ScanFixture { ScanFixture(scanWidth: width, scanHeight: height, peaks: peaks) }
}

private func dumpedScan(_ j: JSON) -> DumpedScan {
    let peaks = arr(j["peaks"]).map { position in
        arr(position).map { p -> SIMD2<Double> in
            let a = arr(p)
            return SIMD2(Double(flt(a[2])), Double(flt(a[3])))
        }
    }
    return DumpedScan(name: str(j["name"]), width: jint(j["scanWidth"]), height: jint(j["scanHeight"]),
                      originX: flt(j["originX"]), originY: flt(j["originY"]), scale: dbl(j["invAngstromPerPixel"]),
                      peaks: peaks)
}

private func entryHash(_ e: PhaseOrientationReference) -> UInt64 {
    // The harness's FNV-1a 64 (main.swift, entryHash) over little-endian 64-bit words: phaseIndex, u, v, w, the rotation's
    // Double bits, the vector count, then per vector h, k, l (Int64 bit patterns), q.x, q.y, length, relativeIntensity.
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

/// The dump, parsed once, and everything built from it that several tests share.
private enum G1 {
    nonisolated(unsafe) static let dump: JSON = {
        let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .appendingPathComponent("Fixtures").appendingPathComponent("slice-4-apps.json")
        let data = try! Data(contentsOf: url)
        return try! JSONSerialization.jsonObject(with: data) as! JSON
    }()

    nonisolated(unsafe) static let reflectionBlocks: [JSON] = objs(dump["reflections"])
    nonisolated(unsafe) static let libraryBlocks: [JSON] = objs(dump["libraries"])
    nonisolated(unsafe) static let runBlocks: [JSON] = objs(dump["runs"])

    static func libraryBlock(_ name: String) -> JSON { libraryBlocks.first { str($0["name"]) == name }! }
    static func runBlock(_ id: String) -> JSON { runBlocks.first { str($0["id"]) == id }! }

    static let scans: [String: DumpedScan] = {
        var out: [String: DumpedScan] = [:]
        for j in objs(dump["scans"]) { let s = dumpedScan(j); out[s.name] = s }
        return out
    }()
    static func scan(_ name: String) -> DumpedScan { scans[name]! }

    static func structure(named name: String) -> CrystalStructure {
        makeStructure(reflectionBlocks.first { str($0["phase"]) == name }!["crystal"] as! JSON)
    }

    static func build(_ name: String) -> PhaseReferenceLibrary {
        let b = libraryBlock(name)
        let settings = isNull(b["referenceSettings"]) ? PhaseReferenceSettings()
            : referenceSettings(b["referenceSettings"] as! JSON)
        return try! PhaseReferenceLibrary.build(phases: definitions(objs(b["phases"])), settings: settings)
    }
    /// `main` (every entry dumped); `main-or` carries no entries of its own (the harness records they equal main's).
    static let main: PhaseReferenceLibrary = build("main")
    static let mainOR: PhaseReferenceLibrary = {
        // main-or has no referenceSettings key: it is built at the defaults, like main.
        let b = libraryBlock("main-or")
        return try! PhaseReferenceLibrary.build(phases: definitions(objs(b["phases"])), settings: PhaseReferenceSettings())
    }()
}

struct S4G1Case: Sendable, CustomTestStringConvertible {
    let index: Int
    let name: String
    var testDescription: String { name }
}

// MARK: - The tests

@Suite(.serialized) struct Slice4AppsParity {

    // MARK: reflections (the list a library is made from, in the port's order)

    private static let reflectionCases: [S4G1Case] =
        G1.reflectionBlocks.enumerated().map { S4G1Case(index: $0.offset, name: str($0.element["phase"])) }

    @Test("reflections: order, hkl, g, |g|, Re F, Im F, |F|^2 bit for bit", arguments: reflectionCases)
    func reflections(_ c: S4G1Case) throws {
        let block = G1.reflectionBlocks[c.index]
        var t = Tally("reflections")
        let s = makeStructure(block["crystal"] as! JSON)
        let kit = try PhaseReferenceLibrary.reflections(of: s, kMaxInvAngstrom: dbl(block["kMax"]))
        t.int(kit.count, block["count"], "\(c.name) count")
        let hkl = arr(block["hkl"]), gLength = arr(block["gLength"]), g = arr(block["g"])
        let re = arr(block["re"]), im = arr(block["im"]), intensity = arr(block["intensity"])
        #expect(kit.count == hkl.count)
        for i in 0..<min(kit.count, hkl.count) {
            let r = kit[i]
            let p = "\(c.name) reflection \(i) \(r.reflection.hkl)"
            let want = ints3(hkl[i])
            t.eq(r.reflection.hkl.h, want.x, "\(p) h"); t.eq(r.reflection.hkl.k, want.y, "\(p) k")
            t.eq(r.reflection.hkl.l, want.z, "\(p) l")
            let gw = arr(g[i])
            t.d(r.reflection.g.x, gw[0], "\(p) g.x"); t.d(r.reflection.g.y, gw[1], "\(p) g.y")
            t.d(r.reflection.g.z, gw[2], "\(p) g.z")
            t.d(r.reflection.gLength, gLength[i], "\(p) |g|")
            t.d(r.structureFactor.real, re[i], "\(p) Re F")
            t.d(r.structureFactor.imaginary, im[i], "\(p) Im F")
            t.d(r.structureFactor.intensity, intensity[i], "\(p) |F|^2")
        }
        t.printed(c.name)
        #expect(t.compared == 1 + 10 * hkl.count)
        #expect(t.failures == 0, t.report)
    }

    /// The lattice rows the reflections' g and the zone-axis directions are built from.
    @Test("reflections: the cells' real and reciprocal lattice rows equal the apps' latReal and latInv", arguments: reflectionCases)
    func latticeRows(_ c: S4G1Case) {
        let crystal = G1.reflectionBlocks[c.index]["crystal"] as! JSON
        let s = makeStructure(crystal)
        var t = Tally("latticeRows")
        let real = arr(crystal["latReal"]), inv = arr(crystal["latInv"])
        for i in 0..<3 {
            let r = arr(real[i]), q = arr(inv[i])
            t.d(s.cell.realLattice[i].x, r[0], "\(c.name) latReal[\(i)].x"); t.d(s.cell.realLattice[i].y, r[1], "\(c.name) latReal[\(i)].y")
            t.d(s.cell.realLattice[i].z, r[2], "\(c.name) latReal[\(i)].z")
            t.d(s.cell.reciprocalLattice[i].x, q[0], "\(c.name) latInv[\(i)].x"); t.d(s.cell.reciprocalLattice[i].y, q[1], "\(c.name) latInv[\(i)].y")
            t.d(s.cell.reciprocalLattice[i].z, q[2], "\(c.name) latInv[\(i)].z")
        }
        t.printed(c.name)
        #expect(t.compared == 18)
        #expect(t.failures == 0, t.report)
    }

    // MARK: libraries

    private static let libraryNames = ["main", "main-or", "al-all-axes", "oblique"]

    @Test("libraries: every entry (phase, axis, rotation, chance fraction, hash, vectors) bit for bit", arguments: libraryNames)
    func libraries(_ name: String) {
        // main-or has no entries of its own: the harness records they equal main's, and they are compared with main's dump.
        let dumped = G1.libraryBlock(name == "main-or" ? "main" : name)
        let library: PhaseReferenceLibrary = name == "main" ? G1.main
            : name == "main-or" ? G1.mainOR : G1.build(name)
        var t = Tally("libraries")
        if name == "main-or" { t.eq(bool(G1.libraryBlock("main-or")["sameEntriesAsMain"]), true, "main-or sameEntriesAsMain") }
        let entries = objs(dumped["entries"])
        t.int(library.entries.count, dumped["entryCount"], "\(name) entryCount")
        t.int(library.matrixPhaseIndex, dumped["matrixPhaseIndex"], "\(name) matrixPhaseIndex")
        t.eq(library.matrixEntryIndices, arr(dumped["matrixEntryIndices"]).map { jint($0) }, "\(name) matrixEntryIndices")
        t.int(library.candidateEntryIndices.count, dumped["candidateEntryCount"], "\(name) candidateEntryCount")
        #expect(library.entries.count == entries.count)
        var counter: [SIMD4<Int>: Int] = [:]
        var fullEntries = 0, vectorValues = 0
        for (i, e) in library.entries.enumerated() where i < entries.count {
            let w = entries[i]
            let key = SIMD4(e.phaseIndex, e.zoneAxis.x, e.zoneAxis.y, e.zoneAxis.z)
            let k = counter[key, default: 0]
            counter[key] = k + 1
            let p = "\(name) entry \(i) [\(e.zoneAxis.x)\(e.zoneAxis.y)\(e.zoneAxis.z)] phase \(e.phaseIndex) k=\(k)"
            t.int(e.phaseIndex, w["p"], "\(p) phaseIndex")
            t.eq(e.zoneAxis, ints3(w["u"]), "\(p) zoneAxis")
            t.int(k, w["k"], "\(p) rotation index")
            t.d(e.inPlaneRotationRad, w["r"], "\(p) inPlaneRotationRad")
            t.int(e.vectors.count, w["n"], "\(p) vector count")
            t.bits(entryHash(e), u64(w["hash"]), "\(p) FNV-1a hash")
            t.d(e.chanceMatchFraction(pairRadius: 0.02, accessibleRadius: 0.9), w["chance"], "\(p) chanceMatchFraction")
            if let v = w["v"] {
                fullEntries += 1
                let rows = arr(v)
                t.eq(e.vectors.count, rows.count, "\(p) dumped vector rows")
                for (j, vec) in e.vectors.enumerated() where j < rows.count {
                    let r = arr(rows[j])
                    t.int(vec.h, r[0], "\(p) vector \(j) h"); t.int(vec.k, r[1], "\(p) vector \(j) k")
                    t.int(vec.l, r[2], "\(p) vector \(j) l")
                    t.d(vec.q.x, r[3], "\(p) vector \(j) q.x"); t.d(vec.q.y, r[4], "\(p) vector \(j) q.y")
                    t.d(vec.length, r[5], "\(p) vector \(j) length"); t.d(vec.relativeIntensity, r[6], "\(p) vector \(j) relativeIntensity")
                    vectorValues += 7
                }
            }
        }
        t.printed("\(name) (\(library.entries.count) entries, \(fullEntries) with every vector dumped, \(vectorValues) vector values)")
        #expect(fullEntries == (name == "main" || name == "main-or" ? 2520 : 5 * (library.entries.count / 180)))
        #expect(t.failures == 0, t.report)
    }

    @Test("libraries: the 49 low-index zone axes in the apps' order")
    func lowIndexAxes() {
        var t = Tally("lowIndexAxes")
        let want = arr(G1.dump["lowIndexZoneAxes"]).map(ints3)
        t.eq(PhaseReferenceLibrary.lowIndexZoneAxes.count, want.count, "count")
        for (i, a) in PhaseReferenceLibrary.lowIndexZoneAxes.enumerated() where i < want.count {
            t.eq(a, want[i], "axis \(i)")
        }
        t.printed("49 axes")
        #expect(want.count == 49)
        #expect(t.failures == 0, t.report)
    }

    @Test("libraries: inPlaneSteps for 15 step values")
    func inPlaneSteps() {
        var t = Tally("inPlaneSteps")
        for row in objs(G1.dump["inPlaneSteps"]) {
            var s = PhaseReferenceSettings()
            s.inPlaneStepDeg = dbl(row["step"])
            let kit = PhaseReferenceLibrary.inPlaneSteps(s)
            let want = arr(row["radians"])
            let p = "step \(s.inPlaneStepDeg)"
            t.int(kit.count, row["count"], "\(p) count")
            t.eq(kit.count, want.count, "\(p) rows")
            for (i, r) in kit.enumerated() where i < want.count { t.d(r, want[i], "\(p) k=\(i)") }
        }
        t.printed("15 steps")
        #expect(objs(G1.dump["inPlaneSteps"]).count == 15)
        #expect(t.failures == 0, t.report)
    }

    // MARK: the experimental vectors and the cartesian zone axis

    @Test("scans: the experimental vectors of every position, two radius / reach pairs", arguments: ["a-matrix-with-precipitates",
        "b-random-vectors", "c-friedel-pairs", "d-special-positions", "e-malformed-count", "f-zero-width"])
    func experimentalVectors(_ name: String) {
        let j = objs(G1.dump["scans"]).first { str($0["name"]) == name }!
        let scan = G1.scan(name)
        var t = Tally("experimentalVectors")
        let defaults = PhaseVectorSettings()
        let pairs: [(key: String, radius: Double, reach: Double)] = [
            ("vectorsDefault", defaults.directBeamRadiusInvAngstrom, defaults.maximumVectorInvAngstrom),
            ("vectorsRadius0Max0p9", 0, 0.9),
        ]
        t.int(scan.width, j["scanWidth"], "\(name) width")
        for pair in pairs {
            let want = arr(j[pair.key])
            t.eq(scan.peaks.count == 0 ? 0 : want.count, scan.peaks.count == 0 ? 0 : scan.peaks.count, "\(name) \(pair.key) positions")
            for (i, peaks) in scan.peaks.enumerated() where i < want.count {
                // both the detector-pixel form with the origin as given and the PeakCoordinates form
                let kit = PhaseVectorMatcher.experimentalVectors(
                    peaks: peaks, originX: scan.originX, originY: scan.originY, invAngstromPerPixel: scan.scale,
                    directBeamRadiusInvAngstrom: pair.radius, maximumVectorInvAngstrom: pair.reach)
                let viaCoordinates = PhaseVectorMatcher.experimentalVectors(
                    peaks: peaks, coordinates: scan.coordinates,
                    directBeamRadiusInvAngstrom: pair.radius, maximumVectorInvAngstrom: pair.reach)
                let w = arr(want[i])
                let p = "\(name) \(pair.key) position \(i)"
                t.eq(kit.count, w.count, "\(p) count")
                t.eq(viaCoordinates.count, w.count, "\(p) count (coordinates form)")
                for (n, q) in kit.enumerated() where n < w.count {
                    let a = arr(w[n])
                    t.d(q.x, a[0], "\(p) vector \(n) x"); t.d(q.y, a[1], "\(p) vector \(n) y")
                }
                for (n, q) in viaCoordinates.enumerated() where n < w.count {
                    let a = arr(w[n])
                    t.d(q.x, a[0], "\(p) vector \(n) x (coordinates form)"); t.d(q.y, a[1], "\(p) vector \(n) y (coordinates form)")
                }
            }
        }
        t.printed(name)
        #expect(t.failures == 0, t.report)
    }

    // MARK: map

    private static let runCases: [S4G1Case] =
        G1.runBlocks.enumerated().map { S4G1Case(index: $0.offset, name: str($0.element["id"])) }

    /// One `map` result list against the dump. Returns nothing; differences go in the tally.
    private func compare(_ map: PhaseMap?, with run: JSON, workers: Int?, tally t: inout Tally) {
        let id = str(run["id"]) + (workers.map { " [workers \($0)]" } ?? "")
        t.eq(map == nil, bool(run["returnedNil"]), "\(id) returnedNil")
        guard let map, !bool(run["returnedNil"]) else { return }
        t.int(map.matrixEntryIndex, run["matrixEntryIndex"], "\(id) matrixEntryIndex")
        t.int(map.matrixPhaseIndex, run["matrixPhaseIndex"], "\(id) matrixPhaseIndex")
        t.eq(map.phaseNames, arr(run["phaseNames"]).map { str($0) }, "\(id) phaseNames")
        t.eq(map.phaseCounts, arr(run["phaseCounts"]).map { jint($0) }, "\(id) phaseCounts")
        let counts = PhaseVerdict.allCases.map { map.count(of: $0) }
        t.eq(counts, arr(run["verdictCounts"]).map { jint($0) }, "\(id) verdictCounts")
        let rows = arr(run["results"])
        t.eq(map.results.count, rows.count, "\(id) result count")
        for (i, r) in map.results.enumerated() where i < rows.count {
            let w = arr(rows[i])
            let p = "\(id) position \(i)"
            t.int(Int(r.verdict.rawValue), w[0], "\(p) verdict")
            t.int(Int(r.phaseIndex), w[1], "\(p) phaseIndex")
            t.int(Int(r.entryIndex), w[2], "\(p) entryIndex")
            t.f(r.score, w[3], "\(p) score")
            t.int(Int(r.matchedCount), w[5], "\(p) matchedCount")
            t.int(Int(r.survivingCount), w[6], "\(p) survivingCount")
            t.int(Int(r.removedCount), w[7], "\(p) removedCount")
            t.int(Int(r.runnerUpPhaseIndex), w[8], "\(p) runnerUpPhaseIndex")
            t.f(r.runnerUpScore, w[9], "\(p) runnerUpScore")
        }
    }

    private func runLibrary(_ run: JSON) -> PhaseReferenceLibrary {
        str(run["library"]) == "main-or" ? G1.mainOR : G1.main
    }

    private func execute(_ run: JSON, workers: Int?) -> PhaseMap? {
        let scan = G1.scan(str(run["scan"]))
        let matrix: Int? = isNull(run["matrixArgument"]) ? nil : jint(run["matrixArgument"])
        return PhaseVectorMatcher.map(
            scan: scan.peaks, scanWidth: scan.width, scanHeight: scan.height, library: runLibrary(run),
            settings: vectorSettings(run["settings"] as! JSON), coordinates: scan.coordinates,
            matrixEntryIndex: matrix, maximumWorkers: workers)
    }

    @Test("map: every PhaseVectorResult field of every position, the fitted matrix entry and the counts", arguments: runCases)
    func map(_ c: S4G1Case) {
        let run = G1.runBlocks[c.index]
        var t = Tally("runs")
        compare(execute(run, workers: 4), with: run, workers: 4, tally: &t)
        t.printed(c.name)
        #expect(t.failures == 0, t.report)
    }

    @Test("map: the malformed scans (count not width x height; zero width) give nil, as the apps' map does")
    func malformedScans() {
        var t = Tally("malformed")
        for id in ["e-malformed-count/main/default", "f-zero-width/main/default"] {
            let run = G1.runBlock(id)
            t.eq(bool(run["returnedNil"]), true, "\(id) dump returnedNil")
            t.eq(execute(run, workers: 2) == nil, true, "\(id) kit result is nil")
            // and through the pixel-form overload an application would call
            let scan = G1.scan(str(run["scan"]))
            let viaPixels = PhaseVectorMatcher.map(
                scan: scan.fixture, library: G1.main, settings: vectorSettings(run["settings"] as! JSON),
                originX: scan.originX, originY: scan.originY, invAngstromPerPixel: scan.scale, maximumWorkers: 2)
            t.eq(viaPixels == nil, true, "\(id) pixel-form overload is nil")
        }
        t.printed("2 malformed scans")
        #expect(t.failures == 0, t.report)
    }

    private static let workerRunIDs = [
        "a-matrix-with-precipitates/main/default",
        "a-matrix-with-precipitates/main/knownVariants",
        "b-random-vectors/main/default",
        "c-friedel-pairs/main/default",
        "d-special-positions/main/default/fit",
    ]

    @Test("map: the default runs at maximumWorkers 1 and at 8 equal the dump (the worker count moves no number)",
          arguments: workerRunIDs)
    func mapAtOneAndEightWorkers(_ id: String) {
        let run = G1.runBlock(id)
        var t = Tally("workers")
        compare(execute(run, workers: 1), with: run, workers: 1, tally: &t)
        compare(execute(run, workers: 8), with: run, workers: 8, tally: &t)
        t.printed(id)
        #expect(t.failures == 0, t.report)
    }

    // MARK: fits and claims

    @Test("matrixFits: fitMatrixOrientation on every scan (entry index, matched, mean distance)")
    func matrixFits() {
        var t = Tally("matrixFits")
        for row in objs(G1.dump["matrixFits"]) {
            let name = str(row["scan"])
            let scan = G1.scan(name)
            let kit = PhaseVectorMatcher.fitMatrixOrientation(
                scan: scan.peaks, scanWidth: scan.width, scanHeight: scan.height, library: G1.main,
                settings: PhaseVectorSettings(), coordinates: scan.coordinates)
            if name == "e-malformed-count" {
                // The source does not check the count and fits 8 positions of a 3 x 3 scan; the port refuses (registration
                // amendment 2026-10-15 16:08, item 2): a recorded refusal of an input the source does not define.
                t.eq(kit == nil, true, "\(name) refused by the port")
                t.eq(isNull(row["fit"]), false, "\(name) the source fitted it")
                continue
            }
            if isNull(row["fit"]) {
                t.eq(kit == nil, true, "\(name) nil")
                continue
            }
            let w = row["fit"] as! JSON
            guard let kit else { t.note("\(name): kit returned nil, dump has a fit"); continue }
            t.int(kit.entryIndex, w["entryIndex"], "\(name) entryIndex")
            t.int(kit.matched, w["matched"], "\(name) matched")
            t.d(kit.meanDistance, w["meanDistance"], "\(name) meanDistance")
        }
        t.printed("6 scans")
        #expect(t.failures == 0, t.report)
    }

    @Test("zoneAxisFits: fitZoneAxis, three sweeps, every fit and every field")
    func zoneAxisFits() throws {
        var t = Tally("zoneAxisFits")
        let rows = objs(G1.dump["zoneAxisFits"])
        #expect(rows.count == 3)
        let al = G1.structure(named: "Al"), beta = G1.structure(named: "beta-double-prime")
        let a = G1.scan("a-matrix-with-precipitates"), b = G1.scan("b-random-vectors")
        let z001 = SIMD3(0, 0, 1), z011 = SIMD3(0, 1, 1), z111 = SIMD3(1, 1, 1)
        let kits: [[PhaseVectorMatcher.ZoneAxisFit]] = [
            try PhaseVectorMatcher.fitZoneAxis(
                scan: a.peaks, scanWidth: a.width, scanHeight: a.height, crystal: al,
                referenceSettings: PhaseReferenceSettings(), settings: PhaseVectorSettings(), coordinates: a.coordinates),
            try PhaseVectorMatcher.fitZoneAxis(
                scan: a.peaks, scanWidth: a.width, scanHeight: a.height, crystal: al,
                referenceSettings: PhaseReferenceSettings(), settings: PhaseVectorSettings(), coordinates: a.coordinates,
                candidateAxes: [z001, z011, z111, SIMD3(1, 1, 2)], inPlaneStepDeg: 10, sampleLimit: 100),
            try PhaseVectorMatcher.fitZoneAxis(
                scan: b.peaks, scanWidth: b.width, scanHeight: b.height, crystal: beta,
                referenceSettings: PhaseReferenceSettings(), settings: PhaseVectorSettings(), coordinates: b.coordinates,
                candidateAxes: Array(PhaseReferenceLibrary.lowIndexZoneAxes.prefix(20)), inPlaneStepDeg: 15, sampleLimit: 50),
        ]
        for (s, row) in rows.enumerated() {
            let fits = objs(row["fits"])
            let p = "sweep \(s) (\(str(row["arguments"])))"
            t.int(kits[s].count, row["count"], "\(p) count")
            t.eq(kits[s].count, fits.count, "\(p) fits")
            for (i, f) in kits[s].enumerated() where i < fits.count {
                let w = fits[i]
                let q = "\(p) fit \(i)"
                t.eq(f.zoneAxis, ints3(w["zoneAxis"]), "\(q) zoneAxis")
                t.d(f.inPlaneRotationRad, w["inPlaneRotationRad"], "\(q) inPlaneRotationRad")
                t.int(f.matchedVectors, w["matchedVectors"], "\(q) matchedVectors")
                t.int(f.totalVectors, w["totalVectors"], "\(q) totalVectors")
                t.d(f.meanDistance, w["meanDistance"], "\(q) meanDistance")
                t.d(f.chanceMatchedVectors, w["chanceMatchedVectors"], "\(q) chanceMatchedVectors")
                t.d(f.sweepMedianFraction, w["sweepMedianFraction"], "\(q) sweepMedianFraction")
            }
        }
        // the port's refusal of a malformed scan (not in the source): an empty list
        let e = G1.scan("e-malformed-count")
        let refused = try PhaseVectorMatcher.fitZoneAxis(
            scan: e.peaks, scanWidth: e.width, scanHeight: e.height, crystal: al,
            referenceSettings: PhaseReferenceSettings(), settings: PhaseVectorSettings(), coordinates: e.coordinates)
        t.eq(refused.isEmpty, true, "malformed scan gives no fits")
        t.printed("3 sweeps")
        #expect(t.failures == 0, t.report)
    }

    @Test("claims: PhaseVectorMatcher.claims for every position of scans a-d (dumped results) and the custom-radius run")
    func claims() {
        var t = Tally("claims")
        func code(_ c: PhaseDiskClaim) -> Int {
            switch c {
            case .notMatched: return -1
            case .matrix: return -2
            case .unexplained: return -3
            case .phase(let p): return p
            }
        }
        for block in objs(G1.dump["claims"]) {
            let runID = str(block["run"])
            let run = G1.runBlock(runID)
            let scan = G1.scan(str(block["scan"]))
            let settings = block["settings"].map { vectorSettings($0 as! JSON) } ?? PhaseVectorSettings()
            let rows = arr(run["results"])
            let matrixIndex = jint(run["matrixEntryIndex"])
            let matrixEntry: PhaseOrientationReference? = matrixIndex >= 0 ? G1.main.entries[matrixIndex] : nil
            let want = arr(block["claims"])
            t.eq(scan.peaks.count, want.count, "\(runID) positions")
            for (i, peaks) in scan.peaks.enumerated() where i < want.count && i < rows.count {
                let w = arr(rows[i])
                var result = PhaseVectorResult()
                result.verdict = PhaseVerdict(rawValue: UInt8(jint(w[0])))!
                result.phaseIndex = Int32(jint(w[1])); result.entryIndex = Int32(jint(w[2]))
                result.score = flt(w[3])
                result.matchedCount = Int32(jint(w[5])); result.survivingCount = Int32(jint(w[6]))
                result.removedCount = Int32(jint(w[7]))
                result.runnerUpPhaseIndex = Int32(jint(w[8])); result.runnerUpScore = flt(w[9])
                let kit = PhaseVectorMatcher.claims(
                    peaks: peaks, coordinates: scan.coordinates, settings: settings, matrixEntry: matrixEntry,
                    result: result, library: G1.main).map(code)
                let viaPixels = PhaseVectorMatcher.claims(
                    peaks: peaks, originX: scan.originX, originY: scan.originY, invAngstromPerPixel: scan.scale,
                    settings: settings, matrixEntry: matrixEntry, result: result, library: G1.main).map(code)
                let expected = arr(want[i]).map { jint($0) }
                t.eq(kit, expected, "\(runID) position \(i) claims")
                t.eq(viaPixels, expected, "\(runID) position \(i) claims (pixel form)")
            }
        }
        t.printed("5 blocks")
        #expect(objs(G1.dump["claims"]).count == 5)
        #expect(t.failures == 0, t.report)
    }

    // MARK: resolution, nearest, score

    @Test("resolution: PhaseVectorResolution and scaledToDetector at 6 scales")
    func resolution() {
        var t = Tally("resolution")
        for row in objs(G1.dump["resolution"]) {
            let scale = dbl(row["invAngstromPerPixel"])
            let r = PhaseVectorResolution(settings: PhaseVectorSettings(), invAngstromPerPixel: scale)
            let p = "scale \(scale)"
            t.d(r.invAngstromPerPixel, row["invAngstromPerPixel"], "\(p) invAngstromPerPixel")
            t.d(r.pairRadiusPixels, row["pairRadiusPixels"], "\(p) pairRadiusPixels")
            t.d(r.matrixRemovalPixels, row["matrixRemovalPixels"], "\(p) matrixRemovalPixels")
            t.d(r.notIndexedAbovePixels, row["notIndexedAbovePixels"], "\(p) notIndexedAbovePixels")
            t.eq(r.advice, isNull(row["advice"]) ? nil : str(row["advice"]), "\(p) advice")
            let s = r.scaledToDetector(PhaseVectorSettings())
            let w = row["scaledToDetector"] as! JSON
            t.d(s.directBeamRadiusInvAngstrom, w["directBeamRadiusInvAngstrom"], "\(p) scaled directBeamRadius")
            t.d(s.maximumVectorInvAngstrom, w["maximumVectorInvAngstrom"], "\(p) scaled maximumVector")
            t.d(s.matrixToleranceInvAngstrom, w["matrixToleranceInvAngstrom"], "\(p) scaled matrixTolerance")
            t.d(s.pairRadiusInvAngstrom, w["pairRadiusInvAngstrom"], "\(p) scaled pairRadius")
            t.int(s.minimumVectors, w["minimumVectors"], "\(p) scaled minimumVectors")
            t.d(s.matrixFallbackExplainedFraction, w["matrixFallbackExplainedFraction"], "\(p) scaled matrixFallback")
            t.int(s.minimumMatchedVectors, w["minimumMatchedVectors"], "\(p) scaled minimumMatchedVectors")
            t.int(s.friedelPairMinimumMatchedVectors, w["friedelPairMinimumMatchedVectors"], "\(p) scaled friedelPairMinimumMatchedVectors")
            t.d(s.chanceMatchMultiple, w["chanceMatchMultiple"], "\(p) scaled chanceMatchMultiple")
            t.d(s.notIndexedAboveInvAngstrom, w["notIndexedAboveInvAngstrom"], "\(p) scaled notIndexedAbove")
            t.d(s.minimumPhaseContrastInvAngstrom, w["minimumPhaseContrastInvAngstrom"], "\(p) scaled minimumPhaseContrast")
            t.d(s.orientationRelationshipToleranceDeg, w["orientationRelationshipToleranceDeg"], "\(p) scaled orientationRelationshipTolerance")
            t.eq(s.completenessAwareCrossPhaseRanking, bool(w["completenessAwareCrossPhaseRanking"]), "\(p) scaled completenessAware")
            t.eq(s.classificationRule.rawValue, str(w["classificationRule"]), "\(p) scaled classificationRule")
            t.int(s.directMatrixMaximumVectors, w["directMatrixMaximumVectors"], "\(p) scaled directMatrixMaximumVectors")
            t.d(s.residualCutoffInvAngstrom, w["residualCutoffInvAngstrom"], "\(p) scaled residualCutoff")
            t.int(s.knownVariantsMinimumSpecificReflections, w["knownVariantsMinimumSpecificReflections"], "\(p) scaled knownVariantsMinimumSpecificReflections")
        }
        t.printed("6 scales")
        #expect(objs(G1.dump["resolution"]).count == 6)
        #expect(t.failures == 0, t.report)
    }

    @Test("nearest, score and containsFriedelPair on 30 library entries")
    func nearestAndScore() {
        var t = Tally("nearest")
        let scratch = PhaseVectorMatcher.Scratch(capacity: 64)
        for row in objs(G1.dump["nearest"]) {
            let e = jint(row["entryIndex"])
            let entry = G1.main.entries[e]
            let us = arr(row["us"]).map { v -> SIMD2<Double> in let a = arr(v); return SIMD2(dbl(a[0]), dbl(a[1])) }
            for (key, radius) in [("nearest0p02", 0.02), ("nearest0p05", 0.05)] {
                let want = arr(row[key])
                for (k, u) in us.enumerated() {
                    let p = "entry \(e) \(key) vector \(k)"
                    let hit = PhaseVectorMatcher.nearest(u, in: entry.vectors, radius: radius)
                    if isNull(want[k]) { t.eq(hit == nil, true, "\(p) nil"); continue }
                    guard let hit else { t.note("\(p): kit nil, dump has a hit"); continue }
                    let w = arr(want[k])
                    t.int(hit.index, w[0], "\(p) index")
                    t.d(hit.distance, w[1], "\(p) distance")
                }
            }
            let sc = PhaseVectorMatcher.score(vectors: us, against: entry, pairRadius: 0.02, scratch: scratch)
            if isNull(row["score0p02"]) {
                t.eq(sc == nil, true, "entry \(e) score nil")
            } else if let sc {
                let w = row["score0p02"] as! JSON
                t.d(sc.score, w["score"], "entry \(e) score")
                t.int(sc.matched, w["matched"], "entry \(e) matched")
                t.int(sc.uniqueReferences, w["uniqueReferences"], "entry \(e) uniqueReferences")
            } else {
                t.note("entry \(e): kit score nil, dump has a score")
            }
            t.eq(PhaseVectorMatcher.containsFriedelPair(us, radius: 0.02), bool(row["friedel0p02"]), "entry \(e) friedel")
        }
        t.printed("30 entries")
        #expect(objs(G1.dump["nearest"]).count == 30)
        #expect(t.failures == 0, t.report)
    }

    // MARK: Euler and the detector frame

    @Test("euler: 270 cases, the angles, crystalToLabMatrix, the py4DSTEM matrix, degrees and the round trip")
    func euler() {
        var t = Tally("euler")
        let cases = objs(G1.dump["euler"])
        for c in cases {
            let label = str(c["label"])
            let e = EulerAngles(py4DSTEMOrientationMatrix: matrix3(c["matrix"]))
            t.f(e.phi1, c["phi1"], "\(label) phi1"); t.f(e.Phi, c["Phi"], "\(label) Phi"); t.f(e.phi2, c["phi2"], "\(label) phi2")
            let lab = flat(e.crystalToLabMatrix), labWant = matrixBits(c["crystalToLab"])
            for i in 0..<9 { t.bits(lab[i].bitPattern, labWant[i], "\(label) crystalToLab[\(i)]") }
            let py = flat(e.py4DSTEMOrientationMatrix), pyWant = matrixBits(c["py4DSTEMMatrix"])
            for i in 0..<9 { t.bits(py[i].bitPattern, pyWant[i], "\(label) py4DSTEMMatrix[\(i)]") }
            let deg = e.degrees, degWant = arr(c["degrees"])
            t.f(deg.0, degWant[0], "\(label) degrees.0"); t.f(deg.1, degWant[1], "\(label) degrees.1"); t.f(deg.2, degWant[2], "\(label) degrees.2")
            let again = EulerAngles(py4DSTEMOrientationMatrix: e.py4DSTEMOrientationMatrix)
            let trip = arr(c["roundTrip"])
            t.f(again.phi1, trip[0], "\(label) roundTrip phi1"); t.f(again.Phi, trip[1], "\(label) roundTrip Phi")
            t.f(again.phi2, trip[2], "\(label) roundTrip phi2")
        }
        t.printed("\(cases.count) cases")
        #expect(cases.count == 270)
        #expect(t.failures == 0, t.report)
    }

    @Test("detectorBases: ZoneAxisFrame.detectorBasis and matrix(detectorBasis:inPlaneAngle:) on every dumped axis, three crystals",
          arguments: ["Al", "triclinic", "beta-double-prime"])
    func detectorBases(_ crystalName: String) {
        var t = Tally("detectorBases")
        let block = objs(G1.dump["detectorBases"]).first { str($0["crystal"]) == crystalName }!
        let axes = objs(block["axes"])
        for a in axes {
            let axis = ints3(a["axis"])
            let p = "\(crystalName) [\(axis.x)\(axis.y)\(axis.z)]"
            let nw = arr(a["n"])
            let n = SIMD3(dbl(nw[0]), dbl(nw[1]), dbl(nw[2]))
            let basis = ZoneAxisFrame.detectorBasis(zoneAxis: n)
            let basisKit = flat(basis), basisWant = matrixBits(a["basis"])
            for i in 0..<9 { t.bits(basisKit[i].bitPattern, basisWant[i], "\(p) basis[\(i)]") }
            for m in objs(a["matrices"]) {
                let angle = flt(m["angle"])
                let kit = flat(ZoneAxisFrame.matrix(detectorBasis: basis, inPlaneAngle: angle))
                let want = matrixBits(m["matrix"])
                for i in 0..<9 { t.bits(kit[i].bitPattern, want[i], "\(p) matrix(angle \(angle))[\(i)]") }
            }
        }
        t.printed("\(crystalName), \(axes.count) axes")
        #expect(axes.count == 49)
        #expect(t.failures == 0, t.report)
    }

    @Test("detectorBases: cartesianZoneAxis from the kit's own lattice equals the apps' unit vector for every dumped axis",
          arguments: ["Al", "triclinic", "beta-double-prime"])
    func cartesianZoneAxis(_ crystalName: String) {
        var t = Tally("cartesianAxis")
        let s = G1.structure(named: crystalName)
        let block = objs(G1.dump["detectorBases"]).first { str($0["crystal"]) == crystalName }!
        for a in objs(block["axes"]) {
            let axis = ints3(a["axis"])
            let p = "\(crystalName) [\(axis.x)\(axis.y)\(axis.z)]"
            let nw = arr(a["n"])
            guard let n = PhaseReferenceLibrary.cartesianZoneAxis(axis, crystal: s) else { t.note("\(p): nil"); continue }
            t.d(n.x, nw[0], "\(p) n.x"); t.d(n.y, nw[1], "\(p) n.y"); t.d(n.z, nw[2], "\(p) n.z")
            let basis = flat(ZoneAxisFrame.detectorBasis(zoneAxis: n)), want = matrixBits(a["basis"])
            for i in 0..<9 { t.bits(basis[i].bitPattern, want[i], "\(p) basis from kit axis [\(i)]") }
        }
        t.printed(crystalName)
        #expect(t.failures == 0, t.report)
    }
}
