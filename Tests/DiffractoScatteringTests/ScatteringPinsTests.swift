//
//  ScatteringPinsTests.swift
//  The kit against the pins written BEFORE its code (registration docs/registration-slice-2.md §4.1–4.4 and amendments 1–5):
//    slice-2-pins.json  tools/slice-2/pins.py: NumPy from the definitions, pymatgen's Laue group and X-ray pattern
//    slice-2-apps.json  tools/slice-2/apps-parity: mac4DSTEM 7d6fdbd1's Crystal.reflections, compiled -O
//  Public API only (DiffractoScattering and DiffractoCore). A failure is a failure of the port, never a tolerance to widen.
//

import CryptoKit
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

// MARK: - The pins, decoded

private struct PinCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct PinSite: Decodable { let z: Int; let frac: [Double]; let occupancy: Double }
private struct PinBuild: Decodable { let how: String; let cifFile: String?; let cifSha256: String?; let cifText: String? }
private struct PinFamily: Decodable {
    let g: Double
    let members: [[Int]]
    let multiplicity: Int
    let orbits: [[Int]]
    let name: String
    let sumF2: Double
}
private struct PinReflections: Decodable { let items: [[Double]] }    // h, k, l, |g|, Re F, Im F, |F| for |F| > 1e-6
private struct PinExtra: Decodable { let keptCount: Int; let rings: Int }
private struct PinStructure: Decodable {
    let name: String
    let build: PinBuild
    let cell: PinCell
    let sites: [PinSite]
    let reflections: PinReflections
    let keptCount: Int
    let families: [PinFamily]
    let extraKMax: [String: PinExtra]?
}
private struct PinCoefficients: Decodable { let z: Int; let coefficients: [[Double]] }
private struct PinPeak: Decodable { let members: [[Int]]; let sumF2Pymatgen: Double }
private struct PinXray: Decodable {
    let name: String
    let elementCoefficients: [String: PinCoefficients]
    let peaks: [PinPeak]
}
private struct PinCutEdge: Decodable { let cell: [Double]; let volume: Double; let reflectionsWithinKMax: Int }
private struct PinFE: Decodable { let g2: [Double]; let rows: [[Double]] }
private struct Pins: Decodable {
    let tableSha256: String
    let table: [[Double]]
    let feSamples: PinFE
    let structures: [PinStructure]
    let xray: [PinXray]
    let cutEdge: PinCutEdge
}

private struct AppsCell: Decodable { let a, b, c, alpha, beta, gamma: Double }
private struct AppsReflection: Decodable { let h, k, l: Int; let g, re, im, i: String }
private struct AppsCase: Decodable {
    let preset: String
    let kMax: Double
    let count: Int
    let cell: AppsCell
    let volumeBits: String
    let sites: [[Double]]          // z, x, y, z, occupancy
    let reflections: [AppsReflection]
}
private struct Apps: Decodable { let cases: [AppsCase] }

enum ScatteringFixtures {
    static let directory: URL = URL(fileURLWithPath: #filePath)
        .deletingLastPathComponent().appendingPathComponent("Fixtures")
    /// Slice 1's CIF fixtures (read by file name, as the pins name them).
    static let coreDirectory: URL = URL(fileURLWithPath: #filePath)
        .deletingLastPathComponent().deletingLastPathComponent()
        .appendingPathComponent("DiffractoCoreTests").appendingPathComponent("Fixtures")

    fileprivate static func decode<T: Decodable>(_ type: T.Type, _ name: String) throws -> T {
        let data = try Data(contentsOf: directory.appendingPathComponent(name))
        let d = JSONDecoder()
        d.keyDecodingStrategy = .convertFromSnakeCase
        return try d.decode(type, from: data)
    }

    fileprivate static let pinsResult: Result<Pins, any Error> = Result { try decode(Pins.self, "slice-2-pins.json") }
    fileprivate static let appsResult: Result<Apps, any Error> = Result { try decode(Apps.self, "slice-2-apps.json") }
    fileprivate static func pins() throws -> Pins { try pinsResult.get() }
    fileprivate static func apps() throws -> Apps { try appsResult.get() }

    fileprivate static func pin(_ name: String) throws -> PinStructure {
        try #require(try pins().structures.first { $0.name == name }, "no pinned structure \(name)")
    }

    fileprivate static func cifText(_ p: PinStructure) throws -> String {
        if let t = p.build.cifText { return t }
        let file = try #require(p.build.cifFile)
        return try String(contentsOf: coreDirectory.appendingPathComponent(file), encoding: .utf8)
    }

    /// Presets: the pinned cell and the pinned sites in the pinned order, symmetry detected. CIF structures: CIFImport of
    /// the fixture's text (β″: the generator's text carried by the pins), which keeps the coordinates as written.
    fileprivate static func build(_ p: PinStructure) throws -> CrystalStructure {
        switch p.build.how {
        case "preset":
            return CrystalStructure.detectingSymmetry(cell: cell(p.cell), sites: p.sites.map(site))
        case "cif":
            return try CIFImport.read(try cifText(p), fileBaseName: p.name).structure
        default:
            Issue.record("unknown build \(p.build.how)")
            throw PinsError.unknownBuild
        }
    }
    fileprivate static func cell(_ c: PinCell) -> UnitCell {
        UnitCell(a: c.a, b: c.b, c: c.c, alphaDeg: c.alpha, betaDeg: c.beta, gammaDeg: c.gamma)
    }
    fileprivate static func site(_ s: PinSite) -> AtomSite {
        AtomSite(z: s.z, fractional: SIMD3(s.frac[0], s.frac[1], s.frac[2]), occupancy: s.occupancy)
    }
    enum PinsError: Error { case unknownBuild }

    static let structureNames = [
        "Au", "Al", "Si", "Mg", "WS2", "Cu3Au-occ", "Weak-cut", "GaAs-zb", "ideal_hcp", "beta_double_prime",
        "gold_P1", "silicon_Fd-3m", "magnesium_P63mmc", "rocksalt_Fm-3m_192ops", "alpha_quartz_P3221",
        "calcite_R-3c_hex", "fluorapatite_P63m", "pyrite_Pa-3", "ws2_2H", "cdi2_P-3m1",
    ]
    static let xrayNames = ["Au", "Si", "Mg", "beta_double_prime", "Cu3Au-occ", "WS2"]
    static let appsCaseIndices = Array(0..<10)
}

private func hex(_ s: String) -> Double { Double(bitPattern: UInt64(s, radix: 16)!) }

private func relErr(_ x: Double, _ y: Double, floor: Double) -> Double { abs(x - y) / max(abs(y), floor) }

private func descendingMembers(_ ms: [MillerIndex]) -> [[Int]] {
    ms.sorted { ($0.h, $0.k, $0.l) > ($1.h, $1.k, $1.l) }.map { [$0.h, $0.k, $0.l] }
}

private func key(_ h: Int, _ k: Int, _ l: Int) -> MillerIndex { MillerIndex(h, k, l) }

struct ScatteringPinsTests {

    // MARK: 4.2 the table

    /// All 103 rows × 10 parameters equal to the fetched file bit for bit; f_e at g² ∈ {0, 0.25, 1, 4} for every Z within
    /// 1e-15 relative of NumPy's value. (Also holds the JSON decoder to exact doubles: a lossy decode fails here first.)
    @Test func testTableBitForBitAndSamples() throws {
        let pins = try ScatteringFixtures.pins()
        #expect(pins.tableSha256 == "be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4")
        #expect(pins.table.count == 103 && ElectronScatteringFactors.parameters.count == 103)
        var bad: [String] = []
        for z in 1...103 {
            let want = pins.table[z - 1], got = ElectronScatteringFactors.parameters[z - 1]
            if want.count != 10 || got.count != 10 { bad.append("Z=\(z) width \(got.count)"); continue }
            for j in 0..<10 where want[j].bitPattern != got[j].bitPattern { bad.append("Z=\(z) p\(j) \(got[j]) vs \(want[j])") }
        }
        #expect(bad.isEmpty, "table differs from the file: \(bad.prefix(5))")
        var worst = 0.0
        var worstAt = ""
        for z in 1...103 {
            for (i, g2) in pins.feSamples.g2.enumerated() {
                let want = pins.feSamples.rows[z - 1][i]
                let got = try #require(ElectronScatteringFactors.f(z: z, gSquared: g2))
                let e = abs(got - want) / abs(want)
                if e > worst { worst = e; worstAt = "Z=\(z) g2=\(g2)" }
            }
        }
        #expect(worst <= 1e-15, "f_e differs from NumPy by \(worst) at \(worstAt)")
        #expect(ElectronScatteringFactors.f(z: 104, gSquared: 0) == nil && ElectronScatteringFactors.f(z: 0, gSquared: 0) == nil)
    }

    /// The CIF texts the pins were made from are the committed fixtures (sha256), so the pins describe these files.
    @Test func testFixtureTextsMatchThePinnedSha256() throws {
        let pins = try ScatteringFixtures.pins()
        #expect(Set(pins.structures.map(\.name)) == Set(ScatteringFixtures.structureNames))
        for p in pins.structures where p.build.how == "cif" {
            let text = try ScatteringFixtures.cifText(p)
            let sha = SHA256.hash(data: Data(text.utf8)).map { String(format: "%02x", $0) }.joined()
            #expect(sha == p.build.cifSha256, "\(p.name): fixture text sha256 \(sha)")
        }
    }

    // MARK: 4.1 NumPy

    /// Every structure at k_max 1.2: Re F, Im F, |F| within 1e-9 relative of max(|F|, 1e-4) on every pinned reflection; the
    /// kept set equal (a reflection within 1e-12 relative of the cut may fall either side); the unlisted ones (|F| ≤ 1e-6
    /// in NumPy) stay below 1e-6; the families equal (count, g 1e-12, members, multiplicity, orbits, name, Σ|F|² 1e-9).
    @Test(arguments: ScatteringFixtures.structureNames)
    func testStructureMatchesNumPyPins(_ name: String) throws {
        let p = try ScatteringFixtures.pin(name)
        let s = try ScatteringFixtures.build(p)
        let list = try s.kinematicReflections(kMax: 1.2)
        var kit: [MillerIndex: KinematicReflection] = [:]
        for r in list.all { kit[r.reflection.hkl] = r }
        var fails: [String] = []

        var pinned = Set<MillerIndex>()
        for it in p.reflections.items where it[3] <= 1.2 {
            let hkl = key(Int(it[0]), Int(it[1]), Int(it[2]))
            pinned.insert(hkl)
            guard let k = kit[hkl] else { fails.append("missing \(hkl)"); continue }
            let f = k.structureFactor
            let floor = max(it[6], 1e-4)
            if abs(f.gLength - it[3]) > 1e-12 * it[3] { fails.append("\(hkl) |g| \(f.gLength) vs \(it[3])") }
            if abs(f.real - it[4]) > 1e-9 * floor { fails.append("\(hkl) Re \(f.real) vs \(it[4])") }
            if abs(f.imaginary - it[5]) > 1e-9 * floor { fails.append("\(hkl) Im \(f.imaginary) vs \(it[5])") }
            if abs(f.magnitude - it[6]) > 1e-9 * floor { fails.append("\(hkl) |F| \(f.magnitude) vs \(it[6])") }
            if it[3] <= 1.2, abs(it[6] - 1e-4) > 1e-12 * 1e-4 {
                let pinKept = it[6] > 1e-4
                let kitKept = f.magnitude > 1e-4
                if pinKept != kitKept { fails.append("\(hkl) cut side: |F| \(f.magnitude) vs \(it[6])") }
            }
        }
        for (hkl, k) in kit where !pinned.contains(hkl) && k.structureFactor.magnitude > 1e-6 + 1e-13 {
            fails.append("unpinned \(hkl) has |F| \(k.structureFactor.magnitude) > 1e-6")
        }
        #expect(list.all.count == kit.count)

        // families
        if list.families.count != p.families.count {
            fails.append("family count \(list.families.count) vs \(p.families.count)")
        } else {
            for (i, (a, b)) in zip(list.families, p.families).enumerated() {
                let fam = a.family
                var why: [String] = []
                if abs(fam.g - b.g) > 1e-12 * b.g { why.append("g \(fam.g) vs \(b.g)") }
                if descendingMembers(fam.members) != b.members { why.append("members") }
                if fam.multiplicity != b.multiplicity { why.append("multiplicity \(fam.multiplicity) vs \(b.multiplicity)") }
                if fam.orbits.map({ [$0.h, $0.k, $0.l] }) != b.orbits { why.append("orbits \(fam.orbits) vs \(b.orbits)") }
                if fam.name != b.name { why.append("name \(fam.name) vs \(b.name)") }
                if abs(a.intensity - b.sumF2) > 1e-9 * b.sumF2 { why.append("sumF2 \(a.intensity) vs \(b.sumF2)") }
                if !why.isEmpty { fails.append("family \(i) \(b.name): " + why.joined(separator: "; ")) }
            }
        }
        let nearCut = p.reflections.items.contains { abs($0[6] - 1e-4) <= 1e-12 * 1e-4 }
        if !nearCut { #expect(list.kept.count == p.keptCount, "\(name): kept \(list.kept.count) vs pinned \(p.keptCount)") }
        #expect(fails.isEmpty, "\(name): \(fails.count) differences, first: \(fails.prefix(6))")
    }

    /// β″ at k_max 1.6 and 2.0: kept count and ring count equal the pins (953 / 1805 rings, 3420 / 6596 reflections).
    @Test func testBetaDoublePrimeLargerKMax() throws {
        let p = try ScatteringFixtures.pin("beta_double_prime")
        let s = try ScatteringFixtures.build(p)
        let extra = try #require(p.extraKMax)
        #expect(Set(extra.keys) == ["1.6", "2.0"])
        for (k, want) in extra {
            let list = try s.kinematicReflections(kMax: Double(k)!)
            #expect(list.kept.count == want.keptCount, "kMax \(k): kept \(list.kept.count) vs \(want.keptCount)")
            #expect(list.families.count == want.rings, "kMax \(k): rings \(list.families.count) vs \(want.rings)")
        }
    }

    // MARK: amendment 4 the cut's strictness

    /// a = b = 100, c = 1 Å, a unit form factor: every reflection has |F| = 1/V = 1e-4, so `>` keeps none (`>=` would keep
    /// all of them).
    @Test func testCutEdgeKeepsNothing() throws {
        let ce = try ScatteringFixtures.pins().cutEdge
        let c = ce.cell
        let s = CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: c[0], b: c[1], c: c[2], alphaDeg: c[3], betaDeg: c[4], gammaDeg: c[5]),
            sites: [AtomSite(z: 1, fractional: [0, 0, 0])])
        #expect(s.cell.volume == ce.volume)
        let list = try s.kinematicReflections(kMax: 1.2, formFactor: { _, _ in 1 })
        // Reflections with |g| = 1.2 exactly in real arithmetic ((h, k, 0), h² + k² = 14400, 1/100 Å⁻¹ steps) are decided by
        // rounding: NumPy's norm and the kit's length may put one on either side. Everything else must agree.
        var onTheEdge = 0
        for h in -120...120 { for k in -120...120 where h * h + k * k == 14400 { onTheEdge += 1 } }
        #expect(abs(list.all.count - ce.reflectionsWithinKMax) <= onTheEdge,
                "kit \(list.all.count) vs pinned \(ce.reflectionsWithinKMax) reflections, \(onTheEdge) lie on |g| = 1.2")
        #expect(list.kept.isEmpty, "\(list.kept.count) kept at |F| = 1e-4")
        let worst = list.all.map { abs($0.structureFactor.magnitude - 1e-4) }.max() ?? 1
        #expect(worst <= 1e-16, "worst | |F| - 1e-4 | = \(worst)")
        // and one ulp-scale step below the edge keeps every one
        let below = try s.kinematicReflections(kMax: 1.2, cut: 1e-4 * (1 - 1e-9), formFactor: { _, _ in 1 })
        #expect(below.kept.count == list.all.count)
        #expect(list.all.count >= ce.reflectionsWithinKMax)
    }

    // MARK: amendment 2 the sign of Im F

    /// Zinc-blende GaAs, non-centrosymmetric about the origin: F = 4 (f_Ga + i^(h+k+l) f_As) / V for unmixed h, k, l,
    /// absent otherwise (1e-13 relative); Im F(1 1 1) < 0 and Im F(−1 −1 −1) > 0. A flipped phase sign fails both.
    @Test func testZincBlendeImaginaryPartSign() throws {
        let s = try ScatteringFixtures.build(try ScatteringFixtures.pin("GaAs-zb"))
        let list = try s.kinematicReflections(kMax: 1.2)
        var fails: [String] = []
        for r in list.all {
            let h = r.reflection.hkl, f = r.structureFactor
            let g2 = f.gLength * f.gLength
            let ga = try #require(ElectronScatteringFactors.f(z: 31, gSquared: g2))
            let asf = try #require(ElectronScatteringFactors.f(z: 33, gSquared: g2))
            let parity = [h.h, h.k, h.l].map { (($0 % 2) + 2) % 2 }.reduce(0, +)
            if parity == 1 || parity == 2 {
                if f.magnitude > 1e-12 { fails.append("\(h) mixed parity |F| \(f.magnitude)") }
                continue
            }
            let n = (((h.h + h.k + h.l) % 4) + 4) % 4
            let (ir, ii): (Double, Double) = [(1, 0), (0, 1), (-1, 0), (0, -1)][n]
            let wantRe = 4 * (ga + ir * asf) / s.cell.volume
            let wantIm = 4 * (ii * asf) / s.cell.volume
            let scale = 4 * (ga + asf) / s.cell.volume
            if abs(f.real - wantRe) > 1e-13 * scale || abs(f.imaginary - wantIm) > 1e-13 * scale {
                fails.append("\(h) F (\(f.real), \(f.imaginary)) vs (\(wantRe), \(wantIm))")
            }
        }
        #expect(fails.isEmpty, "GaAs: \(fails.prefix(5))")
        let plus = try #require(StructureFactors.compute(MillerIndex(1, 1, 1), in: s))
        let minus = try #require(StructureFactors.compute(MillerIndex(-1, -1, -1), in: s))
        #expect(plus.imaginary < 0 && minus.imaginary > 0, "Im F(111) \(plus.imaginary), Im F(-1-1-1) \(minus.imaginary)")
    }

    // MARK: 4.3 pymatgen X-ray

    /// For each pymatgen peak: Σ|F_X|² over its half-space members through the kit's phase sum with pymatgen's X-ray form
    /// factor f(s) = Z − 41.78214 s² Σ a exp(−b s²), s = |g|/2, equals pymatgen's intensity / (2 LP) within 1e-6 relative.
    @Test(arguments: ScatteringFixtures.xrayNames)
    func testXRayPhaseSumsMatchPymatgen(_ name: String) throws {
        let pins = try ScatteringFixtures.pins()
        let block = try #require(pins.xray.first { $0.name == name })
        let s = try ScatteringFixtures.build(try ScatteringFixtures.pin(name))
        var table: [Int: [[Double]]] = [:]
        for (_, e) in block.elementCoefficients { table[e.z] = e.coefficients }
        let xray: FormFactor = { z, g2 in
            guard let coeff = table[z] else { return nil }
            let s2 = g2 / 4
            var sum = 0.0
            for c in coeff { sum += c[0] * exp(-c[1] * s2) }
            return Double(z) - 41.78214 * s2 * sum
        }
        #expect(!block.peaks.isEmpty)
        var worst = 0.0
        var worstAt = ""
        for peak in block.peaks {
            var total = 0.0
            for m in peak.members {
                let f = try #require(StructureFactors.compute(MillerIndex(m[0], m[1], m[2]), in: s, formFactor: xray))
                total += f.intensity * s.cell.volume * s.cell.volume     // the kit's F carries the 1/V; pymatgen's does not
            }
            let e = abs(total - peak.sumF2Pymatgen) / peak.sumF2Pymatgen
            if e > worst { worst = e; worstAt = "\(peak.members.first ?? [])" }
        }
        #expect(worst <= 1e-6, "\(name): worst relative difference \(worst) at \(worstAt)")
    }

    // MARK: 4.4 the apps, bit for bit (amendment 1)

    /// The apps' `Crystal.reflections(kMax:tolerance: 1e-4)` (mac4DSTEM 7d6fdbd1, -O): the kept sets equal, and per kept
    /// reflection the Double bit patterns of |g|, Re F and |F|² equal; Im F bit-equal where |Im F| ≥ 1e-15 and within 1e-15
    /// below that (amendment 1). β″ is built from the CIF here, so its site order differs from the apps' hand list and a
    /// summation-order difference is expected: for it, Re F and Im F are held to 1e-12 relative of max(|F|, 1e-4) (Im F also
    /// passes below 1e-15 within 1e-15), |F|² to 2e-12 relative, |g| bit for bit (registration 4.4 as amended 2026-10-15: the measured worst is 2.7e-13; bit identity is the apps-site-list test).
    @Test(arguments: ScatteringFixtures.appsCaseIndices)
    func testAppsBitForBit(_ index: Int) throws {
        let c = try ScatteringFixtures.apps().cases[index]
        let isBeta = c.preset == "betaDoublePrime"
        let s = isBeta ? try ScatteringFixtures.build(try ScatteringFixtures.pin("beta_double_prime")) : appsStructure(c)
        let (fails, stats) = try compareToApps(c, s, summationOrderTolerance: isBeta)
        #expect(fails.isEmpty, "\(c.preset)@\(c.kMax): \(stats); first: \(fails.prefix(4))")
    }

    /// The cause of the β″ differences above: built from the apps' own site list (same order, same coordinates) the kit is
    /// bit identical to the apps like every other preset.
    @Test(arguments: [8, 9])
    func testAppsBetaDoublePrimeFromTheAppsSiteListIsBitIdentical(_ index: Int) throws {
        let c = try ScatteringFixtures.apps().cases[index]
        #expect(c.preset == "betaDoublePrime")
        let (fails, _) = try compareToApps(c, appsStructure(c), summationOrderTolerance: false)
        #expect(fails.isEmpty, "\(c.preset)@\(c.kMax): \(fails.count) differences, first: \(fails.prefix(6))")
    }
}

private func appsStructure(_ c: AppsCase) -> CrystalStructure {
    let cell = UnitCell(a: c.cell.a, b: c.cell.b, c: c.cell.c, alphaDeg: c.cell.alpha, betaDeg: c.cell.beta,
                        gammaDeg: c.cell.gamma)
    let sites = c.sites.map { AtomSite(z: Int($0[0]), fractional: SIMD3($0[1], $0[2], $0[3]), occupancy: $0[4]) }
    return CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
}

private func compareToApps(_ c: AppsCase, _ s: CrystalStructure, summationOrderTolerance: Bool) throws
    -> (fails: [String], stats: String) {
    var fails: [String] = []
    if s.cell.volume.bitPattern != UInt64(c.volumeBits, radix: 16) { fails.append("volume \(s.cell.volume)") }
    let list = try s.kinematicReflections(kMax: c.kMax)
    var kept: [MillerIndex: StructureFactor] = [:]
    for k in list.kept { kept[k.reflection.hkl] = k.structureFactor }
    var worstRe = 0.0, worstIm = 0.0, worstI = 0.0, worstG = 0.0
    if kept.count != c.count { fails.append("kept \(kept.count) vs the apps' \(c.count)") }
    var appsSet = Set<MillerIndex>()
    for r in c.reflections {
        let hkl = MillerIndex(r.h, r.k, r.l)
        appsSet.insert(hkl)
        guard let f = kept[hkl] else { fails.append("\(hkl) not kept by the kit"); continue }
        let g = hex(r.g), re = hex(r.re), im = hex(r.im), i = hex(r.i)
        if f.gLength.bitPattern != g.bitPattern { fails.append("\(hkl) |g| \(f.gLength) vs \(g)") }
        if summationOrderTolerance {
            let floor = max(i.squareRoot(), 1e-4)
            worstRe = max(worstRe, abs(f.real - re) / floor)
            if abs(im) >= 1e-15 { worstIm = max(worstIm, abs(f.imaginary - im) / floor) }
            worstI = max(worstI, abs(f.intensity - i) / i)
            worstG = max(worstG, abs(f.gLength - g) / g)
            if abs(f.real - re) > 1e-12 * floor { fails.append("\(hkl) Re \(f.real) vs \(re) (rel \(abs(f.real - re) / floor))") }
            let imOK = abs(f.imaginary - im) <= 1e-12 * floor || (abs(im) < 1e-15 && abs(f.imaginary - im) <= 1e-15)
            if !imOK { fails.append("\(hkl) Im \(f.imaginary) vs \(im) (rel \(abs(f.imaginary - im) / floor))") }
            if abs(f.intensity - i) > 2e-12 * i { fails.append("\(hkl) |F|² \(f.intensity) vs \(i)") }
        } else {
            if f.real.bitPattern != re.bitPattern { fails.append("\(hkl) Re \(f.real) vs \(re)") }
            if f.intensity.bitPattern != i.bitPattern { fails.append("\(hkl) |F|² \(f.intensity) vs \(i)") }
            if abs(im) >= 1e-15 {
                if f.imaginary.bitPattern != im.bitPattern { fails.append("\(hkl) Im \(f.imaginary) vs \(im)") }
            } else if abs(f.imaginary - im) > 1e-15 {
                fails.append("\(hkl) Im (noise) \(f.imaginary) vs \(im)")
            }
        }
    }
    for hkl in kept.keys where !appsSet.contains(hkl) { fails.append("\(hkl) kept by the kit, not by the apps") }
    let stats = "worst over \(c.count): Re/floor \(worstRe), Im/floor (|Im| >= 1e-15) \(worstIm), |F|² relative \(worstI), |g| relative \(worstG); failing checks \(fails.count)"
    return (fails, stats)
}
