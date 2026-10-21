//
//  FastDetectionTests.swift
//  Slice 1b, the Cartesian tolerance and the cell-list permutation test (registration
//  docs/registration-slice-1b.md sections 1 and 3): the cell list returns slice 1's operations bit for bit on every
//  slice-1 and slice-2 fixture, a 500-atom symmetric P1 cell is detected in under 2 s, and `.cartesian` tells apart a
//  pseudo-symmetry that `.fractional` cannot in a large cell.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct FastDetectionTests {

    /// slice-2-pins.json's preset structures (cell and sites), the slice-2 fixtures.
    private struct Slice2: Decodable {
        struct Build: Decodable { let how: String }
        struct S: Decodable { let name: String; let build: Build; let cell: CellParams; let sites: [SitePin] }
        let structures: [S]
    }

    private static let slice2Presets: [(String, UnitCell, [AtomSite])] = {
        let url = Fixtures.directory.deletingLastPathComponent().deletingLastPathComponent()
            .appendingPathComponent("DiffractoScatteringTests/Fixtures/slice-2-pins.json")
        guard let data = try? Data(contentsOf: url) else { fatalError("slice-2-pins.json not found at \(url.path)") }
        let decoder = JSONDecoder()
        decoder.keyDecodingStrategy = .convertFromSnakeCase
        guard let parsed = try? decoder.decode(Slice2.self, from: data) else { fatalError("slice-2-pins.json does not decode") }
        return parsed.structures.filter { $0.build.how == "preset" }.map { s in
            (s.name, s.cell.unitCell, s.sites.map {
                AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]), occupancy: $0.occupancy)
            })
        }
    }()

    /// Every fixture the kit has as (name, cell, sites): slice 1b's pins (slice-1 fixtures and the named cells), slice 2's
    /// presets, and the CIFs of the fixtures directory.
    private static var everyFixture: [(String, UnitCell, [AtomSite])] {
        var out = Pins1b.shared.structures.map { ($0.name, $0.cell.unitCell, $0.atomSites) }
        out += slice2Presets
        // near-coincident rows (1.3e-3 apart, both inside the tolerance) and their inversion about a general point: the
        // target of an image is ambiguous, so the choice of slice 1 (the lowest unclaimed index) decides the translation
        let u: [SIMD3<Double>] = [[0.10, 0.20, 0.30], [0.1013, 0.20, 0.30], [0.62, 0.31, 0.12], [0.62, 0.3113, 0.12]]
        let centre = SIMD3<Double>(0.35, 0.6, 0.8)
        let cluster = (u + u.map { centre - $0 }).map { AtomSite(z: 6, fractional: wrapFractional($0)) }
        out.append(("clustered rows", .triclinic(a: 7, b: 8, c: 9, alphaDeg: 80, betaDeg: 85, gammaDeg: 95), cluster))
        let files = (try? FileManager.default.contentsOfDirectory(atPath: Fixtures.directory.path)) ?? []
        for file in files.sorted() where file.hasSuffix(".cif") {
            if let text = try? Fixtures.text(file), let doc = try? CIFImport.read(text, fileBaseName: file) {
                out.append((file, doc.structure.cell, doc.structure.sites))
            }
        }
        return out
    }

    @Test func testCellListEqualsTheReferenceScanBitForBitOnEveryFixture() {
        let fixtures = Self.everyFixture
        #expect(fixtures.count >= 60, "\(fixtures.count) fixtures")
        for (name, cell, sites) in fixtures {
            let fast = SymmetryDetection.detect(sites: sites, cell: cell, search: .cellList)
            let slow = SymmetryDetection.detect(sites: sites, cell: cell, search: .bruteForce)
            #expect(fast == slow, "\(name): cell list \(fast.operations.count) operations, reference \(slow.operations.count)")
            // bit for bit: the translations are the same doubles, not merely close
            for (x, y) in zip(fast.operations, slow.operations) {
                #expect(x.rotation == y.rotation && x.translation == y.translation, "\(name): \(x) vs \(y)")
            }
            // and the whole structure, both settings
            let a = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, search: .cellList)
            let b = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, search: .bruteForce)
            #expect(a == b, "\(name): structure")
        }
    }

    @Test func testCellListEqualsTheReferenceOnTranslationRichCellsAndInCartesianMode() {
        // a 3×3×3 supercell of gold: 108 atoms, 27 × 4 translations × 48 rotations
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        var sites: [AtomSite] = []
        for i in 0..<3 { for j in 0..<3 { for k in 0..<3 { for p in fcc {
            sites.append(AtomSite(z: 79, fractional: (p + SIMD3(Double(i), Double(j), Double(k))) / 3))
        } } } }
        let cell = UnitCell.cubic(a: 3 * 4.0782)
        for tolerance in [SymmetryDetection.Tolerance.default, .cartesian(1e-3), .cartesian(0.05)] {
            let fast = SymmetryDetection.detect(sites: sites, cell: cell, tolerance: tolerance, search: .cellList)
            let slow = SymmetryDetection.detect(sites: sites, cell: cell, tolerance: tolerance, search: .bruteForce)
            #expect(fast == slow, "\(tolerance)")
            #expect(fast.operations.count == 27 * 4 * 48, "\(tolerance): \(fast.operations.count)")
        }
    }

    /// A body-centring that holds to 2e-3 on the first leg and 4e-3 on the second (the translation is applied twice), inside
    /// the fractional 5e-3 and outside its half: both searches find it, so a cell list whose window is narrower than the
    /// tolerance cannot pass; a tolerance of 3e-3 does not find it.
    @Test func testAnOffsetJustInsideTheToleranceIsFoundByBothSearches() {
        let sites = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.5 + 2e-3])]
        let cell = UnitCell.cubic(a: 5)
        let fast = SymmetryDetection.detect(sites: sites, cell: cell, search: .cellList)
        let slow = SymmetryDetection.detect(sites: sites, cell: cell, search: .bruteForce)
        #expect(fast == slow)
        #expect(fast.operations.contains { $0.rotation.isIdentity && !$0.isIdentity }, "the pseudo-centring is a symmetry at 5e-3")
        // the same, with the target across a bin boundary of the cell list (two sites: bins at 0 and 0.5 along each axis)
        let across = [AtomSite(z: 26, fractional: [0, 0, 0.4995]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.0015])]
        let f2 = SymmetryDetection.detect(sites: across, cell: cell, search: .cellList)
        let s2 = SymmetryDetection.detect(sites: across, cell: cell, search: .bruteForce)
        #expect(f2 == s2)
        #expect(s2.operations.contains { $0.rotation.isIdentity && !$0.isIdentity }, "reference finds the centring across the boundary")
        let tight = SymmetryDetection.detect(sites: sites, cell: cell, tolerance: .fractional(3e-3))
        #expect(!tight.operations.contains { $0.rotation.isIdentity && !$0.isIdentity })
    }

    // MARK: 500 atoms

    /// 500 atoms in a 20 Å P1 cell with the point group m-3m: ten generic orbits of 48, and the special positions.
    static func symmetric500() -> (UnitCell, [AtomSite]) {
        let group = PointSymmetry.holohedry(of: .cubic)
        var generator = SystemRandomNumberGeneratorSeeded(seed: 20261015)
        var sites: [AtomSite] = []
        func orbit(_ p: SIMD3<Double>, z: Int) {
            var seen: [SIMD3<Double>] = []
            for r in group.rotations {
                let q = wrapFractional(r.apply(p))
                if !seen.contains(where: { fractionalAxisDelta($0.x, q.x) < 1e-9 && fractionalAxisDelta($0.y, q.y) < 1e-9 && fractionalAxisDelta($0.z, q.z) < 1e-9 }) { seen.append(q) }
            }
            sites += seen.map { AtomSite(z: z, fractional: $0) }
        }
        var generic: [SIMD3<Double>] = []
        while generic.count < 10 {
            let p = SIMD3(Double.random(in: 0.05..<0.45, using: &generator), Double.random(in: 0.05..<0.25, using: &generator),
                          Double.random(in: 0.03..<0.15, using: &generator))
            if p.x > p.y + 0.03, p.y > p.z + 0.03, !generic.contains(where: { length($0 - p) < 0.06 }) { generic.append(p) }
        }
        for (i, p) in generic.enumerated() { orbit(p, z: 6 + i % 3) }
        orbit([0, 0, 0], z: 79); orbit([0.5, 0.5, 0.5], z: 80); orbit([0.5, 0, 0], z: 81); orbit([0, 0.5, 0.5], z: 82)
        orbit([0.1, 0, 0], z: 83); orbit([0.2, 0, 0], z: 84)
        return (.cubic(a: 20), sites)
    }

    @Test func testFiveHundredAtomSymmetricP1CellIsDetectedInUnderTwoSeconds() {
        let (cell, sites) = Self.symmetric500()
        #expect(sites.count == 500, "\(sites.count)")
        let clock = ContinuousClock()
        var detected: CrystalStructure?
        let elapsed = clock.measure { detected = CrystalStructure.detectingSymmetry(cell: cell, sites: sites) }
        let s = detected!
        print("slice 1b: 500-atom symmetric P1 cell, \(s.operations.count) operations, \(elapsed)")
        #expect(s.operations.count == 48 && s.laueClass == .mThreeBarM)
        #expect(elapsed < .seconds(2), "\(elapsed)")
        // the reference path returns the same operations
        let slow = SymmetryDetection.detect(sites: sites, cell: cell, search: .bruteForce)
        #expect(slow.operations == s.operations)
    }

    @Test func testFiveHundredAtomTranslationRichCell() {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        var sites: [AtomSite] = []
        for i in 0..<5 { for j in 0..<5 { for k in 0..<5 { for p in fcc {
            sites.append(AtomSite(z: 79, fractional: (p + SIMD3(Double(i), Double(j), Double(k))) / 5))
        } } } }
        let clock = ContinuousClock()
        var detection: SymmetryDetection.Result?
        let elapsed = clock.measure { detection = SymmetryDetection.detect(sites: sites, cell: .cubic(a: 5 * 4.0782)) }
        print("slice 1b: 500-atom 5×5×5 gold supercell, \(detection!.operations.count) operations, \(elapsed)")
        #expect(detection!.operations.count == 125 * 4 * 48)
        #expect(elapsed < .seconds(30), "\(elapsed)")      // slice 1: about 30 s for 500 atoms at N³
    }

    // MARK: The Cartesian tolerance

    /// A pseudo-centring of 0.013 Å in a 50 Å cell: a fractional 5e-3 is 0.25 Å and reads it as symmetry, a Cartesian 1e-3 Å
    /// does not (Gate B F6).
    @Test func testCartesianToleranceSeparatesAPseudoCentringFractionalCannotSee() {
        let shift = 0.013 / 50
        let sites = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5, 0.5, 0.5 + shift])]
        let cell = UnitCell.cubic(a: 50)
        func hasBodyCentring(_ s: CrystalStructure) -> Bool {
            s.operations.contains { $0.rotation.isIdentity && !$0.isIdentity && $0.isEquivalent(to: SymmetryOperation(rotation: .identity, translation: [0.5, 0.5, 0.5]), translationTolerance: 1e-3) }
        }
        let fractional = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        let cartesian = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, tolerance: .cartesian(1e-3))
        let loose = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, tolerance: .cartesian(0.05))
        #expect(hasBodyCentring(fractional))
        #expect(!hasBodyCentring(cartesian))
        #expect(hasBodyCentring(loose))
        #expect(cartesian.operations.count < fractional.operations.count, "\(cartesian.operations.count) vs \(fractional.operations.count)")
        // a sphere, not a box: an offset of 0.9e-5 along all three axes (1.56e-3 Å over two legs) is inside the per-axis
        // window of 1e-3 Å and outside the Euclidean distance
        let diagonal = [AtomSite(z: 26, fractional: [0, 0, 0]), AtomSite(z: 26, fractional: [0.5 + 0.9e-5, 0.5 + 0.9e-5, 0.5 + 0.9e-5])]
        #expect(!hasBodyCentring(CrystalStructure.detectingSymmetry(cell: cell, sites: diagonal, tolerance: .cartesian(1e-3))))
        #expect(hasBodyCentring(CrystalStructure.detectingSymmetry(cell: cell, sites: diagonal, tolerance: .cartesian(1e-2))))
        #expect(SymmetryDetection.Tolerance.default == .fractional(5e-3))
        // on an exactly symmetric structure the two senses agree
        let au = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: [0, 1, 2, 3].map { i in
            AtomSite(z: 79, fractional: [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]][i])
        })
        let auCart = CrystalStructure.detectingSymmetry(cell: .cubic(a: 4.0782), sites: au.sites, tolerance: .cartesian(1e-3))
        #expect(au.operations == auCart.operations && au.operations.count == 192)
    }
}

/// A small deterministic generator (SplitMix64), so the 500-atom cell is the same on every run.
struct SystemRandomNumberGeneratorSeeded: RandomNumberGenerator {
    var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9E3779B97F4A7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
        return z ^ (z >> 31)
    }
}
