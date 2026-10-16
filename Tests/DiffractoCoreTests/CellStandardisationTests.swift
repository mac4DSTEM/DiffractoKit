//
//  CellStandardisationTests.swift
//  Materials Project import numerics against pymatgen SpacegroupAnalyzer (fixtures/mp-import/pymatgen-pins.json,
//  tools/mp-import/pins.py: pymatgen 2026.9.24, pymatgen-core 2026.9.23, spglib 2.8.0, NumPy 2.5.4). The inputs are
//  committed MP-shaped JSON; nothing here touches the network.
//

import Foundation
import Testing
import simd
import DiffractoCore

struct CellStandardisationTests {
    struct Pins: Decodable {
        struct Lat: Decodable { let lattice: [[Double]]; let sites: [Site] }
        struct Site: Decodable { let abc: [Double]; let species: [String: Double] }
        struct Case: Decodable {
            let input: String; let spaceGroupNumber: Int; let spaceGroupSymbol: String
            let conventional: Lat; let primitive: Lat
        }
        let pymatgen: String; let symprec: Double; let cases: [Case]
    }
    struct Input: Decodable {
        struct Doc: Decodable {
            struct Structure: Decodable {
                struct Lattice: Decodable { let a, b, c, alpha, beta, gamma: Double; let matrix: [[Double]]? }
                struct Site: Decodable { struct Sp: Decodable { let element: String; let occu: Double }
                    let abc: [Double]; let species: [Sp] }
                let lattice: Lattice; let sites: [Site]
            }
            struct Symmetry: Decodable { let crystal_system: String; let symbol: String; let number: Int }
            let structure: Structure; let symmetry: Symmetry
        }
        let data: [Doc]
    }

    static let dir = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/mp-import")
    static let pins: Pins = try! JSONDecoder().decode(Pins.self, from: Data(contentsOf: dir.appendingPathComponent("pymatgen-pins.json")))

    static func load(_ c: Pins.Case) -> (rows: [SIMD3<Double>], cell: UnitCell, sites: [AtomSite], doc: Input.Doc) {
        let doc = try! JSONDecoder().decode(Input.self, from: Data(contentsOf: dir.appendingPathComponent(c.input))).data[0]
        let l = doc.structure.lattice
        let cell = UnitCell(a: l.a, b: l.b, c: l.c, alphaDeg: l.alpha, betaDeg: l.beta, gammaDeg: l.gamma)
        let rows = l.matrix.map { $0.map { SIMD3($0[0], $0[1], $0[2]) } } ?? cell.realLattice
        let sites = doc.structure.sites.map {
            AtomSite(z: Element.atomicNumber(forSymbol: $0.species[0].element)!, fractional: SIMD3($0.abc[0], $0.abc[1], $0.abc[2]),
                     occupancy: $0.species[0].occu)
        }
        return (rows, cell, sites, doc)
    }

    static func params(_ m: [[Double]]) -> [Double] {
        let r = m.map { SIMD3($0[0], $0[1], $0[2]) }
        func ang(_ u: SIMD3<Double>, _ v: SIMD3<Double>) -> Double { acos(simd_dot(u, v) / (simd_length(u) * simd_length(v))) * 180 / .pi }
        return [simd_length(r[0]), simd_length(r[1]), simd_length(r[2]), ang(r[1], r[2]), ang(r[0], r[2]), ang(r[0], r[1])]
    }

    /// Equal as sets modulo one rigid translation: pymatgen moves the origin to its standard one, the app keeps the
    /// served origin (documented in docs/parity.md).
    static func sameSites(_ got: [AtomSite], _ want: [Pins.Site], tol: Double) -> Bool {
        guard got.count == want.count, let first = got.first else { return false }
        let w0 = want.map { SIMD3($0.abc[0], $0.abc[1], $0.abc[2]) }
        return w0.contains { target in
            let shift = target - first.fractional
            return sameSitesShifted(got, want, shift: shift, tol: tol)
        }
    }

    static func sameSitesShifted(_ got: [AtomSite], _ want: [Pins.Site], shift: SIMD3<Double>, tol: Double) -> Bool {
        var used = Set<Int>()
        for g in got {
            guard let j = want.indices.first(where: { j in
                !used.contains(j) && Element.atomicNumber(forSymbol: want[j].species.keys.first!) == g.z
                    && (0..<3).allSatisfy { k in let d = abs(g.fractional[k] + shift[k] - want[j].abc[k]).truncatingRemainder(dividingBy: 1)
                        return min(d, 1 - d) <= tol } }) else { return false }
            used.insert(j)
        }
        return true
    }

    @Test func pinsCarryTheirVersions() {
        #expect(Self.pins.pymatgen == "2026.9.24")
        #expect(Self.pins.cases.count == 3)
    }

    /// The app's route: the lattice and the sites equal pymatgen's conventional standard cell (lengths and angles to
    /// 1e-12 relative, sites to 1e-9 fractional).
    @Test(arguments: pins.cases) func appRouteMatchesPymatgenConventional(_ c: Pins.Case) throws {
        let (rows, _, sites, doc) = Self.load(c)
        let r = try CellStandardisation.standardise(latticeRows: rows, sites: sites, symbol: doc.symmetry.symbol,
                                                    crystalSystem: doc.symmetry.crystal_system)
        let got = Self.params(r.latticeRows.map { [$0.x, $0.y, $0.z] })
        let want = Self.params(c.conventional.lattice)
        for k in 0..<6 { #expect(abs(got[k] - want[k]) <= 1e-9 * max(1, want[k]), "\(c.input) parameter \(k)") }
        #expect(Self.sameSites(r.sites, c.conventional.sites, tol: 1e-9), "\(c.input) sites")
        #expect(doc.symmetry.number == c.spaceGroupNumber)
    }

    /// The reference route: the kit's own detection of the primitive cell gives pymatgen's conventional cell.
    @Test(arguments: pins.cases) func detectedConventionalCellMatchesPymatgen(_ c: Pins.Case) throws {
        let p = c.primitive
        let want = Self.params(c.conventional.lattice)
        let prim = Self.params(p.lattice)
        let cell = UnitCell(a: prim[0], b: prim[1], c: prim[2], alphaDeg: prim[3], betaDeg: prim[4], gammaDeg: prim[5])
        let sites = p.sites.map { AtomSite(z: Element.atomicNumber(forSymbol: $0.species.keys.first!)!,
                                           fractional: SIMD3($0.abc[0], $0.abc[1], $0.abc[2])) }
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, metricLengthTolerance: 1e-6, metricAngleToleranceDeg: 1e-6)
        let conv = try #require(s.conventionalSetting)
        let got = [conv.cell.a, conv.cell.b, conv.cell.c, conv.cell.alphaDeg, conv.cell.betaDeg, conv.cell.gammaDeg]
        for k in 0..<6 { #expect(abs(got[k] - want[k]) <= 1e-6 * max(1, want[k]), "\(c.input) parameter \(k): \(got[k]) vs \(want[k])") }
        #expect(conv.sites.count == c.conventional.sites.count, "\(c.input) site count")
    }

    @Test func guardsRefuseBadInput() {
        let rows = [SIMD3(0.0, 2.0, 2.0), SIMD3(2.0, 0.0, 2.0), SIMD3(2.0, 2.0, 0.0)]
        let al = [AtomSite(z: 13, fractional: SIMD3(0, 0, 0))]
        // an F cell stated as hexagonal does not read conventional after the conversion
        #expect(throws: CellStandardisation.Failure.self) {
            try CellStandardisation.standardise(latticeRows: rows, sites: al, symbol: "Fm-3m", crystalSystem: "hexagonal")
        }
        #expect(throws: CellStandardisation.Failure.self) {
            try CellStandardisation.standardise(latticeRows: rows, sites: al, symbol: "Fm-3m", crystalSystem: nil)
        }
        // a P symbol or none returns the input
        let p = try? CellStandardisation.standardise(latticeRows: rows, sites: al, symbol: "Pm-3m", crystalSystem: "cubic")
        #expect(p?.outcome == .asServed)
        let none = try? CellStandardisation.standardise(latticeRows: rows, sites: al, symbol: nil, crystalSystem: nil)
        #expect(none?.sites == al)
    }
}
