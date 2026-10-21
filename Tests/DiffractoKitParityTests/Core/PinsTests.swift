//
//  PinsTests.swift
//  Agreement with pymatgen (registration section 4), read from slice-1-pins.json (tools/slice-1/pins.py, run in
//  ../References/python/orix before any Swift existed). Agreement is agreement, not truth; a difference is a failure of
//  the port, never a tolerance to widen.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import CryptoKit
import simd
import DiffractoCore

struct PinsTests {

    static let pins: Pins = {
        do { return try Pins.load() } catch { fatalError("slice-1-pins.json: \(error)") }
    }()

    @Test func testPinsCarryTheirVersions() {
        let v = Self.pins.versions
        #expect(v["pymatgen_core"] == "2026.9.23")
        #expect(v["pymatgen"] == "2026.9.24")
        #expect(v["spglib"] == "2.8.0")
        #expect(v["numpy"] == "2.5.4")
        #expect(v["python"]?.hasPrefix("3.14.7") ?? false, "\(String(describing: v["python"]))")
        #expect(!Self.pins.generatedAt.isEmpty)
        #expect(Self.pins.structures.count >= 14)
        #expect(Self.pins.cells.count == 7)
    }

    // MARK: Cells (5.1 parameters; the lattice as pymatgen builds it with vesta=True)

    private func assertMatrix(_ got: [SIMD3<Double>], _ want: [[Double]], scale: Double, _ name: String) {
        for i in 0..<3 { for j in 0..<3 {
            // D17: numpy assert_allclose defaults (rtol 1e-7), with an absolute floor of 1e-12 of the matrix scale for exact zeros.
            #expect(abs(got[i][j] - want[i][j]) <= 1e-12 * scale + 1e-7 * abs(want[i][j]), "\(name) [\(i)][\(j)]")
        } }
    }

    private func assertCell(_ cell: UnitCell, matrix: [[Double]], metric: [[Double]], volume: Double,
                            reciprocal: [[Double]], _ name: String) {
        let scale = max(cell.a, cell.b, cell.c)
        assertMatrix(cell.realLattice, matrix, scale: scale, "\(name) matrix")
        assertMatrix(cell.metricTensor, metric, scale: scale * scale, "\(name) metric")
        assertClose(cell.volume, volume, relative: 1e-7, "\(name) volume")
        assertMatrix(cell.reciprocalLattice, reciprocal, scale: 1 / min(cell.a, cell.b, cell.c), "\(name) reciprocal")
    }

    @Test func testCellsAgreeWithPymatgen() {
        for pin in Self.pins.cells {
            let cell = pin.unitCell
            assertCell(cell, matrix: pin.matrix, metric: pin.metric, volume: pin.volume, reciprocal: pin.reciprocal, pin.name)
            for d in pin.d {
                let hkl = MillerIndex(d.hkl)
                assertClose(cell.d(hkl), d.d, relative: 1e-7, "\(pin.name) d\(hkl)")
                assertClose(cell.gLength(hkl), d.g, relative: 1e-7, "\(pin.name) g\(hkl)")
            }
            #expect(pin.d.count == 7 * 7 * 7 - 1, "\(pin.name)")
        }
    }

    // MARK: Structures

    /// β″ comes from the Swift-built text (the generator's file is gitignored); everything else from its fixture file.
    private func structure(for pin: StructurePin) throws -> CrystalStructure {
        switch pin.name {
        case "beta_double_prime":
            return try CIFImport.read(Fixtures.betaDoublePrimeCIF, fileBaseName: pin.name).structure
        case "beta_double_prime_P1":
            let expanded = try CIFImport.read(Fixtures.betaDoublePrimeCIF).structure
            return CrystalStructure.detectingSymmetry(cell: expanded.cell, sites: expanded.sites)
        default:
            guard let file = pin.cifFile else { throw CheckFailure("no way to build \(pin.name)") }
            return try CIFImport.read(try Fixtures.text(file), fileBaseName: file).structure
        }
    }

    struct CheckFailure: Error { let message: String; init(_ m: String) { message = m } }

    /// The β″ CIF rebuilt in Swift is the generator's byte for byte (the pin's sha256 of the file it wrote).
    @Test func testBetaDoublePrimeCIFTextMatchesTheGenerator() throws {
        let pin = try Self.pins.structure("beta_double_prime")
        let sha = SHA256.hash(data: Data(Fixtures.betaDoublePrimeCIF.utf8)).map { String(format: "%02x", $0) }.joined()
        #expect(sha == pin.cifSha256)
        if let text = pin.cifText { #expect(text == Fixtures.betaDoublePrimeCIF) }
    }

    @Test func testEveryStructureAgreesWithPymatgen() throws {
        for pin in Self.pins.structures {
            let s = try structure(for: pin)
            let name = pin.name
            // The cell.
            #expect(s.cell.a == pin.cell.a, "\(name)"); #expect(s.cell.b == pin.cell.b, "\(name)")
            #expect(s.cell.c == pin.cell.c, "\(name)")
            #expect(s.cell.alphaDeg == pin.cell.alpha, "\(name)"); #expect(s.cell.betaDeg == pin.cell.beta, "\(name)")
            #expect(s.cell.gammaDeg == pin.cell.gamma, "\(name)")
            assertCell(s.cell, matrix: pin.matrix, metric: pin.metric, volume: pin.volume, reciprocal: pin.reciprocal, name)
            // The expanded sites, as a multiset (against the written coordinates where pymatgen snapped ⅓ and ⅔).
            #expect(s.sites.count == pin.expectedSites.count, "\(name) site count")
            var unmatched = Array(s.sites.indices)
            for site in pin.expectedSites {
                let found = unmatched.firstIndex { i in
                    let mine = s.sites[i]
                    return mine.z == site.z && abs(mine.occupancy - site.occupancy) <= 1e-7
                        && fractionalAxisDelta(mine.fractional.x, site.frac[0]) <= 1e-7
                        && fractionalAxisDelta(mine.fractional.y, site.frac[1]) <= 1e-7
                        && fractionalAxisDelta(mine.fractional.z, site.frac[2]) <= 1e-7
                }
                if let found { unmatched.remove(at: found) } else { Issue.record("\(name): no site for Z\(site.z) \(site.frac)") }
            }
            #expect(unmatched.isEmpty, "\(name): unmatched sites \(unmatched)")
            // The symmetry.
            #expect(s.spaceGroupNumber == pin.declaredSpaceGroupNumber, "\(name) declared number")
            #expect(s.laueClass?.symbol == pin.laueClass, "\(name) Laue class")
            #expect(LaueClass(pointGroupSymbol: pin.pointGroup) == s.laueClass, "\(name) point group \(pin.pointGroup)")
            #expect(s.operations.count == pin.operations.count, "\(name) operation count: mine \(s.operations.map(\.xyzString))")
            for op in pin.operations {
                #expect(s.operations.contains { $0.isEquivalent(to: op.operation, translationTolerance: 0.01 / max(s.cell.a, s.cell.b, s.cell.c)) }, "\(name): missing \(op.xyz)")
            }
            #expect(Set(s.laueGroup.rotations) == Set(pin.laueRotations.map { IntegerMatrix($0) }), "\(name) Laue rotations")
            // The file's loop against the detected set (amendment 2).
            #expect(s.fileOperations?.count == pin.fileOperations?.count, "\(name) file operations")
            #expect(s.fileOperationsAreSymmetries == pin.fileOperations?.subsetOfDetected, "\(name) file operations are symmetries")
            if let file = pin.fileOperations {
                #expect(file.equalsDetected == (s.operations.count == file.count), "\(name) loop equals detected")
            }
            // The reflections.
            let list = try s.reflections(kMax: pin.reflections.kMax)
            #expect(list.reflections.count == pin.reflections.count, "\(name) reflection count")
            #expect(list.allowed.count == pin.reflections.allowedCount, "\(name) allowed count")
            for (mine, want) in zip(list.reflections, pin.reflections.items) {
                #expect(mine.hkl == MillerIndex(want.hkl), "\(name) order")
                guard mine.hkl == MillerIndex(want.hkl) else { break }
                assertClose(mine.gLength, want.g, relative: 1e-7, "\(name) g\(mine.hkl)")
                assertClose(mine.d, want.d, relative: 1e-7, "\(name) d\(mine.hkl)")
                #expect(mine.isAllowed == want.allowed, "\(name) allowed \(mine.hkl)")
                if let orbit = want.orbit {
                    #expect(s.laueGroup.orbit(of: mine.hkl) == Set(orbit.map { MillerIndex($0) }), "\(name) orbit \(mine.hkl)")
                }
                #expect(s.laueGroup.canonical(mine.hkl) == MillerIndex(want.canonical), "\(name) canonical \(mine.hkl)")
            }
            // The families.
            let families = list.families()
            #expect(families.count == pin.families.count, "\(name) family count")
            for (mine, want) in zip(families, pin.families) {
                assertClose(mine.g, want.g, relative: 1e-7, "\(name) family g")
                #expect(mine.members == want.members.map { MillerIndex($0) }, "\(name) members of \(want.name)")
                #expect(mine.multiplicity == want.multiplicity, "\(name) multiplicity of \(want.name)")
                #expect(mine.orbits == want.orbits.map { MillerIndex($0) }, "\(name) orbits of \(want.name)")
                #expect(mine.name == want.name, "\(name) name")
            }
        }
    }
}
