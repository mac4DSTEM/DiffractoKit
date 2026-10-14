//
//  PinsTests.swift
//  Agreement with pymatgen (registration section 4), read from slice-1-pins.json (tools/slice-1/pins.py, run in
//  References/pymatgen-venv before any Swift existed). Agreement is agreement, not truth; a difference is a failure of
//  the port, never a tolerance to widen.
//

import XCTest
import CryptoKit
import simd
import DiffractoCore

final class PinsTests: XCTestCase {

    static let pins: Pins = {
        do { return try Pins.load() } catch { fatalError("slice-1-pins.json: \(error)") }
    }()

    func testPinsCarryTheirVersions() {
        let v = Self.pins.versions
        XCTAssertEqual(v["pymatgen_core"], "2026.9.23")
        XCTAssertEqual(v["pymatgen"], "2026.9.24")
        XCTAssertEqual(v["spglib"], "2.8.0")
        XCTAssertEqual(v["numpy"], "2.5.4")
        XCTAssertTrue(v["python"]?.hasPrefix("3.14.7") ?? false, "\(String(describing: v["python"]))")
        XCTAssertFalse(Self.pins.generatedAt.isEmpty)
        XCTAssertGreaterThanOrEqual(Self.pins.structures.count, 14)
        XCTAssertEqual(Self.pins.cells.count, 7)
    }

    // MARK: Cells (5.1 parameters; the lattice as pymatgen builds it with vesta=True)

    private func assertMatrix(_ got: [SIMD3<Double>], _ want: [[Double]], scale: Double, _ name: String) {
        for i in 0..<3 { for j in 0..<3 {
            XCTAssertEqual(got[i][j], want[i][j], accuracy: 1e-12 * scale, "\(name) [\(i)][\(j)]")
        } }
    }

    private func assertCell(_ cell: UnitCell, matrix: [[Double]], metric: [[Double]], volume: Double,
                            reciprocal: [[Double]], _ name: String) {
        let scale = max(cell.a, cell.b, cell.c)
        assertMatrix(cell.realLattice, matrix, scale: scale, "\(name) matrix")
        assertMatrix(cell.metricTensor, metric, scale: scale * scale, "\(name) metric")
        assertClose(cell.volume, volume, relative: 1e-12, "\(name) volume")
        assertMatrix(cell.reciprocalLattice, reciprocal, scale: 1 / min(cell.a, cell.b, cell.c), "\(name) reciprocal")
    }

    func testCellsAgreeWithPymatgen() {
        for pin in Self.pins.cells {
            let cell = pin.unitCell
            assertCell(cell, matrix: pin.matrix, metric: pin.metric, volume: pin.volume, reciprocal: pin.reciprocal, pin.name)
            for d in pin.d {
                let hkl = MillerIndex(d.hkl)
                assertClose(cell.d(hkl), d.d, relative: 1e-12, "\(pin.name) d\(hkl)")
                assertClose(cell.gLength(hkl), d.g, relative: 1e-12, "\(pin.name) g\(hkl)")
            }
            XCTAssertEqual(pin.d.count, 7 * 7 * 7 - 1, pin.name)
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
    func testBetaDoublePrimeCIFTextMatchesTheGenerator() throws {
        let pin = try Self.pins.structure("beta_double_prime")
        let sha = SHA256.hash(data: Data(Fixtures.betaDoublePrimeCIF.utf8)).map { String(format: "%02x", $0) }.joined()
        XCTAssertEqual(sha, pin.cifSha256)
        if let text = pin.cifText { XCTAssertEqual(text, Fixtures.betaDoublePrimeCIF) }
    }

    func testEveryStructureAgreesWithPymatgen() throws {
        for pin in Self.pins.structures {
            let s = try structure(for: pin)
            let name = pin.name
            // The cell.
            XCTAssertEqual(s.cell.a, pin.cell.a, name); XCTAssertEqual(s.cell.b, pin.cell.b, name)
            XCTAssertEqual(s.cell.c, pin.cell.c, name)
            XCTAssertEqual(s.cell.alphaDeg, pin.cell.alpha, name); XCTAssertEqual(s.cell.betaDeg, pin.cell.beta, name)
            XCTAssertEqual(s.cell.gammaDeg, pin.cell.gamma, name)
            assertCell(s.cell, matrix: pin.matrix, metric: pin.metric, volume: pin.volume, reciprocal: pin.reciprocal, name)
            // The expanded sites, as a multiset (against the written coordinates where pymatgen snapped ⅓ and ⅔).
            XCTAssertEqual(s.sites.count, pin.expectedSites.count, "\(name) site count")
            var unmatched = Array(s.sites.indices)
            for site in pin.expectedSites {
                let found = unmatched.firstIndex { i in
                    let mine = s.sites[i]
                    return mine.z == site.z && abs(mine.occupancy - site.occupancy) <= 1e-9
                        && fractionalAxisDelta(mine.fractional.x, site.frac[0]) <= 1e-9
                        && fractionalAxisDelta(mine.fractional.y, site.frac[1]) <= 1e-9
                        && fractionalAxisDelta(mine.fractional.z, site.frac[2]) <= 1e-9
                }
                if let found { unmatched.remove(at: found) } else { XCTFail("\(name): no site for Z\(site.z) \(site.frac)") }
            }
            XCTAssertTrue(unmatched.isEmpty, "\(name): unmatched sites \(unmatched)")
            // The symmetry.
            XCTAssertEqual(s.spaceGroupNumber, pin.declaredSpaceGroupNumber, "\(name) declared number")
            XCTAssertEqual(s.laueClass?.symbol, pin.laueClass, "\(name) Laue class")
            XCTAssertEqual(LaueClass(pointGroupSymbol: pin.pointGroup), s.laueClass, "\(name) point group \(pin.pointGroup)")
            XCTAssertEqual(s.operations.count, pin.operations.count, "\(name) operation count: mine \(s.operations.map(\.xyzString))")
            for op in pin.operations {
                XCTAssertTrue(s.operations.contains { $0.isEquivalent(to: op.operation, translationTolerance: 1e-6) },
                              "\(name): missing \(op.xyz)")
            }
            XCTAssertEqual(Set(s.laueGroup.rotations), Set(pin.laueRotations.map { IntegerMatrix($0) }), "\(name) Laue rotations")
            // The file's loop against the detected set (amendment 2).
            XCTAssertEqual(s.fileOperations?.count, pin.fileOperations?.count, "\(name) file operations")
            XCTAssertEqual(s.fileOperationsAreSymmetries, pin.fileOperations?.subsetOfDetected, "\(name) file operations are symmetries")
            if let file = pin.fileOperations {
                XCTAssertEqual(file.equalsDetected, s.operations.count == file.count, "\(name) loop equals detected")
            }
            // The reflections.
            let list = try s.reflections(kMax: pin.reflections.kMax)
            XCTAssertEqual(list.reflections.count, pin.reflections.count, "\(name) reflection count")
            XCTAssertEqual(list.allowed.count, pin.reflections.allowedCount, "\(name) allowed count")
            for (mine, want) in zip(list.reflections, pin.reflections.items) {
                XCTAssertEqual(mine.hkl, MillerIndex(want.hkl), "\(name) order")
                guard mine.hkl == MillerIndex(want.hkl) else { break }
                assertClose(mine.gLength, want.g, relative: 1e-12, "\(name) g\(mine.hkl)")
                assertClose(mine.d, want.d, relative: 1e-12, "\(name) d\(mine.hkl)")
                XCTAssertEqual(mine.isAllowed, want.allowed, "\(name) allowed \(mine.hkl)")
                if let orbit = want.orbit {
                    XCTAssertEqual(s.laueGroup.orbit(of: mine.hkl), Set(orbit.map { MillerIndex($0) }), "\(name) orbit \(mine.hkl)")
                }
                XCTAssertEqual(s.laueGroup.canonical(mine.hkl), MillerIndex(want.canonical), "\(name) canonical \(mine.hkl)")
            }
            // The families.
            let families = list.families()
            XCTAssertEqual(families.count, pin.families.count, "\(name) family count")
            for (mine, want) in zip(families, pin.families) {
                assertClose(mine.g, want.g, relative: 1e-12, "\(name) family g")
                XCTAssertEqual(mine.members, want.members.map { MillerIndex($0) }, "\(name) members of \(want.name)")
                XCTAssertEqual(mine.multiplicity, want.multiplicity, "\(name) multiplicity of \(want.name)")
                XCTAssertEqual(mine.orbits, want.orbits.map { MillerIndex($0) }, "\(name) orbits of \(want.name)")
                XCTAssertEqual(mine.name, want.name, "\(name) name")
            }
        }
    }
}
