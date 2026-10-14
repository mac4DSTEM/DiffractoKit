//
//  UnitCellTests.swift
//  Planted truth for the cell (registration 5.1): the closed-form d of each system on rational cells, the volume by
//  both formulas, d by the metric and by the vector, the reciprocal rows dual to the real ones, and the metric
//  family of each constructor. No reference code is read.
//

import XCTest
import simd
import DiffractoCore

final class UnitCellTests: XCTestCase {

    static let cells: [(String, UnitCell)] = [
        ("monoclinic", .monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)),
        ("cubic", .cubic(a: 4)),
        ("tetragonal", .tetragonal(a: 3, c: 5)),
        ("orthorhombic", .orthorhombic(a: 3, b: 4, c: 5)),
        ("hexagonal", .hexagonal(a: 3, c: 5)),
        ("rhombohedral", .rhombohedral(a: 4, alphaDeg: 70)),
        ("triclinic", .triclinic(a: 4, b: 5, c: 6, alphaDeg: 80, betaDeg: 85, gammaDeg: 95)),
    ]

    static var indices: [MillerIndex] {
        var out: [MillerIndex] = []
        for h in -3...3 { for k in -3...3 { for l in -3...3 where !(h == 0 && k == 0 && l == 0) {
            out.append(MillerIndex(h, k, l))
        } } }
        return out
    }

    /// V = |det A| equals abc √(1 − cos²α − cos²β − cos²γ + 2 cos α cos β cos γ). Mutant: the volume from the
    /// metric's determinant without the square root.
    func testVolumeByBothFormulas() {
        for (name, cell) in Self.cells {
            let ca = cos(cell.alphaDeg * .pi / 180), cb = cos(cell.betaDeg * .pi / 180), cg = cos(cell.gammaDeg * .pi / 180)
            let closed = cell.a * cell.b * cell.c * (1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg).squareRoot()
            assertClose(cell.volume, closed, relative: 1e-12, name)
        }
    }

    /// aᵢ · bⱼ = δᵢⱼ. Mutant: the reciprocal lattice as G A instead of G⁻¹ A.
    func testReciprocalRowsAreDualToTheRealOnes() {
        for (name, cell) in Self.cells {
            for i in 0..<3 { for j in 0..<3 {
                let want: Double = i == j ? 1 : 0
                XCTAssertEqual(dot(cell.realLattice[i], cell.reciprocalLattice[j]), want, accuracy: 1e-12, "\(name) \(i)\(j)")
            } }
        }
    }

    /// d from G* equals 1/|g| from the Cartesian vector. Mutant: G* as G.
    func testDFromMetricEqualsDFromVector() {
        for (name, cell) in Self.cells {
            for hkl in Self.indices {
                assertClose(cell.d(hkl), 1 / cell.gLength(hkl), relative: 1e-12, "\(name) \(hkl)")
            }
        }
    }

    /// The textbook formulas, one per system (5.1), to 1e-12 relative. Mutants: β read from α; the per-axis tile.
    func testClosedFormDPerSystem() {
        func check(_ cell: UnitCell, _ name: String, _ oneOverDSquared: (Double, Double, Double) -> Double) {
            for hkl in Self.indices {
                let (h, k, l) = (Double(hkl.h), Double(hkl.k), Double(hkl.l))
                let want = 1 / oneOverDSquared(h, k, l).squareRoot()
                assertClose(cell.d(hkl), want, relative: 1e-12, "\(name) \(hkl)")
            }
        }
        // Monoclinic first (β″'s cell).
        do {
            let (a, b, c, beta) = (15.16, 4.05, 6.74, 105.3 * Double.pi / 180)
            check(.monoclinic(a: a, b: b, c: c, betaDeg: 105.3), "monoclinic") { h, k, l in
                (h * h / (a * a) + l * l / (c * c) - 2 * h * l * cos(beta) / (a * c)) / (sin(beta) * sin(beta))
                    + k * k / (b * b)
            }
        }
        check(.cubic(a: 4), "cubic") { h, k, l in (h * h + k * k + l * l) / 16 }
        check(.tetragonal(a: 3, c: 5), "tetragonal") { h, k, l in (h * h + k * k) / 9 + l * l / 25 }
        check(.orthorhombic(a: 3, b: 4, c: 5), "orthorhombic") { h, k, l in h * h / 9 + k * k / 16 + l * l / 25 }
        check(.hexagonal(a: 3, c: 5), "hexagonal") { h, k, l in (4.0 / 3) * (h * h + h * k + k * k) / 9 + l * l / 25 }
        do {
            let (a, alpha) = (4.0, 70 * Double.pi / 180)
            let ca = cos(alpha)
            check(.rhombohedral(a: 4, alphaDeg: 70), "rhombohedral") { h, k, l in
                ((h * h + k * k + l * l) * sin(alpha) * sin(alpha) + 2 * (h * k + k * l + h * l) * (ca * ca - ca))
                    / (a * a * (1 - 3 * ca * ca + 2 * ca * ca * ca))
            }
        }
        do {
            let (a, b, c) = (4.0, 5.0, 6.0)
            let (al, be, ga) = (80 * Double.pi / 180, 85 * Double.pi / 180, 95 * Double.pi / 180)
            let (ca, cb, cg) = (cos(al), cos(be), cos(ga))
            let v = a * b * c * (1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg).squareRoot()
            let s11 = b * b * c * c * sin(al) * sin(al), s22 = a * a * c * c * sin(be) * sin(be)
            let s33 = a * a * b * b * sin(ga) * sin(ga)
            let s12 = a * b * c * c * (ca * cb - cg), s23 = a * a * b * c * (cb * cg - ca), s13 = a * b * b * c * (ca * cg - cb)
            check(.triclinic(a: a, b: b, c: c, alphaDeg: 80, betaDeg: 85, gammaDeg: 95), "triclinic") { h, k, l in
                (s11 * h * h + s22 * k * k + s33 * l * l + 2 * s12 * h * k + 2 * s23 * k * l + 2 * s13 * h * l) / (v * v)
            }
        }
    }

    /// macTEM's S2c hand anchors for β″ (d in Å). Mutant: β read from α (d(200) becomes a/2).
    func testBetaDoublePrimeHandAnchors() {
        let cell = Fixtures.betaDoublePrimeCell
        let anchors: [(MillerIndex, Double)] = [
            (MillerIndex(2, 0, 0), 7.31135), (MillerIndex(0, 0, 1), 6.50112), (MillerIndex(0, 2, 0), 2.025),
            (MillerIndex(1, 1, 0), 3.90306), (MillerIndex(2, 0, 1), 4.32457), (MillerIndex(2, 0, -1), 5.65554),
        ]
        for (hkl, d) in anchors {
            XCTAssertEqual(cell.d(hkl), d, accuracy: 6e-6, "\(hkl)")
        }
        XCTAssertEqual(cell.d(MillerIndex(2, 0, 0)), cell.a * sin(105.3 * .pi / 180) / 2, accuracy: 1e-12)
        XCTAssertEqual(cell.d(MillerIndex(0, 0, 1)), cell.c * sin(105.3 * .pi / 180), accuracy: 1e-12)
        XCTAssertNotEqual(cell.d(MillerIndex(2, 0, 1)), cell.d(MillerIndex(2, 0, -1)), "for β ≠ 90° these are two rings")
    }

    /// Each constructor's metric family, first that applies; a trigonal structure's hexagonal cell is `.hexagonal` here.
    func testMetricFamilyOfEachConstructor() {
        XCTAssertEqual(UnitCell.cubic(a: 4).metricFamily(), .cubic)
        XCTAssertEqual(UnitCell.hexagonal(a: 3, c: 5).metricFamily(), .hexagonal)
        XCTAssertEqual(UnitCell.tetragonal(a: 3, c: 5).metricFamily(), .tetragonal)
        XCTAssertEqual(UnitCell.rhombohedral(a: 4, alphaDeg: 70).metricFamily(), .rhombohedral)
        XCTAssertEqual(UnitCell.orthorhombic(a: 3, b: 4, c: 5).metricFamily(), .orthorhombic)
        XCTAssertEqual(UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3).metricFamily(),
                       .monoclinic(uniqueAxis: .b))
        XCTAssertEqual(UnitCell(a: 3, b: 4, c: 5, alphaDeg: 100, betaDeg: 90, gammaDeg: 90).metricFamily(),
                       .monoclinic(uniqueAxis: .a))
        XCTAssertEqual(UnitCell(a: 3, b: 4, c: 5, alphaDeg: 90, betaDeg: 90, gammaDeg: 100).metricFamily(),
                       .monoclinic(uniqueAxis: .c))
        XCTAssertEqual(UnitCell.triclinic(a: 4, b: 5, c: 6, alphaDeg: 80, betaDeg: 85, gammaDeg: 95).metricFamily(), .triclinic)
        // A tetragonal cell whose c equals a is cubic by metric; the atoms decide (SymmetryDetectionTests).
        XCTAssertEqual(UnitCell.tetragonal(a: 4, c: 4).metricFamily(), .cubic)
        // The tolerance is the caller's: β = 90° + 1e-8 is monoclinic at 1e-9° and orthorhombic at 0.05°.
        let nearly = UnitCell(a: 3, b: 4, c: 5, alphaDeg: 90, betaDeg: 90 + 1e-8, gammaDeg: 90)
        XCTAssertEqual(nearly.metricFamily(), .monoclinic(uniqueAxis: .b))
        XCTAssertEqual(nearly.metricFamily(lengthTolerance: 1e-3, angleToleranceDeg: 0.05), .orthorhombic)
    }

    func testIssuesNameWhatIsWrong() {
        XCTAssertEqual(UnitCell.cubic(a: 4).issues, [])
        XCTAssertEqual(UnitCell(a: 0, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90).issues.first, .invalidLengths)
        XCTAssertEqual(UnitCell(a: 4, b: 4, c: 4, alphaDeg: 180, betaDeg: 90, gammaDeg: 90).issues.first, .invalidAngles)
        // α = β = γ = 130° gives a negative radicand: no finite volume (120° is a flat cell of volume ~1e-6 by rounding).
        let flat = UnitCell(a: 4, b: 4, c: 4, alphaDeg: 130, betaDeg: 130, gammaDeg: 130)
        XCTAssertTrue(flat.issues.contains(.degenerateVolume))
        XCTAssertFalse(flat.isValid)
    }

    func testInterplanarAngles() {
        let cubic = UnitCell.cubic(a: 4)
        XCTAssertEqual(cubic.interplanarAngleDeg(MillerIndex(1, 0, 0), MillerIndex(0, 1, 0)), 90, accuracy: 1e-12)
        XCTAssertEqual(cubic.interplanarAngleDeg(MillerIndex(1, 0, 0), MillerIndex(1, 1, 0)), 45, accuracy: 1e-12)
        XCTAssertEqual(cubic.interplanarAngleDeg(MillerIndex(1, 1, 1), MillerIndex(1, 1, 1)), 0, accuracy: 1e-6)
        // In a hexagonal cell (100) and (010) are 60° apart, not 90°.
        XCTAssertEqual(UnitCell.hexagonal(a: 3, c: 5).interplanarAngleDeg(MillerIndex(1, 0, 0), MillerIndex(0, 1, 0)),
                       60, accuracy: 1e-12)
    }
}
