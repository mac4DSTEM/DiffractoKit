//
//  ReductionTests.swift
//  The lattice reduction: the Niggli-reduced
//  cell against spglib and pymatgen on every pinned cell and 2800 random cells, and against the Niggli conditions
//  written out from the text; the transformation is unimodular and reproduces the reduced metric.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct ReductionTests {

    static let pins = ConventionalCellPins.shared

    private struct Case { let name: String; let cell: UnitCell; let spglib: [Double]; let pymatgen: [Double] }

    private static var cases: [Case] {
        pins.structures.map { Case(name: $0.name, cell: $0.cell.unitCell, spglib: $0.niggliSpglib, pymatgen: $0.niggliPymatgen) }
            + pins.randomCells.enumerated().map {
                Case(name: "random \($1.intended) #\($0)", cell: $1.cell.unitCell, spglib: $1.niggliSpglib, pymatgen: $1.niggliPymatgen)
            }
    }

    @Test func testPinsCarryTheirVersionsAndSizes() {
        let v = Self.pins.versions
        #expect(v["spglib"] == "2.8.0")
        #expect(v["pymatgen_core"] == "2026.10.16")
        #expect(v["numpy"] == "2.5.4")
        #expect(Self.pins.structures.count == 27)
        #expect(Self.pins.randomCells.count == 14 * 200)
        for name in BravaisLattice.allCases {
            let ofThisLattice: Int = Self.pins.randomCells.filter { $0.intended == name.rawValue }.count
            #expect(ofThisLattice == 200, "\(name)")
        }
    }

    /// The reduced cell's six parameters equal spglib's and pymatgen's (two independent implementations of the same
    /// algorithm) to 1e-5, the functions' own tol (the reference's tolerance; was 1e-10; relative on the lengths, degrees on the angles).
    @Test func testNiggliReducedCellAgreesWithSpglibAndPymatgen() {
        var worstSpglib = 0.0, worstPymatgen = 0.0
        for c in Self.cases {
            let reduced = c.cell.niggliReduced()
            #expect(reduced.converged, "\(c.name)")
            let mine = parameters(reduced.cell)
            let ds = parameterDeviation(mine, c.spglib), dp = parameterDeviation(mine, c.pymatgen)
            worstSpglib = max(worstSpglib, ds); worstPymatgen = max(worstPymatgen, dp)
            #expect(ds <= 1e-5, "\(c.name) vs spglib: \(mine) \(c.spglib)")
            #expect(dp <= 1e-5, "\(c.name) vs pymatgen: \(mine) \(c.pymatgen)")
        }
        print("reduction: \(Self.cases.count) cells, worst deviation from spglib \(worstSpglib), from pymatgen \(worstPymatgen)")
    }

    /// The Niggli conditions, checked directly on the reduced basis (not through the oracle); P is an integer matrix of
    /// determinant +1 that reproduces the reduced metric and the volume.
    @Test func testReducedBasisSatisfiesTheNiggliConditionsAndPIsUnimodular() {
        for c in Self.cases {
            let reduced = c.cell.niggliReduced()
            let p = reduced.transformation
            #expect(p.determinant == 1, "\(c.name): det P = \(p.determinant)")
            let v = transformedVectors(c.cell, p)
            let violations = niggliViolations(v, eps: reduced.epsilon)
            #expect(violations.isEmpty, "\(c.name): \(violations)")
            let g = metric(of: v)
            let want = reduced.cell.metricTensor
            let scale = max(want[0].x, want[1].y, want[2].z)
            for i in 0..<3 { for j in 0..<3 {
                #expect(abs(g[i][j] - want[i][j]) <= 1e-12 * scale, "\(c.name) Pᵀ G P [\(i)][\(j)]")
            } }
            #expect(abs(reduced.cell.volume - c.cell.volume) <= 1e-12 * c.cell.volume, "\(c.name) volume")
            #expect(reduced.epsilon == UnitCell.niggliRelativeEpsilon * pow(c.cell.volume, 2.0 / 3.0), "\(c.name) epsilon")
        }
    }

    /// The checker itself is not vacuous: a skewed cell fails it, its reduction passes.
    @Test func testTheConditionCheckerRejectsAnUnreducedBasis() {
        let skew = UnitCell(a: 3, b: 9, c: 30, alphaDeg: 95, betaDeg: 70, gammaDeg: 100)
        let eps = UnitCell.niggliRelativeEpsilon * pow(skew.volume, 2.0 / 3.0)
        #expect(!niggliViolations(skew.realLattice, eps: eps).isEmpty)
        let reduced = skew.niggliReduced()
        #expect(niggliViolations(transformedVectors(skew, reduced.transformation), eps: reduced.epsilon).isEmpty)
        #expect(reduced.cell.volume > 0 && abs(reduced.cell.volume - skew.volume) < 1e-9 * skew.volume)
    }

    /// Near a boundary between two reduced forms (cells within 1e-6 of a tie) the reduction still stops, at a basis that
    /// satisfies the conditions within ε, and a cell is reduced once: reducing the reduced cell changes nothing.
    @Test func testReductionStopsAtBoundariesAndIsIdempotent() {
        for delta in [0.0, 1e-12, 1e-9, 1e-6, -1e-6] {
            for base in [UnitCell(a: 4, b: 4 + delta, c: 4, alphaDeg: 60, betaDeg: 60, gammaDeg: 60),
                         UnitCell(a: 4, b: 4, c: 6 + delta, alphaDeg: 90, betaDeg: 90, gammaDeg: 120),
                         UnitCell(a: 5, b: 5, c: 5, alphaDeg: 90 + delta, betaDeg: 90, gammaDeg: 90),
                         UnitCell(a: 3, b: 4, c: 3 + 2 * delta, alphaDeg: 75, betaDeg: 109.4712206 + delta, gammaDeg: 75)] {
                let r = base.niggliReduced()
                #expect(r.converged, "\(parameters(base))")
                #expect(niggliViolations(transformedVectors(base, r.transformation), eps: r.epsilon).isEmpty)
                let again = r.cell.niggliReduced()
                #expect(parameterDeviation(parameters(again.cell), parameters(r.cell)) <= 1e-9, "\(parameters(base))")
            }
        }
    }

    @Test func testInvalidCellsAreReturnedNotTrapped() {
        let bad = UnitCell(a: 0, b: 1, c: 1, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        let r = bad.niggliReduced()
        #expect(!r.converged)
        #expect(r.transformation == .identity)
        let flat = UnitCell(a: 5, b: 5, c: 5, alphaDeg: 120, betaDeg: 120, gammaDeg: 120)
        #expect(!flat.niggliReduced().converged)
    }
}
