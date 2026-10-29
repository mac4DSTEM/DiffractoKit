//
//  ConventionalCellTests.swift
//  The conventional cell: the Bravais lattice
//  name against spglib's lattice-only dataset and against the lattice each random cell was BUILT as, the cell's volume
//  and lengths, the transformation, and closed-form cells of the cubic lattices.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct ConventionalCellTests {

    static let pins = Pins1b.shared

    /// The name is the lattice each cell was built as (truth) and spglib's (oracle), on 2800 random cells under random
    /// unimodular transformations, 200 per lattice.
    @Test func testBravaisNameOfRandomCellsAgreesWithConstructionAndSpglib() {
        var tally: [String: Int] = [:]
        for (i, pin) in Self.pins.randomCells.enumerated() {
            let conv = pin.cell.unitCell.conventionalCell()
            tally[conv.bravais.rawValue, default: 0] += 1
            #expect(conv.bravais.rawValue == pin.intended, "random #\(i): built as \(pin.intended), got \(conv.bravais) (\(conv.setting))")
            #expect(conv.bravais.rawValue == pin.lattice.bravais, "random #\(i) vs spglib \(pin.lattice.bravais)")
        }
        for name in BravaisLattice.allCases { #expect(tally[name.rawValue] == 200, "\(name): \(tally[name.rawValue] ?? 0)") }
    }

    /// Volume, lengths and the transformation on every pinned cell and random cell: the conventional cell has spglib's
    /// volume, the same lengths (the unique b for a monoclinic lattice, whose a and c are a choice), a right-handed P of
    /// determinant equal to the centring multiplicity whose columns reproduce the cell's own metric.
    @Test func testConventionalCellAgreesWithSpglibAndTheTransformationReproducesIt() {
        struct Case { let name: String; let cell: UnitCell; let lattice: LatticePin }
        let cases = Self.pins.structures.map { Case(name: $0.name, cell: $0.cell.unitCell, lattice: $0.lattice) }
            + Self.pins.randomCells.enumerated().map { Case(name: "random #\($0)", cell: $1.cell.unitCell, lattice: $1.lattice) }
        var worstVolume = 0.0, worstLength = 0.0, worstMetric = 0.0
        for c in cases {
            let conv = c.cell.conventionalCell()
            let p = conv.transformation
            #expect(conv.bravais.rawValue == c.lattice.bravais, "\(c.name)")
            #expect(p.determinant == conv.bravais.centringMultiplicity, "\(c.name): det P \(p.determinant), \(conv.bravais)")
            #expect(conv.centringTranslations.count == conv.bravais.centringMultiplicity, "\(c.name)")
            #expect(conv.centringTranslations.first == SIMD3(0, 0, 0), "\(c.name)")
            let volumeDeviation = abs(conv.cell.volume - c.lattice.convVolume) / c.lattice.convVolume
            worstVolume = max(worstVolume, volumeDeviation)
            #expect(volumeDeviation <= 1e-7, "\(c.name) volume \(conv.cell.volume) vs \(c.lattice.convVolume)")
            // the transformation reproduces the conventional metric and volume from the INPUT cell
            let v = transformedVectors(c.cell, p)
            let g = metric(of: v)
            let want = conv.cell.metricTensor
            let scale = max(want[0].x, want[1].y, want[2].z)
            for i in 0..<3 { for j in 0..<3 {
                let deviation = abs(g[i][j] - want[i][j]) / scale
                worstMetric = max(worstMetric, deviation)
                #expect(deviation <= 1e-7, "\(c.name) Pᵀ G P [\(i)][\(j)] \(g[i][j]) vs \(want[i][j])")
            } }
            #expect(abs(c.cell.volume * Double(p.determinant) - conv.cell.volume) <= 1e-7 * conv.cell.volume, "\(c.name)")
            // lengths
            switch conv.bravais {
            case .mP, .mC:
                let d = abs(conv.cell.b - c.lattice.convUniqueB) / c.lattice.convUniqueB
                worstLength = max(worstLength, d)
                #expect(d <= 1e-10, "\(c.name) unique b")
                #expect(conv.cell.alphaDeg == 90 && conv.cell.gammaDeg == 90 && conv.cell.betaDeg >= 90 - 1e-9, "\(c.name)")
                #expect(conv.uniqueAxis == .b)
            case .aP: break
            default:
                let mine = [conv.cell.a, conv.cell.b, conv.cell.c].sorted()
                for k in 0..<3 {
                    let d = abs(mine[k] - c.lattice.convLengthsSorted[k]) / mine[k]
                    worstLength = max(worstLength, d)
                    #expect(d <= 1e-10, "\(c.name) length \(k): \(mine) vs \(c.lattice.convLengthsSorted)")
                }
            }
            // the conventional cell is idealised to its family exactly
            switch conv.bravais {
            case .cP, .cI, .cF: #expect(conv.cell.a == conv.cell.b && conv.cell.b == conv.cell.c && conv.cell.alphaDeg == 90, "\(c.name)")
            case .hP, .hR: #expect(conv.cell.a == conv.cell.b && conv.cell.gammaDeg == 120, "\(c.name)")
            case .tP, .tI: #expect(conv.cell.a == conv.cell.b && conv.cell.gammaDeg == 90, "\(c.name)")
            default: break
            }
        }
        print("conventional cell: \(cases.count) cells, worst volume \(worstVolume), length \(worstLength), metric \(worstMetric)")
    }

    /// Closed forms: the primitive cells of the cubic lattices have the conventional cell edge a, F and I with 4 and 2
    /// lattice points; a hexagonal-axes rhombohedral lattice has a = 2 a_r sin(α/2), c = a_r √(3 + 6 cos α).
    @Test func testClosedFormConventionalCells() {
        let a = 4.0782
        let fcc = UnitCell.rhombohedral(a: a / 2.0.squareRoot(), alphaDeg: 60).conventionalCell()
        #expect(fcc.bravais == .cF)
        #expect(abs(fcc.cell.a - a) <= 1e-12 * a && fcc.centringMultiplicity == 4)
        #expect(fcc.transformation.determinant == 4)
        let cosBcc = -1.0 / 3.0
        let bcc = UnitCell.rhombohedral(a: 2.8665 * 3.0.squareRoot() / 2, alphaDeg: acos(cosBcc) * 180 / .pi).conventionalCell()
        #expect(bcc.bravais == .cI)
        #expect(abs(bcc.cell.a - 2.8665) <= 1e-9 && bcc.centringMultiplicity == 2)
        let simple = UnitCell.rhombohedral(a: 3, alphaDeg: 90).conventionalCell()
        #expect(simple.bravais == .cP && simple.transformation == .identity)
        let ar = 5.1283, alpha = 55.28
        let rh = UnitCell.rhombohedral(a: ar, alphaDeg: alpha).conventionalCell()
        #expect(rh.bravais == .hR)
        let ca = cos(alpha * .pi / 180)
        #expect(abs(rh.cell.a - 2 * ar * sin(alpha * .pi / 360)) <= 1e-9, "\(rh.cell.a)")
        #expect(abs(rh.cell.c - ar * (3 + 6 * ca).squareRoot()) <= 1e-9, "\(rh.cell.c)")
        #expect(rh.cell.gammaDeg == 120 && rh.centringMultiplicity == 3)
        // the obverse setting: the lattice points are (2/3, 1/3, 1/3) and (1/3, 2/3, 2/3)
        let points = Set(rh.centringTranslations.map { [Int(($0.x * 3).rounded()), Int(($0.y * 3).rounded()), Int(($0.z * 3).rounded())] })
        #expect(points == [[0, 0, 0], [2, 1, 1], [1, 2, 2]], "\(points)")
        // a C-centred monoclinic lattice: the primitive cell of a = 15.16, b = 4.05, c = 6.74, β = 105.3 (β'')
        let beta = UnitCell.monoclinic(a: 15.16, b: 4.05, c: 6.74, betaDeg: 105.3)
        let cPrimitive = Self.primitiveOfC(beta)
        let mc = cPrimitive.conventionalCell()
        #expect(mc.bravais == .mC)
        #expect(abs(mc.cell.a - 15.16) <= 1e-9 && abs(mc.cell.b - 4.05) <= 1e-9 && abs(mc.cell.c - 6.74) <= 1e-9
                && abs(mc.cell.betaDeg - 105.3) <= 1e-8, "\(parameters(mc.cell))")
    }

    /// The primitive cell of the C-centred cell: (a + b)/2, (−a + b)/2, c.
    static func primitiveOfC(_ conv: UnitCell) -> UnitCell {
        let a = conv.realLattice[0], b = conv.realLattice[1], c = conv.realLattice[2]
        let v1 = 0.5 * (a + b), v2 = 0.5 * (b - a)
        func angle(_ x: SIMD3<Double>, _ y: SIMD3<Double>) -> Double {
            acos(dot(x, y) / (length(x) * length(y))) * 180 / .pi
        }
        return UnitCell(a: length(v1), b: length(v2), c: length(c), alphaDeg: angle(v2, c), betaDeg: angle(v1, c),
                        gammaDeg: angle(v1, v2))
    }

    /// A cell that is already conventional maps to itself (P = I): the file's own setting is not relabelled.
    @Test func testConventionalCellsMapToThemselves() {
        let cells: [UnitCell] = [.cubic(a: 4.0782), .hexagonal(a: 3.2094, c: 5.2108), .tetragonal(a: 3, c: 5),
                                 .orthorhombic(a: 3, b: 4, c: 5), .orthorhombic(a: 5, b: 4, c: 3),
                                 .monoclinic(a: 5, b: 4, c: 6, betaDeg: 100)]
        for cell in cells {
            let conv = cell.conventionalCell()
            #expect(conv.transformation == .identity, "\(parameters(cell)): \(conv.transformation) \(conv.setting)")
            #expect(conv.cell == cell)
        }
    }

    @Test func testTheTolerancePrincipleAndDegenerateInput() {
        // a cell 1e-7 off cubic is not cubic at 1e-9 and is at 1e-6
        let near = UnitCell(a: 4, b: 4, c: 4 * (1 + 1e-7), alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        #expect(near.conventionalCell().bravais == .tP)
        #expect(near.conventionalCell(tolerance: 1e-6).bravais == .cP)
        // an invalid cell is returned, not trapped
        let bad = UnitCell(a: -1, b: 1, c: 1, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        #expect(bad.conventionalCell().bravais == .aP)
    }
}
