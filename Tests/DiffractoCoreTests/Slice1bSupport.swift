//
//  Slice1bSupport.swift
//  Role: the pins of tools/slice-1b/pins.py (spglib 2.8.0 and pymatgen as oracles for the reduced cell, the lattice's
//        Bravais name and conventional cell, the operation counts and the rings in the conventional setting), and the
//        Niggli conditions written out from the text, so the tests check the reduced cell directly and not only the
//        oracle (registration docs/registration-slice-1b.md sections 2 and 3).
//

import Foundation
import Testing
import simd
import DiffractoCore

struct Pins1b: Decodable {
    let versions: [String: String]
    let generatedAt: String
    let randomSummary: RandomSummary
    let structures: [StructurePin1b]
    let randomCells: [RandomCellPin]

    static let fileName = "slice-1b-pins.json"
    static let shared: Pins1b = {
        do {
            let data = try Data(contentsOf: Fixtures.directory.appendingPathComponent(fileName))
            let decoder = JSONDecoder()
            decoder.keyDecodingStrategy = .convertFromSnakeCase
            return try decoder.decode(Pins1b.self, from: data)
        } catch { fatalError("slice-1b-pins.json: \(error)") }
    }()

    struct RandomSummary: Decodable { let rejectedDraws: Int; let worstOracleDifference: Double }
}

struct LatticePin: Decodable {
    let bravais: String
    let holohedrySpaceGroup: Int
    let holohedryOrder: Int
    let convVolume: Double
    let convParams: [Double]
    let convLengthsSorted: [Double]
    let convUniqueB: Double
}

struct RingPin: Decodable {
    let g: Double
    let multiplicity: Int
    let name: String
}

struct RingsPin: Decodable {
    let kMax: Double
    let count: Int
    let allowedCount: Int
    let families: [RingPin]
}

struct StructurePin1b: Decodable {
    let name: String
    let source: String
    let cell: CellParams
    let sites: [SitePin]
    let niggliSpglib: [Double]
    let niggliPymatgen: [Double]
    let lattice: LatticePin
    let laueClass: String
    let crystalLaueClass: String
    let pointGroup: String
    let spaceGroupNumber: Int
    let fileOperationsCount: Int
    let conventionalOperationsCount: Int
    let convStdParams: [Double]
    let rings: RingsPin

    var atomSites: [AtomSite] {
        sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.frac[0], $0.frac[1], $0.frac[2]), occupancy: $0.occupancy) }
    }
}

struct RandomCellPin: Decodable {
    let intended: String
    let cell: CellParams
    let niggliSpglib: [Double]
    let niggliPymatgen: [Double]
    let lattice: LatticePin
    let unimodular: [[Int]]
}

/// Cell parameters as the six numbers [a, b, c, α, β, γ].
func parameters(_ cell: UnitCell) -> [Double] { [cell.a, cell.b, cell.c, cell.alphaDeg, cell.betaDeg, cell.gammaDeg] }

/// The largest deviation of two parameter sets: relative on the lengths, absolute in degrees on the angles.
func parameterDeviation(_ x: [Double], _ y: [Double]) -> Double {
    var worst = 0.0
    for i in 0..<3 { worst = max(worst, abs(x[i] - y[i]) / max(abs(y[i]), 1)) }
    for i in 3..<6 { worst = max(worst, abs(x[i] - y[i])) }
    return worst
}

/// The cell's own basis combined by the columns of `p` (Å, Cartesian).
func transformedVectors(_ cell: UnitCell, _ p: IntegerMatrix) -> [SIMD3<Double>] {
    (0..<3).map { j in
        Double(p[0, j]) * cell.realLattice[0] + Double(p[1, j]) * cell.realLattice[1] + Double(p[2, j]) * cell.realLattice[2]
    }
}

/// The Niggli main conditions (Křivý & Gruber 1976; International Tables A, section 9.2) on the reduced basis, written
/// out from the text, each within `eps` (Å²). Returns the conditions that fail.
func niggliViolations(_ v: [SIMD3<Double>], eps: Double) -> [String] {
    let a = dot(v[0], v[0]), b = dot(v[1], v[1]), c = dot(v[2], v[2])
    let xi = 2 * dot(v[1], v[2]), eta = 2 * dot(v[0], v[2]), zeta = 2 * dot(v[0], v[1])
    func equal(_ x: Double, _ y: Double) -> Bool { abs(x - y) <= eps }
    var out: [String] = []
    if a > b + eps { out.append("A ≤ B") }
    if b > c + eps { out.append("B ≤ C") }
    if equal(a, b), abs(xi) > abs(eta) + eps { out.append("A = B ⇒ |ξ| ≤ |η|") }
    if equal(b, c), abs(eta) > abs(zeta) + eps { out.append("B = C ⇒ |η| ≤ |ζ|") }
    let typeI = xi > -eps && eta > -eps && zeta > -eps
    let typeII = xi < eps && eta < eps && zeta < eps
    if !typeI && !typeII { out.append("ξ, η, ζ all positive or all non-positive") }
    if abs(xi) > b + eps { out.append("|ξ| ≤ B") }
    if abs(eta) > a + eps { out.append("|η| ≤ A") }
    if abs(zeta) > a + eps { out.append("|ζ| ≤ A") }
    let s = xi + eta + zeta + a + b
    if s < -eps { out.append("ξ + η + ζ + A + B ≥ 0") }
    if equal(xi, b), zeta > 2 * eta + eps { out.append("ξ = B ⇒ ζ ≤ 2η") }
    if equal(eta, a), zeta > 2 * xi + eps { out.append("η = A ⇒ ζ ≤ 2ξ") }
    if equal(zeta, a), eta > 2 * xi + eps { out.append("ζ = A ⇒ η ≤ 2ξ") }
    if equal(xi, -b), zeta < -eps { out.append("ξ = −B ⇒ ζ = 0") }
    if equal(eta, -a), zeta < -eps { out.append("η = −A ⇒ ζ = 0") }
    if equal(zeta, -a), eta < -eps { out.append("ζ = −A ⇒ η = 0") }
    if equal(s, 0), 2 * (a + eta) + zeta > eps { out.append("sum = 0 ⇒ 2A + 2η + ζ ≤ 0") }
    return out
}

/// Pᵀ G P from public numbers only.
func metric(of vectors: [SIMD3<Double>]) -> [[Double]] {
    (0..<3).map { i in (0..<3).map { j in dot(vectors[i], vectors[j]) } }
}

func isUnimodularWithPositiveDeterminant(_ p: IntegerMatrix) -> Bool { p.determinant == 1 }
