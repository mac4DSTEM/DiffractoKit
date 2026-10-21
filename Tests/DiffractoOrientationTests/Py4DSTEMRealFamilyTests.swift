import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/symmetry/family-real-pins.py. Real py4DSTEM `Crystal`s (0.14.19 at f050d207, pymatgen 2025.1.9 and spglib 2.8.0 in
/// its venv) built from pymatgen-expanded structures of centrosymmetric space groups; `orientation_plan('auto')` takes the plan sector
/// from the detected point group and `symmetry_operators` from `SpacegroupAnalyzer.get_point_group_operations()`, so the operator
/// order is pymatgen's. Written tolerance: 1e-7 absolute (D17; was 1e-12) on family entries (float64 sums of three products of order 1). Columns
/// within 1e-12 of a sector wall (`boundary`) are decided by BLAS rounding in NumPy and are not compared.
private struct RealPins: Decodable {
    struct Lattice: Decodable { let a, b, c, alpha, beta, gamma: Double }
    struct Case: Decodable { let matrix, family: [[Double]]; let hits: [Int]; let boundary: [Bool] }
    struct Crystal: Decodable {
        let `class`: String; let spaceGroup: Int; let symbol: String; let primary: Bool
        let orderFractional: [[Int]]; let operators: [[[Double]]]; let cases: [Case]
    }
    let cases: [Crystal]
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }

struct Py4DSTEMRealFamilyTests {
    private func pins() throws -> [RealPins.Crystal] {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .appendingPathComponent("Fixtures/py4dstem-family-real-pins.json")
        return try JSONDecoder().decode(RealPins.self, from: Data(contentsOf: file)).cases
    }

    private func gap(_ s: CrystalOrientationSymmetry, _ k: RealPins.Case, _ order: Py4DSTEMOperatorOrder) throws -> Double {
        let f = try s.py4DSTEMDirectionFamily(mat(k.matrix), order: order)
        let e = mat(k.family)
        var g = 0.0
        for col in 0..<3 where !k.boundary[col] {
            for row in 0..<3 { g = max(g, abs(f[col][row] - e[col][row])) }
        }
        return g
    }

    @Test func directionFamilyMatchesRealCrystalsForEverySpaceGroup() throws {
        var worst = 0.0, matrices = 0, crystals = 0
        for crystal in try pins() {
            let s = try #require(pinnedSymmetry(crystal.class))
            for k in crystal.cases {
                let g = try gap(s, k, .pymatgen)
                worst = max(worst, g); matrices += 1
                #expect(g < 1e-7, "\(crystal.class) \(crystal.symbol) \(k.matrix): gap \(g)")
            }
            crystals += 1
        }
        #expect(crystals == 57)
        print("PY4DSTEM real crystals: \(crystals) crystals, \(matrices) matrices, worst family gap \(worst)")
    }

    @Test func theHarnessOrderDisagreesWithRealCrystalsWhereTheFallbackOperatorDiffers() throws {
        // The sorted harness list puts a twofold first, a real crystal's list puts the identity first: for 2/m and -3m (-31m), whose
        // plan sector is not a fundamental domain, a direction no operator maps into the sector keeps the first operator.
        var bad: [String: Int] = [:]
        for crystal in try pins() where crystal.primary {
            let s = try #require(pinnedSymmetry(crystal.class))
            for k in crystal.cases where try gap(s, k, .sortedFractional) >= 1e-12 { bad[crystal.class, default: 0] += 1 }
        }
        #expect(bad["2/m", default: 0] == 139 && bad["-3m", default: 0] == 136, "\(bad)")
        #expect(Set(bad.keys) == ["2/m", "-3m"], "every other class agrees under both orders: \(bad)")
        print("PY4DSTEM harness order versus real crystals, mismatching matrices of 151: \(bad.sorted { $0.key < $1.key })")
    }
}
