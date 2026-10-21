import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/symmetry/family-pins.py (py4DSTEM f050d207 `symmetry_reduce_directions`, 0.14.19, NumPy 2.5.0). The
/// written tolerance is 1e-7 absolute (D17; was 1e-12) on family entries (the operator products are float64 sums of three terms, entries of
/// order 1; NumPy's BLAS and the Swift loops differ at 1e-16). Where more than one operator satisfies the in-range test
/// the pin records it (`hits`); those cases are compared too because both implementations take the first operator in the
/// same order. Columns on a sector wall (`boundary`) are not compared: whether the exact `>= 0` test passes is rounding noise.
private struct Pins: Decodable {
    struct Case: Decodable { let matrix, family: [[Double]]; let hits: [Int]; let boundary: [Bool] }
    struct Class: Decodable { let symbol: String; let operatorCount: Int; let cases: [Case] }
    let classes: [Class]
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }

struct Py4DSTEMFamilyTests {
    private func pins() throws -> Pins {
        let file = FixtureDirs.target("DiffractoOrientationTests")
            .appendingPathComponent("Fixtures/py4dstem-family-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
    }

    @Test func directionFamilyMatchesPy4DSTEMForEveryClass() throws {
        var worst = 0.0, count = 0, boundaryColumns = 0
        for c in try pins().classes {
            let s = CrystalOrientationSymmetry(laueClass: try #require(LaueClass(rawValue: c.symbol)))
            for k in c.cases {
                let f = try s.py4DSTEMDirectionFamily(mat(k.matrix), order: .sortedFractional)
                let e = mat(k.family)
                var g = 0.0
                for col in 0..<3 {
                    // A column within 1e-12 of a sector wall is decided by BLAS rounding in NumPy: skipped here, counted below.
                    if k.boundary[col] { boundaryColumns += 1; continue }
                    for row in 0..<3 { g = max(g, abs(f[col][row] - e[col][row])) }
                }
                worst = max(worst, g)
                count += 1
                #expect(g < 1e-7, "\(c.symbol) \(k.hits) \(k.matrix) got \(f) expected \(k.family)")
            }
        }
        print("PY4DSTEM family: \(count) matrices, worst entry gap \(worst), \(boundaryColumns) boundary columns not compared")
    }

    @Test func familyIsNotTheMinimumAngleRotationAndAcceptsNonRotations() throws {
        let s = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        let m = simd_double3x3(simd_quatd(angle: 1.1, axis: simd_normalize(SIMD3(1, 2, 3))))
        let family = try s.py4DSTEMDirectionFamily(m)
        // Columns are folded independently, so the family is generally not orthonormal.
        let gram = family.transpose * family
        #expect(abs(gram[0][1]) > 1e-6 || abs(gram[0][2]) > 1e-6 || abs(gram[1][2]) > 1e-6)
        #expect(throws: OrientationSymmetryError.invalidRotation) {
            try s.py4DSTEMDirectionFamily(.init(.nan))
        }
    }
}
