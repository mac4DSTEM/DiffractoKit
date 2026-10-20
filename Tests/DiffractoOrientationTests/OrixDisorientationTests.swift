import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/orix/disorientation-pins.py (orix 0.15.0). The kit's py4DSTEM matrix is orix's Orientation matrix (symmetry on
/// the left). Tolerances, written per quantity: angles 1e-9 rad (Double arithmetic against orix's float64; arccos near a
/// half-turn is the only conditioning loss), reduced matrix entries 1e-9 where the minimum is unique (next-best angle more
/// than 1e-6 rad away), the Float diagnostic angle of `reduce` 2e-7 rad (Float rounding of an angle below pi).
private struct Pins: Decodable {
    struct Single: Decodable { let matrix, reduced: [[Double]]; let angle: Double; let unique: Bool }
    struct Pair: Decodable { let a, b: [[Double]]; let angle: Double }
    struct Class: Decodable { let symbol: String; let singles: [Single]; let pairs: [Pair] }
    let classes: [Class]
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }

struct OrixDisorientationTests {
    private func pins() throws -> Pins {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .appendingPathComponent("Fixtures/orix-disorientation-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
    }

    @Test func minimumAngleReductionMatchesOrixForEveryClass() throws {
        var worstAngle = 0.0, worstMatrix = 0.0
        for c in try pins().classes {
            let s = try #require(pinnedSymmetry(c.symbol))
            for p in c.singles {
                let r = try s.reduce(mat(p.matrix))
                let da = abs(Double(r.disorientationRad) - p.angle)
                worstAngle = max(worstAngle, da)
                #expect(da < 2e-7, "\(c.symbol) angle")
                if p.unique {
                    let e = mat(p.reduced)
                    var g = 0.0
                    for col in 0..<3 { for row in 0..<3 { g = max(g, abs(r.matrix[col][row] - e[col][row])) } }
                    worstMatrix = max(worstMatrix, g)
                    #expect(g < 1e-9, "\(c.symbol) matrix")
                }
            }
        }
        print("ORIX reduce worst: Float angle \(worstAngle) rad, matrix \(worstMatrix)")
    }

    @Test func disorientationOfPairsMatchesOrixForEveryClass() throws {
        var worst = 0.0
        for c in try pins().classes {
            let s = try #require(pinnedSymmetry(c.symbol))
            for p in c.pairs {
                let d = try s.disorientation(from: mat(p.a), to: mat(p.b))
                let gap = abs(d.angleRad - p.angle)
                worst = max(worst, gap)
                #expect(gap < 1e-9, "\(c.symbol)")
            }
        }
        print("ORIX disorientation worst angle gap: \(worst) rad")
    }

    @Test func disorientationIsSymmetricAndZeroForEquivalentOrientations() throws {
        let all = LaueClass.allCases.map { CrystalOrientationSymmetry(laueClass: $0) }
            + [CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarMOne)]
        for s in all {
            let laue = s.laueClass
            let a = simd_double3x3(simd_quatd(angle: 0.9, axis: simd_normalize(SIMD3(1, 2, 3))))
            for op in s.operators {
                let d = try s.disorientation(from: a, to: op * a)
                #expect(d.angleRad < 1e-7, "\(laue)")
                let trace = d.matrix[0][0] + d.matrix[1][1] + d.matrix[2][2]
                #expect(abs(acos(min(1, max(-1, (trace - 1) / 2))) - d.angleRad) < 1e-6, "\(laue) matrix angle")
            }
            let b = simd_double3x3(simd_quatd(angle: 2.1, axis: simd_normalize(SIMD3(-1, 0.5, 2))))
            let ab = try s.disorientation(from: a, to: b).angleRad, ba = try s.disorientation(from: b, to: a).angleRad
            #expect(abs(ab - ba) < 1e-9)
        }
        let s = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        #expect(throws: OrientationSymmetryError.invalidRotation) {
            try s.disorientation(from: .init(0), to: matrix_identity_double3x3)
        }
    }
}
