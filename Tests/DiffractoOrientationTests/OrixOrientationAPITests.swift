import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/pins/orix/orientation-api-pins.py (orix 0.15.0). Tolerances of orix test_orientation.py: `from_euler` atol 1e-4
/// (`test_from_euler_symmetry`), the rest np.allclose defaults (rtol 1e-5, atol 1e-8) element by element; the measured worst gaps
/// are printed (below 1e-15 when written).
private struct Pins: Decodable {
    struct EulerItem: Decodable { let euler: [Double]; let lab2crystal, crystal2lab: [[Double]] }
    struct Set: Decodable {
        let symbol: String
        let A, B: [[[Double]]]
        let dot: [Double]
        let dotOuter, distanceMatrix: [[Double]]
    }
    let euler: [EulerItem]
    let sets: [Set]
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    var g = 0.0
    for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
    return g
}
/// np.allclose defaults: |x - want| <= 1e-8 + 1e-5 |want|.
private func close(_ x: Double, _ want: Double) -> Bool { abs(x - want) <= 1e-8 + 1e-5 * abs(want) }
private func symmetry(_ symbol: String) -> CrystalOrientationSymmetry {
    if symbol == "-31m" { return CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarOneM) }
    return pinnedSymmetry(symbol)!
}

struct OrixOrientationAPITests {
    private func pins() throws -> Pins {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/orix-orientation-api-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
    }

    @Test func fromEulerMatchesOrixInBothDirections() throws {
        let cubic = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        var worst = 0.0
        for item in try pins().euler {
            let e = SIMD3(item.euler[0], item.euler[1], item.euler[2])
            let lab = try Orientation(euler: e, symmetry: cubic)
            let crystal = try Orientation(euler: e, symmetry: cubic, direction: .crystalToLab)
            #expect(lab.symmetry.laueClass == .mThreeBarM && crystal.symmetry.laueClass == .mThreeBarM)
            let g = max(gap(lab.matrix, mat(item.lab2crystal)), gap(crystal.matrix, mat(item.crystal2lab)))
            worst = max(worst, g)
            #expect(g < 1e-4, Comment(rawValue: "\(item.euler)"))     // orix test_from_euler_symmetry: atol 1e-4
        }
        print("ORIX from_euler worst: \(worst)")
    }

    @Test func fromEulerRefusesNonFiniteAngles() {
        let cubic = CrystalOrientationSymmetry(laueClass: .mThreeBarM)
        #expect(throws: OrientationSymmetryError.invalidRotation) { try Orientation(euler: SIMD3(.nan, 0, 0), symmetry: cubic) }
        #expect(throws: OrientationSymmetryError.invalidRotation) { try Orientation(euler: SIMD3(0, .infinity, 0), symmetry: cubic) }
    }

    @Test func dotDotOuterAndDistanceMatrixMatchOrixForEveryClass() throws {
        var worstDot = 0.0, worstDistance = 0.0
        for set in try pins().sets {
            let s = symmetry(set.symbol)
            let a = try set.A.map { try Orientation(mat($0), symmetry: s) }
            let b = try set.B.map { try Orientation(mat($0), symmetry: s) }
            for (i, want) in set.dot.enumerated() {
                let d = try a[i].dot(b[i])
                worstDot = max(worstDot, abs(d - want))
                #expect(close(d, want), Comment(rawValue: "\(set.symbol) dot \(i)"))
            }
            let outer = try Orientation.dotOuter(a, b)
            #expect(outer.count == a.count && outer.allSatisfy { $0.count == b.count })
            for i in a.indices { for j in b.indices {
                worstDot = max(worstDot, abs(outer[i][j] - set.dotOuter[i][j]))
                #expect(close(outer[i][j], set.dotOuter[i][j]), Comment(rawValue: "\(set.symbol) dotOuter \(i) \(j)"))
            } }
            let distance = try Orientation.distanceMatrix(a)
            #expect(distance.count == a.count && distance.allSatisfy { $0.count == a.count })
            for i in a.indices { for j in a.indices {
                worstDistance = max(worstDistance, abs(distance[i][j] - set.distanceMatrix[i][j]))
                #expect(close(distance[i][j], set.distanceMatrix[i][j]), Comment(rawValue: "\(set.symbol) distance \(i) \(j)"))
            } }
            // A[1] and B[0] are symmetry equivalents of A[0]: the dot is 1 and the rounded, clamped distance exactly 0.
            #expect(distance[0][0] == 0 && distance[0][1] == 0 && distance[1][0] == 0, Comment(rawValue: "\(set.symbol) equivalent distance"))
            let selfDot = try a[0].dot(a[1])
            #expect(abs(outer[0][0] - 1) < 1e-12 && abs(selfDot - 1) < 1e-12, Comment(rawValue: "\(set.symbol) equivalent dot"))
            // The distance is the disorientation angle of `angle(with:)`.
            #expect(close(distance[1][2], try a[1].angle(with: a[2])))
        }
        print("ORIX dot / dot_outer worst: \(worstDot), distance matrix worst: \(worstDistance)")
    }

    @Test func dotRefusesDifferentSymmetries() throws {
        let cubic = try Orientation(matrix_identity_double3x3, symmetry: CrystalOrientationSymmetry(laueClass: .mThreeBarM))
        let hexagonal = try Orientation(matrix_identity_double3x3, symmetry: CrystalOrientationSymmetry(laueClass: .sixOverMmm))
        #expect(throws: OrientationSymmetryError.differentSymmetry) { try cubic.dot(hexagonal) }
        #expect(throws: OrientationSymmetryError.differentSymmetry) { try Orientation.dotOuter([cubic], [hexagonal]) }
        #expect(throws: OrientationSymmetryError.differentSymmetry) { try Orientation.distanceMatrix([cubic, hexagonal]) }
        #expect(try Orientation.distanceMatrix([]).isEmpty)
    }
}
