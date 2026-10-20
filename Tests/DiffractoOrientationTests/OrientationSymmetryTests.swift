import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

private struct SymmetryPins: Decodable {
    struct Row: Decodable {
        let matrix, reduced: [[Double]]
        let angle: Double
    }
    struct Group: Decodable {
        let symbol: String
        let operators: [[[Double]]]
        let cases: [Row]
        let foldedDirections: [[Double]]
        let samples: [Native.Sample]
    }
    struct Native: Decodable {
        struct Sample: Decodable { let count: Int; let directions: [[Double]] }
        let symbol: String
        let colors: [[UInt32]]
        let directions: [[Double]]
        let samples: [Sample]
    }
    struct Python: Decodable {
        struct Record: Decodable {
            let symbol: String
            let directions: [[Double]]
            let bits: [[UInt32]]
            let rgb: [[Double]]
        }
        let records: [Record]
    }
    let groups: [Group]
    let native: [Native]
    let nativeDirections: [[Double]]
    let py4dstem: Python
}

private func matrix(_ r: [[Double]]) -> simd_double3x3 {
    .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) })
}
private func vector(_ a: [Double]) -> SIMD3<Double> { SIMD3(a[0], a[1], a[2]) }
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    (0..<3).flatMap { c in (0..<3).map { abs(a[c][$0] - b[c][$0]) } }.max()!
}
private func pins() throws -> SymmetryPins {
    let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/symmetry-pins.json")
    return try JSONDecoder().decode(SymmetryPins.self, from: Data(contentsOf: file))
}

struct OrientationSymmetryTests {
    @Test func allElevenProperGroupsMatchPymatgenAndAreClosed() throws {
        let fixture = try pins()
        #expect(fixture.groups.count == LaueClass.allCases.count)
        var worst = 0.0
        for p in fixture.groups {
            let laue = try #require(LaueClass(rawValue: p.symbol))
            let ops = CrystalOrientationSymmetry(laueClass: laue).operators
            #expect(ops.count == laue.order / 2)
            #expect(p.operators.count == ops.count, "\(p.symbol): pymatgen proper-operation count")
            for op in ops {
                #expect(abs(simd_determinant(op) - 1) < 2e-14)
                #expect(gap(op.transpose * op, matrix_identity_double3x3) < 2e-14)
                let g = p.operators.map { gap(op, matrix($0)) }.min()!
                worst = max(worst,g)
                #expect(g < 2e-14)
                for other in ops { #expect(ops.contains { gap(op * other, $0) < 2e-14 }) }
            }
        }
        print("SYMMETRY group worst absolute gap: \(worst)")
    }

    @Test func seededRotationsMatchMinimumAngleOracleForEveryLaueClass() throws {
        var worstMatrix = 0.0, worstAngle = 0.0
        for p in try pins().groups {
            let s = CrystalOrientationSymmetry(laueClass: try #require(LaueClass(rawValue:p.symbol)))
            for row in p.cases {
                let input = matrix(row.matrix)
                let result = try s.reduce(input)
                let m = gap(result.matrix, matrix(row.reduced))
                let a = abs(Double(result.disorientationRad) - row.angle)
                worstMatrix = max(worstMatrix,m); worstAngle = max(worstAngle,a)
                #expect(m < 2e-12)
                // Float angle storage; 2e-7 rad covers a Float ulp at pi.
                #expect(a < 2e-7)
                #expect(gap(try s.reduce(result.matrix).matrix, result.matrix) < 2e-12)
                for op in s.operators {
                    let equivalent = try s.reduce(op * input)
                    #expect(gap(equivalent.matrix,result.matrix) < 2e-12)
                }
                let theta = Double(result.disorientationRad)
                for op in s.operators {
                    let q = simd_quatd((op * input).transpose)
                    let alternative = 2 * acos(min(1,abs(q.real)))
                    #expect(theta <= alternative + 2e-7)
                }
            }
        }
        print("SYMMETRY reduction matrix worst: \(worstMatrix); angle worst rad: \(worstAngle)")
    }

    @Test func nativeColorsAndDirectionFoldingArePinnedToOriginalSource() throws {
        let p = try pins()
        var worst = 0.0
        for row in p.native {
            let s = CrystalOrientationSymmetry(laueClass: try #require(LaueClass(rawValue:row.symbol)))
            for (i,d) in p.nativeDirections.enumerated() {
                let c = try s.ipfColor(direction:vector(d),key:.native)
                #expect([c.x.bitPattern,c.y.bitPattern,c.z.bitPattern] == row.colors[i])
                let reduced = try s.reduceDirection(vector(d))
                let g = simd_length(reduced-vector(row.directions[i])); worst = max(worst,g)
                #expect(g < 2e-14)
            }
            for sample in row.samples {
                let actual = try s.sampleFundamentalZone(count:sample.count)
                #expect(actual.count == sample.directions.count)
                for (v,e) in zip(actual,sample.directions) {
                    #expect([v.x.bitPattern,v.y.bitPattern,v.z.bitPattern] == e.map(\.bitPattern))
                }
            }
        }
        print("SYMMETRY native direction worst: \(worst)")
    }

    @Test func pythonColorsEqualPinnedPlotOrientationMapsAtFloatPrecisionBitForBit() throws {
        var worst = 0.0
        for p in try pins().py4dstem.records {
            let s = CrystalOrientationSymmetry(laueClass:try #require(LaueClass(rawValue:p.symbol)))
            for (i,d) in p.directions.enumerated() {
                let c = try s.ipfColor(direction:vector(d),key:.py4DSTEM)
                #expect([c.x.bitPattern,c.y.bitPattern,c.z.bitPattern] == p.bits[i], "\(p.symbol), row \(i)")
                worst = max(worst,zip([Double(c.x),Double(c.y),Double(c.z)],p.rgb[i]).map { abs($0-$1) }.max()!)
            }
        }
        print("SYMMETRY Python Double to public Float RGB worst: \(worst)")
    }

    @Test func samplingCoversDistinctDirectionOrbitsForAllClasses() throws {
        for laue in LaueClass.allCases {
            let s = CrystalOrientationSymmetry(laueClass:laue)
            let sample = try s.sampleFundamentalZone(count:16)
            #expect(sample.count == 16)
            #expect(sample == (try s.sampleFundamentalZone(count:16)))
            for (i,d) in sample.enumerated() {
                #expect(abs(simd_length(d)-1) < 2e-14)
                #expect(simd_length(try s.reduceDirection(d)-d) < 2e-12)
                for op in s.operators {
                    #expect(simd_length(try s.reduceDirection(op*d)-d) < 2e-12)
                    #expect(simd_length(try s.reduceDirection(-op*d)-d) < 2e-12)
                    for other in sample.prefix(i) {
                        #expect(abs(simd_dot(op*d,other)) < 1-1e-12)
                    }
                }
            }
        }
    }

    @Test func remainingEightSectorsAndSamplesArePinnedUsingPymatgenOrbits() throws {
        let p = try pins()
        var worst = 0.0
        for group in p.groups where !group.samples.isEmpty {
            let symmetry = CrystalOrientationSymmetry(laueClass:try #require(LaueClass(rawValue:group.symbol)))
            for (d,e) in zip(p.nativeDirections,group.foldedDirections) {
                let g = simd_length(try symmetry.reduceDirection(vector(d))-vector(e))
                worst = max(worst,g); #expect(g < 2e-12)
            }
            for sample in group.samples {
                let actual = try symmetry.sampleFundamentalZone(count:sample.count)
                #expect(actual.count == sample.directions.count)
                for (v,e) in zip(actual,sample.directions) {
                    let g = simd_length(v-vector(e)); worst = max(worst,g)
                    #expect(g < 2e-12)
                }
            }
        }
        print("SYMMETRY eight new sectors/samples worst: \(worst)")
    }

    @Test func invalidInputsAreExplicitAndEveryClassHasANativeKey() throws {
        let s = CrystalOrientationSymmetry(laueClass:.mThreeBarM)
        #expect(throws:OrientationSymmetryError.invalidRotation) { try s.reduce(.init(0)) }
        #expect(throws:OrientationSymmetryError.invalidRotation) { try s.reduce(-matrix_identity_double3x3) }
        #expect(throws:OrientationSymmetryError.invalidDirection) { try s.reduceDirection(.zero) }
        #expect(throws:OrientationSymmetryError.invalidDirection) { try s.reduceDirection([.nan,0,1]) }
        #expect(throws:OrientationSymmetryError.invalidSampleCount) { try s.sampleFundamentalZone(count:Int.max) }
        #expect(throws:OrientationSymmetryError.invalidSampleCount) { try s.sampleFundamentalZone(count:-1) }
        #expect(try s.sampleFundamentalZone(count:0).isEmpty)
        // 0.6.0: the other eight classes have a native key too (the orix key); nothing throws nativeKeyUnavailable.
        for laue in LaueClass.allCases where ![.mThreeBarM,.sixOverMmm,.twoOverM].contains(laue) {
            let s = CrystalOrientationSymmetry(laueClass:laue)
            #expect(try s.ipfColor(direction:[1,0.3,0.7],key:.native) == s.ipfColor(direction:[1,0.3,0.7],key:.orix))
        }
    }
}
