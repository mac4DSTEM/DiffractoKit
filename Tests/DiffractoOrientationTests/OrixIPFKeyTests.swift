import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/orix/ipf-pins.py (orix 0.15.0, DirectionColorKeyTSL). Tolerance: the Swift key is Double arithmetic
/// transcribed from orix; the public colour is Float, so the written tolerance is 1e-5 absolute (orix's allclose; was 1e-7) on RGB in [0, 1]
/// (Float rounding is 6e-8; the rest is the order of floating-point operations in the 1000-point azimuth table).
private struct OrixPins: Decodable {
    struct Class: Decodable {
        let symbol: String
        let directions: [[Double]]
        let rgb: [[Double]]
    }
    let classes: [Class]
}

struct OrixIPFKeyTests {
    private func pins() throws -> OrixPins {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/orix-ipf-pins.json")
        return try JSONDecoder().decode(OrixPins.self, from: Data(contentsOf: file))
    }

    @Test func orixKeyMatchesOrixForAllElevenLaueClassesAndBothMinus3mSettings() throws {
        let fixture = try pins()
        #expect(fixture.classes.count == 12 && fixture.classes.last?.symbol == "-3m1")
        var worst = 0.0
        for c in fixture.classes {
            let s = try #require(pinnedSymmetry(c.symbol))
            for (d, e) in zip(c.directions, c.rgb) {
                let rgb = try s.ipfColor(direction: SIMD3(d[0], d[1], d[2]), key: .orix)
                let gap = max(abs(Double(rgb.x) - e[0]), abs(Double(rgb.y) - e[1]), abs(Double(rgb.z) - e[2]))
                worst = max(worst, gap)
                #expect(gap < 1e-5, "\(c.symbol) \(d): \(rgb) vs \(e)")  // orix's np.allclose defaults (rtol 1e-5, atol 1e-8); was 1e-7
            }
        }
        print("ORIX IPF worst RGB gap: \(worst)")
    }

    @Test func nativeKeyOfTheEightOtherClassesIsTheOrixKey() throws {
        for c in try pins().classes {
            let s = pinnedSymmetry(c.symbol)!
            guard ![.mThreeBarM, .sixOverMmm, .twoOverM].contains(s.laueClass) else { continue }
            for d in c.directions.prefix(20) {
                let v = SIMD3(d[0], d[1], d[2])
                #expect(try s.ipfColor(direction: v, key: .native) == s.ipfColor(direction: v, key: .orix))
            }
        }
    }

    @Test func symmetricallyEquivalentPolesGetTheSameColour() throws {
        let all = LaueClass.allCases.map { CrystalOrientationSymmetry(laueClass: $0) }
            + [CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarMOne)]
        for s in all {
            let laue = s.laueClass
            let d = SIMD3<Double>(0.31, 0.52, 0.79)
            let reference = try s.ipfColor(direction: d, key: .orix)
            for op in s.operators {
                for sign in [-1.0, 1.0] {
                    let c = try s.ipfColor(direction: sign * (op * d), key: .orix)
                    #expect(simd_length(SIMD3<Double>(c) - SIMD3<Double>(reference)) < 1e-5, "\(laue)")
                }
            }
        }
    }
}

extension OrixIPFKeyTests {
    /// Records how far the app-derived native keys of m-3m, 6/mmm and 2/m are from orix's key (m-3m is documented as
    /// "aligned"; the other two are different conventions). Only the measured gaps are written down.
    @Test func nativeKeysOfTheThreeOriginalClassesVersusOrix() throws {
        var out: [String] = []
        for c in try pins().classes {
            let s = pinnedSymmetry(c.symbol)!
            let laue = s.laueClass
            guard [.mThreeBarM, .sixOverMmm, .twoOverM].contains(laue) else { continue }
            var worst: Float = 0
            for d in c.directions {
                let v = SIMD3(d[0], d[1], d[2])
                let a = try s.ipfColor(direction: v, key: .native), b = try s.ipfColor(direction: v, key: .orix)
                worst = max(worst, abs(a.x - b.x), abs(a.y - b.y), abs(a.z - b.z))
            }
            out.append("\(c.symbol): \(worst)")
            // Upper bound for the construction shared with orix, lower bound for the two different conventions.
            if laue == .mThreeBarM { #expect(worst < 1e-5) } else { #expect(worst > 0.1, "\(c.symbol)") }
        }
        print("NATIVE versus ORIX worst RGB gaps: \(out.joined(separator: "; "))")
    }
}
