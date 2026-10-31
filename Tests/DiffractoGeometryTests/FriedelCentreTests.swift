import DiffractoKitTestSupport
import Foundation
import Testing
@testable import DiffractoGeometry

/// The centre of a pattern whose beam is blocked (a beam stop with and without its arm, a saturated centre), against py4DSTEM
/// f050d207 `get_origin_friedel` (upsample_factor 1) at assert_allclose float64 defaults, rtol 1e-7 (the reference states no
/// tolerance and has no test of the function), and against the planted truth. Fixture: tools/pins/beam-centre/friedel-pins.py.
struct FriedelCentreTests {
    struct Case: Decodable {
        let name: String
        let rows: Int, columns: Int
        let truthRow: Double, truthColumn: Double
        let pixels: [Double]
        let pixelSum: Double
        let py4dstemRow: Double, py4dstemColumn: Double
        let mask: [Int]?

        var pattern: DiffractionPattern { DiffractionPattern(width: columns, height: rows, pixels: pixels.map { Float($0) }) }
        var boolMask: [Bool]? { mask?.map { $0 != 0 } }
    }
    struct Pins: Decodable { let cases: [Case] }

    static let pins: Pins = {
        let url = FixtureDirs.fixtures(of: "DiffractoGeometryTests").appendingPathComponent("py4dstem-friedel-centre-pins.json")
        // swiftlint:disable:next force_try
        return try! JSONDecoder().decode(Pins.self, from: Data(contentsOf: url))
    }()

    static func close(_ a: Double, _ b: Double, rtol: Double = 1e-7) -> Bool { abs(a - b) <= rtol * abs(b) }

    @Test func theFixtureIsReadExactly() {
        #expect(Self.pins.cases.count == 6)
        for c in Self.pins.cases {
            #expect(c.pixels.count == c.rows * c.columns && (c.mask?.count ?? c.pixels.count) == c.pixels.count)
            // the exact sum of the Float32 values (the generator proves a float64 sum of them is exact in any order)
            let sum = c.pattern.pixels.reduce(0.0) { $0 + Double($1) }
            #expect(sum == c.pixelSum, Comment(rawValue: c.name))
        }
    }

    /// With the beam blocked the brightest feature is a Bragg spot (stop and arm, 7.5 px off) or the beam stop rim (stop disc,
    /// about 6 px off), not the centre: the beam finder refuses there rather than answer. An
    /// open beam is the brightest region, centred on the pattern, and is found. The saturated centre was refused before the Friedel
    /// check too, and may be: what is never allowed is a wrong number.
    @Test func theBeamFinderRefusesABlockedBeam() throws {
        for c in Self.pins.cases {
            let fit = PatternCentre.find(pattern: c.pattern)
            if c.name.hasPrefix("stop-") {
                #expect(fit == nil, Comment(rawValue: "\(c.name): \(String(describing: fit)) instead of a refusal"))
            } else if let found = fit {
                let miss = hypot(found.row - c.truthRow, found.column - c.truthColumn)
                #expect(miss < 0.5, Comment(rawValue: "\(c.name): \(found) is \(miss) px from the truth"))
            }
        }
        let open = try #require(Self.pins.cases.first { $0.name == "open-beam" })
        #expect(PatternCentre.find(pattern: open.pattern) != nil)
    }

    @Test func theFriedelCentreMatchesThePinnedReference() throws {
        for c in Self.pins.cases {
            let fit = try #require(FriedelCentre.find(pattern: c.pattern, mask: c.boolMask), Comment(rawValue: c.name))
            #expect(Self.close(fit.row, c.py4dstemRow) && Self.close(fit.column, c.py4dstemColumn),
                    Comment(rawValue: "\(c.name): \(fit.row), \(fit.column) vs py4DSTEM \(c.py4dstemRow), \(c.py4dstemColumn)"))
            #expect(fit.validation == "none")
        }
    }

    /// Agreement is not truth: on these plants the reference itself lands within 0.1 px of the planted centre (0.02 px with
    /// the mask), which is what makes the method worth porting.
    @Test func theFriedelCentreRecoversThePlantedCentre() throws {
        for c in Self.pins.cases {
            let fit = try #require(FriedelCentre.find(pattern: c.pattern, mask: c.boolMask))
            let miss = hypot(fit.row - c.truthRow, fit.column - c.truthColumn)
            #expect(miss < (c.mask == nil ? 0.1 : 0.02), Comment(rawValue: "\(c.name): \(miss) px"))
        }
    }

    @Test func badInputIsRefused() throws {
        let c = try #require(Self.pins.cases.first)
        let p = c.pattern
        #expect(FriedelCentre.find(pattern: DiffractionPattern(width: 1, height: 1, pixels: [1])) == nil)
        #expect(FriedelCentre.find(pattern: DiffractionPattern(width: 8, height: 8, pixels: [Float](repeating: 1, count: 63))) == nil)
        #expect(FriedelCentre.find(pattern: p, mask: [Bool](repeating: true, count: 3)) == nil)
        #expect(FriedelCentre.find(pattern: p, mask: [Bool](repeating: false, count: p.pixels.count)) == nil)
        var nan = p.pixels
        nan[100] = .nan
        #expect(FriedelCentre.find(pattern: DiffractionPattern(width: p.width, height: p.height, pixels: nan)) == nil)
        // a single bright pixel at (0, 0): the peak sits at index 0, whose neighbour -1 wraps to the last row and column as in NumPy
        var corner = [Float](repeating: 0, count: 16 * 16)
        corner[0] = 100
        let c0 = try #require(FriedelCentre.find(pattern: DiffractionPattern(width: 16, height: 16, pixels: corner)))
        #expect(abs(c0.row) < 1e-9 && abs(c0.column) < 1e-9, Comment(rawValue: "\(c0)"))
        // an all-true mask is the unmasked correlation up to rounding, on an open pattern
        let open = try #require(Self.pins.cases.first { $0.name == "open-beam" })
        let plain = try #require(FriedelCentre.find(pattern: open.pattern))
        let masked = try #require(FriedelCentre.find(pattern: open.pattern, mask: [Bool](repeating: true, count: open.pixels.count)))
        #expect(abs(plain.row - masked.row) < 0.05 && abs(plain.column - masked.column) < 0.05, Comment(rawValue: "\(plain) \(masked)"))
    }
}
