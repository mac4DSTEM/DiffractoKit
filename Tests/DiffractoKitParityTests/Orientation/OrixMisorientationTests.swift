import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/pins/orix/misorientation-pins.py (orix 0.15.0). Tolerances: angles atol 1e-4 (orix's test_misorientation tolerance; measured 1e-15),
/// equivalent sets and reduced matrices atol 1e-9 where the smallest angle is unique, the Euler angles atol 1e-9. Every draw is compared with
/// orix's `reduce` loop run with its zone built as `from_symmetry(end, start)` (`orderedZone*`, the cause in Misorientation.swift); where orix's
/// shipped `reduce` is not the smallest angle, the test asserts that it is larger than the minimum.
private struct Pins: Decodable {
    struct Item: Decodable {
        let matrix, reduced: [[Double]]
        let equivalentCount: Int
        let angle, minAngle, orderedZoneAngle: Double
        let orderedZoneReduced: [[Double]]
        let orixIsMinimal, unique: Bool
        let equivalents: [[[Double]]]?
        let grainExchangeCount: Int?
    }
    struct Pair: Decodable { let source, target: String; let sourceOrder, targetOrder: Int; let items: [Item] }
    struct EulerItem: Decodable { let matrix: [[Double]]; let euler: [Double] }
    struct Euler: Decodable { let symbol, orixProperName: String; let region: [Int]; let items: [EulerItem] }
    struct Beta: Decodable { let R0: [[Double]]; let equivalentCount: Int; let reducedAngle, reducedAngleClosedForm: Double }
    let pairs: [Pair]
    let euler: [Euler]
    let beta: Beta
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    var g = 0.0
    for c in 0 ..< 3 { for r in 0 ..< 3 { g = max(g, abs(a[c][r] - b[c][r])) } }
    return g
}
private func symmetry(_ symbol: String) -> CrystalOrientationSymmetry {
    if symbol == "-31m" { return CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarOneM) }
    return pinnedSymmetry(symbol)!
}

struct OrixMisorientationTests {
    private func pins() throws -> Pins {
        let file = FixtureDirs.target("DiffractoOrientationTests").appendingPathComponent("Fixtures/orix-misorientation-pins.json")
        return try JSONDecoder().decode(Pins.self, from: Data(contentsOf: file))
    }

    @Test func equivalentAndReduceMatchOrixForTwelveGroupPairs() throws {
        var worstAngle = 0.0, worstSet = 0.0, worstMatrix = 0.0, orixNotMinimal = 0
        for p in try pins().pairs {
            let source = symmetry(p.source), target = symmetry(p.target)
            #expect(source.operators.count == p.sourceOrder && target.operators.count == p.targetOrder, Comment(rawValue: "\(p.source) \(p.target)"))
            for item in p.items {
                let m = try Misorientation(mat(item.matrix), source: source, target: target)
                let eq = m.equivalent()
                #expect(eq.count == source.operators.count * target.operators.count)
                var distinct: [simd_double3x3] = []
                for e in eq where !distinct.contains(where: { gap($0, e) < 1e-9 }) { distinct.append(e) }
                #expect(distinct.count == item.equivalentCount, Comment(rawValue: "\(p.source)>\(p.target) distinct equivalents"))
                if let pinned = item.equivalents {
                    for e in pinned {
                        let g = distinct.map { gap($0, mat(e)) }.min() ?? 1
                        worstSet = max(worstSet, g)
                        #expect(g < 1e-9, Comment(rawValue: "\(p.source)>\(p.target) equivalent missing"))
                    }
                }
                let r = m.reduce()
                worstAngle = max(worstAngle, abs(r.angleRad - item.minAngle))
                #expect(abs(r.angleRad - item.minAngle) < 1e-9, Comment(rawValue: "\(p.source)>\(p.target) minimum angle"))
                #expect(abs(r.angleRad - item.orderedZoneAngle) < 1e-4, Comment(rawValue: "\(p.source)>\(p.target) orix ordered-zone angle"))
                if item.unique {
                    let g = gap(r.matrix, mat(item.orderedZoneReduced))
                    worstMatrix = max(worstMatrix, g)
                    #expect(g < 1e-9, Comment(rawValue: "\(p.source)>\(p.target) orix ordered-zone reduced matrix"))
                }
                if item.orixIsMinimal {
                    #expect(abs(r.angleRad - item.angle) < 1e-4, Comment(rawValue: "\(p.source)>\(p.target) orix angle"))
                    if item.unique {
                        let g = gap(r.matrix, mat(item.reduced))
                        worstMatrix = max(worstMatrix, g)
                        #expect(g < 1e-9, Comment(rawValue: "\(p.source)>\(p.target) reduced matrix"))
                    }
                } else {
                    orixNotMinimal += 1
                    #expect(item.angle > item.minAngle)        // orix's own reduce is larger than the minimum: the recorded DEVIATION
                }
                if let ge = item.grainExchangeCount {
                    var d: [simd_double3x3] = []
                    for e in m.equivalent(grainExchange: true) where !d.contains(where: { gap($0, e) < 1e-9 }) { d.append(e) }
                    #expect(d.count == ge, Comment(rawValue: "\(p.source) grain exchange"))
                }
            }
        }
        print("ORIX misorientation worst: angle \(worstAngle), equivalent set \(worstSet), reduced matrix \(worstMatrix); orix reduce above the minimum on \(orixNotMinimal) draws")
    }

    @Test func eulerFundamentalRegionMatchesOrixForEveryClass() throws {
        var worst = 0.0
        for c in try pins().euler {
            let s = symmetry(c.symbol)
            for item in c.items {
                let o = try Orientation(mat(item.matrix), symmetry: s)
                let e = o.inEulerFundamentalRegion()
                let want = SIMD3(item.euler[0], item.euler[1], item.euler[2])
                // angles near 0 and 2 pi are the same angle: compare on the circle
                var g = 0.0
                for i in 0 ..< 3 {
                    var d = abs(e[i] - want[i])
                    d = min(d, abs(d - 2 * .pi))
                    g = max(g, d)
                }
                worst = max(worst, g)
                #expect(g < 1e-8, Comment(rawValue: "\(c.symbol) \(e) vs \(want)"))
            }
        }
        print("ORIX Euler region worst: \(worst)")
    }

    @Test func betaDoublePrimeInAluminium() throws {
        let b = try pins().beta
        let m = try Misorientation(mat(b.R0), source: CrystalOrientationSymmetry(laueClass: .twoOverM),
                                   target: CrystalOrientationSymmetry(laueClass: .mThreeBarM))
        var distinct: [simd_double3x3] = []
        for e in m.equivalent() where !distinct.contains(where: { gap($0, e) < 1e-9 }) { distinct.append(e) }
        #expect(distinct.count == b.equivalentCount && b.equivalentCount == 24)
        #expect(abs(m.reduce().angleRad - b.reducedAngle) < 1e-4)
        #expect(abs(b.reducedAngle - b.reducedAngleClosedForm) < 1e-9)
    }
}
