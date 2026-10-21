import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// Pins: tools/orix/direction-pins.py (orix 0.15.0 `Vector3d.in_fundamental_sector`, orix's Laue groups in the kit frame).
/// The kit's zone is not orix's (largest-(z,x,y) orbit member, or the apps' native sector for m-3m, 6/mmm, 2/m), so the
/// comparison is of orbits under orix's own group. Written tolerances: 1e-5 (orix's allclose; was 1e-9) on unit-vector components for orbit membership against orix, 1e-9 for the kit's own invariance check; formerly one 1e-9 for orbit
/// membership and for the reduced pole being the same for every orbit member (the kit's own ties are compared at 1e-13;
/// the random poles are generic, the designed poles lie on mirror planes and axes, where the orbit members coincide).
private struct DirectionPinsFlat: Decodable {
    struct Class: Decodable { let symbol: String; let group: [[[Double]]]; let poles, orixFolded, probes: [[Double]] }
    let classes: [Class]
}
private func v(_ a: [Double]) -> SIMD3<Double> { SIMD3(a[0], a[1], a[2]) }

struct OrixDirectionTests {
    private func flat() throws -> [DirectionPinsFlat.Class] {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/orix-direction-pins.json")
        return try JSONDecoder().decode(DirectionPinsFlat.self, from: Data(contentsOf: file)).classes
    }
    private func group(_ c: DirectionPinsFlat.Class) -> [simd_double3x3] {
        c.group.map { simd_double3x3(rows: $0.map { SIMD3($0[0], $0[1], $0[2]) }) }
    }
    /// Smallest distance from `a` to the orbit of `b` under `g`.
    private func orbitGap(_ a: SIMD3<Double>, _ b: SIMD3<Double>, _ g: [simd_double3x3]) -> Double {
        g.map { simd_length($0 * b - a) }.min()!
    }

    @Test func reduceDirectionIsAnOrbitRepresentativeAgreeingWithOrix() throws {
        var worstMember = 0.0, worstOrix = 0.0, worstInvariance = 0.0, count = 0
        for c in try flat() {
            let s = try #require(pinnedSymmetry(c.symbol))
            let g = group(c)
            #expect(g.count == s.laueClass.order, "\(c.symbol)")
            for (p, folded) in zip(c.poles, c.orixFolded) {
                let d = simd_normalize(v(p))
                let r = try s.reduceDirection(d)
                // 1. the reduced pole is a member of the pole's orbit and of orix's folded pole's orbit
                let m = orbitGap(r, d, g), o = orbitGap(r, v(folded), g)
                worstMember = max(worstMember, m); worstOrix = max(worstOrix, o)
                #expect(m < 1e-5 && o < 1e-5, "\(c.symbol) \(p): member \(m), orix \(o)")
                // 2. every orbit member reduces to the same pole (the zone holds one pole per orbit)
                for gi in g {
                    let r2 = try s.reduceDirection(gi * d)
                    let gap = simd_length(r2 - r)
                    worstInvariance = max(worstInvariance, gap)
                    #expect(gap < 1e-9, "\(c.symbol) \(p): \(r) versus \(r2)")
                }
                count += 1
            }
        }
        print("ORIX reduceDirection: \(count) poles; worst orbit-membership gap \(worstMember), versus orix's fold \(worstOrix), invariance \(worstInvariance)")
    }

    @Test func sampledPolesLieInTheZoneAreInequivalentAndCoverIt() throws {
        var table: [String] = []
        for c in try flat() {
            let s = try #require(pinnedSymmetry(c.symbol))
            let g = group(c)
            // Coverage as the covering radius over the probe poles (angle to the nearest orbit member of the nearest sample).
            var radii: [Int: Double] = [:]
            for n in [1, 3, 8, 16, 64] {
                let samples = try s.sampleFundamentalZone(count: n)
                #expect(samples.count == n, "\(c.symbol) \(n)")
                for x in samples {
                    #expect(simd_length(try s.reduceDirection(x) - x) < 1e-9, "\(c.symbol): sample not in the zone")
                    #expect(abs(simd_length(x) - 1) < 1e-12)
                }
                // No two samples are equivalent under orix's group.
                for i in 0..<samples.count { for j in (i + 1)..<samples.count {
                    #expect(orbitGap(samples[i], samples[j], g) > 1e-6, "\(c.symbol) \(n): samples \(i), \(j) are equivalent")
                } }
                let worst = c.probes.map { p -> Double in
                    let probe = v(p)
                    let best = samples.map { sample in g.map { simd_dot(probe, $0 * sample) }.max()! }.max()!
                    return acos(min(1, best))
                }.max()!
                radii[n] = worst
            }
            // Regression bound, not an optimality claim: more samples never cover worse than 1 sample, and 64 samples cover
            // the zone to within the coarse bound written below (measured radii are printed and recorded in the docs).
            #expect(radii[64]! < radii[8]! && radii[8]! <= radii[1]!, "\(c.symbol): \(radii)")
            if c.symbol != "-1" { #expect(radii[64]! < 0.35, "\(c.symbol) 64 samples: \(radii[64]!) rad") }
            table.append("\(c.symbol): " + [1, 3, 8, 16, 64].map { String(format: "%.3f", radii[$0]!) }.joined(separator: "/"))
        }
        print("ORIX sampling covering radius (rad, probe poles, counts 1/3/8/16/64): " + table.joined(separator: "; "))
    }
}
