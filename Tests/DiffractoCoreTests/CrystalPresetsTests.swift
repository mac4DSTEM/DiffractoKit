import Foundation
import Testing
import simd
import DiffractoCore

/// PresetPins: tools/presets/pins.py (pymatgen 2026.9.24, spglib 2.8.0). The pins build each structure independently, from the
/// published cell and Wyckoff positions with `Structure.from_spacegroup`. Tolerances, written per quantity: lattice
/// parameters, volume and d(hkl) 1e-12 relative (closed forms in float64 on both sides); fractional site positions 1e-9
/// (the kit's hand lists carry 2/3 as a Double, pymatgen's orbit expansion computes it); counts, Z, space-group number
/// and operation count exact.
private struct PresetPins: Decodable {
    struct Site: Decodable { let z: Int; let abc: [Double] }
    struct Lattice: Decodable { let a, b, c, alpha, beta, gamma: Double }
    struct Preset: Decodable {
        let lattice: Lattice
        let volume: Double
        let spaceGroupNumber: Int
        let pointGroup: String
        let operationCount: Int
        let siteCount: Int
        let sites: [Site]
        let d: [Double]
    }
    let hkl: [[Int]]
    let presets: [String: Preset]
}

private func kitStructure(_ id: String) throws -> CrystalStructure {
    switch id {
    case "au_fcc": return CrystalPresets.gold
    case "al_fcc": return CrystalPresets.aluminum
    case "ni_fcc": return CrystalPresets.nickel
    case "cu_fcc": return CrystalPresets.copper
    case "fe_bcc": return CrystalPresets.iron
    case "si_diamond": return CrystalPresets.silicon
    case "sc_cubic": return CrystalPresets.simpleCubic(a: 3.0, z: 84)
    case "mg_hcp": return CrystalPresets.magnesium
    case "ws2_2h": return CrystalPresets.tungstenDisulfide
    case "beta_double_prime": return CrystalPresets.betaDoublePrime
    default: throw CocoaError(.fileReadUnknown)
    }
}

struct CrystalPresetsTests {
    private func pins() throws -> PresetPins {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/presets-pins.json")
        return try JSONDecoder().decode(PresetPins.self, from: Data(contentsOf: file))
    }

    private func wrapped(_ x: Double) -> Double { let m = x - x.rounded(.down); return m > 1 - 1e-9 ? 0 : m }

    @Test func everyPresetMatchesPymatgenAndSpglib() throws {
        let fixture = try pins()
        #expect(fixture.presets.count == 10)
        var worstD = 0.0, worstSite = 0.0, worstVolume = 0.0
        for (id, p) in fixture.presets.sorted(by: { $0.key < $1.key }) {
            let s = try kitStructure(id)
            // Cell and volume.
            let l = p.lattice
            for (name, kit, ref) in [("a", s.cell.a, l.a), ("b", s.cell.b, l.b), ("c", s.cell.c, l.c),
                                     ("alpha", s.cell.alphaDeg, l.alpha), ("beta", s.cell.betaDeg, l.beta),
                                     ("gamma", s.cell.gammaDeg, l.gamma)] {
                #expect(abs(kit - ref) <= 1e-12 * ref, "\(id) \(name)")
            }
            worstVolume = max(worstVolume, abs(s.cell.volume - p.volume) / p.volume)
            #expect(abs(s.cell.volume - p.volume) <= 1e-12 * p.volume, "\(id) volume")
            // Sites: the same multiset of (Z, fractional position mod 1).
            #expect(s.sites.count == p.siteCount, "\(id) site count")
            var unmatched = p.sites
            for site in s.sites {
                let f = SIMD3(wrapped(site.fractional.x), wrapped(site.fractional.y), wrapped(site.fractional.z))
                func gap(_ e: PresetPins.Site) -> Double {
                    var g = 0.0
                    for k in 0..<3 { let d = abs(f[k] - e.abc[k]); g = max(g, min(d, 1 - d)) }
                    return e.z == site.z ? g : .infinity
                }
                guard let i = unmatched.indices.min(by: { gap(unmatched[$0]) < gap(unmatched[$1]) }) else { Issue.record("\(id) extra site"); continue }
                worstSite = max(worstSite, gap(unmatched[i]))
                #expect(gap(unmatched[i]) < 1e-9, "\(id) site \(f)")
                unmatched.remove(at: i)
            }
            #expect(unmatched.isEmpty, "\(id) missing sites")
            // d spacings.
            for (h, e) in zip(fixture.hkl, p.d) {
                let d = s.cell.d(MillerIndex(h[0], h[1], h[2]))
                worstD = max(worstD, abs(d - e) / e)
                #expect(abs(d - e) <= 1e-12 * e, "\(id) d\(h)")
            }
            // Symmetry: the Laue class from the point group, and the operation count of the cell.
            let expected: LaueClass = ["m-3m": .mThreeBarM, "6/mmm": .sixOverMmm, "2/m": .twoOverM][p.pointGroup]!
            #expect(s.laueClass == expected, "\(id) Laue class")
            #expect(s.operations.count == p.operationCount, "\(id) operations \(s.operations.count) vs \(p.operationCount)")
        }
        print("PRESETS worst: volume \(worstVolume), d \(worstD), site \(worstSite)")
    }

    @Test func cubicFamiliesAndTheLiteratureCellContents() {
        #expect(CrystalPresets.cubic(.fcc, a: 4, z: 29).sites.count == 4)
        #expect(CrystalPresets.cubic(.bcc, a: 3, z: 26).sites.count == 2)
        #expect(CrystalPresets.cubic(.simpleCubic, a: 3, z: 84).sites.count == 1)
        #expect(CrystalPresets.cubic(.diamond, a: 5, z: 14).sites.count == 8)
        #expect(CrystalPresets.CubicFamily.allCases.map(\.provenanceID) == ["fcc", "bcc", "simple_cubic", "diamond"])
        // Mg10 Si12: the cell content Andersen et al. state.
        let beta = CrystalPresets.betaDoublePrime
        #expect(beta.sites.count == 22)
        #expect(beta.sites.filter { $0.z == 12 }.count == 10 && beta.sites.filter { $0.z == 14 }.count == 12)
        #expect(beta.isUsable && CrystalPresets.tungstenDisulfide.isUsable)
    }
}
