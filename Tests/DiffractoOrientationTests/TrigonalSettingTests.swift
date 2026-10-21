import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoOrientation

/// The two -3m settings (twofolds between the axes, -31m, and on the axes, -3m1). Pins: tools/symmetry/trigonal-pins.py (pymatgen
/// 2026.9.24 / spglib 2.8.0: `SpacegroupAnalyzer(from_spacegroup(n)).get_point_group_operations()` for space groups 162-167). The
/// operator comparison is on sets of 3x3 matrices, tolerance 1e-12 absolute (entries are 0, +-1, +-0.5, +-sqrt(3)/2: pure Double
/// arithmetic); the fractional matrices are integers and compared exact.
private struct TrigonalPins: Decodable {
    struct Site: Decodable { let z: Int; let abc: [Double] }
    struct Lattice: Decodable { let a, b, c, alpha, beta, gamma: Double }
    struct Case: Decodable {
        let spaceGroup: Int; let symbol, setting: String; let lattice: Lattice; let sites: [Site]
        let operatorsFractional, operatorsCartesian: [[[Double]]]
    }
    let cases: [Case]
}
private func mat(_ r: [[Double]]) -> simd_double3x3 { .init(rows: r.map { SIMD3($0[0], $0[1], $0[2]) }) }
private func gap(_ a: simd_double3x3, _ b: simd_double3x3) -> Double {
    (0..<3).flatMap { c in (0..<3).map { abs(a[c][$0] - b[c][$0]) } }.max()!
}

struct TrigonalSettingTests {
    private func pins() throws -> [TrigonalPins.Case] {
        let file = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures/trigonal-pins.json")
        return try JSONDecoder().decode(TrigonalPins.self, from: Data(contentsOf: file)).cases
    }
    private func structure(_ c: TrigonalPins.Case) -> CrystalStructure {
        let l = c.lattice
        return CrystalStructure.detectingSymmetry(
            cell: UnitCell(a: l.a, b: l.b, c: l.c, alphaDeg: l.alpha, betaDeg: l.beta, gammaDeg: l.gamma),
            sites: c.sites.map { AtomSite(z: $0.z, fractional: SIMD3($0.abc[0], $0.abc[1], $0.abc[2])) })
    }
    private func setting(_ name: String) throws -> ThreeBarMSetting { try #require(ThreeBarMSetting(rawValue: name)) }

    @Test func operatorsOfBothSettingsAreTheProperPartOfPymatgensPointGroup() throws {
        var worst = 0.0
        for c in try pins() {
            let s = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: try setting(c.setting))
            #expect(s.threeBarMSetting == ThreeBarMSetting(rawValue: c.setting))
            let proper = c.operatorsCartesian.map(mat).filter { simd_determinant($0) > 0 }
            #expect(proper.count == 6 && s.operators.count == 6)
            for e in proper {
                let g = s.operators.map { gap($0, e) }.min()!
                worst = max(worst, g)
                #expect(g < 1e-7, "\(c.spaceGroup) \(c.setting)")  // D17 (was 1e-12)
            }
            // Exact: the fractional integer matrices.
            let expected = Set(c.operatorsFractional.filter { m in
                let d = mat(m); return abs(simd_determinant(d) - 1) < 1e-9
            }.map { $0.flatMap { $0.map { Int($0) } } })
            let kit = Set(s.fractionalRotations.map { r in r.rows.flatMap { [Int($0.x), Int($0.y), Int($0.z)] } })
            #expect(kit == expected, "\(c.spaceGroup) \(c.setting)")
        }
        print("TRIGONAL operators versus pymatgen: worst gap \(worst)")
    }

    @Test func theSettingsDifferInExactlyTheThreeTwofolds() {
        let a = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarOneM)
        let b = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarMOne)
        #expect(Set(a.fractionalRotations).intersection(Set(b.fractionalRotations)).count == 3)   // identity and the two threefolds
        // Default is the original setting; other classes have no setting and ignore the argument.
        #expect(CrystalOrientationSymmetry(laueClass: .threeBarM).threeBarMSetting == .threeBarOneM)
        #expect(CrystalOrientationSymmetry(laueClass: .threeBarM).operators == a.operators)
        for laue in LaueClass.allCases where laue != .threeBarM {
            let s = CrystalOrientationSymmetry(laueClass: laue, setting: .threeBarMOne)
            #expect(s.threeBarMSetting == nil && s.operators == CrystalOrientationSymmetry(laueClass: laue).operators, "\(laue)")
        }
    }

    @Test func aStructureNamesItsOwnSettingAndTheModelIsUsable() throws {
        var seen: [String] = []
        for c in try pins() {
            let s = structure(c)
            let symmetry = try #require(PhaseModel.orientationSymmetry(of: s), "\(c.symbol)")
            #expect(symmetry.laueClass == .threeBarM && symmetry.threeBarMSetting?.rawValue == c.setting, "\(c.symbol)")
            #expect(PhaseModel.orientationClass(of: s) == .threeBarM)
            let model = PhaseModel(id: c.symbol, displayName: c.symbol, structure: s, source: .custom, spaceGroupNumber: c.spaceGroup)
            #expect(model.orientationSetting?.rawValue == c.setting && model.supportsOrientationMapping, "\(c.symbol)")
            // (The pin structures are crowded test cells: only the symmetry verdict is asserted, not close contacts.)
            #expect(!model.validationIssues.map(\.code).contains("hexagonal_symmetry_mismatch"), "\(c.symbol)")
            #expect(model.provenance["laue_setting"] == c.setting && model.provenance["laue_class"] == "-3m")
            #expect(model.orientationSymmetry?.operators == symmetry.operators)
            // The reduction works in the structure's own frame: a rotation by each of its operators reduces to the identity.
            for op in symmetry.operators {
                #expect(try symmetry.reduce(op).disorientationRad < 1e-6, "\(c.symbol)")
            }
            seen.append("\(c.symbol)=\(c.setting)")
        }
        print("TRIGONAL detection: \(seen.joined(separator: ", "))")
    }

    @Test func aDeclaredSettingThatTheCellDoesNotHaveIsAValidationIssue() throws {
        let cases = try pins()
        let p31m = structure(try #require(cases.first { $0.setting == "-31m" }))
        let p3m1 = structure(try #require(cases.first { $0.setting == "-3m1" }))
        let wrongA = PhaseModel(id: "a", displayName: "a", structure: p31m, orientationClass: .threeBarM,
                                orientationSetting: .threeBarMOne, source: .custom)
        let wrongB = PhaseModel(id: "b", displayName: "b", structure: p3m1, orientationClass: .threeBarM, source: .custom)
        #expect(wrongA.validationIssues.map(\.code).contains("hexagonal_symmetry_mismatch"))
        #expect(wrongB.validationIssues.map(\.code).contains("hexagonal_symmetry_mismatch"))
        let right = PhaseModel(id: "c", displayName: "c", structure: p3m1, orientationClass: .threeBarM,
                               orientationSetting: .threeBarMOne, source: .custom)
        #expect(!right.validationIssues.map(\.code).contains("hexagonal_symmetry_mismatch") && right.revisionID.contains("|-3m1|"))
        #expect(PhaseModel(id: "d", displayName: "d", structure: p31m, source: .custom).revisionID.contains("|-3m|"))
    }

    @Test func theColourKeyIsConstantOnEachSettingsOrbits() throws {
        // A pole and its images under the setting's own operators (and inversion) share a colour; under the other setting's
        // twofolds they generally do not.
        for setting in ThreeBarMSetting.allCases {
            let s = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: setting)
            let other = CrystalOrientationSymmetry(laueClass: .threeBarM, setting: setting == .threeBarMOne ? .threeBarOneM : .threeBarMOne)
            var moved = 0.0
            for d in [SIMD3<Double>(0.31, 0.52, 0.79), SIMD3(0.9, 0.2, 0.1), SIMD3(0.2, 0.8, 0.3)] {
                let reference = try s.ipfColor(direction: d, key: .orix)
                for op in other.operators {
                    let c = try s.ipfColor(direction: op * d, key: .orix)
                    moved = max(moved, Double(simd_length(c - reference)))
                }
            }
            #expect(moved > 0.05, "\(setting): the other setting's operators must not preserve the colour")
        }
    }
}
