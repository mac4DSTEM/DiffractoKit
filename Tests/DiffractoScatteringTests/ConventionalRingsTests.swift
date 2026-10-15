//
//  ConventionalRingsTests.swift
//  Slice 1b truth that needs slice 2's cut (registration docs/registration-slice-1b.md section 3): the rings after the
//  |F| cut of a PRIMITIVE cell, taken in its conventional setting, are the conventional cell's rings exactly (names,
//  g to 1e-12, multiplicities): Au 8, Si 10, and β″'s 417 rings from 1442 kept reflections.
//

import Foundation
import Testing
import simd
import DiffractoCore
import DiffractoScattering

struct ConventionalRingsTests {

    private static func primitiveFCC(a: Double) -> UnitCell { .rhombohedral(a: a / 2.0.squareRoot(), alphaDeg: 60) }

    private static func sameRings(_ x: KinematicReflectionList, _ y: KinematicReflectionList, _ label: String) {
        #expect(x.kept.count == y.kept.count, "\(label) kept")
        #expect(x.families.count == y.families.count, "\(label) rings")
        for (p, q) in zip(x.families, y.families) {
            #expect(p.family.name == q.family.name, "\(label) name \(p.family.name) vs \(q.family.name)")
            #expect(p.family.multiplicity == q.family.multiplicity, "\(label) multiplicity of \(p.family.name)")
            #expect(abs(p.family.g - q.family.g) <= 1e-12 * q.family.g, "\(label) g of \(p.family.name)")
            #expect(abs(p.intensity - q.intensity) <= 1e-9 * q.intensity, "\(label) intensity of \(p.family.name)")
        }
    }

    @Test func testGoldAndSiliconPrimitiveCellsGiveTheConventionalCellsRings() throws {
        let au = CrystalStructure.detectingSymmetry(cell: Self.primitiveFCC(a: 4.0782), sites: [AtomSite(z: 79, fractional: [0, 0, 0])])
        let auConv = try #require(au.conventionalSetting).structure
        let reference = try Presets.gold.kinematicReflections(kMax: 1.2)
        let got = try auConv.kinematicReflections(kMax: 1.2)
        #expect(got.families.count == 8)
        Self.sameRings(got, reference, "Au")

        let si = CrystalStructure.detectingSymmetry(cell: Self.primitiveFCC(a: 5.4309), sites: [
            AtomSite(z: 14, fractional: [0, 0, 0]), AtomSite(z: 14, fractional: [0.25, 0.25, 0.25])])
        let siConv = try #require(si.conventionalSetting).structure
        let siGot = try siConv.kinematicReflections(kMax: 1.2)
        #expect(siGot.families.count == 10)
        Self.sameRings(siGot, try Presets.silicon.kinematicReflections(kMax: 1.2), "Si")
    }

    @Test func testBetaDoublePrimePrimitiveCellGivesTheConventionalRings() throws {
        let conv = Presets.betaDoublePrime
        let a = conv.cell.realLattice[0], b = conv.cell.realLattice[1], c = conv.cell.realLattice[2]
        let v1 = 0.5 * (a + b), v2 = 0.5 * (b - a)
        func angle(_ x: SIMD3<Double>, _ y: SIMD3<Double>) -> Double { acos(dot(x, y) / (length(x) * length(y))) * 180 / .pi }
        let cell = UnitCell(a: length(v1), b: length(v2), c: length(c), alphaDeg: angle(v2, c), betaDeg: angle(v1, c), gammaDeg: angle(v1, v2))
        var primitive: [AtomSite] = []
        for s in conv.sites {
            let p = wrapFractional(SIMD3(s.fractional.x + s.fractional.y, s.fractional.y - s.fractional.x, s.fractional.z))
            if !primitive.contains(where: { $0.z == s.z && abs($0.occupancy - s.occupancy) < 1e-9 && fractionalAxisDelta($0.fractional.x, p.x) < 1e-9
                    && fractionalAxisDelta($0.fractional.y, p.y) < 1e-9 && fractionalAxisDelta($0.fractional.z, p.z) < 1e-9 }) {
                primitive.append(AtomSite(z: s.z, fractional: p, occupancy: s.occupancy))
            }
        }
        #expect(primitive.count == 11)
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: primitive)
        let setting = try #require(s.conventionalSetting)
        #expect(setting.bravais == .mC && setting.operations.count == 8 && setting.laueClass == .twoOverM)
        let got = try setting.structure.kinematicReflections(kMax: 1.2)
        #expect(got.kept.count == 1442)
        #expect(got.families.count == 417)
        Self.sameRings(got, try conv.kinematicReflections(kMax: 1.2), "β″")
    }
}
