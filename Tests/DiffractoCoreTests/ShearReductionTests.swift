//
//  ShearReductionTests.swift
//  Slice 1b, Gate B finding F1 (registration docs/registration-slice-1b.md section 6, "Gate B on 4876b18"): a cell sheared
//  by hundreds or thousands of lattice units read aP / -1 with no warning, because the reduction shaved it back by unit
//  shears under a 2000-step cap. The reduction shears by integer multiples now, so any shear reduces in a few steps; when
//  the reduction still does not settle (rounding in the six parameters at a shear of 1e6) the result says so
//  (`reductionConverged`) and the conventional route falls back to the file setting's detection.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct ShearReductionTests {

    /// The cell whose basis vectors are the integer combinations `columns` of the basis `base`.
    private static func sheared(_ base: [SIMD3<Double>], _ columns: [SIMD3<Int>]) -> UnitCell {
        var v: [SIMD3<Double>] = []
        for c in columns {
            let x: Double = Double(c.x), y: Double = Double(c.y), z: Double = Double(c.z)
            v.append(x * base[0] + y * base[1] + z * base[2])
        }
        func angle(_ i: Int, _ j: Int) -> Double {
            // atan2 of |cross| and dot: accurate also for the nearly parallel vectors of a strongly sheared cell
            let d: Double = simd_dot(v[i], v[j])
            let x: Double = simd_length(simd_cross(v[i], v[j]))
            return atan2(x, d) * 180 / .pi
        }
        return UnitCell(a: simd_length(v[0]), b: simd_length(v[1]), c: simd_length(v[2]),
                        alphaDeg: angle(1, 2), betaDeg: angle(0, 2), gammaDeg: angle(0, 1))
    }

    private static let cubic4: [SIMD3<Double>] = [SIMD3(4, 0, 0), SIMD3(0, 4, 0), SIMD3(0, 0, 4)]
    /// The primitive cell of NaCl's face-centred lattice, a = 5.64 Å.
    private static let fccPrimitive: [SIMD3<Double>] = [SIMD3(0, 2.82, 2.82), SIMD3(2.82, 0, 2.82), SIMD3(2.82, 2.82, 0)]

    /// The shears of the finding and its neighbours: every one of the three shear steps (A5, A6, A7) is needed.
    private static func shears(_ k: Int) -> [(String, [SIMD3<Int>])] {
        [("c + k(a + b)", [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(k, k, 1)]),
         ("b + k a", [SIMD3(1, 0, 0), SIMD3(k, 1, 0), SIMD3(0, 0, 1)]),
         ("c + k b", [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(0, k, 1)]),
         ("c + k a", [SIMD3(1, 0, 0), SIMD3(0, 1, 0), SIMD3(k, 0, 1)])]
    }

    /// The six parameters of a cell sheared by k carry a relative rounding of about 1e-16 k²; the symmetry test of the
    /// conventional cell takes a relative tolerance of 1e-9 by default, which the rounding passes up to k = 1200 and
    /// not beyond (the lattice is then exactly cubic only to the rounding): a larger shear is read at 1e-5.
    private static func tolerance(_ k: Int) -> Double { k <= 1200 ? UnitCell.conventionalCellTolerance : 1e-5 }

    /// Cubic a = 4 under the finding's shears, k from 1 to 10^5: the reduction converges in a few dozen steps (the unit
    /// shear took k of them: 1200 steps and beyond its cap of 2000), the reduced cell is the cube, and the Bravais
    /// lattice is cP.
    @Test func testShearedCubicCellReducesToTheCube() {
        for k in [1, 300, 1000, 1200, 1500, 3000, 10_000, 100_000] {
            for (name, columns) in Self.shears(k) {
                let cell = Self.sheared(Self.cubic4, columns)
                let r = cell.niggliReduced()
                #expect(r.converged, "k = \(k) \(name): not converged after \(r.steps) steps")
                #expect(r.steps <= 100, "k = \(k) \(name): \(r.steps) steps")
                #expect(parameterDeviation(parameters(r.cell), [4, 4, 4, 90, 90, 90]) <= 1e-5,
                        "k = \(k) \(name): \(parameters(r.cell))")
                let conventional = cell.conventionalCell(tolerance: Self.tolerance(k))
                #expect(conventional.reductionConverged, "k = \(k) \(name)")
                #expect(conventional.bravais == .cP, "k = \(k) \(name): \(conventional.bravais) \(conventional.setting)")
            }
        }
    }

    /// NaCl's primitive cell with b + 1200 a (the finding), and larger shears: cF, a = 5.64 Å.
    @Test func testShearedPrimitiveFaceCentredCellIsFaceCentredCubic() {
        for k in [1200, 5000, 10_000, 100_000] {
            let cell = Self.sheared(Self.fccPrimitive, [SIMD3(1, 0, 0), SIMD3(k, 1, 0), SIMD3(0, 0, 1)])
            let conventional = cell.conventionalCell(tolerance: Self.tolerance(k))
            #expect(conventional.reductionConverged, "k = \(k)")
            #expect(conventional.bravais == .cF, "k = \(k): \(conventional.bravais) \(conventional.setting)")
            #expect(abs(conventional.cell.a - 5.64) <= 1e-3, "k = \(k): \(conventional.cell.a)")
        }
    }

    /// At a shear of 10^6 lattice units the cell is flatter than `UnitCell.isValid` allows (V < 1e-6 a b c), and a
    /// zigzag shear of 10^3 or more makes the six parameters too coarse to settle: the result is then an `aP` that SAYS
    /// the reduction did not converge and that it is the cell as given, never an unflagged aP. A shear that still
    /// reduces gives the right Bravais lattice.
    @Test func testAnExtremeShearIsNeverASilentAP() {
        var flagged = 0
        for k in [100_000, 1_000_000] {
            var all = Self.shears(k)
            all.append(("zigzag", [SIMD3(1, k, k * 7 % 97), SIMD3(0, 1, k / 2), SIMD3(0, 0, 1)]))
            for (name, columns) in all {
                let cell = Self.sheared(Self.cubic4, columns)
                let conventional = cell.conventionalCell(tolerance: 1e-5)
                if conventional.reductionConverged {
                    #expect(conventional.bravais != .aP, "k = \(k) \(name): a silent aP (\(conventional.setting))")
                } else {
                    flagged += 1
                    #expect(conventional.bravais == .aP && conventional.transformation == .identity, "k = \(k) \(name)")
                    #expect(conventional.setting.contains("not reducible"), "k = \(k) \(name): \(conventional.setting)")
                }
            }
        }
        #expect(flagged >= 2, "the cells beyond the validity rule are flagged: \(flagged)")
    }

    /// An invalid cell is flagged, the good cell is not; the structure route falls back to the file setting's detection
    /// and carries the flag to `ConventionalSetting`.
    @Test func testANonConvergedReductionFallsBackToTheFileSettingAndSaysSo() {
        let flat = UnitCell(a: 4, b: 4, c: 4, alphaDeg: 120, betaDeg: 120, gammaDeg: 120)   // zero volume
        let conventional = flat.conventionalCell()
        #expect(!conventional.reductionConverged)
        #expect(conventional.transformation == .identity)
        #expect(conventional.setting.contains("not reducible"))

        let cube = UnitCell(a: 4, b: 4, c: 4, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        #expect(cube.conventionalCell().reductionConverged)

        let site = AtomSite(z: 11, fractional: SIMD3(0, 0, 0), occupancy: 1)
        let structure = CrystalStructure.detectingSymmetry(cell: flat, sites: [site])
        let fileSetting = CrystalStructure.detectingSymmetry(cell: flat, sites: [site], conventional: false)
        let setting = structure.conventionalSetting
        #expect(setting != nil && setting?.reductionConverged == false)
        #expect(structure.operations == fileSetting.operations)
    }
}
