//
//  ToleranceFixesTests.swift
//  DiffractoKit 0.2.0, the Gate B findings on tolerances (tier 1 for F5, F6 and F8, whose outcomes a scientist reads;
//  tier 2 for F3): the position tolerance is refused where it cannot work and no structure is ever returned without the
//  identity operation; the metric tolerances read NaN, negative and infinite as 0; an absurd cell and an invalid Niggli
//  epsilon never trap and are never reported as reduced. Every test names under "RED UNDER" the production change that
//  turns it red.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore

struct ToleranceFixesTests {
    typealias Tolerance = SymmetryDetection.Tolerance

    static func silicon() throws -> String { try Fixtures.text("silicon_Fd-3m.cif") }

    static func read(_ text: String, _ configure: (inout CIFImport.Policy) -> Void) -> Result<CIFDocument, CIFImportError> {
        var policy = CIFImport.Policy.default
        configure(&policy)
        do { return .success(try CIFImport.read(text, policy: policy)) } catch { return .failure(error) }
    }

    /// The tolerances that cannot work: not a number, zero or negative, infinite, below the 1e-9 floor (down to the
    /// denormals), and a window of half a cell or more.
    static let refused: [Tolerance] = [
        .fractional(.nan), .fractional(0), .fractional(-0.0), .fractional(-1), .fractional(-1e308), .fractional(-.infinity),
        .fractional(.infinity), .fractional(1e-300), .fractional(5e-324), .fractional(Double.leastNormalMagnitude),
        .fractional(1e-10), .fractional(0.5), .fractional(0.7), .fractional(1e300),
        .cartesian(.nan), .cartesian(0), .cartesian(-0.0), .cartesian(-1), .cartesian(-.infinity), .cartesian(.infinity),
        .cartesian(1e-300), .cartesian(5e-324), .cartesian(1e-10), .cartesian(1e300), .cartesian(5), .cartesian(1000),
    ]

    // MARK: F5 and F6

    /// `CIFImport.read` refuses every such position tolerance with `invalidStructure` and a message naming the tolerance;
    /// none traps and none imports a structure.
    /// RED UNDER: the tolerance check removed from `SymmetryDetection.detect` / `detectingInConventionalCell` (the first
    /// negative or denormal value traps in the cell list, NaN and 0 import a structure with no operation).
    @Test func cifImportRefusesAnUnusablePositionTolerance() throws {
        let text = try Self.silicon()
        for tolerance in Self.refused {
            switch Self.read(text, { $0.positionTolerance = tolerance }) {
            case .success: Issue.record("imported under \(tolerance)")
            case .failure(let error):
                guard case .invalidStructure(let messages) = error else { Issue.record("\(error) under \(tolerance)"); continue }
                #expect(messages.contains { $0.contains("position tolerance") }, "\(messages)")
                #expect("\(error)".contains("position tolerance"))
            }
        }
    }

    /// The same tolerances through `detectingSymmetry` (both routes) and `SymmetryDetection.detect`: no trap, the identity
    /// alone, `symmetryRefusal` / `refusal` naming the tolerance, and the structure's own refusal path
    /// (`validationIssues`, `isUsable`) reports it.
    /// RED UNDER: as above.
    @Test func detectionRefusesAnUnusablePositionTolerance() {
        let cell = UnitCell.cubic(a: 4)
        let sites = [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.5, 0.5, 0]),
                     AtomSite(z: 13, fractional: [0.5, 0, 0.5]), AtomSite(z: 13, fractional: [0, 0.5, 0.5])]
        for tolerance in Self.refused {
            for conventional in [true, false] {
                let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, tolerance: tolerance,
                                                                   conventional: conventional)
                #expect(structure.operations == [.identity], "\(tolerance)")
                #expect(structure.symmetryRefusal?.contains("position tolerance") == true)
                #expect(!structure.isUsable)
                #expect(structure.validationIssues.contains { $0.code == "invalid_symmetry_tolerance" })
                #expect(Set(structure.laueGroup.rotations) == [.identity, .inversion])      // {±R} of the identity
                #expect(structure.conventionalSetting == nil)
            }
            let result = SymmetryDetection.detect(sites: sites, cell: cell, tolerance: tolerance)
            #expect(result.operations == [.identity])
            #expect(result.refusal?.contains("position tolerance") == true)
        }
    }

    /// A tolerance in its domain gives no refusal, and the structure is the ordinary one. The domain is a window of
    /// 1e-9 up to (not including) 0.5 of a cell axis, as a fraction (`.fractional`) or as the length times the reciprocal
    /// axis (`.cartesian`).
    /// RED UNDER: the check refusing a usable tolerance (the floor raised, the cap lowered).
    @Test func aUsableToleranceIsNotRefused() throws {
        let cell = UnitCell.cubic(a: 4)
        let sites = [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.5, 0.5, 0]),
                     AtomSite(z: 13, fractional: [0.5, 0, 0.5]), AtomSite(z: 13, fractional: [0, 0.5, 0.5])]
        for tolerance in [Tolerance.fractional(1e-9), .fractional(1e-6), .fractional(5e-3), .fractional(0.25), .fractional(0.499),
                          .cartesian(1e-8), .cartesian(1e-3), .cartesian(0.5), .cartesian(1.9)] {
            #expect(tolerance.refusal(for: cell) == nil, "\(tolerance)")
            #expect(tolerance.isValid(for: cell))
            let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, tolerance: tolerance)
            #expect(structure.symmetryRefusal == nil)
            #expect(structure.operations.first == .identity)
            #expect(structure.isUsable)
        }
        for tolerance in Self.refused { #expect(tolerance.refusal(for: cell) != nil && !tolerance.isValid(for: cell), "\(tolerance)") }
        // the default is unchanged
        #expect(Tolerance.default.refusal(for: cell) == nil)
        let silicon = try CIFImport.read(try Self.silicon())
        #expect(silicon.structure.operations.count == 192)
        // a Cartesian window is tested against the cell it will be used in: 1 Å is half a 2 Å cell
        #expect(Tolerance.cartesian(1).refusal(for: .cubic(a: 2)) != nil)
        #expect(Tolerance.cartesian(1).refusal(for: .cubic(a: 20)) == nil)
    }

    /// No structure is returned without the identity operation, whatever the usable tolerance and however close to a cell
    /// edge the coordinates lie (the wrap snaps within 1e-9 of 0 or 1, which a tolerance at the floor must survive).
    /// RED UNDER: the refused result built with no operations (the identity dropped), or a detection that loses the
    /// identity of a site at a cell edge (the end-of-`detect` guard is the second line of defence behind the floor, and
    /// cannot be reached from a tolerance in the domain).
    @Test func theIdentityIsAlwaysThere() {
        let cell = UnitCell.cubic(a: 4)
        let edges: [Double] = [0, 1e-12, 5e-10, 1e-9, 1.5e-9, 1 - 1e-9, 1 - 5e-10, 1 - 1.5e-9, 0.999999999, 0.5, 0.25 + 1e-10]
        for x in edges {
            let sites = [AtomSite(z: 13, fractional: [x, 0, 1 - x]), AtomSite(z: 8, fractional: [0.5, x, 0.5])]
            for t in [1e-9, 1.0000001e-9, 3e-9, 1e-8, 1e-6, 5e-3, 0.1, 0.4999] {
                for tolerance in [Tolerance.fractional(t), .cartesian(t * 4)] {
                    for conventional in [true, false] {
                        for search in [SymmetryDetection.Search.cellList, .bruteForce] {
                            let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, tolerance: tolerance,
                                                                               conventional: conventional, search: search)
                            #expect(structure.operations.contains { $0.isIdentity }, "x \(x) t \(t) \(tolerance)")
                            #expect(structure.laueGroup.rotations.contains(.identity))
                        }
                    }
                }
            }
        }
    }

    // MARK: F3

    /// A metric tolerance that is NaN, negative or infinite is read as 0 (every length and angle exactly equal): the
    /// silicon cell is still read cubic at 0 because its lengths are written equal; a cell written with a 1e-12 difference
    /// is not. Through the CIF policy, `detectingSymmetry` and the metric family.
    /// RED UNDER: the sanitising of the two metric tolerances removed (NaN and negative read as "nothing equal", infinite
    /// as "everything equal").
    @Test func metricTolerancesReadNaNNegativeAndInfiniteAsZero() throws {
        let text = try Self.silicon()
        let zeroLength = try CIFImport.read(text, policy: { var p = CIFImport.Policy.default; p.metricLengthTolerance = 0; return p }())
        let zeroAngle = try CIFImport.read(text, policy: { var p = CIFImport.Policy.default; p.metricAngleToleranceDeg = 0; return p }())
        for bad in [Double.nan, -1, -1e-9, -.infinity, .infinity] {
            let length = Self.read(text) { $0.metricLengthTolerance = bad }
            let angle = Self.read(text) { $0.metricAngleToleranceDeg = bad }
            #expect(try length.get() == zeroLength, "\(bad)")
            #expect(try angle.get() == zeroAngle, "\(bad)")
        }
        let skew = UnitCell(a: 5, b: 5 * (1 + 1e-12), c: 5, alphaDeg: 90, betaDeg: 90, gammaDeg: 90)
        for bad in [Double.nan, -1, .infinity] {
            #expect(skew.metricFamily(lengthTolerance: bad) == skew.metricFamily(lengthTolerance: 0))
            #expect(skew.metricFamily(angleToleranceDeg: bad) == skew.metricFamily(angleToleranceDeg: 0))
            #expect(UnitCell.cubic(a: 5).metricFamily(lengthTolerance: bad, angleToleranceDeg: bad) == .cubic)
            let sites = [AtomSite(z: 13, fractional: [0, 0, 0])]
            #expect(CrystalStructure.detectingSymmetry(cell: skew, sites: sites, metricLengthTolerance: bad, conventional: false)
                    == CrystalStructure.detectingSymmetry(cell: skew, sites: sites, metricLengthTolerance: 0, conventional: false))
            #expect(CrystalStructure.detectingSymmetry(cell: skew, sites: sites, metricLengthTolerance: bad)
                    == CrystalStructure.detectingSymmetry(cell: skew, sites: sites, metricLengthTolerance: 0))
        }
        #expect(skew.metricFamily(lengthTolerance: 0) != .cubic)
    }
}

// MARK: - F8: absurd cells and the Niggli epsilon

struct AbsurdCellTests {
    /// The fluorapatite file with its first cell length replaced.
    static func fluorapatite(a: String) throws -> String {
        let text = try Fixtures.text("fluorapatite_P63m.cif")
        return text.replacingOccurrences(of: "_cell_length_a 9.3973", with: "_cell_length_a \(a)")
    }

    /// A cell axis a billion times shorter than the others is refused by the import (two sites a few thousandths of an
    /// Å apart) or flagged as not reduced; it never traps (the integer matrix products of the conventional route
    /// overflowed). One that is a trillion times longer imports and says its reduction did not converge.
    /// RED UNDER: the checked conjugation by the transformation (`IntegerMatrix.conjugating`, `conjugatedBy`) replaced by the
    /// trapping `multiplied(by:)` (the first case traps there: a transformation with entries of 4.7e9 meets the file's
    /// operations). The other checked products of the route (the candidate loop, the plane reduction, the adjugate, the
    /// shears) are defence in depth for cells that no test case reaches.
    @Test func anAbsurdCifNeverTraps() throws {
        for a in ["1e-9", "1e-7", "1e-5", "1e-3", "1e-12", "1e-15", "1e-300"] {
            do {
                let document = try CIFImport.read(try Self.fluorapatite(a: a))
                #expect(document.structure.conventionalSetting?.reductionConverged != true, "\(a)")
            } catch let error as CIFImportError {
                guard case .invalidStructure(let messages) = error else { Issue.record("\(error)"); continue }
                #expect(!messages.isEmpty)
            }
        }
        let long = try CIFImport.read(try Self.fluorapatite(a: "1e12"))
        #expect(long.structure.conventionalSetting?.reductionConverged == false)
        #expect(long.structure.conventionalSetting?.bravais == .aP)
    }

    /// The same cells without a file: the conventional cell is either reduced or flagged (`reductionConverged == false`
    /// with the cell as given, `aP` and the identity), `centringMultiplicity` can be read, and the detection runs.
    /// RED UNDER: as above.
    @Test func anAbsurdCellIsFlaggedNotTrapped() {
        for (a, b, c) in [(1e-9, 9.4, 6.9), (1e-12, 9.4, 9.4), (9.4, 1e-9, 6.9), (6.9, 9.4, 1e-9), (1e9, 1.0, 1.0),
                          (1e-9, 1e-9, 1e9), (1e-6, 5.0, 7.0), (3.0, 3.0, 1e-10)] {
            for gamma in [90.0, 120.0, 97.3] {
                let cell = UnitCell(a: a, b: b, c: c, alphaDeg: 90, betaDeg: 90, gammaDeg: gamma)
                let conventional = cell.conventionalCell()
                if !conventional.reductionConverged {
                    #expect(conventional.bravais == .aP && conventional.transformation == .identity, "\(a) \(b) \(c)")
                    #expect(conventional.setting.contains("not reducible"))
                }
                #expect(conventional.centringMultiplicity >= 1, "\(a) \(b) \(c)")
                let structure = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 13, fractional: [0, 0, 0])])
                #expect(structure.operations.first == .identity)
            }
        }
    }

    /// A Niggli epsilon that is NaN, negative or infinite is refused like an invalid cell: unreduced and not converged,
    /// never "converged" on a cell that needs reducing. A sheared cell reduces at the default and at 0.
    /// RED UNDER: the epsilon guard of `niggliReduced` removed (NaN reports converged with the identity on the sheared
    /// cell); the tolerance sanitising of `conventionalCell(tolerance:)` removed (NaN reads the cubic cell as aP).
    @Test func anInvalidNiggliEpsilonIsNotReportedAsReduced() {
        let hex = UnitCell(a: 4, b: 4, c: 7, alphaDeg: 90, betaDeg: 90, gammaDeg: 120)
        // a cell that needs a reduction: b' = b + 3a of the hexagonal cell, a' = a
        let bx = -2.0 + 3 * 4, by = 4 * sin(2 * Double.pi / 3)
        let needs = UnitCell(a: 4, b: (bx * bx + by * by).squareRoot(), c: 7, alphaDeg: 90, betaDeg: 90,
                             gammaDeg: atan2(by, bx) * 180 / .pi)
        #expect(needs.niggliReduced().converged)
        #expect(needs.niggliReduced().transformation != .identity)
        for bad in [Double.nan, -1, -1e-9, -.infinity, .infinity] {
            let reduced = needs.niggliReduced(relativeEpsilon: bad)
            #expect(!reduced.converged, "\(bad)")
            #expect(reduced.transformation == .identity && reduced.steps == 0 && reduced.epsilon == 0)
            #expect(reduced.cell == needs)
        }
        #expect(needs.niggliReduced(relativeEpsilon: 1e-5).converged)
        // the tolerance of the conventional route: NaN, negative and infinite read as 0
        // (at 0 the cubic cell reads hR because cos 90° is not exactly 0; nonsense is "every matrix is a symmetry": aP)
        #expect(UnitCell.cubic(a: 4).conventionalCell(tolerance: 0).bravais == .hR)
        #expect(UnitCell.tetragonal(a: 3, c: 5).conventionalCell(tolerance: 0).bravais == .mC)
        for bad in [Double.nan, -1, .infinity] {
            #expect(needs.conventionalCell(tolerance: bad) == needs.conventionalCell(tolerance: 0), "\(bad)")
            #expect(hex.conventionalCell(tolerance: bad) == hex.conventionalCell(tolerance: 0), "\(bad)")
            #expect(UnitCell.cubic(a: 4).conventionalCell(tolerance: bad) == UnitCell.cubic(a: 4).conventionalCell(tolerance: 0), "\(bad)")
            #expect(UnitCell.tetragonal(a: 3, c: 5).conventionalCell(tolerance: bad).bravais == .mC, "\(bad)")
        }
    }
}
