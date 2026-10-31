//
//  PolicyTests.swift
//  DiffractoKit 0.2.0: the policies of DiffractoCore are parameters with the suite's values
//  as defaults. The default path is held by the existing suite; here a non-default value moves the outcome in the
//  stated direction and an invalid one is refused or taken as documented. Every test names under "RED UNDER" the production
//  mutation that turns a default-path test red.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore

// MARK: - ReflectionList.Policy

struct ReflectionListPolicyTests {
    typealias Policy = ReflectionList.Policy

    static func aluminium(_ cell: UnitCell = .cubic(a: 3)) -> CrystalStructure {
        CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 13, fractional: [0, 0, 0])])
    }

    /// The documented defaults are the constants of the list.
    /// RED UNDER: any default of `ReflectionList.Policy` changed.
    @Test func theDefaultsAreTheConstants() {
        let p: Policy = .default
        #expect(p.maximumKMax == 4)
        #expect(p.maximumReflections == 10_000_000)
        #expect(p.maximumEnumeration == 100_000_000)
        #expect(p.familyTolerance == 1e-9)
        #expect(p.absenceTolerance == 1e-9)
        #expect(p.maximumKMax == ReflectionList.maximumKMax)
        #expect(p.familyTolerance == ReflectionList.familyTolerance)
        #expect(p.maximumReflections == ReflectionList.maximumReflections)
    }

    /// The k_max cap: 5 Å⁻¹ is refused at the default cap of 4 and built at a cap of 6; 2 Å⁻¹ is built at the default and
    /// refused at a cap of 1. A cap that is not a positive number refuses every k_max.
    /// RED UNDER: the policy's cap ignored (the 4 of `maximumKMax` read instead).
    @Test func theKMaxCapMovesWhatIsBuilt() throws {
        let structure = Self.aluminium(.cubic(a: 1))
        var refused = false
        do { _ = try structure.reflections(kMax: 5) } catch { guard case .invalidKMax = error else { throw error }; refused = true }
        #expect(refused)
        let raised = try structure.reflections(kMax: 5, policy: Policy(maximumKMax: 6))
        #expect(raised.reflections.count > 0)
        _ = try structure.reflections(kMax: 2)
        var lowered = false
        do { _ = try structure.reflections(kMax: 2, policy: Policy(maximumKMax: 1)) } catch { guard case .invalidKMax = error else { throw error }; lowered = true }
        #expect(lowered)
        for bad in [Double.nan, 0, -4] {
            var rejected = false
            do { _ = try structure.reflections(kMax: 1, policy: Policy(maximumKMax: bad)) } catch { guard case .invalidKMax = error else { throw error }; rejected = true }
            #expect(rejected)
        }
    }

    /// The count caps: a 3 Å cell at 1 Å⁻¹ holds about 113 reflections in a 343-index box; a cap of 50 reflections or of 100
    /// indices refuses it as `cellTooLargeOrInvalid`.
    /// RED UNDER: the reflection cap or the enumeration cap of the policy ignored.
    @Test func theCountCapsMoveWhatIsBuilt() throws {
        let structure = Self.aluminium()
        let built = try structure.reflections(kMax: 1)
        #expect(built.reflections.count > 50)
        var tooMany = false
        do { _ = try structure.reflections(kMax: 1, policy: Policy(maximumReflections: 50)) }
        catch { guard case .cellTooLargeOrInvalid = error else { throw error }; tooMany = true }
        #expect(tooMany)
        var tooWide = false
        do { _ = try structure.reflections(kMax: 1, policy: Policy(maximumEnumeration: 100)) }
        catch { guard case .cellTooLargeOrInvalid = error else { throw error }; tooWide = true }
        #expect(tooWide)
        var negative = false
        do { _ = try structure.reflections(kMax: 1, policy: Policy(maximumReflections: -1)) }
        catch { guard case .cellTooLargeOrInvalid = error else { throw error }; negative = true }
        #expect(negative)
    }

    /// The ring rule: a tetragonal cell with c/a = 1 + 1e-7 has (1 0 0) and (0 0 1) in two rings at the default tolerance and
    /// in one ring at 1e-6; a NaN or negative tolerance is taken as 0 (no ring shares two |g|).
    /// RED UNDER: the policy's family tolerance ignored; its default changed.
    @Test func theFamilyToleranceGroupsTheRings() throws {
        let structure = Self.aluminium(.tetragonal(a: 1, c: 1 + 1e-7))
        let separate = try structure.reflections(kMax: 1.5)
        let joined = try structure.reflections(kMax: 1.5, policy: Policy(familyTolerance: 1e-6))
        #expect(joined.families().count < separate.families().count)
        let none = try structure.reflections(kMax: 1.5, policy: Policy(familyTolerance: 0))
        let nan = try structure.reflections(kMax: 1.5, policy: Policy(familyTolerance: .nan))
        let negative = try structure.reflections(kMax: 1.5, policy: Policy(familyTolerance: -1))
        #expect(nan.families() == none.families())
        #expect(negative.families() == none.families())
        #expect(none.families().count >= separate.families().count)
        let direct = ReflectionList.families(of: separate.reflections, laueGroup: separate.laueGroup, millerBravais: false,
                                             familyTolerance: 1e-6)
        #expect(direct.count == joined.families().count)
    }

    /// The absence tolerance: (2 0 0) under a translation of 1/2 + 1e-7 has a phase 2e-7 from an integer: absent at 1e-9,
    /// present at 1e-6. Through the list and through `isAllowed`.
    /// RED UNDER: the absence tolerance ignored (the 1e-9 literal read instead).
    @Test func theAbsenceToleranceDecidesWhatIsForbidden() throws {
        let translation = SymmetryOperation(rotation: .identity, translation: SIMD3(0.5 + 1e-7, 0, 0))
        let hkl = MillerIndex(2, 0, 0)
        #expect(!ReflectionList.isAllowed(hkl, by: [translation]))
        #expect(ReflectionList.isAllowed(hkl, by: [translation], tolerance: 1e-6))
        #expect(!ReflectionList.isAllowed(hkl, by: [translation], tolerance: .nan))
        #expect(!ReflectionList.isAllowed(hkl, by: [translation], tolerance: -1))
        let structure = CrystalStructure(cell: .cubic(a: 3), sites: [AtomSite(z: 13, fractional: [0, 0, 0])],
                                         operations: [.identity, translation])
        let strict = try structure.reflections(kMax: 1)
        let tolerant = try structure.reflections(kMax: 1, policy: Policy(absenceTolerance: 1e-6))
        let a = try #require(strict.reflections.first { $0.hkl == hkl })
        let b = try #require(tolerant.reflections.first { $0.hkl == hkl })
        #expect(!a.isAllowed)
        #expect(b.isAllowed)
    }
}

// MARK: - CrystalStructure.ValidationPolicy

struct ValidationPolicyTests {
    typealias Policy = CrystalStructure.ValidationPolicy

    static func pair(_ first: Double, _ second: Double, gap: Double = 0.075) -> CrystalStructure {
        CrystalStructure(cell: .cubic(a: 4),
                         sites: [AtomSite(z: 13, fractional: [0, 0, 0], occupancy: first),
                                 AtomSite(z: 13, fractional: [gap, 0, 0], occupancy: second)],
                         operations: [.identity])
    }

    /// The documented defaults.
    /// RED UNDER: either default of `ValidationPolicy` changed.
    @Test func theDefaultsAreTheDocumentedValues() {
        let p: Policy = .default
        #expect(p.closeContactLimit == 0.5)
        #expect(p.splitSiteOccupancyTolerance == 0.02)
        #expect(p.closeContactLimit == CrystalStructure.closeContactLimit)
    }

    /// The close-contact limit: two full sites 0.3 Å apart are a contact at the default 0.5 Å and none at 0.2; two sites
    /// 0.3 Å apart are a contact at a limit of 0.35. A limit that is not a positive number finds no contact.
    /// RED UNDER: the policy's close-contact limit ignored.
    @Test func theCloseContactLimitMovesTheIssue() throws {
        let structure = Self.pair(1, 1)
        let closest = try #require(structure.shortestCloseContact)
        #expect(abs(closest - 0.3) < 1e-9)
        #expect(structure.validationIssues.map(\.code) == ["atomic_sites_too_close"])
        let relaxed = Policy(closeContactLimit: 0.2)
        #expect(structure.shortestCloseContact(using: relaxed) == nil)
        #expect(structure.validationIssues(using: relaxed).isEmpty)
        #expect(structure.isUsable(under: relaxed))
        #expect(!structure.isUsable)
        let strict = Policy(closeContactLimit: 0.35)
        #expect(structure.shortestCloseContact(using: strict) != nil)
        let wide = Self.pair(1, 1, gap: 0.2)
        #expect(wide.shortestCloseContact == nil)
        #expect(wide.shortestCloseContact(using: Policy(closeContactLimit: 1.0)) != nil)
        for bad in [Double.nan, 0, -1] {
            #expect(structure.shortestCloseContact(using: Policy(closeContactLimit: bad)) == nil)
        }
    }

    /// The split-site rule: occupancies 0.5 and 0.51 (sum 1.01) are a split site by default and a contact at a tolerance of 0;
    /// 0.51 and 0.52 (sum 1.03) are a contact by default and a split site at 0.05. A NaN or negative tolerance is taken as 0.
    /// RED UNDER: the policy's split-site tolerance ignored.
    @Test func theSplitSiteToleranceMovesTheIssue() {
        let split = Self.pair(0.5, 0.51)
        #expect(split.shortestCloseContact == nil)
        #expect(split.shortestCloseContact(using: Policy(splitSiteOccupancyTolerance: 0)) != nil)
        #expect(split.shortestCloseContact(using: Policy(splitSiteOccupancyTolerance: .nan)) != nil)
        #expect(split.shortestCloseContact(using: Policy(splitSiteOccupancyTolerance: -1)) != nil)
        let over = Self.pair(0.51, 0.52)
        #expect(over.shortestCloseContact != nil)
        #expect(over.shortestCloseContact(using: Policy(splitSiteOccupancyTolerance: 0.05)) == nil)
    }
}

// MARK: - SymmetryDetection.Policy and SymmetryOperation's bound

struct DetectionPolicyTests {
    typealias Policy = SymmetryDetection.Policy

    /// The documented defaults.
    /// RED UNDER: any default of `SymmetryDetection.Policy` changed.
    @Test func theDefaultsAreTheDocumentedValues() {
        let p: Policy = .default
        #expect(p.occupancyTolerance == 0.02)
        #expect(p.translationGrid == 1.0 / 24)
        #expect(p.gridSnapFactor == 3)
        #expect(p.precisionLimitFactor == 10)
        #expect(p.occupancyTolerance == SymmetryDetection.occupancyTolerance)
        #expect(p.translationGrid == SymmetryDetection.translationGrid)
    }

    static func bodyCentred(second: Double) -> [AtomSite] {
        [AtomSite(z: 55, fractional: [0, 0, 0], occupancy: 1), AtomSite(z: 55, fractional: [0.5, 0.5, 0.5], occupancy: second)]
    }

    /// The occupancy tolerance: sites of occupancy 1 and 0.97 are one species at a tolerance of 0.05 (the body centring is
    /// found: 96 operations) and two at the default 0.02 (48). A NaN or negative tolerance is taken as 0.
    /// RED UNDER: the policy's occupancy tolerance ignored.
    @Test func theOccupancyToleranceDecidesTheCentring() {
        let cell = UnitCell.cubic(a: 4)
        let sites = Self.bodyCentred(second: 0.97)
        let strict = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        let tolerant = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, policy: Policy(occupancyTolerance: 0.05))
        #expect(tolerant.operations.count > strict.operations.count)
        let equal = CrystalStructure.detectingSymmetry(cell: cell, sites: Self.bodyCentred(second: 1))
        #expect(equal.operations.count == tolerant.operations.count)
        for bad in [Double.nan, -1] {
            let taken = CrystalStructure.detectingSymmetry(cell: cell, sites: Self.bodyCentred(second: 1),
                                                           policy: Policy(occupancyTolerance: bad))
            #expect(taken.operations.count == equal.operations.count)
            let differing = CrystalStructure.detectingSymmetry(cell: cell, sites: sites, policy: Policy(occupancyTolerance: bad))
            #expect(differing.operations.count == strict.operations.count)
        }
    }

    /// The translation grid: four sites on a line at x = 0.1, 0.2335, 0.2, 0.1331 (written to three decimals) make the mirror
    /// x -> -x + t with t = 0.3333 by least squares (the pairs say 0.3335 and 0.3331); it snaps to 1/3 within three of its
    /// uncertainties at the default grid, and stays 0.3333 at a grid of 1/5, a snap factor of 0, or an invalid grid.
    /// RED UNDER: the policy's grid or snap factor ignored.
    @Test func theTranslationGridSnapsATranslation() {
        let xs: [Double] = [0.1, 0.2335, 0.2, 0.1331]
        let sites: [AtomSite] = xs.map { AtomSite(z: 13, fractional: SIMD3($0, 0, 0)) }
        func translations(_ policy: Policy) -> [Double] {
            let found = SymmetryDetection.detect(sites: sites, cell: .orthorhombic(a: 6, b: 4, c: 5), policy: policy)
            let mirrors = found.operations.filter { $0.rotation.rows[0].x == -1 && $0.rotation.rows[1].y == 1 }
            return mirrors.map { $0.translation.x }
        }
        func snapped(_ values: [Double]) -> Bool { values.contains { abs($0 - 1.0 / 3) < 1e-12 } }
        func raw(_ values: [Double]) -> Bool { values.contains { abs($0 - 0.3333) < 1e-9 } }
        let atDefault = translations(.default)
        #expect(snapped(atDefault))
        #expect(!raw(atDefault))
        let coarse = translations(Policy(translationGrid: 0.2))
        #expect(raw(coarse))
        #expect(!snapped(coarse))
        let unsnapped = translations(Policy(gridSnapFactor: 0))
        #expect(raw(unsnapped))
        #expect(!snapped(unsnapped))
        for bad in [0.0, -1.0, Double.nan] {
            let invalid = translations(Policy(translationGrid: bad))
            #expect(raw(invalid))
            #expect(!snapped(invalid))
        }
    }

    /// The precision-limit factor: coordinates written to 0.05 report the rotations that fail the position tolerance but
    /// pass at ten half-steps as unresolved; a factor of 0 reports none, and 20 reports at least as many.
    /// RED UNDER: the policy's factor ignored (the 10 literal read instead).
    @Test func thePrecisionLimitFactorReportsUnresolvedRotations() {
        let sites: [AtomSite] = [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 8, fractional: [0.2, 0.3, 0.1])]
        func unresolved(_ policy: Policy) -> Int {
            SymmetryDetection.detect(sites: sites, cell: .cubic(a: 4), coarsestHalfStep: 0.05, policy: policy)
                .unresolvedRotations.count
        }
        let atDefault = unresolved(.default)
        #expect(atDefault > 0)
        #expect(unresolved(Policy(precisionLimitFactor: 0)) == 0)
        #expect(unresolved(Policy(precisionLimitFactor: .nan)) == 0)
        #expect(unresolved(Policy(precisionLimitFactor: 20)) >= atDefault)
    }

    /// The coefficient bound of a symmetry-operation text: 1500 is unparsable at the default 1000 and parsed at 2000; 0.5 is
    /// unparsable at a bound of 0.4. A NaN or negative bound is taken as 0.
    /// RED UNDER: the bound ignored (the 1000 literal read instead).
    @Test func theOperationMagnitudeBoundIsAParameter() throws {
        #expect(SymmetryOperation.defaultMaximumMagnitude == 1000)
        var refused = false
        do { _ = try SymmetryOperation(xyz: "1500+x,y,z") } catch { guard case .unparsable = error else { throw error }; refused = true }
        #expect(refused)
        let parsed = try SymmetryOperation(xyz: "1500+x,y,z", maximumMagnitude: 2000)
        #expect(parsed.translation.x == 1500)
        _ = try SymmetryOperation(xyz: "1/2+x,y,z")
        for bound in [0.4, Double.nan, -1] {
            var tight = false
            do { _ = try SymmetryOperation(xyz: "1/2+x,y,z", maximumMagnitude: bound) } catch { guard case .unparsable = error else { throw error }; tight = true }
            #expect(tight)
        }
        _ = try SymmetryOperation(xyz: "x,y,z", maximumMagnitude: 0)
    }
}

// MARK: - CIFImport.Policy

struct CIFImportPolicyTests {
    typealias Policy = CIFImport.Policy

    /// A 10 Å cubic cell with one site written to one decimal and a mirror at x = 0.304 that moves it by 0.008 (0.08 Å; a
    /// closed list, as the import requires): the two images are the same site within the rounding error and 0.10 Å.
    static let nearlyCoincident = """
        data_near
        _cell_length_a 10.0
        _cell_length_b 10.0
        _cell_length_c 10.0
        _cell_angle_alpha 90.0
        _cell_angle_beta 90.0
        _cell_angle_gamma 90.0
        loop_
        _symmetry_equiv_pos_as_xyz
          'x,y,z'
          '-x+0.608,y,z'
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
          Fe1 Fe 0.3 0.3 0.3

        """

    static let bodyCentredPartial = """
        data_cs
        _cell_length_a 4.0
        _cell_length_b 4.0
        _cell_length_c 4.0
        _cell_angle_alpha 90.0
        _cell_angle_beta 90.0
        _cell_angle_gamma 90.0
        loop_
        _atom_site_label
        _atom_site_type_symbol
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        _atom_site_occupancy
          Cs1 Cs 0.0 0.0 0.0 1.0
          Cs2 Cs 0.5 0.5 0.5 0.97

        """

    /// The documented defaults.
    /// RED UNDER: any default of `CIFImport.Policy` changed.
    @Test func theDefaultsAreTheDocumentedValues() {
        let p: Policy = .default
        #expect(p.maximumMergeDistance == 0.10)
        #expect(p.minimumMergeTolerance == 1e-6)
        #expect(p.maximumCoordinateHalfStep == 1e-2)
        #expect(p.rightAngleToleranceDeg == 0.05)
        #expect(p.maximumOperationMagnitude == 1000)
        #expect(p.operationMetricTolerance == 0.01)
        #expect(p.positionTolerance == .default)
        #expect(p.metricLengthTolerance == 1e-9)
        #expect(p.metricAngleToleranceDeg == 1e-9)
        #expect(p.symmetry == .default)
        #expect(p.validation == .default)
    }

    /// A read given `.default` equals a read given none.
    @Test func aReadGivenTheDefaultEqualsAReadGivenNone() throws {
        let plain = try CIFImport.read(Self.nearlyCoincident)
        let withPolicy = try CIFImport.read(Self.nearlyCoincident, policy: .default)
        #expect(plain == withPolicy)
        let beta = try CIFImport.read(Fixtures.betaDoublePrimeCIF)
        let betaPolicy = try CIFImport.read(Fixtures.betaDoublePrimeCIF, policy: .default)
        #expect(beta == betaPolicy)
    }

    /// The merge distance: the two images 0.08 Å apart merge at the default 0.10 Å (one site) and stay apart at 0.05 Å, where
    /// the close-contact rule then refuses the structure; with the contact limit lowered too, both sites are kept. A NaN or
    /// negative distance merges only coincident images.
    /// RED UNDER: the policy's merge distance ignored (the 0.10 literal read instead).
    @Test func theMergeDistanceDecidesWhichImagesAreOneSite() throws {
        let merged = try CIFImport.read(Self.nearlyCoincident)
        #expect(merged.structure.sites.count == 1)
        var refused = false
        do { _ = try CIFImport.read(Self.nearlyCoincident, policy: Policy(maximumMergeDistance: 0.05)) }
        catch { guard case .invalidStructure = error else { throw error }; refused = true }
        #expect(refused)
        var both = Policy(maximumMergeDistance: 0.05)
        both.validation = CrystalStructure.ValidationPolicy(closeContactLimit: 0.05)
        let kept = try CIFImport.read(Self.nearlyCoincident, policy: both)
        #expect(kept.structure.sites.count == 2)
        for bad in [Double.nan, -1] {
            var rejected = false
            do { _ = try CIFImport.read(Self.nearlyCoincident, policy: Policy(maximumMergeDistance: bad)) }
            catch { guard case .invalidStructure = error else { throw error }; rejected = true }
            #expect(rejected)
        }
    }

    /// The operation magnitude bound reaches the loop: a bound of 0.001 refuses the 0.008 translation as unsupported.
    /// RED UNDER: the policy's bound not passed to the operation parser.
    @Test func theOperationBoundReachesTheLoop() {
        var refused = false
        do { _ = try CIFImport.read(Self.nearlyCoincident, policy: Policy(maximumOperationMagnitude: 0.001)) }
        catch { if case .unsupportedSymmetryOperation = error { refused = true } }
        #expect(refused)
    }

    /// The symmetry policy reaches the detection: sites of occupancy 1 and 0.97 are centred at an occupancy tolerance of 0.05.
    /// RED UNDER: the policy's symmetry policy not passed to `detectingSymmetry`.
    @Test func theSymmetryPolicyReachesTheDetection() throws {
        let plain = try CIFImport.read(Self.bodyCentredPartial)
        var tolerant = Policy()
        tolerant.symmetry = SymmetryDetection.Policy(occupancyTolerance: 0.05)
        let centred = try CIFImport.read(Self.bodyCentredPartial, policy: tolerant)
        #expect(centred.structure.operations.count > plain.structure.operations.count)
    }

    /// The half-step cap: 0.3 (one decimal) is uncertain by the cap, 1e-2 by default; 0.3333 by 1e-4 under any larger cap
    /// and by the cap under a smaller one; exact dyadic values (0.25) by 0 whatever the cap; a NaN or negative cap is 0.
    /// RED UNDER: the cap argument ignored (the 1e-2 literal read instead).
    @Test func theHalfStepCapBoundsTheWrittenUncertainty() {
        #expect(CIFImport.coordinateHalfStep("0.3") == 1e-2)
        #expect(CIFImport.coordinateHalfStep("0.3", cap: 1e-3) == 1e-3)
        #expect(CIFImport.coordinateHalfStep("0.3333") == 1e-4)
        #expect(CIFImport.coordinateHalfStep("0.3333", cap: 5e-2) == 1e-4)
        #expect(CIFImport.coordinateHalfStep("0.3333", cap: 1e-5) == 1e-5)
        #expect(CIFImport.coordinateHalfStep("0.25", cap: 0.5) == 0)
        #expect(CIFImport.coordinateHalfStep("0.3", cap: .nan) == 0)
        #expect(CIFImport.coordinateHalfStep("0.3", cap: -1) == 0)
    }

    /// The right-angle tolerance of the centring judgement: space group 5 (C2) in a cell with alpha = gamma = 90.04 degrees
    /// lists 4 operations when those angles count as right (0.05) and 2 (primitive) at a tolerance of 0.01. A NaN tolerance is 0.
    /// RED UNDER: the tolerance argument ignored (the 0.05 literal read instead).
    @Test func theRightAngleToleranceJudgesTheCentring() {
        let centred = CIFImport.requiredOperationCount(spaceGroup: 5, alphaDeg: 90.04, betaDeg: 100, gammaDeg: 90.04)
        #expect(centred == 4)
        let strict = CIFImport.requiredOperationCount(spaceGroup: 5, alphaDeg: 90.04, betaDeg: 100, gammaDeg: 90.04,
                                                      rightAngleToleranceDeg: 0.01)
        #expect(strict == 2)
        let nan = CIFImport.requiredOperationCount(spaceGroup: 5, alphaDeg: 90.04, betaDeg: 100, gammaDeg: 90.04,
                                                   rightAngleToleranceDeg: .nan)
        #expect(nan == 2)
        let loose = CIFImport.requiredOperationCount(spaceGroup: 5, alphaDeg: 90.4, betaDeg: 100, gammaDeg: 90.4,
                                                     rightAngleToleranceDeg: 0.5)
        #expect(loose == 4)
    }
}
