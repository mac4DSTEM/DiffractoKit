//
//  ConventionalSymmetryTests.swift
//  The symmetry of any cell setting: the
//  operations, Laue class, Bravais lattice and rings of every pinned structure in the file's setting and in the
//  conventional setting, against spglib; and the planted truth (primitive fcc, bcc and β″, rutile with its 4-fold along
//  a, magnesium with γ = 60°, corundum in rhombohedral axes, a skewed hexagonal lattice).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct ConventionalSymmetryTests {

    static let pins = ConventionalCellPins.shared

    /// Every pinned structure, detected once (default `conventional: true`).
    static let detected: [(pin: ConventionalStructurePin, structure: CrystalStructure)] = pins.structures.map {
        ($0, CrystalStructure.detectingSymmetry(cell: $0.cell.unitCell, sites: $0.atomSites))
    }

    static func named(_ name: String) -> CrystalStructure {
        detected.first { $0.pin.name == name }!.structure
    }

    // MARK: Both settings against spglib

    @Test func testOperationsLaueClassAndBravaisInBothSettingsAgreeWithSpglib() throws {
        for (pin, s) in Self.detected {
            let name = pin.name
            #expect(s.setting == .file, "\(name)")
            let conv = try #require(s.conventionalSetting, "\(name): no conventional setting")
            #expect(conv.bravais.rawValue == pin.lattice.bravais, "\(name) Bravais")
            #expect(conv.structure.setting == .conventional, "\(name)")
            // the file's setting: the operations spglib finds in the cell as given
            #expect(s.operations.count == pin.fileOperationsCount, "\(name) file-setting operations: \(s.operations.count) vs \(pin.fileOperationsCount)")
            #expect(s.laueClass?.symbol == pin.laueClass, "\(name) file-setting Laue class \(String(describing: s.laueClass))")
            // the conventional setting: the file's count times the lattice's centring multiplicity (as spglib's lattice-only cell)
            #expect(conv.operations.count == pin.conventionalOperationsCount, "\(name) conventional operations \(conv.operations.count) vs \(pin.conventionalOperationsCount)")
            #expect(conv.laueClass?.symbol == pin.laueClass, "\(name) conventional Laue class")
            #expect(conv.sites.count == s.sites.count * conv.transformation.determinant, "\(name) site count")
            // the operations of the two settings are the same group seen in two bases
            for op in conv.operations {
                let back = try #require(conv.fileOperation(of: op), "\(name): \(op) has no integer image")
                #expect(s.operations.contains { $0.isEquivalent(to: back, translationTolerance: 1e-6) }, "\(name): \(back) missing in the file's setting")
            }
            #expect(Set(s.laueGroup.rotations) == Set(PointSymmetry.laueGroup(of: s.operations).rotations))
        }
    }

    /// The rings are setting-independent in g and multiplicity: the file's setting, the conventional setting and
    /// spglib's standardised cell agree; the names agree where the conventional cell is spglib's (same parameters).
    @Test func testRingsInBothSettingsAgreeWithSpglibsStandardisedCell() throws {
        var namesCompared = 0, namesSkipped = 0, rings = 0
        var differentCells: [String] = []
        for (pin, s) in Self.detected {
            let name = pin.name
            let conv = try #require(s.conventionalSetting).structure
            let sameCell = parameterDeviation(parameters(conv.cell), pin.convStdParams) <= 1e-8
            for (label, structure) in [("file", s), ("conventional", conv)] {
                let list = try structure.reflections(kMax: pin.rings.kMax)
                #expect(list.allowed.count == pin.rings.allowedCount || label == "file" && !(structure.cell.volume > 0),
                        "\(name) \(label) allowed \(list.allowed.count) vs \(pin.rings.allowedCount)")
                let families = list.families()
                #expect(families.count == pin.rings.families.count, "\(name) \(label) rings \(families.count) vs \(pin.rings.families.count)")
                for (mine, want) in zip(families, pin.rings.families) {
                    rings += 1
                    assertClose(mine.g, want.g, relative: 1e-7, "\(name) \(label) ring g")
                    #expect(mine.multiplicity == want.multiplicity, "\(name) \(label) multiplicity of \(want.name)")
                    if label == "conventional" {
                        if sameCell { #expect(mine.name == want.name, "\(name) conventional ring name"); namesCompared += 1 } else { namesSkipped += 1; if !differentCells.contains(name) { differentCells.append(name) } }
                    }
                }
            }
        }
        print("conventional rings: \(rings) rings compared in g and multiplicity; names compared \(namesCompared), skipped \(namesSkipped) (conventional cell differs from spglib's in parameters: \(differentCells))")
        #expect(namesCompared >= 1200)
    }

    // MARK: The transformation between the settings

    @Test func testIndicesAndOperationsTransformBothWays() throws {
        for (pin, s) in Self.detected {
            let conv = try #require(s.conventionalSetting)
            let cell = pin.cell.unitCell
            for hkl in [MillerIndex(1, 0, 0), MillerIndex(1, 2, -1), MillerIndex(-2, 1, 3), MillerIndex(0, 0, 4)] {
                // h' = h P: the same reciprocal vector in the two bases, so the same g
                let hc = conv.conventionalIndices(of: hkl)
                let gFile = cell.gVector(hkl), gConv = conv.cell.gVector(hc)
                #expect(abs(length(gFile) - length(gConv)) <= 1e-9 * length(gFile), "\(pin.name) |g| of \(hkl)")
                #expect(conv.fileIndices(ofConventional: hc) == hkl, "\(pin.name) round trip of \(hkl)")
            }
            // an operation survives the round trip
            for op in s.operations.prefix(6) {
                let there = try #require(conv.conventionalOperation(of: op))
                let back = try #require(conv.fileOperation(of: there))
                #expect(back.isEquivalent(to: op, translationTolerance: 1e-9), "\(pin.name): \(op) → \(there) → \(back)")
            }
        }
    }

    // MARK: Planted truth

    @Test func testPrimitiveCellsOfFccBccAndBetaDoublePrimeDetectTheirConventionalSymmetry() throws {
        let au = Self.named("au_fcc_primitive")
        #expect(au.operations.count == 48 && au.laueClass == .mThreeBarM)
        let auConv = try #require(au.conventionalSetting)
        #expect(auConv.operations.count == 192 && auConv.laueClass == .mThreeBarM && auConv.bravais == .cF)
        #expect(auConv.sites.count == 4 && abs(auConv.cell.a - 4.0782) <= 1e-9)
        // the rings after the cut are the conventional cell's rings exactly: gold's 8 families at 1.2 Å⁻¹, as the
        // conventional P1 list gives
        let gold = Self.named("gold_P1")
        let goldFamilies = try gold.reflections(kMax: 1.2).families()
        let auFamilies = try auConv.structure.reflections(kMax: 1.2).families()
        #expect(goldFamilies.count == 8 && auFamilies.count == 8)
        for (x, y) in zip(goldFamilies, auFamilies) {
            #expect(x.name == y.name && x.multiplicity == y.multiplicity)
            assertClose(x.g, y.g, relative: 1e-12, "Au ring g")
        }
        // the primitive cell's own rings (file setting) are the same rings by g and multiplicity
        let filePrimitive = try au.reflections(kMax: 1.2).families()
        #expect(filePrimitive.map(\.multiplicity) == goldFamilies.map(\.multiplicity))

        let fe = Self.named("fe_bcc_primitive")
        let feConv = try #require(fe.conventionalSetting)
        #expect(fe.operations.count == 48 && feConv.operations.count == 96 && feConv.laueClass == .mThreeBarM && feConv.bravais == .cI)
        #expect(abs(feConv.cell.a - 2.8665) <= 1e-9 && feConv.sites.count == 2)

        // β″: C-centred primitive cell, 2/m, 8 operations in the conventional cell, the predicted 1448 and 419
        let beta = Self.named("beta_double_prime_primitive")
        let betaConv = try #require(beta.conventionalSetting)
        #expect(beta.operations.count == 4 && betaConv.operations.count == 8 && betaConv.laueClass == .twoOverM && betaConv.bravais == .mC)
        #expect(betaConv.sites.count == 22)
        #expect(abs(betaConv.cell.a - 15.16) <= 1e-9 && abs(betaConv.cell.b - 4.05) <= 1e-9 && abs(betaConv.cell.c - 6.74) <= 1e-9
                && abs(betaConv.cell.betaDeg - 105.3) <= 1e-8, "\(parameters(betaConv.cell))")
        let list = try betaConv.structure.reflections(kMax: 1.2)
        #expect(list.allowed.count == 1448)
        let families = list.families()
        #expect(families.count == 419)
        // and they are the β″ file's own rings, names included (the cell is the file's cell)
        let file = try CIFImport.read(Fixtures.betaDoublePrimeCIF, fileBaseName: "beta").structure
        let fileFamilies = try file.reflections(kMax: 1.2).families()
        #expect(fileFamilies.map(\.name) == families.map(\.name))
        #expect(fileFamilies.map(\.multiplicity) == families.map(\.multiplicity))
    }

    @Test func testNonStandardSettingsOfRutileMagnesiumAndCorundum() throws {
        let rutile = Self.named("rutile_unique_axis_a")
        #expect(rutile.operations.count == 16 && rutile.laueClass == .fourOverMmm)
        let rutileConv = try #require(rutile.conventionalSetting)
        #expect(rutileConv.operations.count == 16 && rutileConv.bravais == .tP && rutileConv.cell.a == rutileConv.cell.b)
        #expect(abs(rutileConv.cell.c - 2.9587) <= 1e-9 && abs(rutileConv.cell.a - 4.5937) <= 1e-9, "\(parameters(rutileConv.cell))")

        let mg = Self.named("magnesium_gamma_60")
        #expect(mg.operations.count == 24 && mg.laueClass == .sixOverMmm)
        let mgConv = try #require(mg.conventionalSetting)
        #expect(mgConv.cell.gammaDeg == 120 && abs(mgConv.cell.a - 3.2094) <= 1e-9 && mgConv.bravais == .hP)
        // the rings are the hexagonal cell's: Miller-Bravais names there, three indices in the γ = 60° file cell
        let mgRings = try mgConv.structure.reflections(kMax: 1.2).families()
        #expect(mgRings.count == 25 && mgRings.first?.name.split(separator: " ").count == 4)
        let hex = try Self.named("magnesium_P63mmc").reflections(kMax: 1.2).families()
        #expect(hex.map(\.name) == mgRings.map(\.name) && hex.map(\.multiplicity) == mgRings.map(\.multiplicity))

        let corundum = Self.named("corundum_rhombohedral_axes")
        #expect(corundum.operations.count == 12 && corundum.laueClass == .threeBarM)
        let cConv = try #require(corundum.conventionalSetting)
        #expect(cConv.operations.count == 36 && cConv.laueClass == .threeBarM && cConv.bravais == .hR)
        #expect(cConv.sites.count == 30 && abs(cConv.cell.a - 4.7602) <= 1e-9 && abs(cConv.cell.c - 12.9933) <= 1e-9)
        let cRings = try cConv.structure.reflections(kMax: 1.2).families()
        let cFirstNameWords: Int? = cRings.first.map { (f: ReflectionFamily) -> Int in f.name.split(separator: " ").count }
        let cFirstName: String? = cRings.first?.name
        #expect(cFirstNameWords == 4, "Miller-Bravais names in the hexagonal axes: \(String(describing: cFirstName))")
        // in the rhombohedral file cell the names keep three indices
        let fileRings = try corundum.reflections(kMax: 1.2).families()
        let fileFirstNameWords: Int? = fileRings.first.map { (f: ReflectionFamily) -> Int in f.name.split(separator: " ").count }
        #expect(fileFirstNameWords == 3)
        #expect(fileRings.map(\.multiplicity) == cRings.map(\.multiplicity))

        let skew = Self.named("triclinic_hexagonal_reduced")
        #expect(skew.operations.count == 24 && skew.laueClass == .sixOverMmm)
        #expect(skew.conventionalSetting?.bravais == .hP)
        let super7 = Self.named("supercell_7x")
        #expect(super7.operations.count == 112 && super7.conventionalSetting?.operations.count == 112)
    }

    /// A file in its conventional setting keeps its numbers and its identity transformation; `conventional: false` is
    /// the plain detection in the file's cell.
    @Test func testConventionalFalseIsSliceOnesDetectionAndConventionalFilesAreUnchanged() throws {
        for (pin, s) in Self.detected where pin.source == "slice-1" {
            let plain = CrystalStructure.detectingSymmetry(cell: pin.cell.unitCell, sites: pin.atomSites, conventional: false)
            #expect(plain.conventionalSetting == nil && plain.setting == .file, "\(pin.name)")
            if s.conventionalSetting?.transformation == .identity {
                #expect(plain.operations == s.operations, "\(pin.name): a conventional file's operations differ")
                #expect(plain.laueClass == s.laueClass)
            } else {
                // the file's cell is not the conventional one: the full detection finds at least what the holohedry test did
                #expect(plain.operations.count <= s.operations.count, "\(pin.name)")
                for op in plain.operations {
                    #expect(s.operations.contains { $0.isEquivalent(to: op, translationTolerance: 1e-6) }, "\(pin.name): \(op)")
                }
            }
        }
        // a primitive fcc cell detects only the holohedry of rhombohedral axes without the reduction
        let au = Self.pins.structures.first { $0.name == "au_fcc_primitive" }!
        let plain = CrystalStructure.detectingSymmetry(cell: au.cell.unitCell, sites: au.atomSites, conventional: false)
        #expect(plain.operations.count == 12 && plain.laueClass == .threeBarM, "rhombohedral-axes holohedry only: \(plain.operations.count)")
    }
}

extension ConventionalSymmetryTests {
    /// Two rows of one element on one position (Fe²⁺ and Fe³⁺ at 0.5 each) in a primitive fcc cell stay two sites each
    /// in the conventional cell: 8 sites, not 4 (carried into the conventional setting).
    @Test func testCoSitedRowsAreNotMergedByTheCentringImages() throws {
        let a = 4.0782 / 2.0.squareRoot()
        let cell = UnitCell.rhombohedral(a: a, alphaDeg: 60)
        let sites = [AtomSite(z: 26, fractional: [0, 0, 0], occupancy: 0.5, label: "Fe2"),
                     AtomSite(z: 26, fractional: [0, 0, 0], occupancy: 0.5, label: "Fe3")]
        let s = CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
        let conv = try #require(s.conventionalSetting)
        #expect(conv.bravais == .cF && conv.sites.count == 8, "\(conv.sites.count)")
        #expect(s.sites.count == 2 && conv.operations.count == 192 && s.operations.count == 48)
        // a centred lattice's images are one per centring translation: the formula of the conventional cell is 4 × the file's
        #expect(conv.structure.formula == "Fe4", "\(conv.structure.formula)")
    }
}
