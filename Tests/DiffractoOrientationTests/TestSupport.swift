//
//  TestSupport.swift
//  Role: the structures and scan fixtures the slice 4 tests build. The structures are the cells and sites the planted-truth
//        tests of vector matching use (an fcc cell, the hexagonal and monoclinic cells, the published beta-double-prime
//        refinement), built as `CrystalStructure` with the identity operation only: a vector-matching library reads the
//        kinematic reflections, not the symmetry. The scan fixtures wrap the kit's `[[SIMD2<Double>]]` scan form.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering
import DiffractoOrientation

/// A cell with its sites and the identity operation.
func makeCrystal(a: Double, b: Double, c: Double,
                 alphaDeg: Double = 90, betaDeg: Double = 90, gammaDeg: Double = 90,
                 sites: [AtomSite]) -> CrystalStructure {
    CrystalStructure(cell: UnitCell(a: a, b: b, c: c, alphaDeg: alphaDeg, betaDeg: betaDeg, gammaDeg: gammaDeg),
                     sites: sites, operations: [.identity])
}

extension CrystalStructure {
    /// Face-centred cubic (4 atoms).
    static func fcc(a: Double, z: Int) -> CrystalStructure {
        makeCrystal(a: a, b: a, c: a, sites: [
            AtomSite(z: z, fractional: [0, 0, 0]),
            AtomSite(z: z, fractional: [0.5, 0.5, 0]),
            AtomSite(z: z, fractional: [0.5, 0, 0.5]),
            AtomSite(z: z, fractional: [0, 0.5, 0.5]),
        ])
    }

    /// Simple cubic (1 atom).
    static func sc(a: Double, z: Int) -> CrystalStructure {
        makeCrystal(a: a, b: a, c: a, sites: [AtomSite(z: z, fractional: [0, 0, 0])])
    }

    /// Pure aluminium at room temperature (4.0495 Å): a 4.04 Å cell moves 5 of 29 241 positions of one measured scan,
    /// because the {220} ring straddles k-max = 0.70 Å⁻¹.
    static var aluminum: CrystalStructure { fcc(a: 4.0495, z: 13) }

    /// β″ (Mg₅Si₆), the canonical experimental refinement: S. J. Andersen, H. W. Zandbergen, J. Jansen, C. Traeholt,
    /// U. Tundal, O. Reiso, Acta Materialia 46(9), 3283-3298 (1998), doi:10.1016/S1359-6454(97)00493-X; cell from the
    /// abstract, coordinates from Table 3, set 3 (the C2/m refinement). The six published sites are expanded under C2/m
    /// (unique axis b) so the expansion can be read; it must give Mg₁₀Si₁₂ = 22 atoms, the cell content the paper states.
    static var betaDoublePrime: CrystalStructure {
        func fourI(_ z: Int, _ x: Double, _ zc: Double) -> [AtomSite] {
            [SIMD3(x, 0, zc), SIMD3(-x, 0, -zc),
             SIMD3(x + 0.5, 0.5, zc), SIMD3(-x + 0.5, 0.5, -zc)]
                .map { AtomSite(z: z, fractional: $0) }
        }
        let mg = 12, si = 14
        var sites: [AtomSite] = [
            AtomSite(z: mg, fractional: [0, 0, 0]),
            AtomSite(z: mg, fractional: [0.5, 0.5, 0]),
        ]
        sites += fourI(mg, 0.3459, 0.089)
        sites += fourI(mg, 0.430, 0.652)
        sites += fourI(si, 0.0565, 0.649)
        sites += fourI(si, 0.1885, 0.224)
        sites += fourI(si, 0.2171, 0.617)
        return makeCrystal(a: 15.16, b: 4.05, c: 6.74, alphaDeg: 90, betaDeg: 105.3, gammaDeg: 90, sites: sites)
    }
}

/// The reflections a library is made from (`try!`: the fixtures are valid structures).
func refl(_ structure: CrystalStructure, _ kMax: Double) -> [KinematicReflection] {
    try! PhaseReferenceLibrary.reflections(of: structure, kMaxInvAngstrom: kMax)
}

/// A peak at a detector position, as the pixel form of a scan holds it.
func pixel(_ x: Float, _ y: Float) -> SIMD2<Double> { SIMD2(Double(x), Double(y)) }

/// A scan as the tests build it: row-major positions, each an array of peak positions.
struct ScanFixture {
    var scanWidth: Int
    var scanHeight: Int
    var peaks: [[SIMD2<Double>]]
}

extension PhaseVectorMatcher {
    /// `map` on a fixture, pixel form.
    static func map(scan: ScanFixture, library: PhaseReferenceLibrary, settings: PhaseVectorSettings,
                    originX: Float, originY: Float, invAngstromPerPixel: Double,
                    matrixEntryIndex: Int? = nil, maximumWorkers: Int? = nil,
                    progress: (@Sendable (Double) -> Void)? = nil,
                    isCancelled: (@Sendable () -> Bool)? = nil) -> PhaseMap? {
        map(scan: scan.peaks, scanWidth: scan.scanWidth, scanHeight: scan.scanHeight, library: library,
            settings: settings, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
            matrixEntryIndex: matrixEntryIndex, maximumWorkers: maximumWorkers, progress: progress,
            isCancelled: isCancelled)
    }

    /// `fitZoneAxis` on a fixture, pixel form.
    static func fitZoneAxis(scan: ScanFixture, crystal: CrystalStructure,
                            referenceSettings: PhaseReferenceSettings, settings: PhaseVectorSettings,
                            originX: Float, originY: Float, invAngstromPerPixel: Double,
                            candidateAxes: [SIMD3<Int>] = PhaseReferenceLibrary.lowIndexZoneAxes,
                            inPlaneStepDeg: Double = 5, sampleLimit: Int = 400) -> [ZoneAxisFit] {
        try! fitZoneAxis(scan: scan.peaks, scanWidth: scan.scanWidth, scanHeight: scan.scanHeight, crystal: crystal,
                         referenceSettings: referenceSettings, settings: settings, originX: originX, originY: originY,
                         invAngstromPerPixel: invAngstromPerPixel, candidateAxes: candidateAxes,
                         inPlaneStepDeg: inPlaneStepDeg, sampleLimit: sampleLimit)
    }
}

/// The symmetry a pin file's class symbol names; "-3m1" is `-3m` in the twofolds-on-a setting (the pin files list it last).
func pinnedSymmetry(_ symbol: String) -> CrystalOrientationSymmetry? {
    if symbol == "-3m1" { return CrystalOrientationSymmetry(laueClass: .threeBarM, setting: .threeBarMOne) }
    return LaueClass(rawValue: symbol).map { CrystalOrientationSymmetry(laueClass: $0) }
}
