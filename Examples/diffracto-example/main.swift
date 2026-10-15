//
//  main.swift: DiffractoKit's first ten minutes, as a program (README, "First ten minutes").
//
//  1. Read a CIF (gold, written out below) into a `CrystalStructure`: cell, atoms, detected symmetry.
//  2. List the rings the structure gives at an electron-diffraction scale: d, multiplicity and the summed |F|^2.
//  3. Compute the relativistic electron wavelength at 200 kV.
//  4. Render a synthetic gold {111} ring on a detector, fit it, and calibrate: the pixel size in reciprocal space and
//     the camera length, with the refusals the kit would give for a ring it does not trust.
//
//  Built only on the public API of DiffractoCore, DiffractoScattering and DiffractoGeometry. No file is read: the gold
//  structure is the string below (a, 4.0782 A, is a measurement of gold, cited as such; the four atoms are the fcc basis).
//
//  Run:  swift run -c release diffracto-example
//

import Foundation
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

let goldCIF = """
data_gold
_cell_length_a 4.0782
_cell_length_b 4.0782
_cell_length_c 4.0782
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Au1 Au 0.0 0.0 0.0
Au2 Au 0.5 0.5 0.0
Au3 Au 0.5 0.0 0.5
Au4 Au 0.0 0.5 0.5
"""

func f(_ x: Double, _ digits: Int) -> String { String(format: "%.\(digits)f", x) }

do {
    // 1. The structure. A CIF that cannot be trusted throws a `CIFImportError` that names the tag or value at fault.
    let document = try CIFImport.read(goldCIF)
    let gold = document.structure
    print("\(gold.formula): a = \(f(gold.cell.a, 4)) A, V = \(f(gold.cell.volume, 3)) A^3, "
          + "Laue class \(gold.laueClass?.symbol ?? "?") (\(gold.crystalSystem?.description ?? "?")), "
          + "\(gold.operations.count) operations")

    // 2. The rings to k = 1/d = 1.2 per A, after the |F| cut the kit applies by default (StructureFactors.defaultCut).
    let list = try gold.kinematicReflections(kMax: 1.2)
    print("\(list.families.count) rings, \(list.kept.count) reflections:")
    for ring in list.families.prefix(6) {
        let family = ring.family
        print("  \(family.name.padding(toLength: 9, withPad: " ", startingAt: 0)) d = \(f(family.d, 4)) A  "
              + "multiplicity \(family.multiplicity)  sum|F|^2 = \(String(format: "%.4e", ring.intensity)) A^-4")
    }

    // 3. The wavelength: relativistic, CODATA 2018, in angstroms for a voltage in volts.
    let volts = 200_000.0
    let lambda = ElectronWavelength.angstroms(volts: volts)
    print("lambda(200 kV) = \(f(lambda, 6)) A")

    // 4. A synthetic gold {111} ring: a camera length of 300 mm on 0.05 mm pixels is 6000 px, and the ring sits at
    //    R = L tan(2 theta) with 2 theta = 2 asin(lambda / 2d). The kit has to find R, and from it Δq and L, back.
    let d111 = list.families[0].family.d
    let truthLengthPixels = 6000.0
    let truthRadius = truthLengthPixels * tan(RingCalibration.twoTheta(wavelengthAngstroms: lambda, d: d111))
    let rows = 200, columns = 200
    let centreRow = 100.4, centreColumn = 99.6          // the true centre, not on a pixel
    let sigma = 2.2                                      // the ring's width in pixels
    var pixels = [Float](repeating: 0, count: rows * columns)
    for r in 0..<rows {
        for c in 0..<columns {
            let radius = ((Double(r) - centreRow) * (Double(r) - centreRow)
                          + (Double(c) - centreColumn) * (Double(c) - centreColumn)).squareRoot()
            let x = radius - truthRadius
            pixels[r * columns + c] = Float(exp(-x * x / (2 * sigma * sigma)))
        }
    }
    let pattern = DiffractionPattern(width: columns, height: rows, pixels: pixels)

    // The start centre is "a click away" (1 px off); the annulus brackets the ring. 0.05 mm is the detector's pixel pitch.
    let result = try RingCalibration.fit(pattern: pattern, centerQX: centreRow + 1.0, centerQY: centreColumn - 0.8,
                                         innerRadius: truthRadius - 20, outerRadius: truthRadius + 20,
                                         d: d111, highTensionVolts: volts, pixelPitchMM: 0.05)
    if let fit = result.fit {
        print("fit: \(fit.model.rawValue), centre (\(f(fit.centerQX, 3)), \(f(fit.centerQY, 3))) "
              + "[true (\(centreRow), \(centreColumn))], R = \(f((fit.a + fit.b) / 2, 3)) px [true \(f(truthRadius, 3))]")
    }
    if result.refusals.isEmpty, let calibration = result.calibration {
        print("calibrated: dq = \(String(format: "%.6e", calibration.deltaQ)) per A per px, "
              + "L = \(f(calibration.lengthPixels ?? .nan, 1)) px = \(f(calibration.lengthMillimetres ?? .nan, 2)) mm "
              + "[true \(f(truthLengthPixels, 1)) px]")
    } else {
        print("not calibrated; refusals: \(result.refusals.map(\.name))")
    }
} catch {
    print("error: \(error.localizedDescription)")
    exit(1)
}
