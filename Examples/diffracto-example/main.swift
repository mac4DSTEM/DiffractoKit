// DiffractoKit, the first ten minutes (this is the README's program, character for character).
//
// Gold from a CIF -> its rings after the |F| cut -> a synthetic ring pattern at a known camera length -> the ring fitted
// and the camera calibrated back. The CIF below is written out in the program: the cubic cell of gold, a = 4.0782 A (the
// room-temperature lattice parameter, the value used throughout the kit's tests) and the four atoms of
// the face-centred basis. Run it with:  swift run -c release diffracto-example
import Foundation
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

let cifText: String = """
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

func fixed(_ x: Double, _ digits: Int) -> String { String(format: "%.\(digits)f", x) }
func sci(_ x: Double) -> String { String(format: "%.4e", x) }

do {
    // 1. The structure. A CIF that cannot be trusted throws an error that names the tag or value at fault.
    let gold: CrystalStructure = try CIFImport.read(cifText).structure
    let laue: String = gold.laueClass?.symbol ?? "?"
    print("\(gold.formula): a = \(fixed(gold.cell.a, 4)) A, Laue class \(laue), \(gold.operations.count) operations")

    // 2. The rings to k = 1/d = 1.2 per A after the |F| cut (default 1e-4): d, multiplicity and the summed |F|^2.
    let list: KinematicReflectionList = try gold.kinematicReflections(kMax: 1.2)
    print("\(list.families.count) rings, \(list.kept.count) reflections:")
    for ring in list.families.prefix(6) {
        let name: String = ring.family.name.padding(toLength: 9, withPad: " ", startingAt: 0)
        let d: String = fixed(ring.family.d, 4)
        print("  \(name) d = \(d) A  multiplicity \(ring.family.multiplicity)  sum|F|^2 = \(sci(ring.intensity)) A^-4")
    }

    // 3. A synthetic {111} ring. A camera length of 300 mm on 0.05 mm pixels is 6000 px, the electron wavelength at
    //    200 kV is lambda, and the ring sits at R = L tan(2 theta). The kit has to find R, and from it dq and L, back.
    let volts: Double = 200_000
    let pitchMM: Double = 0.05
    let lambda: Double = ElectronWavelength.angstroms(volts: volts)
    let d111: Double = list.families[0].family.d
    let trueLength: Double = 300.0 / pitchMM
    let twoTheta: Double = RingCalibration.twoTheta(wavelengthAngstroms: lambda, d: d111)
    let trueRadius: Double = trueLength * tan(twoTheta)
    let trueDeltaQ: Double = (1.0 / d111) / trueRadius
    print("lambda(200 kV) = \(fixed(lambda, 6)) A")

    let size: Int = 200
    let trueRow: Double = 100.4
    let trueColumn: Double = 99.6
    let sigma: Double = 2.2
    var pixels: [Float] = [Float](repeating: 0, count: size * size)
    for row in 0..<size {
        for column in 0..<size {
            let dr: Double = Double(row) - trueRow
            let dc: Double = Double(column) - trueColumn
            let x: Double = (dr * dr + dc * dc).squareRoot() - trueRadius
            let value: Double = exp(-x * x / (2.0 * sigma * sigma))
            pixels[row * size + column] = Float(value)
        }
    }
    let pattern: DiffractionPattern = DiffractionPattern(width: size, height: size, pixels: pixels)

    // 4. Fit the ring (start centre one click away; the annulus brackets the ring), then calibrate with its spacing.
    let ring: RingRadiusFit = try RingCalibration.fitRing(pattern: pattern, centerQX: trueRow + 1.0,
                                                          centerQY: trueColumn - 0.8, innerRadius: trueRadius - 20.0,
                                                          outerRadius: trueRadius + 20.0)
    let result: RingFitCalibration = try RingCalibration.calibrating(ring, d: d111, highTensionVolts: volts,
                                                                     pixelPitchMM: pitchMM)
    if let fit = result.fit {
        let radius: Double = (fit.a + fit.b) / 2.0
        print("fit: centre (\(fixed(fit.centerQX, 3)), \(fixed(fit.centerQY, 3))) [true (\(trueRow), \(trueColumn))], "
              + "R = \(fixed(radius, 3)) px [true \(fixed(trueRadius, 3))]")
    }
    if result.refusals.isEmpty, let calibration = result.calibration {
        let length: Double = calibration.lengthPixels ?? .nan
        let millimetres: Double = calibration.lengthMillimetres ?? .nan
        print("dq = \(sci(calibration.deltaQ)) per A per px [true \(sci(trueDeltaQ))]")
        print("L = \(fixed(length, 1)) px = \(fixed(millimetres, 2)) mm [true \(fixed(trueLength, 1)) px = 300.00 mm]")
    } else {
        print("not calibrated; refusals: \(result.refusals.map(\.name))")
    }
} catch {
    print("error: \(error.localizedDescription)")
    exit(1)
}
