# Calibrating a camera from one ring

Fit one diffraction ring of known spacing and read the pixel size in reciprocal space and the camera length from it.

## What you need

A ``DiffractionPattern`` (`width`, `height`, `Float` pixels; rows are `qx`, columns are `qy`), a start centre within a few
pixels of the true one, an annulus (inner and outer radius, in pixels) that brackets the ring, the ring's spacing d in
angstroms, and, for a camera length, the high tension in volts and the detector pixel pitch in millimetres.

## Fit, then calibrate

``RingCalibration/fitRing(pattern:centerQX:centerQY:innerRadius:outerRadius:progress:isCancelled:)`` measures the ring's
radius without a spacing, so a ring can be measured before its family is known. The ring's radius is (a + b) / 2 of the
fitted ellipse. ``RingCalibration/calibrating(_:d:highTensionVolts:pixelPitchMM:)`` adds d, the high tension and the
pitch; ``RingCalibration/fit(pattern:centerQX:centerQY:innerRadius:outerRadius:d:highTensionVolts:pixelPitchMM:progress:isCancelled:)``
does both in one call.

```swift
import DiffractoGeometry

let ring: RingRadiusFit = try RingCalibration.fitRing(pattern: pattern, centerQX: 101.4, centerQY: 98.8,
                                                      innerRadius: 44.0, outerRadius: 84.0)
let result: RingFitCalibration = try RingCalibration.calibrating(ring, d: 2.3545, highTensionVolts: 200_000,
                                                                 pixelPitchMM: 0.05)
if result.refusals.isEmpty, let calibration = result.calibration {
    print(calibration.deltaQ, calibration.lengthPixels ?? .nan, calibration.lengthMillimetres ?? .nan)
}
```

dq = (1/d) / R needs no wavelength; the camera length is L = R / tan(2 theta), with 2 theta = 2 asin(lambda / 2d) and
lambda from ``ElectronWavelength/angstroms(volts:)``. Its pixel value times the pitch is the length in millimetres.

## Refusals

``RingFitCalibration/refusals`` lists, by name, why a ring must not be calibrated from, and ``RingFitCalibration/calibration``
is then nil: a fit that fell back to a conic, too little of the ring covered (``RingCalibration/Refusal``
`insufficientAngularCoverage`), a ring wider than its annulus, an excessive residual or ellipticity, or no ring evidence.
Two refusals withhold only the camera length: no high tension (`noHighTension`) and a ring beyond the Bragg limit. An
impossible d, high tension or pitch throws `invalidInput`.

The fit is synchronous and CPU bound (about 0.1 s at R of 100 px). Run it off the main actor, and pass `progress` or
`isCancelled` (or cancel the task) to stop it.
