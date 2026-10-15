# ``DiffractoGeometry``

Wavelength, ring fits and the calibration of a diffraction pattern, with the refusals that say when not to trust one.

## Overview

DiffractoGeometry connects a crystal's rings to a detector. ``ElectronWavelength`` gives the relativistic wavelength for a
high tension, ``RingCalibration`` turns one ring of known spacing into the pixel size in reciprocal space (dq, in inverse
angstroms per pixel) and the camera length, ``EllipseCalibration`` fits a ring or a detector ellipse in a pattern,
``RingIndexing`` decides which family a measured |q| belongs to, and ``ZoneAxis`` handles beam directions and the
detector basis.

The kit never defaults an unknown: without a high tension you get dq and no camera length. A result that must not be
used comes back with named refusals, not a number alone. Start with <doc:CalibratingACameraFromOneRing>.

## Topics

### Essentials

- <doc:CalibratingACameraFromOneRing>
- ``RingCalibration``
- ``DiffractionPattern``
- ``ElectronWavelength``

### Calibrating from a ring

- ``RingCalibration``
- ``RingRadiusFit``
- ``RingFitCalibration``
- ``RingEvidence``

### Fitting rings and ellipses

- ``EllipseCalibration``
- ``EllipseCalibrationFit``
- ``EllipseCalibrationModel``
- ``EllipseProfileParameters``
- ``EllipseTransform``

### Naming a ring

- ``RingIndexing``

### Beam directions

- ``ZoneAxis``
