# DiffractoKit

Crystallography and electron diffraction for macOS, as a Swift package: the numbers a scientist reads, never UI. DiffractoKit
is the Mac-native counterpart of the crystallography parts of **pymatgen**, **py4DSTEM** and **pyxem**: the unit cell and its
reciprocal lattice, symmetry found from the atoms, CIF import, d-spacings and the list of reflections with systematic
absences, electron scattering factors and kinematic structure factors, rings and zone axes, and the calibration of a
diffraction pattern (wavelength, camera length, reciprocal-space pixel size, the ring and ellipse fit). It is the library
that the macEM suite's apps, mac4DSTEM and macTEM, import by tagged version, and it is built so that other people can build
their own apps on it.

Mac-native: Swift 6 strict concurrency, `Foundation` and `simd` only (no Python, no dependency), macOS 27. Headless, one
SwiftPM product per module. Where a path is accelerated (Accelerate, none yet), it keeps a CPU reference path with a
registered tolerance. **Status: 0.0.x, the API may still move**: pin an exact tag and bump on purpose.

## The modules

One SwiftPM product per module; a module imports only the ones above it in this table.

| Module | Holds | State |
|---|---|---|
| `DiffractoCore` | `UnitCell` (all seven crystal systems: real and reciprocal lattice, metric, volume, d and g of an hkl, the angle between planes), `SymmetryOperation` and `PointSymmetry` (the eleven Laue classes, `LaueClass`, `CrystalSystem`), symmetry found from the atom sites (`SymmetryDetection`, `CrystalStructure`), `CIFImport` (CIF text to `CrystalStructure`, with named refusals), `ReflectionList` (every hkl up to a k_max with its systematic-absence flag), orbits and families (rings) with their names | slice 1, `docs/registration-slice-1.md`; slice 1b (Niggli reduction, the conventional cell, a Cartesian tolerance) is registered, not written |
| `DiffractoScattering` | `ElectronScatteringFactors` (Lobato and Van Dyck, the table of py4DSTEM, Z = 1 to 103), `StructureFactors` (F(hkl) in A^-2), `CrystalStructure.kinematicReflections(kMax:)`: the reflections after an \|F\| cut, and their rings with the summed \|F\|^2 | slice 2, `docs/registration-slice-2.md` |
| `DiffractoGeometry` | `ElectronWavelength` (relativistic, CODATA 2018), `RingCalibration` (camera length, reciprocal-space pixel size from one ring of known d, with the refusals that say when a ring must not be calibrated from; `fit` runs the whole chain on a pattern), `EllipseCalibration` and `EllipseTransform` (the ring and detector-ellipse fit, py4DSTEM's model), `RingIndexing` (which family a measured \|q\| belongs to, and when it must not be named), `ZoneAxis` (beam directions, the detector basis, excitation errors) | slice 3, `docs/registration-slice-3.md`; 0.0.4 amendments in `docs/registration-0.0.4.md` |
| `DiffractoOrientation` | nothing yet: a placeholder for orientation results, symmetry reduction, fundamental zones, IPF keys and Euler conventions, to be checked against pyxem, diffsims and orix. Do not link it | planned |

## Install

```swift
// Package.swift
dependencies: [
    .package(url: "https://github.com/mac4DSTEM/DiffractoKit", exact: "0.1.0"),
],
targets: [
    .target(name: "MyApp", dependencies: [
        .product(name: "DiffractoCore", package: "DiffractoKit"),
        .product(name: "DiffractoScattering", package: "DiffractoKit"),
        .product(name: "DiffractoGeometry", package: "DiffractoKit"),
    ]),
]
```

The package needs macOS 27 and a Swift 6.2 toolchain (`swift-tools-version: 6.2`). It carries no unsafe flags and no
dependency: your build compiles it optimised itself; a Debug build of the reflection loops and of the ring fit is slow.
`DiffractoGeometry` depends on the other two, but each file that names a type of a module should import that module itself
(`MemberImportVisibility` is on in this package).

## First ten minutes

The program below is `Examples/diffracto-example`; run it with `swift run -c release diffracto-example`. It reads gold from a
CIF, lists its rings, computes the wavelength, and fits and calibrates a synthetic ring.

```swift
import DiffractoCore
import DiffractoScattering
import DiffractoGeometry

// 1. A structure: from CIF text (named errors for a file it cannot trust), or built in code.
let gold = try CIFImport.read(cifText).structure
//   let cell = UnitCell.cubic(a: 4.0782)
//   let gold = CrystalStructure.detectingSymmetry(cell: cell, sites: [AtomSite(z: 79, fractional: [0, 0, 0]), ...])
print(gold.formula, gold.laueClass?.symbol ?? "?", gold.cell.d(MillerIndex(1, 1, 1)))   // Au4 m-3m 2.3545...

// 2. The rings to k = 1/d = 1.2 per A, after the |F| cut (default 1e-4): d, multiplicity, summed |F|^2.
let list = try gold.kinematicReflections(kMax: 1.2)
for ring in list.families.prefix(6) {
    print(ring.family.name, ring.family.d, ring.family.multiplicity, ring.intensity)
}

// 3. The relativistic electron wavelength, in angstroms, for volts.
let lambda = ElectronWavelength.angstroms(volts: 200_000)               // 0.025079... A

// 4. Fit a ring in a pattern and calibrate. `pattern` is a DiffractionPattern(width:height:pixels:) of Float32;
//    the start centre is (row, column), the annulus brackets the ring, d is the ring's spacing in A.
let result = try RingCalibration.fit(pattern: pattern, centerQX: 101.4, centerQY: 98.8,
                                     innerRadius: 44, outerRadius: 84,
                                     d: list.families[0].family.d, highTensionVolts: 200_000, pixelPitchMM: 0.05)
if result.refusals.isEmpty, let calibration = result.calibration {
    print(calibration.deltaQ, calibration.lengthPixels ?? .nan, calibration.lengthMillimetres ?? .nan)
} else {
    print(result.refusals.map(\.name))      // e.g. "insufficientAngularCoverage", "excessiveEllipticity"
}
```

Output of the program on the machine it was written on (macOS 27, Apple Swift 6.4, release build):

```
Au4: a = 4.0782 A, V = 67.827 A^3, Laue class m-3m (cubic), 192 operations
8 rings, 112 reflections:
  {1 1 1}   d = 2.3545 A  multiplicity 8  sum|F|^2 = 1.2701e+00 A^-4
  {2 0 0}   d = 2.0391 A  multiplicity 6  sum|F|^2 = 7.8331e-01 A^-4
  {2 2 0}   d = 1.4419 A  multiplicity 12  sum|F|^2 = 8.6897e-01 A^-4
  {3 1 1}   d = 1.2296 A  multiplicity 24  sum|F|^2 = 1.2543e+00 A^-4
  {2 2 2}   d = 1.1773 A  multiplicity 8  sum|F|^2 = 3.8030e-01 A^-4
  {4 0 0}   d = 1.0195 A  multiplicity 6  sum|F|^2 = 2.0531e-01 A^-4
lambda(200 kV) = 0.025079 A
fit: Asymmetric Gaussian ring, centre (100.400, 99.600) [true (100.4, 99.6)], R = 63.911 px [true 63.911]
calibrated: dq = 6.645293e-03 per A per px, L = 6000.0 px = 300.00 mm [true 6000.0 px]
```

Conventions to know before you use a number: lengths in angstroms; reciprocal vectors **without** 2 pi (g = h b1 + k b2 +
l b3, d = 1/\|g\|, \|g\| in A^-1); the real lattice has a along x and b in the xy plane (the py4DSTEM and VESTA convention);
F(hkl) is in A^-2 (f in A over V in A^3); a pattern's rows are `qx` and its columns `qy`, as in py4DSTEM.

## Reading a result as a scientist

- A refusal is a result, not a failure: `RingCalibration` and `RingIndexing` say by name why a ring or an index must not be
  used (a conic fallback, too little of the ring, a ring wider than its annulus, an ellipse too strong to take (a + b) / 2, no
  high tension, beyond the Bragg limit). `CIFImport` throws errors that name the tag or value at fault.
- The kit never defaults an unknown. A missing high tension gives the reciprocal-space pixel size and no camera length.
- Symmetry is found from the atoms, not read from the file's operations; the file's loop is kept as `fileOperations`, and
  `fileOperationsAreSymmetries` says whether it was right.

## Oracles and pins

Every number has a reference, and the references are pinned (`CLAUDE.md`):

- **pymatgen 2026.9.24** with pymatgen-core 2026.9.23, spglib 2.8.0, NumPy 2.5.4, Python 3.14.7 (cells, expanded CIF sites,
  symmetry operations, d-spacings, reflections within a sphere, the X-ray peak sums), SciPy 1.18.1 (the ring fit's
  least-squares optimum), **py4DSTEM** f050d207's `scattering_factors.txt` (sha256 `be845d1d...`, bit for bit), and the suite
  apps' compiled code (mac4DSTEM 7d6fdbd1, macTEM 9482831) where a port must be bit identical.
- The pins are generated before the Swift code (`tools/slice-*/pins.py`, committed as JSON under `Tests/*/Fixtures/`); the
  Swift tests read the pins and a planted truth (a closed form, a group-theory count, a cell content); a difference is a
  failure of the port, never a tolerance to widen. Agreement with a package is agreement, not truth.
- `docs/parity.md` lists, function by function, what is held against pymatgen, py4DSTEM and pyxem / diffsims / orix, and
  how each row is checked. It also says what is **not** held.
- Nothing here is validated on a real diffraction pattern: that validation is each app's to register against a tag
  (`NOTICE`).

## Build and test

```sh
tools/test.sh                      # the gate: swift test -c release, stamps .build/last-green on a clean tree
tools/last-green.sh                # the last green commit and whether HEAD equals it
tools/consumer-check/run.sh        # a separate package using the public API only
swift run -c release diffracto-example
```

Never run `swift build --build-tests -c release` (it fails on this toolchain); `tools/test.sh` is the gate. The tests read
committed pins and need no Python. The oracle venv is not in the repository: `docs/registration-slice-1.md` section 2 and
`CLAUDE.md` say how it is pinned and provisioned (`References/pymatgen-venv`), for regenerating a pin.

## Licence

GPL-3.0 (`LICENSE`, `NOTICE`). The code is ported from mac4DSTEM and macTEM (both GPL-3.0; `github.com/mac4DSTEM`), whose
lattice and reflection conventions follow py4DSTEM's `Crystal` class (GPL-3.0): the real-lattice matrix with a along x and b in
the xy plane, reciprocal vectors without 2 pi, the per-axis index bound. The electron scattering table of
`DiffractoScattering` is py4DSTEM's (Lobato and Van Dyck 2014), also GPL-3.0. **pymatgen** (MIT) is an oracle: cells,
symmetry and d-spacings are checked against it, and no pymatgen code is copied or ported; a part ported from pymatgen in the
future would be MIT on its own, but one package with any GPL part is GPL, so the whole is GPL-3.0 and an app built on it must
be GPL-compatible. The CIF fixtures under `Tests/` come from mac4DSTEM's tests (GPL-3.0) with their provenance in a comment
line; the beta'' cell content (Andersen et al. 1998) is written by a generator into the gitignored `References/`, as the apps
do. Per-source notices are in `NOTICE`.

## Status and documents

0.0.x: the API may move between tags; the apps pin exact tags. `docs/status.md` says what is live and what is unverified.

| Document | What |
|---|---|
| `docs/status.md` | what is live, the hand-off, what is not verified |
| `docs/parity.md` | the checklist against pymatgen, py4DSTEM, pyxem / diffsims / orix |
| `docs/registration-slice-1.md`, `-1b.md`, `-2.md`, `-3.md`, `registration-0.0.4.md` | what each slice makes, the definitions, the pins, the planted truth, the gates, the amendments |
| `docs/consumers-0.0.3.md` | what an app changes to import 0.0.3, and the numbers that move |
| `docs/slice-1-mutants-2026-10-14.md` | the mutation run behind slice 1 |
| `tools/consumer-check/README.md` | the separate package that uses the public API only |
| `CLAUDE.md` | the rules for working in this repository |
