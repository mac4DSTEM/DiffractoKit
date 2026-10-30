# consumer-check

A separate SwiftPM package that depends on DiffractoKit **by path** (`.package(path: "../..")`) and uses, from outside the
package, the public API of slices 1 to 3: `UnitCell`, `MillerIndex`, `IntegerMatrix`, `SymmetryOperation`, `PointSymmetry`,
`LaueClass`, `CrystalSystem`, `Element`, `AtomSite`, `CrystalStructure`, `SymmetryDetection`, `CIFImport`, `CIFDocument`,
`ReflectionList`, `Reflection`, `ReflectionFamily`, `ElectronScatteringFactors`, `StructureFactors`,
`KinematicReflectionList`, and DiffractoGeometry's `ElectronWavelength`, `RingCalibration` (closed forms, `calibrate`,
`FitReport`, `refusals`, `fit`), `RingIndexing` (the list, family and `ReflectionFamily` overloads, every refusal),
`ZoneAxis` (`lowIndex`, `cartesian`, `detectorBasis`, `isUnitDirection`, `project`), `DiffractionPattern`,
`EllipseTransform`, `EllipseCalibration` (`fit1D`, `fitBestAvailable`, `fitAmorphousRing`, the errors), the typed errors.
`package` access does not cross a package boundary, so the target compiles only while each symbol it names is `public`;
the values it checks are ones a reader verifies by hand (β″'s d(200) = a sin β / 2, gold's eight families, λ at 200 kV =
0.025079340450548003 Å, Au {111} at R 150 px and 200 kV: L_px 14082.307097295722 and Δq 0.0028313372405787253 Å⁻¹/px,
a q at Au {111}'s g with σ = 1e-3 q accepted as {1 1 1}, β″'s [0 0 1] = (cos β, 0, sin β), a clean Gaussian ring
rendered by `GeometryChecks.swift` and recovered inside the apps' bands (radius 1 %, b/a 0.01, centre 1 px), and a 2:1
ellipse's offset (4, 3) corrected to (2, 3)).

DiffractoKit 0.1.0 adds `RingFitChecks.swift` and two sections of `main.swift`: `RingCalibration.q` / `radius`, `RingEvidence`
(by hand and on a planted ring), `fitRing` / `calibrating` (equal to `fit`), `requireRingEvidence`, the three new refusals and
`FitReport`'s new fields, `progress:` / `isCancelled:` on every fit, `LocalizedError` on the three error types, and the
tetragonal settings. 397 checks at 0.0.4.

Slice 1b adds `ConventionalChecks.swift`: `UnitCell.niggliReduced` / `ReducedCell` (fcc gold's primitive cell, a skewed basis of
it, an invalid cell), `UnitCell.conventionalCell` / `ConventionalCell` / `BravaisLattice` (gold primitive -> cF with 4 centring
translations, bcc -> cI, a C-centred monoclinic primitive cell -> mC with b the unique axis), `CrystalStructure.setting` /
`conventionalSetting` / `ConventionalSetting` and its four conversions (gold: 192 operations in the conventional cell and 48 in
the primitive file cell; rutile with its 4-fold axis along a: 4/mmm, 16 operations, the file's (1 0 0) the conventional
(0 0 1)), and `detectingSymmetry(conventional:tolerance:search:)` with `SymmetryDetection.Tolerance` / `Search` (a 0.012 Å
displacement in a 30 Å cell: symmetric at the fractional 5e-3, not at the Cartesian 1e-3 Å; the cell list equals the brute-force
scan). 530 checks in all.

DiffractoKit 0.2.0 adds `PolicyChecks.swift`: every `Policy` type (`ReflectionList`, `CrystalStructure.ValidationPolicy`,
`SymmetryDetection`, `CIFImport`, `EllipseCalibration`, `RingEvidence`, `RingCalibration`) with its `.default` read field by
field, one non-default use whose outcome a reader verifies by hand (a 0.3 Å contact usable under a 0.2 Å limit, a 120 degree
gap fitted under a coverage bound of 20 bins, occupancies 0.50 and 0.51 two species at tolerance 0), the refusals of an invalid
value (a position tolerance outside its domain, an `Int.max` bin count, a NaN residual cap), `SymmetryOperation(xyz:maximumMagnitude:)`,
`RingIndexing`'s `minimumMargin:` and `reach`, and `retainingAll:` (`kept` and `families` equal, `all` empty).

```
tools/consumer-check/run.sh
CONSUMER_CHECK_JOBS=2 tools/consumer-check/run.sh   # swift build with -j 2
```

Its sibling `tools/readme-consumer` is a second, much smaller consumer: it compiles exactly the main README's program and
checks the output against the README's (`tools/readme-consumer/run.sh`).
