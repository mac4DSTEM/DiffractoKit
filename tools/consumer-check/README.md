# consumer-check

A separate SwiftPM package that depends on DiffractoKit **by path** (`.package(path: "../..")`) and uses, from outside the
package, the public API: `UnitCell`, `MillerIndex`, `IntegerMatrix`, `SymmetryOperation`, `PointSymmetry`,
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

`RingFitChecks.swift` covers the ring fit (`RingCalibration.q` / `radius`, `RingEvidence`, `fitRing` / `calibrating`, the
refusals, `FitReport`, `progress:` / `isCancelled:`, `LocalizedError` on the error types, the tetragonal settings);
`ConventionalChecks.swift` the Niggli reduction, the conventional cell and `SymmetryDetection` with its tolerance and search;
`PolicyChecks.swift` every `Policy` type read field by field, one non-default use each and the refusals of invalid values.
About 530 checks in all.

```
tools/consumer-check/run.sh
CONSUMER_CHECK_JOBS=2 tools/consumer-check/run.sh   # swift build with -j 2
```

Its sibling `tools/readme-consumer` is a second, much smaller consumer: it compiles exactly the main README's program and
checks the output against the README's (`tools/readme-consumer/run.sh`).
