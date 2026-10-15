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

287 checks (161 on DiffractoGeometry, in `GeometryChecks.swift`). This is still a **skeleton**: it proves the surface is public and reachable. The negative controls (SpectroKit's
`controls.py` pattern: revert one declaration to `package` in a copy, require the build to fail naming the symbol) and the
per-consumer symbol table come with the first tag; the per-consumer table is `docs/consumers-0.0.3.md`.

```
tools/consumer-check/run.sh
```
