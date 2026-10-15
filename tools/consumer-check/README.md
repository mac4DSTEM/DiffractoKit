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

DiffractoKit 0.1.0 adds `Geometry004Checks.swift` and two sections of `main.swift`: `RingCalibration.q` / `radius`, `RingEvidence`
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
scan). 530 checks in all, and 190 negative controls (144 at 0.0.4, 46 for slice 1b; 13 enum-case, conformance and parameter declarations are listed under `SKIPPED`).

Negative controls (SpectroKit's `controls.py` pattern): `controls.py` copies the repository under `.build/consumer-check/mutants/`,
checks that the unedited copy builds, then reverts ONE declaration from `public` to `package` at a time and requires the
consumer build to fail with the compiler's protection-level diagnostic naming that symbol (one warm build directory, an
incremental build per control). The controls are one per new 0.1.0 public symbol plus one per public type the apps use
(`docs/consumers-0.0.3.md`) and members of it the check calls. An enum case and a protocol conformance cannot be reverted to
`package`: `controls.py` lists them under `SKIPPED` with the reason. `controls.py <root> --dry` checks that every edit target
occurs exactly once without building.

```
tools/consumer-check/run.sh              # the plain check
tools/consumer-check/run.sh --controls   # the plain check, then the negative controls
CONSUMER_CHECK_JOBS=2 tools/consumer-check/run.sh --controls   # every swift build with -j 2
```
