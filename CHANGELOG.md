# Changelog

All notable changes to DiffractoKit. Versions are semantic (0.x: a minor release may change the API). The numbers each module
reproduces, and against what, are in `docs/parity.md`.

## Unreleased (proposed 0.6.0)

Orientation symmetry beyond the three app classes, pinned to orix; the py4DSTEM direction families; the preset crystals, the
Materials Project models and the phase model that the two apps carry.

### Added
- `IPFColorKey.orix`: orix 0.15.0's `DirectionColorKeyTSL` (the EDAX/TSL key as in MTEX) for all eleven Laue classes, in the kit's
  standard settings. `IPFColorKey.native` now has a key for every class: the app's own m-3m, 6/mmm and b-unique 2/m keys
  are unchanged (pinned to their Float bits), and the other eight classes use the orix key. Measured against orix on 2 777
  poles: worst RGB gap 2.98e-8 (Float rounding). The app's native keys versus orix: m-3m 3.5e-6, 6/mmm 0.30, 2/m 1.0 (different
  conventions for the last two). 2/m and -3m are mapped from orix's settings by a fixed frame (`docs/orientation-symmetry.md`).
- `CrystalOrientationSymmetry.disorientation(from:to:)`: the smallest angle between two orientations of one crystal (Double
  radians) with the reduced matrix. `reduce(_:)` and the pair angle were cross-checked against orix 0.15.0 for all eleven
  classes on seeded random orientations and pairs: Float angle within 1.19e-7 rad, reduced matrix within 6.7e-16, pair angle within
  1.4e-14 rad. **No change to the reduction was needed.**
- `CrystalOrientationSymmetry.py4DSTEMDirectionFamily(_:)`: py4DSTEM f050d207's `symmetry_reduce_directions` (`Orientation.family`),
  all eleven classes, bit identical to the real function on 7 161 pinned columns (171 columns on a sector wall are not compared). It
  returns families, not a rotation; `reduce` is still the minimum-angle rotation.
- `CrystalPresets` (DiffractoCore): fcc, bcc, simple cubic, diamond, hcp and the cubic families; Au, Al, Ni, Cu, Fe, Si, Mg at the
  apps' lattice constants; beta'' (Mg5Si6, C2/m) and 2H-WS2 from their literature refinements. Pinned to independent pymatgen
  `Structure.from_spacegroup` expansions and spglib (volume and d 4e-16 relative, sites 3e-11, Laue class and operation count exact).
- `MaterialsProjectSummaryResponse` and the other wire types, `MaterialsProjectImport` (`request`, `decode`, `firstDocument`,
  `structure(from:)`, `check`, `parseFormula`, `Failure`), `PhaseExpectation`, `PhaseExpectationVerdict` (DiffractoCore): the
  apps' Materials Project JSON models and verdicts on the kit's `CellStandardisation`; no network access.
- `PhaseModel`, `PhaseModelLibrary`, `PhaseModelSource`, `PhaseModelIssue` (DiffractoOrientation): the apps' `CrystalModel` checks
  (cell, basis, occupancy, close contacts, `unsupported_element`, a declared symmetry the cell does not have), provenance and the
  FNV-1a content fingerprint (recomputed independently), built-in models by id, and `PhaseModel(materialsProject:fetchedAt:)`.
  `PhaseModel.orientationClass(of:)` generalises the apps' cubic / hexagonal / monoclinic test to the eleven classes: a class is
  named only when the structure's own Laue operators in its cell are the class's standard operators (agrees with spglib's Laue class
  on 15 space groups). The -3m1 groups (P-3m1, R-3m in hexagonal axes: twofolds on a) return nil.
- `tools/orix/` (provision script, `ipf-pins.py`, `disorientation-pins.py`), `tools/symmetry/family-pins.py`, `tools/presets/pins.py`;
  consumer-check 693 checks (was 666).

### Changed
- `OrientationSymmetryError.nativeKeyUnavailable` is no longer thrown (kept so `switch` statements compile).

### Migration
- `IPFColorKey` has a new case, `.orix`: an exhaustive `switch` over it needs one more case.
- `ipfColor(direction:key: .native)` for a class other than m-3m, 6/mmm and 2/m no longer throws `nativeKeyUnavailable`; it returns
  the orix colour. A consumer that treated the refusal as "no key" now gets a key.
- Nothing else changes: no existing number moves. mac4DSTEM and macTEM may drop their own `Crystal` presets, `ReferenceStructures`,
  `CrystalModel` checks and `MaterialsProjectImport` wire types for the kit's (their `ACOMCrystalSymmetry` family becomes
  `PhaseModel.orientationClass`; the apps' `.identity` is `nil`, and the kit also reduces tetragonal, orthorhombic, cubic m-3, trigonal
  and 6/m cells that the apps called unreduced).

## 0.5.0 (2026-10-20)

ACOM matching over a whole scan, so the apps can drop their own scan loop (mac4DSTEM `OrientationMatching.matchAll`, CPU path).

- `ACOMMatcher.matchAll(peaks:width:height:originX:originY:invAngstromPerPixel:selection:maxWorkers:isCancelled:progress:)`:
  positions dealt round-robin over worker threads (default: the active processor count), `isCancelled` polled before every
  position (a cancelled run returns nil), `progress` in 1 % batches, non-decreasing, the last call 1. Each position is exactly `ACOMMatcher.match`,
  so no pinned number moves.
- `ACOMScanSelection` (`full`, `preview(maxDimension:)`, `square(centerX:centerY:radius:)`) and `ACOMScanResults` (the grid,
  `expandedPreview(selection:)`), taken from the app's `ACOMScanSelection` and `OrientationMap` container; the map's presentation
  (Euler read-out text, reliability and IPF-Z images) stays in the app.
- `tools/acom-benchmark`: the speed measurement. On the default 600-direction library (32 x 128 bins, 19 peaks a pattern), under a
  machine load of about 9: 1.45 ms a pattern on one worker, 0.24 ms on eight (4 200 patterns a second; 2026-10-20, load 3); a 256 x 256 scan is about
  2 minutes on one worker and 20 seconds on eight. A cancel request is therefore honoured within one pattern (about 2 ms).

Not here: the apps' Metal backend (measured slower than this CPU path in the app, kept in the app as a named option) and strain
mapping (a Bragg-vector lattice fit with no crystal in it; it belongs to CuboKit).

### Migration
Additive: a consumer of 0.4.1 recompiles unchanged. mac4DSTEM replaces its CPU `matchAll` loop and `ACOMScanSelection` with the
kit's; the Metal path and the `OrientationMap` presentation stay.

## 0.4.1 (2026-10-16)

Dependency range only; no API, source or number changes.

- ImagoKit is now `.upToNextMinor(from: "0.2.0")` (was `exact: "0.2.0"`), so an app that pins exact kit tags resolves one ImagoKit.

### Migration
None for consumers of 0.4.0 beyond re-pinning to 0.4.1. A consumer may now use any ImagoKit 0.2.x.

## 0.4.0 (2026-10-16)

Orientation symmetry for all eleven Laue classes, the ACOM matcher, the Materials Project cell import, and helpers for
application front ends (crystal loading, ring overlay, indexing display, ring scale).

### Migration
Additive: no existing signature or number changes, so a consumer of 0.3.1 recompiles unchanged. The one change to the
package graph is a new dependency, ImagoKit 0.2.0 (exact), for `FFT1D` (ImagoFFT), used by `DiffractoOrientation` only;
a consumer that also pins ImagoKit must pin 0.2.0 (0.2.0 is additive over 0.1.1). The 0.3.1 fix (symmetry detection of
non-lattice rotations) is included in this tag.

### Added
- `RingScale` (DiffractoGeometry): |q| of a ring radius and the radius of a |q| under a `RingCalibration` (the exact closed forms with a beam energy, the linear Δq·R without one), the linear helper `RingScale.linear(radius:deltaQ:)` and `fitVersusFile`, ported from an application; its forms agree with `RingCalibration.q`/`radius` (P24). The application's `fileScale` is not ported (the kit has no `CalibrationUnitConversion`).
- `IndexingDisplay.index` and `DisplayedIndexing` (DiffractoGeometry): the numbers a room shows for `RingIndexing`'s decision (sigma_q from a percent tolerance, delta %, |dd|/d, z, margin, odds and the "> 10⁶ : 1" cap text), ported from an application; the decision is untouched (P23).
- `RingOverlay` and `RingEllipse` (DiffractoGeometry): a ring fit as an image-coordinate ellipse (x = column, y = row), the drawing polyline, the predicted ring of another family, and the point-to-ellipse distance, ported from an application (P22).
- `CrystalLoading.load(bytes:fileName:)` and `LoadedCrystal`: a CIF from bytes with a 5 MiB cap, strict UTF-8, `CIFImport` and a lower-case hex SHA-256 of the bytes, ported from an application; refusals `CrystalLoading.Refusal.tooLarge`/`.notUTF8`, importer errors unchanged (P21).
- `CellStandardisation.standardise`: the primitive-to-conventional conversion of a centred cell (fixed ITA matrices for F, I, C, A, B, R) with its guards, its F conversion pinned to pymatgen `SpacegroupAnalyzer` standard cells (P20; I by a closed form). The reference route is `CrystalStructure.detectingSymmetry(...).conventionalSetting`.
- `CrystalOrientationSymmetry`: proper operation tables and minimum-angle rotation reduction for all eleven Laue classes,
  with deterministic quaternion ties, direction folding and fundamental-zone sampling. Standard settings are explicit;
  there is no unknown-class or identity fallback.
- `IPFColorKey.py4DSTEM`: the pinned automatic-sector, full-mask beam key of py4DSTEM f050d207, bit identical after
  converting Python's Double RGB to the public Float RGB type. Pass Python family poles (original poles for 2/m).
- `IPFColorKey.native`: the previously existing m-3m, 6/mmm and b-unique 2/m keys, pinned to its original Float bits. Other
  native keys throw a named refusal. Invalid rotations, poles and sampling counts are refused explicitly.
- Pinned pymatgen point-group operations, SciPy minimum-angle fixtures, py4DSTEM colour fixtures and original-source
  native colours/samples; generator, seed and measured tolerances in `docs/parity.md` (P17–P19).
- ACOM matcher (`DiffractoOrientation`): `ACOMTemplateLibrary` (template library of a simplified polar-correlation
  route: fundamental-zone beam directions of any Laue class, polar images with radial kernel, azimuthal blur, per-ring mean
  subtraction and unit norm, per-ring azimuthal FFT) and `ACOMMatcher` (best template, in-plane angle, score, runner-up
  reliability and the mirror pass, result as `OrientationResult`), with `ACOMTemplateSettings`, `ACOMReflection`, `ACOMPeak`,
  `ACOMError`. `ACOMExcitationModel` chooses the excitation error: `.flat`, `.appCurvature` (the sign used by the applications it was ported from) or `.py4DSTEM`
  (py4DSTEM's `excitation_errors`). This is not py4DSTEM's SO(3) correlogram; docs/parity.md says what is compared. Pinned (P25)
  to a NumPy transcription of the original application source (m-3m, 6/mmm, 2/m) and, for the spots and planted orientations, to py4DSTEM
  f050d207 run (spots 5e-15 relative). Measured: with `.py4DSTEM` the route is within 4 degrees of the truth on 17 of 18 planted
  patterns; one 2/m pattern returns a different beam direction (127 degrees).
- Dependency: ImagoKit 0.2.0 (exact), for `FFT1D` (ImagoFFT), used by DiffractoOrientation only.

### Conventions
- py4DSTEM's direction-family reduction is not ported (returns families, not a rotation). Neither the minimum-angle
  rotation nor the native direction fold substitutes for Python's family array. See `docs/orientation-symmetry.md`.
- Additive API; no orientation-map container added.

## 0.3.1 (2026-10-15)

### Fixed
- Symmetry detection: a CIF whose operation loop lists rotations that are not symmetries of the cell's lattice (the 48 signed
  permutations of x, y, z written for a 60° rhombohedral primitive cell of MgO) gave `-1` or `-3m` instead of `m-3m`. Such a
  rotation was tested only against the atoms, passed as an operation, and the closure of the set then dropped the others. A
  file's rotation now joins the candidates only if it preserves the cell's metric (to a relative 1e-3); the lattice bounds the
  class and the atoms decide it. Reflections are unchanged on every pinned structure; `SymmetryDetection.Result.candidates` no longer lists such a rotation.
- Not a change: hcp written in a cell with γ = 60° is read as 6/mmm (it is hcp); 0.1.0's 2/m was the metric not being recognised.

## 0.3.0

One release for two steps (there is no separate 0.2.0 tag): the generality pass, where every policy is a parameter with the
previous value as its documented default, the platform floor is lowered and the docs stand alone (the default path is bit
identical to 0.1.1 on 31,725 outputs); and vector-matched phase mapping with the value types of an orientation result in
`DiffractoOrientation`.

### Added: phase mapping and orientation results (`DiffractoOrientation`)
- `PhaseReferenceLibrary` (`build(phases:settings:)`, `Failure`, `lowIndexZoneAxes` (49 directions), `inPlaneSteps`,
  `projectedVectors`, `rotate`, `reflections(of:kMaxInvAngstrom:)`, `unsupportedElements(of:)`), `PhaseDefinition`,
  `ReferenceVector`, `PhaseOrientationReference` (`chanceMatchFraction`), `LatticeVector`, `OrientationRelationship`,
  `PhaseReferenceSettings` (7 settings).
- `PhaseVectorMatcher` (`map`, `classify`, `classifyKnownVariants`, `fitMatrixOrientation`, `fitZoneAxis`, `claims`,
  `experimentalVectors`, `nearest`, `score`, `containsFriedelPair`, `ZoneAxisFit`), `PhaseVectorSettings` (17 settings,
  `ClassificationRule`), `PhaseVectorResolution`, `PhaseVerdict`, `PhaseVectorResult`, `PhaseMap`, `PhaseDiskClaim`,
  `PeakCoordinates` (detector pixels with the origin subtracted in `Float`, or calibrated Å⁻¹).
- `EulerAngles` (Bunge z-x-z as SciPy, from and to py4DSTEM's orientation matrix), `ZoneAxisFrame` (`detectorBasis(zoneAxis:)`,
  `matrix(detectorBasis:inPlaneAngle:)`), `OrientationResult`.
- Every setting is a parameter with a documented default (the values tuned on measured SPED scans and on Thronsen et al.'s dataset); scans
  are `[[SIMD2<Double>]]`; progress and cancellation are closures; `map(maximumWorkers:)`.
- Regression tests: py4DSTEM f050d207 (`generate_diffraction_pattern` at λ → 0: hkl sets, |q| within 1.1e-15, the in-plane
  frame rotation per zone axis, intensities within py4DSTEM's complex64 precision), SciPy 1.18.1 Euler angles on 1000 rotations
  and the gimbal cases, NumPy closed forms, and the original Swift implementation bit for bit (about 1.68 million values).

### Added: the generality pass
- Policy types, each `Sendable, Equatable` with a defaulted memberwise `init` and `static let default`:
  `ReflectionList.Policy` (k_max cap 4 Å⁻¹, reflection and enumeration caps 1e7 and 1e8, the ring rule 1e-9, absence 1e-9),
  `CrystalStructure.ValidationPolicy` (close contact 0.5 Å, split-site occupancy 0.02), `SymmetryDetection.Policy`,
  `SymmetryOperation(xyz:maximumMagnitude:)` (1000), `CIFImport.Policy`, `EllipseCalibration.Policy` (18 fields plus
  `evidence`), `RingEvidence.Policy` (8), `RingCalibration.Policy`, `RingIndexing` `minimumMargin:` (ln 19). The |F| cut
  (`cut:`, 1e-4 Å⁻²) is documented as a policy. Where each is passed: the doc comment of each type.
- `kinematicReflections(..., retainingAll: Bool = true)`: with `false` the list from before the |F| cut is never held (`all` is
  empty; `kept` and `families` identical). A 30 Å cell at k_max 2.5 (1,767,026 reflections, 441,710 kept): 820 MB peak and
  2.3 s by default, 279 MB and 2.5 s with `false`.
- `CrystalStructure.symmetryRefusal`, `SymmetryDetection.Result.refusal`, `SymmetryDetection.Tolerance.refusal(for:)` and
  `isValid(for:)`, `ConventionalCell.reductionConverged`.
- A standalone README with the first ten minutes, `Examples/diffracto-example`, a DocC catalogue per module, CI.

### Changed
- The platform floor is macOS 15 and iOS 18 (was macOS 27); numbers unchanged there.
- Application names were removed from sources, test names and comments; test names changed accordingly.
- `ReflectionListError.cellTooLargeOrInvalid`'s text names the estimated count, the default limit and the policy field.

### Fixed
- Traps: `EllipseCalibration.Policy(angularBinCount:)` with a huge count; a symmetry position tolerance that is negative,
  denormal or tiny (now refused: the window must be 1e-9 up to 0.5 of a cell axis); a CIF with an absurd cell (overflow in the
  conventional-cell route: refused by validation or flagged `reductionConverged == false`); `IntegerMatrix.determinant` wraps.
- A position tolerance of 0 or NaN returned a structure without the identity operation; it is refused, and the identity is
  always present.
- A cell sheared by hundreds of units was silently left unreduced (found in review); a non-converged reduction is flagged.
- NaN, negative or infinite metric tolerances and Niggli epsilons are read as documented (0, or unreduced and unconverged).
- `RingIndexing.reach` answers NaN where `index` refuses the margin.
- Against the original phase-mapping implementation:
  - A rotation about the beam whose R22 rounds to 1 − 2⁻⁵³ (with R20 = R21 = 0) was lost by the Euler decomposition (both
    angles 0); it is folded into φ1.
  - An exact score tie between two phases was broken by the process's hash seed; it goes to the lower phase index.
  - An in-plane step of 1e-300 degrees trapped in `Int(_:)` and 1e-9 degrees asked for 3.6e11 rotations: `build` throws
    `libraryTooLarge` before allocating, `fitZoneAxis` throws it past `PhaseReferenceLibrary.maximumInPlaneStepCount` (360 000),
    and `inPlaneSteps` holds a finer step at that count. Steps within the cap give the same rotations as before.
  - `PhaseMap(width:height:...)` trapped on a width × height that overflows (it now holds no results); `map` returns nil for a
    hand-built library whose entries name a missing phase (`PhaseReferenceLibrary.isConsistent`); `classify` and
    `classifyKnownVariants` give no verdict for candidate indices outside the library.
  - A cancelled `PhaseReferenceLibrary.build` returned a library silently lacking phases; it now throws `CancellationError`.

### Migration
**Policies.** Nothing is required: every new parameter defaults to the previous value. Notes for a caller who sets policies:
`EllipseCalibration.Policy(angularBinCount:)` derives the degeneracy and sparse bounds, while setting `angularBinCount` on
`.default` keeps 30 and 12; through `RingCalibration.fitRing` the calibration policy's evidence replaces the ellipse policy's,
and loosening coverage needs both `RingCalibration.Policy.minimumOccupiedAngularBins` and the ellipse policy's
`degeneracyBoundBins`; removing the `ReflectionList.Policy` caps removes the memory and time protection (a k_max of 100 reached
10 GB). Error texts quote the default bounds. A kinematic pattern cannot tell +c from −c on a polar structure: no consumer may
claim a polarity from a kinematic match. Earlier per-release migration notes are in `docs/consumers-0.0.3.md` to `docs/consumers-0.3.0.md` and in git history.

**Phase mapping.** An application that carried its own vector matching can delete it when it adopts the kit's. For code moving from such an implementation: `Crystal` becomes `CrystalStructure` (`PhaseDefinition.crystal`); a structure the
kit refuses throws `PhaseReferenceLibrary.Failure.structureRefused`; `BraggVectors` becomes `scan:scanWidth:scanHeight:` with
`PeakCoordinates` (pass each peak's `Float` x, y as `SIMD2<Double>`; the pixel form reproduces the old arithmetic bit for bit);
the cancellation token becomes `progress:` / `isCancelled:` closures; `ACOMOrientation` is `ZoneAxisFrame`; `fitZoneAxis`
throws; a scan whose position count is not width × height gives nil (or `[]`) from every entry point, and an out-of-range
`matrixEntryIndex` makes `map` return nil; a cancelled `build` throws `CancellationError`. The decision rule of `.knownVariants` and its
two constants follow Thronsen et al., Ultramicroscopy 255 (2024) 113861 (CC BY 4.0; NOTICE). No other number moves.

## 0.1.1

Niggli reduction, the conventional cell in both settings, the Cartesian symmetry tolerance, detection on a cell list;
the CI type-check fixes.

## 0.1.0

Cells of the seven systems, symmetry from the atoms, CIF import, reflections, absences, orbits and rings;
electron scattering factors and structure factors; the wavelength, ring and ellipse calibration, camera length, indexing and
zone axes; error texts; progress and cancellation of the fits.
