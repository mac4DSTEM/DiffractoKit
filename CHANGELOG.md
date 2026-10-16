# Changelog

All notable changes to DiffractoKit. Versions are semantic (0.x: a minor may change the API). The numbers each module
reproduces, and against what, are in `docs/parity.md`.

## Unreleased

### Added
- `CrystalOrientationSymmetry`: proper operation tables and minimum-angle rotation reduction for all eleven Laue classes,
  with deterministic quaternion ties, direction folding and fundamental-zone sampling. Standard settings are explicit;
  there is no unknown-class or identity fallback.
- `IPFColorKey.py4DSTEM`: the pinned automatic-sector, full-mask beam key of py4DSTEM f050d207, bit identical after
  converting Python's Double RGB to the public Float RGB type. Pass Python family poles (original poles for 2/m).
- `IPFColorKey.native`: the app's existing m-3m, 6/mmm and b-unique 2/m keys, pinned to its original Float bits. Other
  native keys throw a named refusal. Invalid rotations, poles and sampling counts are refused explicitly.
- Pinned pymatgen point-group operations, SciPy minimum-angle fixtures, py4DSTEM colour fixtures and original-source
  native colours/samples; generator, seed and measured tolerances in `docs/parity.md` (P17–P19).
- ACOM matcher (`DiffractoOrientation`): `ACOMTemplateLibrary` (template library of the apps' simplified polar-correlation
  route: fundamental-zone beam directions of any Laue class, polar images with radial kernel, azimuthal blur, per-ring mean
  subtraction and unit norm, per-ring azimuthal FFT) and `ACOMMatcher` (best template, in-plane angle, score, runner-up
  reliability and the mirror pass, result as `OrientationResult`), with `ACOMTemplateSettings`, `ACOMReflection`, `ACOMPeak`,
  `ACOMError`. `ACOMExcitationModel` chooses the excitation error: `.flat`, `.appCurvature` (the apps' sign) or `.py4DSTEM`
  (py4DSTEM's `excitation_errors`). This is not py4DSTEM's SO(3) correlogram; docs/parity.md says what is compared. Pinned (P21)
  to a NumPy transcription of the app source (m-3m, 6/mmm, 2/m) and, for the spots and planted orientations, to py4DSTEM
  f050d207 run (spots 5e-15 relative). Measured: with `.py4DSTEM` the route is within 4 degrees of the truth on 17 of 18 planted
  patterns; one 2/m pattern returns a different beam direction (127 degrees).
- Dependency: ImagoKit 0.1.1 (exact), for `FFT1D` (ImagoFFT), used by DiffractoOrientation only.

### Conventions
- py4DSTEM's direction-family reduction is not ported (returns families, not a rotation). Neither the minimum-angle
  rotation nor the native direction fold substitutes for Python's family array. See `docs/orientation-symmetry.md`.
- Additive API; no package dependency, matcher, UI or orientation-map container added.

## 0.3.1 (unreleased: the owner tags it)

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
- No app or suite name remains in sources, test names or comments; test names changed accordingly.
- `ReflectionListError.cellTooLargeOrInvalid`'s text names the estimated count, the default limit and the policy field.

### Fixed
- Traps: `EllipseCalibration.Policy(angularBinCount:)` with a huge count; a symmetry position tolerance that is negative,
  denormal or tiny (now refused: the window must be 1e-9 up to 0.5 of a cell axis); a CIF with an absurd cell (overflow in the
  conventional-cell route: refused by validation or flagged `reductionConverged == false`); `IntegerMatrix.determinant` wraps.
- A position tolerance of 0 or NaN returned a structure without the identity operation; it is refused, and the identity is
  always present.
- A cell sheared by hundreds of units was silently left unreduced (slice 1b's review); a non-converged reduction is flagged.
- NaN, negative or infinite metric tolerances and Niggli epsilons are read as documented (0, or unreduced and unconverged).
- `RingIndexing.reach` answers NaN where `index` refuses the margin.
- Against the phase-mapping implementation it was ported from:
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
claim a polarity from a kinematic match. The app-specific notes of earlier releases are in git history
(`docs/consumers-0.0.3.md` to `docs/consumers-0.3.0.md`; they leave the repository at the next lean pass and stay in git history).

**Phase mapping.** An app that carried its own vector matching deletes it when it adopts the kit's (mac4DSTEM:
`PhaseVectorMatching.swift`, `PhaseReferenceLibrary.swift`, `EulerAngles` and `ACOMOrientation`, and their tests). For code moving from the app implementation: `Crystal` becomes `CrystalStructure` (`PhaseDefinition.crystal`); a structure the
kit refuses throws `PhaseReferenceLibrary.Failure.structureRefused`; `BraggVectors` becomes `scan:scanWidth:scanHeight:` with
`PeakCoordinates` (pass each peak's `Float` x, y as `SIMD2<Double>`; the pixel form reproduces the old arithmetic bit for bit);
the cancellation token becomes `progress:` / `isCancelled:` closures; `ACOMOrientation` is `ZoneAxisFrame`; `fitZoneAxis`
throws; a scan whose position count is not width × height gives nil (or `[]`) from every entry point, and an out-of-range
`matrixEntryIndex` makes `map` return nil; a cancelled `build` throws `CancellationError`. The decision rule of `.knownVariants` and its
two constants follow Thronsen et al., Ultramicroscopy 255 (2024) 113861 (CC BY 4.0; NOTICE). No other number moves. On
re-pinning to 0.3.0, mac4DSTEM deletes its own `PhaseVectorMatching` and `PhaseReferenceLibrary` copies.

## 0.1.1

Slice 1b: Niggli reduction, the conventional cell in both settings, the Cartesian symmetry tolerance, detection on a cell list;
the CI type-check fixes.

## 0.1.0

Slices 1 to 3: cells of the seven systems, symmetry from the atoms, CIF import, reflections, absences, orbits and rings;
electron scattering factors and structure factors; the wavelength, ring and ellipse calibration, camera length, indexing and
zone axes; error texts; progress and cancellation of the fits.
