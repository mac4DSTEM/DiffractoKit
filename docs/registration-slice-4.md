# Slice 4: vector matching and phase mapping (`DiffractoOrientation`, 0.3.0) — registration

Registered 2026-10-15, 15:45 CEST, before any Swift of slice 4, by DiffractoKit's fourth session (the owner's overrule of
15:20: vector matching and phase mapping come into the kit now, as 0.3.0; suite `ROADMAP.md`, phase 2 step 4). Tier 1: every
label, distance and count it returns is a number a scientist reads. A change after code starts is a dated amendment (§9).

## 1. What slice 4 makes

A port of mac4DSTEM 7d6fdbd1's vector-matched phase mapping into `DiffractoOrientation`, public API only, read with
`git -C ../mac4DSTEM show 7d6fdbd1:<path>`:

| Source (7d6fdbd1) | Ported | Stays in the app |
|---|---|---|
| `mac4DSTEM/Core/Crystal/PhaseReferenceLibrary.swift` (572 lines) | all of it: `ReferenceVector`, `PhaseOrientationReference` (with `chanceMatchFraction`), `LatticeVector`, `OrientationRelationship`, `PhaseDefinition`, `PhaseReferenceSettings`, `PhaseReferenceLibrary` (`build`, `Failure`, `lowIndexZoneAxes`, `inPlaneSteps`, `cartesianZoneAxis`, `projectedVectors`, `rotate`) | nothing |
| `mac4DSTEM/Core/Crystal/PhaseVectorMatching.swift` (1698 lines) | all of it: `PhaseVectorSettings` (with `ClassificationRule`), `PhaseVectorResolution`, `PhaseVerdict`, `PhaseVectorResult`, `PhaseMap`, `PhaseVectorMatcher` (`experimentalVectors`, `nearest`, `score`, `containsFriedelPair`, `ZoneAxisFit`, `fitZoneAxis`, `fitMatrixOrientation`, `crossPhaseWinsOver`, `classify`, `classifyKnownVariants`, `specificReflectionCount`, the matrix challenge, `projectedAzimuth`, `orientationConsistent`, `map`), `PhaseDiskClaim` and `claims` | nothing |
| `mac4DSTEM/Core/Analysis/OrientationResult.swift` (894 lines) | the value types of an orientation result: `EulerAngles` (the py4DSTEM / SciPy Bunge convention, both directions), `ACOMOrientation.detectorBasis(zoneAxis:)` and `matrix(detectorBasis:inPlaneAngle:)` (the frame), `OrientationResult` (one position's result) | the map and image presentation (`FloatImage`, colour keys, read-out text), the backend enum, the cubic symmetry and IPF key (they go with template matching, slice 5) |

`displayName` and `displayName`-like UI strings stay where they are in the source when they are part of the result's meaning
(`PhaseVerdict.displayName`), and leave when they are app wording.

Not in slice 4: `OrientationPlan`, `OrientationMatcher` (polar-correlation template matching, Accelerate and Metal) and the
symmetry reduction: they are slice 5, against pyxem / orix and py4DSTEM's `crystal_ACOM`, registered on their own.

## 2. Substitutions (the only places the port may differ from the source)

1. **`Crystal` → `CrystalStructure`.** The library reads, per phase, the kit's kinematic reflections at `kMaxInvAngstrom`
   (`CrystalStructure.kinematicReflections(kMax:cut:)`'s kept list at the default |F| cut, which slice 2 proved bit identical to
   the app's `Crystal.reflections(kMax:)` on its presets): `g` Cartesian in Å⁻¹ (no 2π), `intensity` = |F|², `hkl`. The order of
   that list must equal the app's order (the density cap's sort and the final length sort break ties on input order; a
   different order is a different library). `latReal` / `latInv` rows are the kit's real and reciprocal lattice vectors in the
   same Cartesian frame (a along +x, slice 1). `Crystal.unsupportedElements` is the kit's scattering-table coverage.
2. **`BraggVectors` → `[[SIMD2<Double>]]`.** `map`, `fitMatrixOrientation` and `fitZoneAxis` take the scan as row-major
   `[[SIMD2<Double>]]` (one array per position) plus `scanWidth` and `scanHeight`, and refuse (nil) a count that does not equal
   width · height, as the source refuses a malformed scan. Two input forms, both public: **detector pixels** with
   `originX`, `originY` (Float) and `invAngstromPerPixel` (Double), converted exactly as the source does (`Float` subtraction,
   then `Double`, then the scale: so an app passing its peaks' `Float` positions as `SIMD2<Double>(Double(x), Double(y))` gets
   bit-identical vectors), and **calibrated Å⁻¹ vectors about the origin** (`invAngstromPerPixel` 1, origin 0). `BraggPeak` in
   `experimentalVectors` and `claims` becomes `SIMD2<Double>` the same way.
3. **`AnalysisCancellationToken` → the kit's convention** (slice 3b's G7): `progress: (@Sendable (Double) -> Void)?` and
   `isCancelled: (@Sendable () -> Bool)?`. A cancelled run returns nil (no partial map), as the source.
4. **Concurrency.** `map` keeps `DispatchQueue.concurrentPerform` over positions and gains `maximumWorkers: Int? = nil`
   (nil = the active processor count, as the source). Each position's result depends only on its own vectors, so the worker
   count moves no number (gate G1 runs 1 and 8). Swift 6 strict concurrency with 0 warnings; the unchecked pointer box may be
   replaced by per-worker buffers merged after, if the result is bit identical.
5. **Access and names.** `package` → `public`; no app type, no app or suite name, no repository path in Sources or test names.
   The source's comments carry the measurements that set each default (dates, datasets, Gate B records): they are kept as
   *the reason for the default* with the dataset named (Thronsen et al. 2024, dataset A; "a 4x-binned cube at 0.045741 Å⁻¹ per
   pixel"), without the app's file paths.

Nothing else changes: every expression on a numeric path is the source's, operand for operand.

## 3. The 24 settings: parameters with the tuned values as documented defaults

They are already parameters in the source; the port keeps every field `public var`, every default unchanged, and adds
`public static let default`. The values are the owner's, tuned on the owner's and Thronsen et al.'s data; each field's doc
comment says what it guards and what moved it (the measurement), so a stranger can tell a tuned policy from physics.

`PhaseVectorSettings` (17): `directBeamRadiusInvAngstrom` 0.15, `maximumVectorInvAngstrom` 0 (no outer reach),
`matrixToleranceInvAngstrom` 0.02, `pairRadiusInvAngstrom` 0.02, `minimumVectors` 2, `matrixFallbackExplainedFraction` 0 (off),
`minimumMatchedVectors` 3, `friedelPairMinimumMatchedVectors` 2, `chanceMatchMultiple` 5, `notIndexedAboveInvAngstrom` 0.015,
`minimumPhaseContrastInvAngstrom` 0 (off), `orientationRelationshipToleranceDeg` 10, `completenessAwareCrossPhaseRanking` false,
`classificationRule` `.search`, `directMatrixMaximumVectors` 1, `residualCutoffInvAngstrom` 0.07,
`knownVariantsMinimumSpecificReflections` 1.

`PhaseReferenceSettings` (7): `kMaxInvAngstrom` 1.6, `excitationSlabInvAngstrom` 0.05, `excitationWidthInvAngstrom` 0.03,
`minimumIntensityFraction` 0.05, `maximumVectorsPerEntry` 48, `inPlaneStepDeg` 2, `maximumEntries` 4000.

Not settings (constants of the method, registered as such): the 1e-9 length floor of a projected vector, the 1e-12 zone-axis
norm floor, the 0.9 switch of the detector basis's reference axis (§4), the low-index set u, v, w in −2…2 reduced and
deduplicated up to sign (62 axes), the `fitZoneAxis` defaults (`inPlaneStepDeg` 5, `sampleLimit` 400: arguments already).

## 4. Definitions (what the numbers are)

- **Reference vector.** For a phase, a zone axis [uvw] (n = (u·a + v·b + w·c)/|…|, the kit's real lattice rows) and the
  reflections g with |g| ≤ kMax: s = g·n (flat Ewald sphere; no curvature); kept when |s| ≤ slab; excited intensity
  I = |F|² · exp(−(s/w)²), w = max(width, 1e-9); in-plane q = (g·e1, g·e2), kept when 1e-9 < |q| ≤ kMax and I > 0; the floor
  keeps I ≥ minimumIntensityFraction · max I; above `maximumVectorsPerEntry` the strongest are kept (ties: shorter first);
  relative intensity I / max I; sorted by |q|. Each in-plane step θ_k = k · 2π/N, N = max(1, round(360/step)) (one step,
  θ = 0, for a step that is not finite, ≤ 0 or > 360) rotates every q by θ_k.
- **The frame (the py4DSTEM frame rotation).** e1 = normalise(n × r), e2 = n × e1, with r = x̂ when |n_x| < 0.9 else ŷ.
  py4DSTEM f050d207 (`Crystal.parse_orientation`, default `proj_x`) builds x' = (z × (ŷ × z)) normalised (x' from ŷ × z, or
  −ẑ near ŷ), y' = z × x'. The two in-plane frames differ by a rotation about n (and possibly a reflection) that depends on
  the axis: the pins record, per zone axis, the angle φ_py and the handedness such that py4DSTEM's (g·x', g·y') equals the kit's
  q rotated by φ_py (or reflected then rotated). The library's own frame is the contract (it is what the apps' maps were made
  in); φ_py is the registered conversion for a reader comparing with py4DSTEM. py4DSTEM's excitation weight is
  exp(−s²/(2σ²)) with a curved Ewald sphere: the kit's equals it at σ = w/√2 and λ → 0.
- **Euler angles.** `EulerAngles(py4DSTEMOrientationMatrix: M)`: R = Mᵀ (crystal → lab), then extrinsic z-x-z as SciPy
  `Rotation.from_matrix(R).as_euler("zxz")` (Bunge φ1, Φ, φ2), with SciPy's gimbal-lock rule (|sin Φ| ≤ 1e-8: φ2 = 0, the
  observable rotation in φ1), angles wrapped to (−π, π], stored as Float. `crystalToLabMatrix` = SciPy `from_euler("zxz")`.
  `ACOMOrientation.matrix(detectorBasis: B, inPlaneAngle: a)` = B · R_z(a) with R_z's columns (c, −s, 0), (s, c, 0), (0, 0, 1).
- **Matching.** Experimental vector u = ((x − x0) · scale, (y − y0) · scale), kept when direct-beam radius < |u| and (outer
  reach ≤ 0 or |u| < reach). `nearest`: the closest reference within the radius (a length-band prune that never changes the
  answer). `score`: the mean, over unique references hit, of their best distance. Chance fraction p = min(1, V · ρ² / R²) for V
  reference vectors, pair radius ρ, accessible radius R. The verdict rules (`.search`: matrix removal at its tolerance, the
  minimum-vectors exclusion, the matched-vector floor with the Friedel-pair exception, the chance guard, the not-indexed
  distance, the contrast margin, the orientation-relationship filter, the matrix challenge, the matrix fallback; `.knownVariants`:
  Thronsen et al.'s argmin with the residual cutoff and the specific-reflection guard) are the source's, read as the definition.

## 5. References and truth

- **Planted truth (the gate):** the apps' tests of the ported types, moved to Swift Testing with every expectation kept:
  mac4DSTEM 7d6fdbd1 `mac4DSTEMTests/PhaseVectorMatchingTests.swift` (classes `PhaseVectorMatchingTests`,
  `PhaseVectorGateBTests`, `PhaseVectorResolutionTests` and the zone-axis tests of the ported types; 2637 lines, 69 tests in all,
  of which the presentation tests of `PhaseMapPresentationTests` stay in the app) and `PhaseReferenceLibrarySlabTests.swift`
  (2 tests); the orientation-result tests of the ported types from `ResultOrientationTests.swift`. Planted scans: a known
  phase at a known zone axis and rotation, labelled that phase; random vectors at chance; a matrix-only scan all matrix.
- **The apps' compiled code, bit for bit** (`tools/slice-4/apps-parity`, as slices 2 and 3b): the two source files and the
  value types of `OrientationResult.swift` from `git show`, with minimal shims for `Crystal` (backed by the same reflection
  list), `BraggVectors`, `BraggPeak` and the cancellation token, compiled at -O and -Onone into `.build/apps-parity-4/`; dumps of
  every library entry (all vectors, bit patterns) and every per-position result on the fixtures of §6.
- **py4DSTEM f050d207** (the venv under `../References/CuboKit/py4dstem-venv`, read-only, through its link; every run under
  `../tools/memcap.sh`): `Crystal.generate_diffraction_pattern` on the same phases and zone axes at σ = w/√2,
  `tol_excitation_error_mult` = slab/σ, `tol_intensity` 0, `k_max` = kMax, and a wavelength driven to λ → 0 (or py4DSTEM's
  structure factors and g vectors with the flat excitation applied in NumPy, if its API will not take λ → 0): the set of hkl, |q|,
  and the relative intensity before the density cap; and `parse_orientation` for φ_py. At the real 200 kV wavelength the
  difference (reflections in or out of the slab because of curvature) is a parity row, not a gate.
- **SciPy 1.18.1** (the pymatgen venv): `Rotation.from_matrix` / `as_euler("zxz")` / `from_euler` for the Euler pins.
- **NumPy from these definitions** (not from any Swift): `chanceMatchFraction`, `score`, `nearest`, `containsFriedelPair`,
  `inPlaneSteps`, the 62 low-index axes, `PhaseVectorResolution`'s pixel numbers and `scaledToDetector`.
- **No reference exists for the verdict rules themselves.** py4DSTEM's vector-based phase tool (`crystal_phase.py`,
  `Crystal_Phase`) fits phase weights by non-negative least squares on orientation-matched Bragg peaks: a different method. It
  is a parity row ("not held: different method"), with an agreement count on the planted scans if time allows, never a gate.
  Thronsen et al.'s repository carries no licence and is not used; the method is from the paper (CC BY 4.0).

## 6. Pins and fixtures (generated before the Swift)

`tools/slice-4/pins.py` (under `../tools/memcap.sh`, deterministic, no clock) writes `Tests/DiffractoOrientationTests/Fixtures/
slice-4-pins.json`: phases Al (Fm-3m, 4.0495 Å), Au, Cu, θ′-Al2Cu, β″-Mg5Si6 (the kit's slice 1 fixture), one hexagonal (Mg) and
one triclinic cell; zone axes [001], [011], [111], [112], [123] and the β″ [010] at the default settings; per entry the hkl set,
|q| and relative intensities (py4DSTEM), φ_py and handedness; 1000 seeded random rotations plus the gimbal cases (Φ = 0, π, and
within 1e-9 of them) through SciPy; the closed-form pieces above. `tools/slice-4/apps-parity/run.sh` writes
`slice-4-apps.json` from the apps' code on the same phases and on planted scans (seeded): a 16 × 16 Al [001] matrix scan with
θ′ and β″ positions planted at known rotations, a random-vector scan, a scan of one Friedel pair per position, and each
setting at one non-default value.

## 7. Gates

- **G1 bit identity with the apps' code:** every library entry (hkl, q, |q|, relative intensity, bit patterns) and every
  `PhaseVectorResult` field on the §6 scans equal `slice-4-apps.json`, at `maximumWorkers` 1 and 8, at -O; -Onone agrees with -O
  or the difference is recorded.
- **G2 py4DSTEM:** the hkl sets equal; |q| within 1e-12 Å⁻¹; relative intensity within 1e-9 relative (the measured λ → 0 floor
  is recorded and, if larger, the tolerance is amended with the reason before the test is written); the kit's q rotated by φ_py
  equals py4DSTEM's within 1e-12.
- **G3 Euler:** φ1, Φ, φ2 equal SciPy's after Float rounding (|Δ| ≤ one Float ulp of the angle + 1e-12) on all 1000 + gimbal
  cases, modulo 2π on the wrap; `crystalToLabMatrix` within 1e-12 of SciPy's matrix; the round trip M → angles → matrix within
  Float precision (2e-7).
- **G4 closed forms:** the NumPy pins within 1e-15 relative (exact where integer).
- **G5 planted truth:** every ported app test passes with its expectations unedited (a change of expectation is a finding).
- **G6 parameters:** each of the 24 settings has a test that a non-default value moves the outcome the documented way, and an
  invalid value (NaN, negative, 0 where meaningless) is handled as the source handles it (documented, never a trap).
- **G7 public API:** `tools/consumer-check` gains the slice 4 symbols (library, matcher, settings, results, Euler), with negative
  controls; Swift 6, 0 warnings; DocC article for `DiffractoOrientation`; the README gains a phase-map example only if it
  compiles in `tools/readme-consumer`.
- Every new test broken once (fail-first log); one Gate B refuter (Sonnet) on the frozen commit.

## 8. Planned mutants (each must turn a test red)

1. excitation weight exp(−(s/w)²) → exp(−s²/(2w²)); 2. slab test |s| ≤ slab → <; 3. the 0.9 basis switch → 0.5; 4. the density
cap's tie-break dropped; 5. `inPlaneSteps` rounding → truncation; 6. `nearest` radius ≤ → <; 7. the chance guard removed; 8. the
Friedel exception removed; 9. `notIndexedAbove` compared with ≥ instead of >; 10. matrix removal at the pair radius instead
of its own tolerance; 11. the known-variants specific-reflection guard off; 12. the Euler transpose dropped; 13. the gimbal
threshold 1e-8 → 1e-4; 14. the direct-beam cut ≥ instead of >; 15. `chanceMatchFraction` without the clamp; 16. the
experimental vector computed in Double before the subtraction (the Float order lost).

## 9. Amendments

### 2026-10-15, 15:55 CEST: from the pins (`tools/slice-4/pins.py`, written before the Swift tests; exit 3 on three registered numbers, none edited)

1. **G2 intensity tolerance.** py4DSTEM f050d207 holds its structure factors in complex64: its relative intensities differ
   from float64 by up to 1.15e-5 (β″'s weak reflections; 3.1e-6 θ′; 1.8-3.6e-7 elsewhere). G2 becomes: against py4DSTEM, hkl
   sets equal, |q| within 1e-12 Å⁻¹, the rotated q within 1e-12, relative intensity within 5e-5 relative; against the NumPy
   float64 reference in the same file (py4DSTEM's convention: Lobato table, no Debye-Waller, occupancy, 1/V, no 2π, the flat
   excitation), relative intensity within 1e-12. Measured λ → 0 through py4DSTEM's own API (`setup_diffraction(1e30)`): 0 on six
   phases, 1.5e-14 on β″ against NumPy-flat; hkl sets identical, positions 2.2e-16.
2. **The low-index set has 49 axes, not 62.** The source's code (u, v, w in −2…2, divided by their gcd, deduplicated up to
   sign) gives 49; "62" came from a comment in the source ("a 62-axis sweep"), which is wrong about its own code. The code is
   the contract (G1); the kit's doc says 49.
3. **The gimbal threshold.** The source switches to the gimbal branch at |sin Φ| ≤ 1e-8; SciPy 1.18.1 locks at Φ below about
   1e-7. For Φ in (1e-8, 1e-7) the two give different, equally valid splits of the same rotation between φ1 and φ2. The kit
   keeps the source's rule (G1, bit identity with the apps); G3 compares the rotation matrix there, not the angles (the pins flag
   those cases `rules_agree: false`); a parity row records the window. Mutant 13 stays caught at Φ = 1e-6 and 1e-5.
4. **`PhaseVectorResolution` defined** (missing in §4): pixels = tolerance / scale (NaN when scale ≤ 0); advice only when
   scale > 0, the pair radius in pixels is finite and below 1; `scaledToDetector` sets pair radius and matrix tolerance to one
   pixel, the not-indexed distance to 0.75 pixel and the direct-beam radius to max(setting, 3 pixels), unchanged when scale ≤ 0.
   Tested against this closed form in the Swift test.
5. **Definitions the pins chose, recorded:** `score` = each experimental vector assigned its nearest reference within ρ, the
   mean over the references hit of their best distance; `containsFriedelPair` = some pair with |vᵢ + vⱼ| ≤ radius; `inPlaneSteps`
   rounds half away from zero (steps 11 and 13 separate rounding from truncation). If the source differs, the source is the
   contract and the pin is corrected with this amendment's reason.
6. **Structures.** θ′-Al2Cu has no definition in the apps or the kit: the pins use a body-centred tetragonal cell, a 4.04 Å,
   c 5.80 Å, Cu at (0,0,0) and Al at (0,½,¼), (½,0,¼), (0,½,¾), (½,0,¾), plus the body centring (Al4Cu2), written out site by
   site (the owner may confirm the structure; it is a test phase, not a claim about θ′). The triclinic phase is slice 1's cell
   (4, 5, 6 Å; 80, 85, 95°) with two sites. β″ is the apps' preset in its site order.
7. **The 200 kV parity record** (py4DSTEM's curved Ewald sphere against the flat one, same σ and slab): 33 of 37 axes keep the
   same hkl set (relative intensities 1-24 % apart); θ′ [112] +6, θ′ [123] +21, β″ [111] +51, β″ [123] −54. The frame rotation
   φ_py is a proper rotation (det +1) on every axis; the per-axis 2 × 2 map is in the pins (`frame_map`).

### 2026-10-15, 16:15 CEST: from the apps-parity harness (`tools/slice-4/apps-parity`, before the G1 tests)

The harness compiles the apps' code from 7d6fdbd1 with no hand-written shim (`BraggPeak`, `BraggVectors`, `EulerAngles` and
`ACOMOrientation` cut out of their files by line range) and dumps 72 `map` runs, four libraries, fits, claims, 270 Euler cases
and the detector bases (`slice-4-apps.json`, 13.5 MB; git stores it compressed, about 2.2 MB). Every `PhaseVectorSettings` field
moves at least one position on the scans.

1. **-O against -Onone.** Every `map` result, fit, claim and basis is byte identical. Three last-bit differences, all from
   code generation, none in a verdict: Im F of some β″ and triclinic reflections (|Δ| ≤ 2e-18; slice 2's amendment 1), the
   in-plane step k = 13 of every library (at -O `rotate` takes cos and sin of one angle as one `__sincos_stret` call, which
   differs from two libm calls at that angle), and `crystalToLabMatrix` in 2 of 270 cases (same cause). G1 is at -O against the
   -O dump, as slices 2 and 3b; the kit's `rotate` keeps cos and sin of one angle in one function so the optimiser pairs them.
2. **What the source refuses, precisely** (the port matches it, and §2.2 is read this way): only `map` refuses a position count
   that is not width × height (and a width or height ≤ 0, and an orientation-relationship filter that leaves no candidate);
   `fitMatrixOrientation` and `fitZoneAxis` do not check the count. On cancellation `map` returns nil, `fitMatrixOrientation`
   its best so far and `fitZoneAxis` a partial list. The port's additions (a malformed scan gives `[]` from `fitZoneAxis`, an
   out-of-range matrix entry gives nil from `map` where the source traps) are refusals of inputs the source does not define,
   recorded in `docs/consumers-0.3.0.md`.
3. **`inPlaneSteps`** is θ_k = Double(k) · (2π / Double(N)) with the quotient computed once; N = round(360/step) rounding half
   away from zero. §4's "k · 2π/N" is read as this expression.
4. **Reflection order** of the source: h, then k, then l ascending over ±ceil(kMax · |aᵢ|), (0,0,0) skipped, |g| ≤ kMax
   inclusive, |F| > 1e-4, then sorted by |g|; in all six dumped lists every exact tie is in (h, k, l) ascending order. The port's
   order (|g|, then hkl ascending, on the kit's kept list) is this order; G1 checks it list by list.
5. **A finding in the source: phase ties depend on the hash seed.** `classify` (`.search`) ranks phases out of a Dictionary,
   so an exact score tie between two phases is broken by the process's random hash seed: the same input can be labelled
   differently in two runs. The harness pins `SWIFT_DETERMINISTIC_HASHING=1`; a random-seed run was byte identical on the
   dumped scans (no exact tie there). Planned fix in the kit (tier 1, a test first): ties broken by the lower phase index, the
   order the library lists the phases; identical to the source wherever the source is deterministic; named in
   `consumers-0.3.0.md` for the app.
6. **Euler details** beyond §4: Φ = acos of the clamped R22; |sin Φ| > 1e-8 for the general branch; wrapped angles with
   |angle| < 1e-12 snapped to 0; angles stored as Float and `crystalToLabMatrix` built from the Float-rounded angles.
