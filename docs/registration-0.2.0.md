# DiffractoKit 0.2.0: the generality pass (registration)

Registered 2026-10-15, about 14:10 CEST, before any code, by DiffractoKit's third session (the owner's 0.2.0 brief, relayed by
the coordinator at 13:52: "every suite policy a parameter with the suite's value as the documented default"). Tier 1 where
a default could move a number (it must not), tier 2 for the parameter paths, tier 3 for wording and docs.

## 1. What becomes a parameter

A **policy** is a value chosen by the suite (an intensity cut, a refusal bound, a ring-grouping tolerance, a fit acceptance
rule, a tolerance tuned on the owner's data) rather than fixed by physics or by the reference implementation's definition.
Each becomes a parameter with the suite's value as its documented default, grouped where several act together:

- `StructureFactors`: the |F| cut (1e-4 Å⁻²) — already `cut:`; documented as a policy.
- Reflection lists: the k_max refusal (4 Å⁻¹), the expected-count cap (1e7), the ring-grouping tolerance (1e-9 relative on
  |g|, the "ring rule").
- `CrystalStructure` validation: the close-contact limit (0.5 Å) and the split-site occupancy sum (1.02).
- `SymmetryDetection`: the position (5e-3 fractional) and occupancy (0.02) tolerances, the metric tolerances (1e-9) — the
  metric ones already parameters; the translation grid (1/24).
- `CIFImport`: the refusals that are policy rather than grammar (the precision-limited rotation test's factor 10, the
  coefficient bound 1000), as a `CIFImport.Policy` with `.default`.
- The ring fit: the coverage bins (36, 30, 12), the residual caps (0.2, 0.4, 0.22), the one-ring ratio 1.10, the ring-width
  floor 0.25 / 0.30 px and the lowered floor 0.01 px, the MAD cut 4 · 1.4826 and 70 % inliers, the axis ratio 0.2–5 —
  as `EllipseCalibration.Policy` with `.default` — and the calibration's refusal bounds (30 bins, residual 0.22,
  ellipticity 0.093) and the ring-evidence bounds (§1.3 of registration-0.0.4) as `RingCalibration.Policy` with
  `.default`.
- `RingIndexing`: the coverage k = 3 (already a parameter) and the ln 19 margin.

Not policies (stay constants): physical constants (CODATA), the reference table (py4DSTEM's), closed forms, the
Bragg-limit bound 1/√2 (geometry), numerical-method constants that change no accepted result (LM damping, step limits).
The lane that ports this lists every constant it classified, with its class, in an amendment.

## 2. Gates

- **The default path is bit identical**: every existing test, pin and parity fixture (35 parity rings at -O, the slice 1–3
  pins, the 0.1.0 tests) passes unchanged with the defaults; `Policy.default` equals the constants they replace.
- **Each parameter path has a tier-2 test**: a non-default value changes the outcome in the stated direction (a lower cut
  keeps more reflections, a looser coverage bound accepts a gapped ring, and so on), and an invalid value (negative,
  NaN, out of range) is refused or clamped as documented.
- Break once: one mutant per policy group (the default changed) must turn a default-path test red.
- consumer-check: each `Policy` type, its `.default` and one non-default use, with negative controls.

## 3. Also in 0.2.0 (tier 3 / tier 2)

No app or suite name in Sources, test comments, Package.swift, README, DocC or the example (`git grep -i` is the proof;
the registrations, consumers notes, status, parity, NOTICE and tools READMEs keep them as the record). A README and a first
example that stand alone; a second consumer package that compiles only the README's code; a DocC catalogue per module with
one article; the platform floor at macOS 15 and iPadOS 18 if build, tests and numbers are unchanged there (the owner's
decision of 13:58), else macOS 27 with the reason; no build products in the repository, every `tools/` file referenced.

## 4. Amendments

### 2026-10-15, 15:30 CEST: the classification of every constant, and the gates' results

Written at the candidate (policies as parameters on top of 6c16c9a), from the porting lane's report and the diff. The groups
and values of §1 stand; what the lane found beyond §1 is added. Classes: **policy now a parameter**, **already a parameter**,
**constant** (physics or reference, numerical, measurement definition).

**Policy, now a parameter** (the default is the old value; `Policy.default`):
- `ReflectionList.Policy`: `maximumKMax` 4 Å⁻¹, `maximumReflections` 1e7, `maximumEnumeration` 1e8, `familyTolerance` 1e-9
  (relative, the ring rule), `absenceTolerance` 1e-9.
- `CrystalStructure.ValidationPolicy`: `closeContactLimit` 0.5 Å, `splitSiteOccupancyTolerance` 0.02 (sum bound 1.02).
- `SymmetryDetection.Policy`: `occupancyTolerance` 0.02, `translationGrid` 1/24, `gridSnapFactor` 3, `precisionLimitFactor` 10;
  `SymmetryOperation(xyz:maximumMagnitude:)` 1000.
- `CIFImport.Policy`: `maximumMergeDistance` 0.10 Å, `minimumMergeTolerance` 1e-6, `maximumCoordinateHalfStep` 1e-2,
  `rightAngleToleranceDeg` 0.05, `maximumOperationMagnitude` 1000, position tolerance, metric tolerances, and the nested
  symmetry and validation policies.
- `EllipseCalibration.Policy` (18 fields: bins 36 / 30 / 12, strong-signal fraction 0.2, one-ring ratio 1.10, sample minimums
  20 and 32, relative signal 1e-6, impulse factor 2, residual caps 0.2 / 0.22 / 0.4, axis ratio 5, MAD cut 4 · 1.4826, inliers
  70 %, width floors 0.25 (+0.05) and 0.01 px) with `RingEvidence.Policy` (8 fields: sector pixels 8, amplitude fraction 0.5,
  significance 3, bins 30, drop 0.25, shift 0.25 px, shift significance 4, secondary peak 0.25).
- `RingCalibration.Policy`: bins 36, occupied 30, residual 0.22, ellipticity 0.093, and the evidence policy.
- `RingIndexing.minimumMargin` ln 19 (a `minimumMargin:` argument on the indexing calls).

**Already a parameter** (kept, documented as policy): the |F| cut (`cut:`, 1e-4 Å⁻²); the indexing coverage k (3);
`SymmetryDetection.Tolerance` (`.fractional(5e-3)` default, `.cartesian`); `metricLengthTolerance`, `metricAngleToleranceDeg`
(1e-9 / default); `search:`; `conventional:`; `UnitCell.niggliReduced(relativeEpsilon:)` and `conventionalCell(tolerance:)`;
the annulus, the maximum residuals and the progress and cancellation arguments of the fits.

**Constant**:
- *Physics or reference:* CODATA-2018 constants and the wavelength; the electron scattering table (py4DSTEM f050d207's); the
  closed forms of the cell, q(R) and R(q); the Bragg-limit bound 1/sqrt(2) (geometry).
- *Numerical:* Levenberg-Marquardt damping and step limits and the like; none changes an accepted result.
- *Measurement definition:* what a sector, a drop, a shift and a profile are in `RingEvidence`; the family keys of a ring.
  (The lane's full per-constant list is in its report; this is the classification by group, not a line-by-line copy.)
- *Policy-like, left constant and listed for a later minor:* the profile-fit plausibility bounds, the cell-volume floor 1e-6,
  the Niggli epsilon default (1e-5, an argument already). Error texts quote the default bounds (except
  `RingCalibration.Refusal.reason(policy:)`). The consumer check has no Policy checks yet.

**Gates' results (to be filled at the freeze):**
- Suite: 352 / 0 merged on 6c16c9a (306 before, no existing expectation edited; the apps-parity fixtures and every pin
  unchanged, so the default path is bit identical); 348 on the lane's base afe6983.
- Mutants: 28 / 28 red (`tools/kit-0.2.0/mutants-policy.py`, one default changed per group, on the lane's base).
- Consumer check: 537 / 537 at the floor, 0 warnings; Policy checks and their negative controls NOT yet added (0.2.x).
- Refuter: not yet run; the 0.2.x session runs the one Gate B refuter on this frozen candidate before the owner tags 0.2.0.
