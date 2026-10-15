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

(none yet)
