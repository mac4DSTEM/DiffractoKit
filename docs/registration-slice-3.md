# DiffractoKit slice 3: DiffractoGeometry, the geometry the apps switch to (registration)

Registered 2026-10-15, about 01:35 CEST, before any code, by the DiffractoKit session (the owner's instruction relayed by
the coordinator at 00:55: "begin to use DiffractoKit's slice 3 in the apps"). Tier 1: wavelength, camera length, Δq, the
accepted (hkl) family and the fitted ring are numbers a scientist reads. Seeds read only by `git show`: mac4DSTEM
`7d6fdbd1` (pushed) and macTEM `9482831` (local, its V1-8 S2c commit). What is written here is fixed before the code; a
change after code starts is a dated amendment with its reason.

## 0. Scope, and what is deliberately not slice 3

Slice 3 is two parts, ported in this order, each with its own pins, tests, break-once proof and one Gate B refuter:

- **3a, closed forms and decisions** (about 700 seed lines): the electron wavelength; the camera length and Δq from one
  ring of known d, with the fit-report refusals (macTEM `RingCalibration.swift`); the indexing decision of a measured
  |q| against a ring list (macTEM `DSpacingIndexing.index`, over slice 2's `KinematicFamily`); zone-axis primitives
  (mac4DSTEM `PhaseReferenceLibrary`'s low-index axes, Cartesian zone axis and projection, `ACOMOrientation.detectorBasis`,
  and the excitation error of `OrientationPlan.project`).
- **3b, the ring and ellipse fit** (about 1100 seed lines): mac4DSTEM `7d6fdbd1` `EllipseCalibration.swift` (conic, angular
  coverage, the one-ring check, the asymmetric-Gaussian profile fit, ADR 069's two starts and the 0.30 px floor guard),
  the three-field `DiffractionPattern`, `Calibration.ellipseTransform` / `ellipseCorrectedOffset`, and macTEM's
  `RingCalibration+Fit` (fit → report → refusals → calibrate).

**Not slice 3** (the coordinator's brief named vector matching; this registration moves it, with the reason, and the owner
can overrule): `PhaseVectorMatching` (1698 lines) with its reference library build (572) is app-coupled (`BraggVectors`,
the ACOM types), carries 24 tuned settings measured on the owner's datasets the kit does not hold, and ships
`validation: "none"`; it is slice 4, registered on its own after 3b, with `BraggVectors` replaced by a plain
`[[SIMD2<Double>]]`. Also not here: the ACOM orientation plan and templates (`DiffractoOrientation`), `QCalibration`
from Bragg-vector maps, `CalibrationUnitConversion` (units metadata, a readers concern), the annulus guard (refuted as
registered, never landed in mac4DSTEM; not ported), the seed check (not landed).

## 1. What 3a makes (module `DiffractoGeometry`, depends on DiffractoCore and DiffractoScattering)

- `ElectronWavelength.angstroms(volts:)`: λ = hc / √(eV (2 m₀c² + eV)), V in volts, CODATA 2018 (h, c, e exact; m₀ =
  9.1093837015e-31 kg), hc/e and m₀c²/e DERIVED in code from them, as macTEM `RingCalibration.wavelengthAngstroms`.
  **The one wavelength of the suite**: mac4DSTEM's `DPC.electronWavelengthAngstrom` (12.2639/√…) is 2.7e-5 below it and
  its `ParallaxPreprocessing` form 1.2e-7 below (rounded constants); stated in `docs/consumers-0.0.3.md` as what a
  mac4DSTEM bump moves.
- `RingCalibration` (macTEM's, renamed nothing): `twoTheta(wavelengthAngstroms:d:)` = 2 arcsin(λ/2d), NaN unless λ/2d <
  1/√2; `cameraLengthPixels(radius:twoTheta:)` = R / tan 2θ, NaN unless R > 0 and 0 < 2θ < π/2; `qPerPixel(d:radius:)` =
  (1/d)/R (no 2π, needs no λ); `cameraLengthMillimetres(pixelPitchMM:cameraLengthPixels:)`; `calibrate(radius:d:
  highTensionVolts:pixelPitchMM:)` returning every quantity with the high tension it used (never defaulted); the
  `FitReport` and `refusals(_:highTensionVolts:)` in macTEM's order (model, coverage 30 of 36 bins or sparse, σ at the
  annulus width, residual > 0.22, ellipticity > 0.093, no high tension; a non-finite input gives `[.invalidInput]` alone),
  `withholdsDeltaQ` false only for no high tension and the Bragg limit. `EllipseCalibrationModel` (the fit's model tag)
  comes with 3a as the small enum it is.
- `RingIndexing` (macTEM's `DSpacingIndexing.index` and its types, generalised): `index(q:sigma:in:kMax:coverage:)` over
  `[KinematicFamily]` (or `[ReflectionFamily]`) built to k_max: z = |q − g|/σ_q, δ = g/q − 1, candidates by (z, g);
  `ambiguous` when two or more have z ≤ k or the margin (z_next² − z_best²)/2 < ln 19 (exactly ln 19 accepted);
  reach r = √(k² + 2 ln 19)·σ_q; `beyondKMax` when q + max(r, |q − g_best|, |q − g_next|) > k_max (or q + max(r, nearest)
  with nothing within); `noneWithinTolerance(nearest)`; `invalidInput`; the verdict order invalidInput, ambiguous,
  beyondKMax, noneWithinTolerance, accepted; `defaultCoverage` 3; `validation: "none: …"`. macTEM's `lattice(of:)`,
  `orbit(of:in:)`, `Family` and `families(of:kMax:)` are NOT ported: slices 1 and 2 replace them for every crystal system
  (so `unsupportedLattice` disappears; a tetragonal, orthorhombic, trigonal or triclinic cell is indexed).
- `ZoneAxis`: the 49 low-index directions (mac4DSTEM `PhaseReferenceLibrary.lowIndexZoneAxes`), `cartesian(_:cell:)`
  (u a₁ + v a₂ + w a₃ normalised), `detectorBasis(_:)` (reference x unless |n.x| ≥ 0.9, then y; e1 = normalize(n × ref),
  e2 = n × e1; columns e1, e2, n), the projection of a reflection list onto a zone axis with the flat sphere (s_g = g·n)
  or the curved one (s_g = g·n + λ|g|²/2, mac4DSTEM `OrientationPlan.project`'s sign, its DEVIATION from py4DSTEM's
  lab-z sign stated), the slab |s_g| ≤ s_max, (x, y) = (g·e1, g·e2).

## 2. What 3b makes

`DiffractionPattern` (width, height, row-major Float pixels; qx = ROW, qy = COLUMN, as the apps), `EllipseCalibration`
(public: `fit1D`, `fitBestAvailable`, `fitAmorphousRing`, the fit, the model enum, the errors) bit for bit with mac4DSTEM
`7d6fdbd1` (every constant of the digest kept: coverage 36/30/12, residual caps 0.2/0.4/0.22, one-ring ratio 1.10, floor
0.30 px, axis ratio 0.2–5, MAD 4·1.4826 and 70 % inliers, the annulus boundaries `inner < r ≤ outer` in `fit1D` and
`inner < r < outer` in `fitAmorphousRing` kept for identity), `Calibration.ellipseTransform` / `ellipseCorrectedOffset`,
and `RingCalibration.fit(pattern:…)` (macTEM's `+Fit`, on the 7d6fdbd1 fit: macTEM's own copy is pre-ADR-069, so its
S3-b numbers move when it imports the kit; stated in the consumers note).

## 3. References and pins

- **3a.** macTEM `tools/ring-calibration-parity/reference.py` (decimal, 60 digits; 320 cases: V ∈ {60, 80, 120, 200, 300}
  kV, Au {111} {200} {220} {311} {422} at a 4.0782 and Al {111} {200} {220} at 4.0495, R ∈ {40, 80, 200, 250} px, p ∈
  {0.014, 0.05} mm, 68 refusal plants, 18 Bragg-limit plants, 10 radius plants), ported to `tools/slice-3/ring-pins.py`
  writing a fixture; band 1e-12 relative. macTEM `tools/indexing-parity/reference.py` (NumPy) ported to
  `tools/slice-3/index-pins.py` on the kit's own ring lists (slice 2's), T1 deterministic plants (every family, u ∈ {1e-4,
  1e-3, 3e-3, 1e-2, 3e-2}, s ∈ {±2.9, ±2, ±1, ±0.5, 0}) and T4 refusal plants; verdict, payload and list equal, δ 1e-14,
  z 1e-10·max(1, z), margin 1e-9·max(1, margin). Zone axes: pymatgen's `Lattice.get_cartesian_coords` for the axes, its
  reciprocal lattice for g · n, to 1e-12.
- **3b.** The apps' own fixtures as planted truth, reproduced bit for bit or to their registered tolerance: mac4DSTEM's
  `EllipseNoisyRingTests` plants (SplitMix64 + Box–Muller renderer, identical Float pixels; its bands 1 % radius, 0.01
  b/a, 1 px centre) and its four `EllipseStopRuleTests` rings (R(1 + 1.75 (s/R)²), a/b = 1 + ε to 1e-8, θ to 1e-6);
  `CalibrationEllipseRow29Tests`; macTEM's `RingFitCalibrationTests` plants (on the 7d6fdbd1 fit: re-measured, the
  differences stated); the clean-sweep rings of mac4DSTEM `tools/ellipse-calibration-test` against SciPy 1.18.1's
  `least_squares` optimum (centre 1e-5 px, axes 1e-7, θ 1e-6 rad; NumPy and SciPy are in `References/pymatgen-venv`);
  and the apps' compiled `EllipseCalibration.swift` at 7d6fdbd1 (by `git show`, as slice 2's apps-parity) on the same
  plants, bit for bit. The owner's Ni65Cu35 32×32 mean-pattern fixture is used only with a `../data/README.md` manifest
  line, read via `$MACEM_DATA`, never copied.
- **py4DSTEM** at f050d207 is named where the apps cite it (fit_ellipse, the excitation error); its file:line citations are
  checked against the pin before they are written into the kit.

## 4. Truth (independent of every oracle)

- λ at 80, 120, 200, 300 kV: 0.041757160772834212, 0.033492152726723738, 0.025079340450548003, 0.019687489006848796 Å
  (60-digit arithmetic, macTEM's literals); λ = λ_classical / √(1 + eV/(2m₀c²)).
- Au {111} (d 2.3546 Å), R 150 px, 200 kV: 2θ 0.010651260937060768, L_px 14082.307097295722, Δq 0.0028313372405787253 Å⁻¹/px,
  L 197.15229936214011 mm at p 0.014 mm; the exact length exceeds the small-angle one by 4e-5 (200 kV Au {111}) and 9.5e-4
  (80 kV Au {422}); the ellipticity bound 0.093 is where (a+b)/2 and √(ab) differ by 9.8864669599e-4.
- Indexing: macTEM's dyadic lists (families {1, 1.5, 2}, σ 0.125, q 1: z 0, next 4, margin 8; z = k exactly accepted;
  margin exactly ln 19 accepted by an ulp scan; q + R = k_max exactly) and its 41 named cases re-expressed on the kit's
  lists; on β″ at σ_q/q = 1e-2 exactly eight families are nameable at s = 0 ({2 0 0}, {0 0 1}, {2 0 -1}, {2 0 1}, {1 1 0},
  {1 1 -1}, {3 1 0}, {3 1 -1}), 23 at 3e-3, 4 at 3e-2 (macTEM's ADR 024).
- Zone axes: β″ [0 0 1] is (cos β, 0, sin β) in the apps' frame; the 49 axes are the distinct primitive [u v w] with |u|,
  |v|, |w| ≤ 2 up to sign... (the exact list is mac4DSTEM's, pinned by value); a reflection with g·n = 0 lies in the zero-
  order Laue zone; the detector basis is orthonormal and right-handed (det +1) for every axis.
- 3b: the apps' planted rings and their bands as above; a perfectly circular clean ring recovers R(1 + 1.75 (s/R)²) on
  the conic path and R on the profile path.

## 5. Gates

1. Pins within their bands; the apps' code bit for bit where the registration says so.
2. Truth exact or to its stated bound.
3. Every new test broken once (`tools/slice-3/mutants.py`, slice 1's pattern touching the module after each edit).
4. One Gate B refuter per part on its frozen commit; Gate D when a cause is unproven.
5. `tools/test.sh` green with the stamp; consumer-check extended; `docs/consumers-0.0.3.md` says what mac4DSTEM and macTEM
   change to import DiffractoKit for these (their own sessions do the switch).

## 6. What is not claimed

No real pattern is calibrated or indexed here; validation on the owner's data is each app's to register. `accepted` means
"19 times as likely as the next listed family at equal priors under a Gaussian σ_q", not "the ring is this family"
(macTEM's T3 and T5 numbers stand). ADR 069's measured rates (2009 synthetic rows: wrong 414 → 79, right 1401 → 1814) are
mac4DSTEM's, on its harness sets, cited, not re-claimed.

## 7. Amendments

(none yet)
