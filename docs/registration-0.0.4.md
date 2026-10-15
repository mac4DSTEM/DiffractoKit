# DiffractoKit 0.0.4: the kit loop's first asks (registration)

Registered 2026-10-15, 11:14 CEST, before any code, by DiffractoKit's third session. Sources of the asks: slice 3b's
Gate B (registration-slice-3 §7, F1–F3, F6, F7), macTEM's Diffraction-room design (its §0 and §3, read-only, G1, G2, G7,
G8) and macTEM's 0.0.3 switch (five misfits), relayed by the coordinator. What is written here is fixed before the code; a
change after code starts is a dated amendment with its reason. Tiers: §1, §2, §3 and §6 are tier 1 (a number a scientist
reads moves or is newly given); §4 and §5 are tier 2 (behaviour, no number moves); §7 is tier 3.

## 0. What deliberately does not change

- The arithmetic of every fit that the apps' 7d6fdbd1 code accepts on a pattern without impulses and with a dynamic
  range of at least 1 stays bit identical (the 35 parity rings at -O and -Onone, `EllipseAppsParityTests`). Every change
  below either turns a thrown error into a fit on inputs the apps refused, removes an isolated impulse pixel, or adds a
  check that only the calibration path (`RingCalibration`) and an opt-in flag apply.
- `EllipseCalibration.fit1D`'s pre-fit coverage count (`occupiedAngularBins`) stays the apps' count. mac4DSTEM's
  detector-ellipse path (b/a and θ, `fitBestAvailable` with default arguments) keeps its acceptance rule: mac4DSTEM
  refuted a mean-radius annulus guard on its own WS2 data (its `ring-fit-annulus-guard-gateD-2026-10-14.md`: correct b/a
  from rings fitted just past the outer radius), so the radius-only refusals below belong to the calibration, where the
  radius is the number read. mac4DSTEM opts in with `requireRingEvidence: true` in its own lane if it wants them.
- Two rings closer than about twice their width are one ring to any single-ring profile model; no refusal is claimed for
  them (§1.6).

## 1. The ring-fit amendment (slice 3b's F1, F2, F3, F6, F7 and macTEM's sigmaAtAnnulusWidth plant)

### 1.1 Isolated impulses (F3: one hot pixel refuses a perfect ring; a strongly negative one moves the fit)

A pixel inside the annulus is an **impulse** when it lies outside the closed range [min, max] of its finite in-pattern
8-neighbours (at least 3 of them) by more than **2 × the annulus's robust range**, the robust range being P99.5 − P0.5 of
the annulus's finite values (nearest-rank percentiles). Impulses are skipped exactly as a non-finite pixel already is, in
`fit1D`'s and `fitAmorphousRing`'s sample sets and in `fitAmorphousRing`'s global maximum. A smooth ring has no pixel
outside its neighbours' range by twice the ring's own contrast, so no parity ring loses a pixel (checked: the parity
fixtures stay bit identical). Gate: a perfect ring with one pixel at 5×, 20× and 1000× its peak, and one at −5× and −50×,
fits with the centre within 0.05 px and R within 0.05 % of the same ring without it; the hot-pixel ring the apps refused as
`insufficientAngularCoverage(1)` is fitted.

### 1.2 Scale (F6: `insufficientSignal` is scale-dependent)

The signal test becomes relative: a fit needs a finite dynamic range above |maximum| × 1e-6 and at least
`Double.leastNormalMagnitude` (was: above max(1, |maximum|) × 1e-6). When the dynamic range is below 1, both fits run on the
samples multiplied by the power of two 2^k that brings the range into [1, 2) and scale I0, I1 and the background back by
2^−k: multiplication by a power of two is exact, so the fit of a ring of amplitude 2^−k·A is the fit of amplitude A with its
intensities scaled, bit for bit, and every absolute floor in the optimisers (`max(1, …)` scales, the 1e-12 damping floor,
the 1e-20 pivot) acts on the normalised problem. Ranges of 1 and above run unchanged. Gate: a ring at amplitudes 2^−40,
1e-6, 1e-3, 1, 1e3 and 1e6 (background in proportion) gives centre, a, b and θ equal to the amplitude-1 fit bit for bit at
powers of two and within 1e-9 relative otherwise.

### 1.3 Ring evidence about the fitted ellipse (F1 and F2)

A new public measurement, `EllipseCalibration.ringEvidence(pattern:fit:innerRadius:outerRadius:) -> RingEvidence`, reads
the data about a fit, model-free in the azimuth:

- Every finite, non-impulse annulus pixel gets its elliptical offset δ = ρR − R about the fitted ellipse (ρ the model's
  normalised elliptical radius, R = (a + b)/2, as the profile model measures it) and its azimuth about the fitted centre,
  in 36 bins as the coverage rule. Only pixels with −(3w_in + 2) ≤ δ ≤ 3w_out + 2 are read (w the fitted ring sigmas; for
  a conic, both are the half-maximum sigma of the seed profile, else (outer − inner)/8).
- **Per sector**, a linear least-squares fit y = α·g(δ) + τ·g′(δ) + β + γ·δ, g the fitted asymmetric-Gaussian ring shape
  of unit peak: α is the sector's ring amplitude (with its standard error from the residuals), t = −τ/α its radial offset
  from the fitted ellipse (to first order), β + γδ a local line under it. A sector needs at least 8 pixels.
- **Amplitude** A = the median α over sectors that have 8 pixels. A sector **holds ring evidence** when α ≥ 0.5·A and
  α ≥ 3 × its standard error. `occupiedBins` counts them; the calibration requires 30 of 36, the same bound as the
  coverage rule (F1: noise in a gap makes pixels bright, not a ring at the fitted radius).
- **Peak**: the profile pooled over all read pixels in 0.5 px bins of δ must fall from its value at δ = 0 (the bin mean
  over |δ| ≤ 0.5) by at least **0.25·A** on each side, the inner minimum taken over −(3w_in + 2) ≤ δ ≤ −w_in and the outer
  over w_out ≤ δ ≤ 3w_out + 2, each inside the annulus (`innerDrop`, `outerDrop` as fractions of A; NaN when a side has no
  bin). A ring is a peak; a disk edge (a plateau on one side) and a blob (monotonic) are not (F2).
- **Geometry**: over the sectors with evidence, a weighted (1/SE²) fit of t = c0 + c1 cos φ + s1 sin φ + c2 cos 2φ + s2 sin 2φ.
  `centreShift` = √(c1² + s1²), `ellipseShift` = √(c2² + s2²) (px), each with its standard error. A fit whose sectors put
  the ring off the fitted centre or ellipse by more than **0.25 px and 4 standard errors** is not the ring the data hold
  (F1's gap noise moves the centre 0.7–2.7 px; F2's cos² texture reads a circle as a 4.6 % ellipse, ±0.69 px at R 30).
- `passes` is all of: occupiedBins ≥ 30, both drops ≥ 0.25, neither shift beyond its bound.

### 1.4 Where the evidence acts

- `EllipseCalibration.fitBestAvailable(…, requireRingEvidence: Bool = false)`: with `true`, the returned fit carries
  `ringEvidence` (a new optional field, nil otherwise), and the profile fit's floor guard changes as §1.5 says. With the
  default, nothing changes.
- `RingCalibration.fit` and `fitRing` (§2) use `true`. `FitReport` gains `ringEvidence` and the annulus's inner and
  outer radii (optional, defaulted; the old initializer still compiles), and `refusals` adds, after `excessiveEllipticity`
  and before `noHighTension`, each Δq-withholding: `ringOutsideAnnulus(radius:inner:outer:)` when R = (a + b)/2 is not
  inside (inner, outer) (F2's blob); `insufficientAngularCoverage(occupiedBins:)` when the evidence holds fewer than 30
  sectors (the case is reused, the count is the evidence's); `notARingPeak(innerDrop:outerDrop:)`;
  `ringGeometryInconsistent(centreShift:ellipseShift:)`. Each refusal appears once (coverage is listed where it is today).

### 1.5 The one-sided floor (F7: one edge sharper than 0.3 px is refused as `ringWidthAtFloor`)

With `requireRingEvidence: true`, a profile fit with exactly one ring sigma at or below 0.30 px is kept when its ring
evidence passes (§1.3), and refused as `ringWidthAtFloor` otherwise; both at the floor is still refused (ADR 069's dead
ring: its R carries no ring). Gate: a ring with σ_in 3 px and σ_out 0.15 and 0.05 px is fitted with R within 0.1 px of the
plant's peak and calibrated; a dead-ring plant (the profile drained, mac4DSTEM's floor-guard family, re-rendered) is still
refused.

### 1.6 What is not claimed

Two rings 1–4 px apart at width 2–3 px are fitted as one broadened ring (slice 3b's F2, +2 to +7 % on the stronger
ring): their sum is a single-peaked profile that falls on both sides, so no test of one ring's shape can refuse it without
refusing real non-Gaussian rings (macTEM's Lorentzian and skewed plants are accepted at 0.0002 % and 0.55 %). Stated in
the consumers note with the bound measured on the plants (|R_fit − R_strong| ≤ the separation) and the remedy that exists:
a crystal's ring list says whether two families lie within the ring's width (`RingIndexing` returns `ambiguous`).

### 1.7 macTEM's `refuse-bound` plant (its ring-fit-parity line 37)

Its plant is a 3-px ring at R 120 on a 1/r³ background of the ring's peak, annulus ±20 px (width 40). The refusal it
expected, `sigmaAtAnnulusWidth`, fires only when a ring sigma reaches the annulus width; the fit reads the ring to 0.014 %
there (macTEM's own log, refusal tier, accepted max). Registered prediction: the kit accepts it (Δq within 0.05 % of the
plant) before and after this amendment, so the plant, not the kit, is wrong, and a plant with a ring sigma at or above the
annulus width (σ 8 px in a 16 px annulus) is refused by name. Both are tests.

### 1.8 Gates (tier 1)

Each fix has a failing test first on the slice 3b plant family it names (re-rendered in Swift with fresh parameters and
seeds, not the refuter's), then the fix; the 35 parity rings at -O and -Onone stay bit identical (no parity ring loses a
pixel, changes scale or meets the floor; where one does, it is named here by amendment before the fixture changes);
SciPy's sweep unchanged; `RingFitCalibrationTests`' eight macTEM plants still calibrated; every new test broken once
(`tools/kit-0.0.4/mutants.py`); one Gate B refuter on the frozen commit.

Refuted by: a parity ring changing a bit; a clean or noisy ring of the existing tests newly refused; any F1/F2 plant of the
new tests calibrated; the gates of §1.1, §1.2 and §1.5 missed.

## 2. G1: a ring-radius fit without a spacing, and the calibration from it

`RingCalibration.fitRing(pattern:centerQX:centerQY:innerRadius:outerRadius:) throws -> RingRadiusFit` (the fit, the report
and the fit's own refusals: everything but `noHighTension`, with nil fit when coverage refused before fitting) and
`RingCalibration.calibrating(_:d:highTensionVolts:pixelPitchMM:) throws -> RingFitCalibration`. `fit` becomes exactly
`calibrating(fitRing(…), d:…)`. Gate: on every ring of `RingFitCalibrationTests`, the parity rings and the §1 plants, the
two calls equal `fit` field for field (`==`), refusals in the same order. Tier 1 (the radius is a number); no number moves.

## 3. G2: exact q from a radius, and its inverse

`RingCalibration.q(radius:cameraLengthPixels:wavelengthAngstroms:)` = 2 sin(½ atan(R/L)) / λ (Å⁻¹, no 2π) and
`RingCalibration.radius(q:cameraLengthPixels:wavelengthAngstroms:)` = L tan(2 asin(qλ/2)). Total: NaN unless every input is
finite and positive and, for the inverse, qλ/2 < 1/√2 (2θ < 90°, as `twoTheta`); an overflow is NaN. Identity: q of the
calibrating ring's own (R, L_px, λ) is 1/d. Pins: `tools/kit-0.0.4/q-pins.py`, the truth in 60-digit `decimal` (as
macTEM's ring-calibration reference) and NumPy 2.5.4's float64 beside it, on V ∈ {60, 80, 120, 200, 300} kV, L_px ∈ {500,
2000, 14082.307097295722, 1e5}, R ∈ {1, 40, 150, 600, 2000} px and the inverse on the q those give plus q at 0.99 of the
Bragg limit; band 1e-12 relative to the decimal truth (the round trip radius(q(R)) within 1e-12 relative of R); the
q(R, L, λ) = 1/d identity on slice 3a's 320 cases to 1e-12.

## 4. G7: cancellation and progress (tier 2: no number moves)

`EllipseCalibration.fit1D`, `fitAmorphousRing`, `fitBestAvailable`, `RingCalibration.fit` and `fitRing` take
`progress: (@Sendable (Double) -> Void)? = nil` (monotone, 0 to 1) and stop with `CancellationError` when
`Task.isCancelled` or an `isCancelled: (@Sendable () -> Bool)? = nil` closure says so, checked once per optimiser
iteration and per Jacobian column. Gate: with both closures given (never cancelling) every output equals the plain call
bit for bit on the parity rings; a cancel at the first check throws `CancellationError` within one column. The profile
fit's Jacobian columns and normal matrix are computed in parallel (each sum in its sequential order, so bit identical,
held by the same parity gate); the measured speedup is reported, not gated.

## 5. G8: error texts (tier 2)

`CIFImportError`, `KinematicError` and `ReflectionListError` conform to `LocalizedError`, `errorDescription` equal to their
`description`. Gate: a test per error type that the localized description equals the description.

## 6. macTEM's switch: tetragonal along a or b (misfit d, tier 1)

`UnitCell.metricFamily` calls a cell tetragonal only when a = b ≠ c (90°); a cell with a = c ≠ b or b = c ≠ a reads
orthorhombic, so its 4-fold is never tested and a tetragonal structure in such a setting reads mmm. Fixed: the metric is
tetragonal when any two lengths are equal (angles 90°) with the third different, and the holohedry 4/mmm is built about the
odd axis. Pins: `tools/kit-0.0.4/tetragonal-pins.py`, pymatgen 2026.9.24's `SpacegroupAnalyzer` (spglib 2.8.0) point
group on three tetragonal structures (rutile TiO2, anatase TiO2, β-Sn) each in its three axis settings (the 4-fold along
a, b, c by permuting the cell and the coordinates) and on two controls (an orthorhombic cell with a = c by accident whose
atoms are not 4-fold, and a cubic one): the kit's Laue class and its operation count equal pymatgen's on all; truth: the
setting changes nothing (the three settings give the same class, the same count of operations and the same ring
multiplicities). Hexagonal settings with the 6-fold along a or b have the same shape of defect; they are recorded, not
fixed here (slice 1b's conventional cell removes both).

## 7. Tier 3 (docs and API hygiene)

(a) The placeholder enums `DiffractoGeometry`, `DiffractoScattering` and `DiffractoOrientation` that shadow their module
names are removed. (b) The overflow rule stays as Gate B on 3a fixed it (λ NaN above about 1.3e154 V; L_px and Δq NaN on
overflow), stated in the consumers note, whose item 15 is corrected. (c) The note states the -Onone cost in an app's
Debug configuration. (e) The note names the close-contact rule behind `invalidStructure`.

## 8. Amendments

### 2026-10-15, 11:30 CEST, after the ring-fit code started: what the first runs changed

Fail first, on the unmodified code (`RingFitAmendmentTests`, 12 tests): F1 (3 of 4 gapped noisy rings calibrated, the
fourth refused for another reason), F2 (the disk edge at blur 2.5 calibrated; all three textured circles read as
ellipses), F3 (5×, 20×, 1000× hot pixels refused), F6 (2^−40 refused as `insufficientSignal`; 2^−3 not bit equal) and F7
(σ_out 0.15 and 0.05 refused) red; the blob and macTEM's refuse-bound plant already as predicted (§1.7 held: calibrated,
Δq within 0.05 %).

- **§1.2, the normalisation.** "A range below 1 into [1, 2)" rescaled the parity rings (a unit ring's Float maximum lies
  just below 1) and moved their profile fits by up to 1.7e6 ulps. Now: a range below 0.5 is brought into [0.5, 1); every
  range of 0.5 or more runs the apps' arithmetic. Scale invariance is exact within the normalised regime (2^−40 equals
  2^−3 bit for bit, intensities scaled by 2^−37 exactly) and within 1e-6 relative across regimes (the gate's band; the
  Float pixels round s·v unlike v by up to 6e-8 relative, and the optimisers' step scales differ by less than 2 between
  the regimes; measured 4.5e-9 rad in θ). The "within 1e-9 relative" of §1.2 was set before that was known.
- **§1.3, the geometry fit's weights**: 1/(SE(t)² + 0.01²) (a noise-free sector has SE 0) and the covariance scaled by
  max(1, χ²/(n − 5)), conservative toward acceptance. `ringEvidence` takes the start centre too (the annulus is about it):
  `ringEvidence(pattern:fit:centerQX:centerQY:innerRadius:outerRadius:)`.
- **§1.5, the one-sided floor.** Kept as registered, the fit of a ring with an edge sharper than the floor read R 0.12 px
  (σ_out 0.15) and 0.24 px (σ_out 0.05) inside the plant: the 0.25 px floor cannot represent the edge and the peak moves
  to compensate. So a one-sided-floor fit whose evidence passes is refined once more with that floor lowered to 0.01 px
  (kept when its residual is not worse and its centre and axes pass the same guards), and its evidence must pass again.
  Measured: R within 0.1 px of the plant for both. The floor of every other fit stays 0.25 px (bit identity).
- **§1.7, the plant that must be refused.** σ 8 px in a 16 px annulus is not σ ≥ W; it is fitted and refused
  (`ringGeometryInconsistent`). The plant is now σ 16 px in a 16 px annulus: refused, measured as
  `ringGeometryInconsistent` too (the annulus shows ±0.5σ, so the fitted σ stays below 16 and `sigmaAtAnnulusWidth` does
  not fire); the test accepts any of `sigmaAtAnnulusWidth`, `notARingPeak`, `ringGeometryInconsistent`.
- **A finding, not a gap of the five:** at Gaussian noise 0.4 of the peak, a whole symmetric 3 px ring at R 42 is
  calibrated with R +0.78 % (the profile trades its radius against its asymmetry: σ_in 3.27, σ_out 2.76 fitted), the same
  before and after this amendment; the evidence cannot see it (the sectors agree with the fitted shape). The test's band
  is 1 % with the reason; a radius standard uncertainty is macTEM's ask G3 (its slice 2), named in the consumers note.
- **§4, cancellation** is checked once per optimiser iteration (both optimisers) and before each refinement, not per
  Jacobian column: the columns now run in parallel. Measured identity: the parity rings with both closures equal the
  apps bit for bit (`EllipseAppsParityTests`' `bestWithClosures`).

### 2026-10-15, about 12:45 CEST, after G2's code: the underflow rule

`q` and `radius` also return NaN when the result underflows to 0 (a positive R or q below about 1e-300 relative to L or
1/λ), not only on overflow: zero is not the q of a ring. The pins plant it. Measured against the 60-digit truth: q 8.2e-16,
radius 3.3e-15 on the grid and 4.5e-14 at 0.99 of the Bragg limit, round trip 5.7e-16, the 1/d identity on slice 3a's 320
cases 5.5e-16 (band 1e-12); the truth computed two ways (algebraic and series) agreeing to 1e-50. The inverse's last check is
`isFinite && != 0` so the Bragg guard is seen by its own test (its mutant stayed green with `> 0`). Mutants 6 of 6 red.
