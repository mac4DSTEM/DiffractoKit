# DiffractoKit 0.1.0: what mac4DSTEM and macTEM change to re-pin

Written 2026-10-15 by DiffractoKit's third session, for mac4DSTEM's and macTEM's sessions. 0.1.0 is the kit loop's first
tag (prepared as 0.0.4; the owner made it 0.1.0 at 13:16, and the apps re-pin once, at 0.1.0, skipping 0.0.4): macTEM's Diffraction-room asks (G1, G2, G7, G8), the five misfits of macTEM's 0.0.3 switch, and the five ring-fit gaps
slice 3b's refuter found in the apps' fit, fixed in the kit. Registered in `docs/registration-0.0.4.md` (with its dated
amendments). The tag is the owner's; until it exists nothing here can be pinned. Read `docs/consumers-0.0.3.md` first: this
note lists only what changed since 0.0.3.

## In one table

| Who | Deletes or changes | Why |
|---|---|---|
| macTEM | its **G1 interim** (the four public calls composed in TEMCore, pinned equal to `RingCalibration.fit`) | `RingCalibration.fitRing(pattern:centerQX:centerQY:innerRadius:outerRadius:)` and `calibrating(_:d:highTensionVolts:pixelPitchMM:)` |
| macTEM | its **G2 interim** (exact q and its inverse in TEMCore, registered R3/R4) | `RingCalibration.q(radius:cameraLengthPixels:wavelengthAngstroms:)` and `radius(q:cameraLengthPixels:wavelengthAngstroms:)`, pinned against 60-digit decimal and NumPy |
| macTEM | its `LocalizedError` wrappers of the kit's errors (`CrystalLoading`) | `CIFImportError`, `KinematicError`, `ReflectionListError` are `LocalizedError` |
| macTEM | "cannot cancel" in the room (Remove drops the answer) | `progress:` and `isCancelled:` on every fit, and `Task` cancellation, throw `CancellationError` |
| both | any `switch` over `RingCalibration.Refusal` | four new cases (below); a `switch` that lists every case no longer compiles |
| both | any reference to the enums `DiffractoGeometry`, `DiffractoScattering`, `DiffractoOrientation` | removed: they shadowed their modules |
| both | any `.tetragonal` pattern or comparison | `UnitCell.MetricFamily.tetragonal(uniqueAxis:)` (below) |
| mac4DSTEM | nothing is required; the ring-fit fixes reach it through the pin (below), and `requireRingEvidence: true` is its own lane's decision | |

## What numbers move

1. **The calibration path refuses more (macTEM's `RingCalibration.fit`, and `fitRing`).** Since 0.1.0 the calibration runs
   the fit with `requireRingEvidence: true` and reads the data sector by sector about the fitted ellipse
   (`EllipseCalibration.ringEvidence`, `RingEvidence`). Five refusals, each withholding Δq:
   - `insufficientAngularCoverage(occupiedBins:)` now also when fewer than 30 of 36 sectors hold the ring AT THE FIT
     (slice 3b's F1: noise in a 90°–120° gap made the gap's bins "occupied"; such rings were calibrated with the centre
     0.7–2.7 px and R −0.9 to −3.5 % off). The payload is the smaller of the two counts.
   - `ringOutsideAnnulus(radius:inner:outer:)` (new): (a + b)/2 not inside the annulus (a central blob).
   - `notARingPeak(innerDrop:outerDrop:)` (new): the profile about the fit does not fall by a quarter of the ring's
     amplitude on both sides (a bright-disk edge, which 0.0.3 calibrated at −10 to −18 %).
   - `secondRingInAnnulus(height:)` (new, after Gate B): the pooled profile holds a second peak of a quarter of the
     ring's amplitude inside the annulus (a neighbouring ring pulled 13 of 14 such fits onto itself or between the two).
   - `ringGeometryInconsistent(centreShift:ellipseShift:)` (new): the sectors put the ring more than 0.25 px and 4
     standard errors off the fitted centre or ellipse (a textured circle read as a 4.6 % ellipse; F1's centre bias).
   Order: model, coverage, sigma, residual, ellipticity, outside the annulus, not a peak, second ring, geometry, high
   tension.
   Measured on the kit's sets: the 35 parity rings, macTEM's eight `RingFitCalibrationTests` plants, clean and noisy
   whole rings (noise to 0.4 of the peak) and macTEM's `refuse-bound` plant are still calibrated; the F1 and F2 plants are
   refused. **macTEM re-runs its `ring-fit-parity` and S3-b harnesses at the re-pin:** a ring 0.0.3 calibrated may now be
   refused by name; each such row is a registered change, not a regression to wave through.
2. **One ring sigma at the 0.30 px floor is no longer refused on the calibration path** when the ring evidence passes:
   the fit is refined once more with that side's floor at 0.01 px (slice 3b's F7). A ring with one edge sharper than
   0.25 px is calibrated with R within 0.1 px (0.0.3 refused it as `ringWidthAtFloor`). Both sigmas at the floor is still
   refused (ADR 069's dead ring). Off the calibration path (`requireRingEvidence: false`, the default) nothing changes.
3. **Isolated impulses are ignored by every fit** (`fit1D`, `fitAmorphousRing`, `fitBestAvailable`, both apps): a pixel
   outside the range of its 8 neighbours by more than twice the annulus's robust range (P99.5 − P0.5) is skipped like a
   NaN. A ring with one hot pixel of 5× to 1000× its peak is fitted (0.0.3 refused it as `insufficientAngularCoverage(1)`)
   and a strongly negative pixel no longer moves it. A pattern without such a pixel fits bit for bit as in 0.0.3 (the
   parity gate). **Real data with hot pixels will fit differently** (better); a golden recorded on such a pattern moves.
4. **Scale.** The signal test is relative (0.0.3 refused a ring of amplitude ≤ 1e-6 as `insufficientSignal`), and a
   dynamic range below 0.5 runs normalised by an exact power of two: a pattern in units where the ring is tiny fits as
   the same pattern scaled up. Ranges of 0.5 and more fit bit for bit as in 0.0.3.
5. **Not refused, and stated:** two rings closer than about twice their width are one broadened ring to the fit (its R lies
   between them); at Gaussian noise 0.4 of the ring's peak a whole 3 px ring at R 42 px reads R +0.78 % (the profile
   trades radius against asymmetry), the same as in 0.0.3. A radius standard uncertainty is macTEM's ask G3, not in 0.1.0.

Nothing else moves: the closed forms, the indexing, the crystal, CIF and scattering numbers are 0.0.3's bit for bit, and
`EllipseCalibration`'s default path is bit identical to mac4DSTEM `7d6fdbd1` on the 35 parity rings at -O and -Onone, with
and without the new closures.

## Gate B on the frozen candidate (f3933d2)

The candidate was frozen at f3933d2. Its refuter (registration-0.0.4 §8, "Gate B on f3933d2") found one wrong acceptance
and one wrong progress report, both fixed before the tag (67f8d22) and in 0.1.0:
a second ring inside the annulus could pull the calibration onto the other ring or between the two (now refused as
`secondRingInAnnulus(height:)`, listed after `notARingPeak`), and progress could step back through the one-sided refit. Also recorded there: a ring on
a steep 1/r² tail or a composite background is refused (conservative); an exponential outer tail reads R up to −2.1 % and
a flat-top ring of width 14 px reads R at its outer edge (accepted: the evidence reads the data at the fitted shape); the
impulse rule acts on a one-pixel spike anywhere in the pattern (a one-pixel central beam moves a fit by ≤ 1.1e-4 px) and
on isolated high counts of a sparse Poisson ring (R by ≤ 0.005 px).

## What is new (additive)

- `RingCalibration.fitRing(...) -> RingRadiusFit` (`fit`, `report`, `refusals`, `radius`, `isUsable`) and
  `RingCalibration.calibrating(_:d:highTensionVolts:pixelPitchMM:)`; `fit` is exactly `calibrating(fitRing(...))`, field for
  field (a test on every plant, with and without V).
- `RingCalibration.q(radius:cameraLengthPixels:wavelengthAngstroms:)` = 2 sin(½ atan(R/L))/λ (Å⁻¹, no 2π) and
  `radius(q:cameraLengthPixels:wavelengthAngstroms:)` = L tan(2 asin(qλ/2)); NaN outside the domain (inverse: qλ/2 below
  1/√2). q of the calibrating ring's own (R, L_px, λ) is 1/d. Pinned against 60-digit decimal truth (two
  independent forms) and NumPy: worst 8.2e-16 (q), 4.5e-14 (radius, at 0.99 of the Bragg limit); NaN on overflow and on
  underflow to 0.
- `EllipseCalibration.ringEvidence(pattern:fit:centerQX:centerQY:innerRadius:outerRadius:) -> RingEvidence` (per-sector
  amplitudes and offsets, `occupiedBins`, `innerDrop`/`outerDrop`, `centreShift`/`ellipseShift` with errors, `isPeak`,
  `geometryHolds`, `passes`, the registered bounds as static constants); `EllipseCalibrationFit.ringEvidence`;
  `fitBestAvailable(..., requireRingEvidence:)`, `fitAmorphousRing(..., requireRingEvidence:)`;
  `FitReport(annulusInnerRadius:annulusOuterRadius:ringEvidence:)` (defaulted; the 0.0.3 initializer still compiles).
- `progress: (@Sendable (Double) -> Void)?` and `isCancelled: (@Sendable () -> Bool)?` on `fit1D`, `fitAmorphousRing`,
  `fitBestAvailable`, `RingCalibration.fit` and `fitRing`; a cancelled `Task` or a `true` from the closure throws
  `CancellationError` at the next optimiser iteration. No speedup in 0.1.0: the fit's cost is 0.0.3's (macTEM's
  66 s at ±90 px on 2048² stands); a parallel Jacobian is a later lane.
- `LocalizedError` on `CIFImportError`, `KinematicError`, `ReflectionListError` (`errorDescription` is the description).
- **Tetragonal along a or b (a number moves).** A cell with a = c ≠ b or b = c ≠ a (90°) was read orthorhombic, so its
  4-fold was never tested and a tetragonal structure in such a setting read mmm with half its operations and split rings.
  `UnitCell.MetricFamily.tetragonal` is now `.tetragonal(uniqueAxis:)` (as `.monoclinic`), and 4/mmm is built about the odd
  axis. Pinned against pymatgen 2026.9.24 / spglib 2.8.0 (symprec 1e-4, 1e-3 and 1e-2 agree): rutile (136, 16 point
  operations), anatase and β-Sn (141, 16 point and 32 space-group operations) in all three settings, and two controls
  (an orthorhombic cell with a = c by accident, mmm; SrTiO3, m-3m); the three settings give the same class, count and ring
  multiplicities. **Source change:** a `case .tetragonal:` in an app's `switch` becomes `case .tetragonal:` with a pattern
  (`case .tetragonal(let axis)` or `case .tetragonal(_)`); a comparison `== .tetragonal` becomes
  `== .tetragonal(uniqueAxis: .c)` where c-unique is meant. The same defect with the 6-fold along a or b is recorded, not
  fixed (slice 1b's conventional cell removes both). What moves in an app: only a CIF or cell in such a setting, which
  now reads 4/mmm (16 operations) and merges its rings as spglib does.

## The consumer check

`tools/consumer-check` (a separate package, the public API only, Debug, 0 warnings) passes 404 checks (288 at 0.0.3),
each new 0.1.0 symbol with a real call and a checked result. Its negative controls (`run.sh --controls`,
`controls.py`, SpectroKit's pattern) revert one declaration at a time from `public` to `package` in a copy and require the
consumer build to fail with the protection-level diagnostic naming it: 147 of 147 caught, covering every new 0.1.0 symbol
and every type (with at least one member) that the 0.0.3 note's mapping tables name. Eight declarations cannot be made
`package` and are listed with the reason in `controls.py`: the four new `Refusal` cases, `MetricFamily.tetragonal`, and
the three `LocalizedError` conformances (their `errorDescription` witnesses are controlled).

## Corrections to the 0.0.3 note

- **Item 11, `invalidStructure`, names its rule:** besides a NaN or infinite form factor, occupancy or coordinate,
  `CrystalStructure.validationIssues` refuses two sites closer than `CrystalStructure.closeContactLimit` = 0.5 Å
  (minimum image over the 27 neighbouring cells; below the shortest contact in any crystal, H–H in solid hydrogen at
  0.74 Å) unless their occupancies sum to at most 1.02 (a split site, the apps' rule); code `atomic_sites_too_close`,
  read through `shortestCloseContact`. A CIF with such a pair imports as `invalidStructure`.
- **The Laue class to ACOM family mapping needs the metric too** (mac4DSTEM's switch, relayed 12:1x): the 0.0.3 table
  mapped m-3m to cubic, 6/mmm to hexagonal and 2/m (b unique) to monoclinic, with the metric condition only for 2/m. Every
  row needs it: an fcc structure written in its primitive 60° rhombohedral cell has Laue class m-3m, but its Cartesian
  frame is not the cubic axes, so a fundamental zone drawn for cubic is wrong. The mapping is m-3m with
  `structure.cell.metricFamily() == .cubic`, 6/mmm with `== .hexagonal`, 2/m with `== .monoclinic(uniqueAxis: .b)`,
  everything else `.identity`. `UnitCell.metricFamily(lengthTolerance:angleToleranceDeg:)` is public (defaults 1e-9 and
  1e-9°); mac4DSTEM's `ACOMCrystalSymmetry.reduction(for:)` at 46759af6 does exactly this. (Slice 1b's conventional cell will
  let a consumer reduce such a cell first.)
- **Two outcomes of the 0.0.3 switch, by design** (mac4DSTEM registered them as unpredicted): the kit refuses split-atom
  gold (two whole atoms 0.016 Å apart) as `invalidStructure` by the close-contact rule below; and a polar P6_3mc structure
  (SiGe, mp-978534), which the apps' old 6/m check refused, now reads Laue 6/mmm and imports as hexagonal: the Laue class
  is the point group with inversion added, which is what a kinematic pattern shows under Friedel's law (a dynamical or
  anomalous-scattering consumer must not read it as the point group).
- **Item 15, the overflow rule:** the closed forms are macTEM's except at overflow, where the kit is total and NaN
  (slice 3a's Gate B, F1 and F2): λ is NaN above about 1.3e154 V (macTEM's gives 0), and L_px and Δq are NaN where
  the quotient overflows (R 1e300, a subnormal d), where macTEM's give +∞. No physical input reaches either; a consumer
  that tested for 0 or ∞ tests `isFinite`. The kit keeps NaN: a quantity that cannot be computed is never a plausible
  number, and NaN is the rule every closed form here follows.
- **The ring fit's cost in an app's Debug build:** Xcode builds a package dependency at -Onone in an app's Debug
  configuration, whatever the kit's own tests use. One `RingCalibration.fit` measured 2.0 s at -Onone against 0.04 s at
  -O (macTEM's probe, on its own ring). An app that runs fits in Debug sets the package's optimisation for that
  configuration (`SWIFT_OPTIMIZATION_LEVEL` per target through an xcconfig, or a Release-built scheme for fit work); the
  kit's manifest can carry no unsafe flags.

## mac4DSTEM

Re-pins for the ring-fit fixes (items 3 and 4 above reach its detector-ellipse path; items 1 and 2 only if it opts in).
Its pinned consumer tests (`EllipseNoisyRingTests`, `EllipseStopRuleTests`, `CalibrationEllipseRow29Tests`) are among the
kit's parity rings and stay bit identical. What is its own decision, in its ring-fit lane: whether the detector-ellipse
path passes `requireRingEvidence: true` and refuses on `fit.ringEvidence?.passes == false`. The kit did not make that the
default because mac4DSTEM refuted a radius-only annulus guard on its WS2 data (its `ring-fit-annulus-guard-gateD-2026-10-14`:
rings fitted just past the outer radius with a correct b/a); the evidence's coverage and geometry checks speak to b/a, the
annulus and peak checks to the radius, and its data decide.

## macTEM

Re-pins and deletes its two interims (G1 and G2) in the same commit, each replaced by the kit call it stood in for; its
pinned tests of the interims become tests of the kit calls with the same numbers. Then the room: fit with `fitRing` (no d
needed for ring k), name the family, `calibrating`; q and the predicted radii from `q` and `radius`; the fit job's Cancel
through `Task` cancellation or `isCancelled`; the error texts from `localizedDescription`. Its `ring-fit-parity` and S3-b
numbers are re-measured (item 1); its `refuse-bound` plant expected the wrong refusal (registration-0.0.4 §1.7: a 3 px ring
in a 40 px annulus cannot meet `sigmaAtAnnulusWidth`; the kit calibrates it within 0.05 %), so that plant's expectation is
corrected in its harness, and a ring with σ at the annulus width is refused.
