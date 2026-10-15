# Status

What is live, one table; then the hand-off. Rules are in `CLAUDE.md`; the registrations in `docs/registration-*.md`.

## Live (2026-10-15 15:30 CEST: the 0.1.1 tag candidate and the 0.2.0 candidate)

| Area | State | Evidence |
|---|---|---|
| Slice 1, `DiffractoCore` | Cells of the seven systems, symmetry from the atoms (tetragonal along a, b or c since 0.1.0), CIF import, reflections, absences, orbits, rings; Gate B closed | `registration-slice-1.md`; pins vs pymatgen on 19 structures + 7 cells; tetragonal pins (spglib, three settings); mutants 29/29 + 12/12 |
| Slice 2, `DiffractoScattering` | Electron scattering factors (py4DSTEM f050d207's table), F(hkl) in Å⁻², the \|F\| > 1e-4 cut, rings after it; Gate B closed | `registration-slice-2.md`; NumPy pins on 20 structures, pymatgen X-ray sums, the apps' compiled code bit for bit; mutants 17/17 |
| Slice 3a, `DiffractoGeometry` | CODATA-2018 wavelength, camera length and Δq from one ring, the exact q(R) and R(q) (0.1.0, G2), indexing for every crystal system, zone axes; Gate B closed | `registration-slice-3.md`, `registration-0.0.4.md` §3; 60-digit decimal pins (q 8.2e-16, radius 4.5e-14 at 0.99 of the Bragg limit; mpmath refuter 2.6e-16) |
| Slice 3b + the 0.1.0 ring-fit amendment | mac4DSTEM 7d6fdbd1's ring and ellipse fit, bit identical on every input without an impulse and with a range ≥ 0.5; impulses ignored, scale-free; the calibration path reads the ring evidence (coverage at the fit, a peak, one ring, the geometry) and refuses by name; a one-sided sharp edge kept with a lowered-floor refit; `fitRing` / `calibrating` (G1); progress and cancellation (G7); Gate B on f3933d2 closed (F1 second ring and F2 progress fixed, F3–F6 recorded) | `registration-0.0.4.md` §1, §2, §4, §8; parity 35 rings at -O and -Onone with and without closures; refuter: 181 ring plants, 836 G1 comparisons; mutants 18/18 |
| Error texts (G8) | `LocalizedError` on `CIFImportError`, `KinematicError`, `ReflectionListError` | `registration-0.0.4.md` §5 |
| Slice 1b | On `main` for 0.1.1: Niggli reduction, conventional cell, both settings, Cartesian tolerance, cell-list detection; Gate B (one refuter) on 4876b18: no wrong number or class on realistic inputs; F1 (multi-unit shears, a non-converged reduction flagged by `reductionConverged`) fixed in 6c16c9a, F2-F4 documented. **0.1.1 tag candidate: 4876b18** (the owner tags) | `registration-slice-1b.md` §6; `consumers-0.1.1.md`; mutants 27/27 + 3 new red; suite 306/306 at 6c16c9a |
| 0.2.0, the generality pass | Candidate: every suite policy a parameter with the old value as default (`Policy` types in Core, Scattering and Geometry; default path bit identical), no app or suite name in Sources, standalone README/example/`tools/readme-consumer`/DocC, the floor macOS 15 and iOS 18 (302/302 and consumer-check 537/537 there, iOS build ok), CI type-check fixes (nothing at or above 40 ms). Policies are staged and uncommitted; the freeze fills the gates (suite count, mutants 28/28 via `tools/kit-0.2.0/mutants-policy.py`, consumer check, refuter). Tag commit: to be filled | `registration-0.2.0.md` §4; `consumers-0.2.0.md` |
| Tests | Swift Testing throughout | `tools/test.sh` (the board row names the last stamp and counts) |
| Consumer check | 404 checks through the public API, Debug, 0 warnings; 147/147 negative controls over every new 0.1.0 symbol and every type the apps use (`run.sh --controls`) | `tools/consumer-check/` |
| Parity | `docs/parity.md` against pymatgen, py4DSTEM and pyxem/diffsims/orix (held / partly / not held per row, each with its check) | |
| Docs and CI | A standalone README with the first ten minutes, `Examples/diffracto-example`, `.spi.yml`, `.github/workflows/ci.yml` (runs on the owner's first push) | |
| Consumers | `docs/consumers-0.2.0.md` and `consumers-0.1.1.md` (since 0.1.0); `docs/consumers-0.1.0.md` (what each app changes to re-pin, the numbers that move, macTEM's two interims to delete); `docs/consumers-0.0.3.md` for the switch itself | |

Validation: none of this is validated on a real pattern; agreement with pymatgen, spglib, NumPy, SciPy, mpmath, py4DSTEM's table
and the apps' code, and recovery of planted truth, are what is claimed.

## Hand-off (0.2.x)

The suite board's `lanes/diffractokit` holds the next prompt. In order:

1. The owner pushes `main` and tags 0.1.1 at 4876b18 (slice 1b with its consumer check and the CI type-check fix; 0.1.0 is
   tagged before it); the CI file runs on the push.
2. 0.2.0: commit the policy change (staged, tier 1), run the freeze (suite count, mutants 28/28, consumer check with Policy
   checks, one refuter), fill `registration-0.2.0.md` §4 and `consumers-0.2.0.md` (tag commit), then the owner tags 0.2.0.
3. **Slice 4 next: vector matching and phase mapping move into DiffractoKit as 0.3.0** (the owner's decision of 2026-10-15
   15:20); `DiffractoOrientation` gets its content. Registered on its own before any code.
4. Policies not yet parameters (profile-fit plausibility bounds, cell volume 1e-6, Niggli epsilon), if an app asks.

References is a link to `../References/DiffractoKit` (the move happens at this closeout; scripts read the same
`References/` paths).

## Unverified / not done

- `testRinglessFitOnARealMeanPatternFallsBackToTheConic` is disabled by name: the owner's Ni65Cu35 mean pattern has no
  `../data/README.md` manifest line yet.
- No database CIF (COD, ICSD) read yet; SKCIFParser and Gemmi are named as parsing comparisons, not run.
- The 0.1.0 ring-fit evidence is measured on synthetic plants only; macTEM's real Au pattern rings are re-measured at its re-pin.
- The parallel Jacobian was withdrawn (Swift 6 warnings); the fit's speed is 0.0.3's.
