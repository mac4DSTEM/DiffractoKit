# Status

What is live, one table; then the hand-off. Rules are in `CLAUDE.md`; the registrations in `docs/registration-*.md`.

## Live (2026-10-15 afternoon, end of the third session: the 0.1.0 candidate)

| Area | State | Evidence |
|---|---|---|
| Slice 1, `DiffractoCore` | Cells of the seven systems, symmetry from the atoms (tetragonal along a, b or c since 0.1.0), CIF import, reflections, absences, orbits, rings; Gate B closed | `registration-slice-1.md`; pins vs pymatgen on 19 structures + 7 cells; tetragonal pins (spglib, three settings); mutants 29/29 + 12/12 |
| Slice 2, `DiffractoScattering` | Electron scattering factors (py4DSTEM f050d207's table), F(hkl) in Å⁻², the \|F\| > 1e-4 cut, rings after it; Gate B closed | `registration-slice-2.md`; NumPy pins on 20 structures, pymatgen X-ray sums, the apps' compiled code bit for bit; mutants 17/17 |
| Slice 3a, `DiffractoGeometry` | CODATA-2018 wavelength, camera length and Δq from one ring, the exact q(R) and R(q) (0.1.0, G2), indexing for every crystal system, zone axes; Gate B closed | `registration-slice-3.md`, `registration-0.0.4.md` §3; 60-digit decimal pins (q 8.2e-16, radius 4.5e-14 at 0.99 of the Bragg limit; mpmath refuter 2.6e-16) |
| Slice 3b + the 0.1.0 ring-fit amendment | mac4DSTEM 7d6fdbd1's ring and ellipse fit, bit identical on every input without an impulse and with a range ≥ 0.5; impulses ignored, scale-free; the calibration path reads the ring evidence (coverage at the fit, a peak, one ring, the geometry) and refuses by name; a one-sided sharp edge kept with a lowered-floor refit; `fitRing` / `calibrating` (G1); progress and cancellation (G7); Gate B on f3933d2 closed (F1 second ring and F2 progress fixed, F3–F6 recorded) | `registration-0.0.4.md` §1, §2, §4, §8; parity 35 rings at -O and -Onone with and without closures; refuter: 181 ring plants, 836 G1 comparisons; mutants 18/18 |
| Error texts (G8) | `LocalizedError` on `CIFImportError`, `KinematicError`, `ReflectionListError` | `registration-0.0.4.md` §5 |
| Slice 1b | Ported and green (bc99fa0: Niggli reduction, conventional cell, both settings, Cartesian tolerance, cell-list detection; 27/27 mutants), then **reverted on `main` before the 0.1.0 tag** pending its Gate B refuter (owner: 1b is 0.1.x) | `registration-slice-1b.md` §6 |
| Tests | Swift Testing throughout | `tools/test.sh` (the board row names the last stamp and counts) |
| Consumer check | 404 checks through the public API, Debug, 0 warnings; 147/147 negative controls negative controls over every new 0.1.0 symbol and every type the apps use (`run.sh --controls`) | `tools/consumer-check/` |
| Parity | `docs/parity.md` against pymatgen, py4DSTEM and pyxem/diffsims/orix (held / partly / not held per row, each with its check) | |
| Docs and CI | A standalone README with the first ten minutes, `Examples/diffracto-example`, `.spi.yml`, `.github/workflows/ci.yml` (runs on the owner's first push) | |
| Consumers | `docs/consumers-0.1.0.md` (what each app changes to re-pin, the numbers that move, macTEM's two interims to delete); `docs/consumers-0.0.3.md` for the switch itself | |

Validation: none of this is validated on a real pattern; agreement with pymatgen, spglib, NumPy, SciPy, mpmath, py4DSTEM's table
and the apps' code, and recovery of planted truth, are what is claimed.

## Hand-off (0.1.x)

The suite board's `lanes/diffractokit` holds the next prompt. In order:

1. The owner pushes `main` and tags 0.1.0 (commands on the row); the CI file runs for the first time then.
2. Slice 1b back on `main` (`git revert` of the revert commit, or `git cherry-pick bc99fa0`), its one Gate B refuter
   (registration-slice-1b §4), the consumer-check extension prepared for it (patch named in the row), parity rows, and
   `docs/consumers-0.1.1.md`.
3. CI type-check fixes if the runner reports any (a local scan found no expression above 200 ms).
4. Slice 4, vector matching, registered on its own (the owner may overrule its move on the row).
5. 0.2.0, the generality pass (owner): suite policies as parameters, no app or suite names in public docs.

## Unverified / not done

- `testRinglessFitOnARealMeanPatternFallsBackToTheConic` is disabled by name: the owner's Ni65Cu35 mean pattern has no
  `../data/README.md` manifest line yet.
- No database CIF (COD, ICSD) read yet; SKCIFParser and Gemmi are named as parsing comparisons, not run.
- The 0.1.0 ring-fit evidence is measured on synthetic plants only; macTEM's real Au pattern rings are re-measured at its re-pin.
- The parallel Jacobian was withdrawn (Swift 6 warnings); the fit's speed is 0.0.3's.
