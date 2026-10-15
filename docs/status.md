# Status

What is live, one table; then the hand-off. Rules are in `CLAUDE.md`; the registrations in `docs/registration-*.md`.

## Live (2026-10-15 morning, end of the second session)

| Area | State | Evidence |
|---|---|---|
| Slice 1, `DiffractoCore` | Cells of the seven systems, symmetry from the atoms, CIF import, reflections, absences, orbits, rings; Gate B closed (four defects fixed) | `registration-slice-1.md`; pins vs pymatgen on 19 structures + 7 cells; mutants 29/29 (`docs/slice-1-mutants-2026-10-14.md`) |
| Slice 2, `DiffractoScattering` | Electron scattering factors (py4DSTEM f050d207's Lobato table), F(hkl) in Å⁻², the \|F\| > 1e-4 cut, rings after the cut; Gate B closed (no wrong number; guards and labels fixed) | `registration-slice-2.md`; NumPy pins on 20 structures, pymatgen X-ray sums, the apps' compiled code bit for bit; β″ 417 rings / 1442 reflections; mutants 17/17 |
| Slice 3a, `DiffractoGeometry` | CODATA-2018 wavelength, camera length and Δq from one ring with refusals, the indexing decision for every crystal system, zone axes and the excitation error; Gate B closed (no number moved; guards, an overload, two wording errors) | `registration-slice-3.md`; macTEM's 60-digit pins (320 cases), its T1 tables reproduced exactly; lane mutants 10/10 and 19/19; Gate B mutants in `tools/slice-3/mutants.py` |
| Slice 3b, `DiffractoGeometry` | mac4DSTEM 7d6fdbd1's ring and ellipse fit (ADR 069), the transform, `RingCalibration.fit`; Gate B: no port defect, two traps fixed, five ring-fit gaps recorded for mac4DSTEM | bit identical to the apps' compiled code (35 parity rings at -O and -Onone; the refuter's ~940 more); SciPy sweep; lane mutants 9/9 |
| Slice 1b | Registered only: Niggli reduction and the conventional cell, a Cartesian tolerance option, fast detection | `registration-slice-1b.md` |
| Tests | Swift Testing throughout (slice 1 converted, its 75 names kept) | `tools/test.sh` (the board row names the last stamp and counts) |
| Gate | `tools/test.sh` stamps only on a clean tree with no untracked file in Sources/Tests (CuboKit's rule) | `.build/last-green` |
| Consumer check | 288 checks through the public API of all three modules, Debug, 0 warnings (no negative controls yet) | `tools/consumer-check/run.sh` |
| Consumers | `docs/consumers-0.0.3.md`: what mac4DSTEM and macTEM change to import slices 1–3, the numbers that move, the ring-fit findings they inherit | |

Validation: none of this is validated on a real pattern; agreement with pymatgen, NumPy, SciPy, py4DSTEM's table and the
apps' code, and recovery of planted truth, are what is claimed.

## Hand-off

The suite board's `lanes/diffractokit` holds the next prompt. In order:

1. The owner tags 0.0.3 (or names the version) after reviewing `docs/consumers-0.0.3.md`; the apps' sessions switch, each
   as a reviewed tier-1 lane with its own registration where a number moves.
2. The consumer-check's negative controls (SpectroKit's `controls.py` pattern) before the first tag the apps pin.
3. Slice 1b (registered): reduction to the conventional cell, the Cartesian tolerance, cell-list detection.
4. Slice 4, vector matching (mac4DSTEM's `PhaseVectorMatching` and its reference library), registered on its own.
5. mac4DSTEM's ring-fit lane: the five upstream findings of slice 3b's Gate B (registration-slice-3 §7).

## Unverified / not done

- `testRinglessFitOnARealMeanPatternFallsBackToTheConic` is disabled by name: the owner's Ni65Cu35 mean pattern has no
  `../data/README.md` manifest line yet.
- No database CIF (COD, ICSD) read yet; SKCIFParser and Gemmi are named as parsing comparisons, not run.
- The Debug run of the 3b parity tests takes about 10 minutes; `tools/test.sh` (release) runs them against the -O fixture.
