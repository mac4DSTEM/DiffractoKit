# Status

What is live, one table; then the hand-off. Rules are in `CLAUDE.md` and the suite's `../CLAUDE.md` ("How we work": the
standard process since 2026-10-15 17:25: tests that fail before and pass after, one review per change, pinned regression
fixtures, `CHANGELOG.md`, CI as the gate). The registrations `docs/registration-*.md` are the record of the first build-out.

## Live (0.5.0 pending release, 2026-10-20; last tag 0.4.1)

| Area | State | Evidence |
|---|---|---|
| Slice 1, `DiffractoCore` | Cells of the seven systems, symmetry from the atoms, CIF import, reflections, absences, orbits, rings | pins vs pymatgen and spglib (`docs/parity.md` P1-P3, P12-P13) |
| Slice 1b | Niggli reduction, conventional cell, both settings, Cartesian tolerance, cell-list detection; tagged in 0.1.1 (4876b18) | P12, P13 |
| Slice 2, `DiffractoScattering` | Electron scattering factors (py4DSTEM f050d207's table), F(hkl), the cut, rings; `retainingAll: false` keeps memory to the kept list (0.2.0) | P4-P7 |
| Slice 3, `DiffractoGeometry` | Wavelength, q(R) and R(q), camera length, indexing, zone axes; the ring and ellipse fit with its evidence and refusals | P8-P11 |
| 0.2.0, the generality pass | Every policy a parameter with the old value as default; floor macOS 15 / iOS 18; standalone README, example, DocC; the 0.2.0 review's 8 findings fixed (traps, the tolerance domain, absurd cells, NaN tolerances) and mac4DSTEM's 3 asks (memory of the kinematic list: 820 MB → 279 MB on a 30 Å cell at k_max 2.5; the refusal text; the polarity note) | `CHANGELOG.md` 0.2.0; default path bit identical to 0.1.1 on 31,725 outputs |
| Slice 4, `DiffractoOrientation` (0.3.0) | Vector-matched phase mapping (`PhaseReferenceLibrary`, `PhaseVectorMatcher`, 24 settings as documented defaults) and the orientation-result value types (`EulerAngles`, `ZoneAxisFrame`, `OrientationResult`), ported from mac4DSTEM 7d6fdbd1; two fixes against the source (a rotation lost at R22 = 1 − 2⁻⁵³; phase ties by hash seed) | `CHANGELOG.md` 0.3.0; P14 py4DSTEM, P15 SciPy and NumPy, P16 the source bit for bit (about 1.68 million values) |
| Orientation symmetry (0.4.0) | `CrystalOrientationSymmetry` supports all 11 proper Laue subgroups, minimum-angle rotation reduction and direction sampling; explicit `.py4DSTEM` / `.native` IPF keys, no identity fallback | P17–P19; `docs/orientation-symmetry.md`; py4DSTEM direction-family reduction not ported (returns families, not a rotation) |
| ACOM matcher (0.4.0) | `ACOMTemplateLibrary` and `ACOMMatcher` in `DiffractoOrientation`: the apps' simplified polar-correlation route (not py4DSTEM's SO(3) correlogram), with `ACOMExcitationModel.py4DSTEM` for py4DSTEM-convention patterns; ImagoKit 0.2.x (`.upToNextMinor(from: "0.2.0")`) for `FFT1D`. Pinned to a NumPy transcription of the app source and to py4DSTEM's spots and planted orientations. Known: with py4DSTEM's patterns 17 of 18 planted orientations are found within 4 degrees; one 2/m pattern is matched to another beam direction | P25; `tools/acom-matcher/` |
| ACOM scan matching (0.5.0) | `ACOMMatcher.matchAll` with `ACOMScanSelection`, `ACOMScanResults`, cancellation and progress; 1.45 ms a pattern on one worker, 0.24 ms on eight | `ACOMScanTests`; `tools/acom-benchmark/run.sh` |
| Tests | Swift Testing; `swift test -c release`: 550 (Core 158, Scattering 38, Geometry 209, Orientation 145; 3 existing known issues and 1 existing data skip) | `tools/test.sh` |
| macTEM pieces (0.4.0) | `CrystalLoading`, `RingOverlay`/`RingEllipse`, `IndexingDisplay`, `RingScale`; `CellStandardisation` (Materials Project import) | P20-P24 |
| Consumer check | 666 checks through the public API (Debug, 0 warnings); 240 / 240 negative controls (nightly material) | `tools/consumer-check/` |
| Parity | py4DSTEM 3 held / 16 partly / 9 not held; pymatgen 13 / 16 / 22; pyxem / diffsims / orix 0 / 9 / 4 | `docs/parity.md` |
| Hook | No sleep or poll loops; no machine-wide kills (killall, pkill and kill by pattern refused unless pinned to `/claude-501/` or a session UUID); proof 127 / 127 | `tools/hooks/` |
| CI | The first runs failed on the runner's Xcode 26.6 compiler crash (catch clauses naming an error case); fixed in e3516e3, unconfirmed until the owner's next push | `.github/workflows/ci.yml` |

Validation: none of this is validated on a real pattern; agreement with pymatgen, spglib, NumPy, SciPy, mpmath, py4DSTEM and the
apps' code, and recovery of planted truth, are what is claimed.

## Hand-off

**One line (2026-10-20):** done: ACOM scan matching committed for 0.5.0 (tests 550 green, consumer-check 666, benchmark run), release-ready, awaiting tag; next: mac4DSTEM switches its CPU `matchAll` and `ACOMScanSelection` onto 0.5.0, then the open kit gaps (`.native` IPF key beyond m-3m, 6/mmm, 2/m; orix/MTEX cross-check of the minimum-angle reduction; py4DSTEM direction-family reduction).

The suite board's `lanes/diffractokit` holds the next prompt. In order:

1. The release shape (asked of the owner through the coordinator, 17:55): (A) no separate 0.2.0, tag `main` as 0.3.0 after
   slice 4's review; or (B) 0.2.0 on a release commit off `main`, then 0.3.0. The tag commands are on the row.
2. Slice 4's one review (running at this closeout; its findings land as focused commits with tests).
3. Slice 5: template matching (`OrientationPlan`, `OrientationMatcher`) remains. Symmetry and IPF primitives are implemented (P17–P19); matcher integration, orientation-family reduction and further orix comparison remain separate work.
4. The lean pass: `docs/consumers-*.md` (replaced by `CHANGELOG.md`), the registrations and `docs/slice-1-mutants-2026-10-14.md`
   move to `../_archive/DiffractoKit/docs/` once the comments and docs that cite them are rewritten (`git grep` proof);
   `tools/kit-0.2.0/mutants-policy.py` and the slice mutants scripts retire with them (mutation runs are retired).
5. Policies not yet parameters (profile-fit plausibility bounds, cell volume 1e-6), if an app asks.

## Unverified / not done

- CI on GitHub after e3516e3 (the runner's compiler).
- `testRinglessFitOnARealMeanPatternFallsBackToTheConic` is disabled by name: the owner's Ni65Cu35 mean pattern has no
  `../data/README.md` manifest line yet.
- No database CIF (COD, ICSD) read yet.
- The 0.1.0 ring-fit evidence is measured on synthetic plants only.
- θ′-Al2Cu in the slice 4 pins is a test structure (Al4Cu2 body-centred, a 4.04 Å, c 5.80 Å), not checked against a
  published CIF.
