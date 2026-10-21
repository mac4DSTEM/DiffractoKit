# tools/kit-0.0.4: pins and break-once runs of DiffractoKit 0.0.4

Registered in the design note (archived in the suite's private folder). Each script works in a copy of the repository, never the working tree.

| Script | What | Result (2026-10-15) |
|---|---|---|
| `q-pins.py` | G2: q and its inverse, 60-digit decimal truth computed two ways (algebraic and series, agreeing to 1e-50) beside NumPy 2.5.4; writes `Tests/DiffractoGeometryTests/Fixtures/q-pins.json` (byte-identical on rerun) | 100 forward, 120 inverse, 320 identity, 51 NaN plants |
| `tetragonal-pins.py` | §6: pymatgen 2026.9.24 / spglib 2.8.0 (symprec 1e-4) on rutile, anatase, β-Sn in three settings and two controls; writes `Tests/DiffractoCoreTests/Fixtures/tetragonal-pins.json` | the kit equals spglib in every setting |
| `mutants-ringfit.py` | §1, §2, §4 and Gate B's F1, F2: the ring-fit amendment, the G1 split, G7, the second ring, monotone progress | 18 of 18 red (first pass 15 of 17, with a parallel-order mutant since removed with the parallel path: the occupancy fraction and the annulus refusal had no test that saw them alone; `aRingWithAWeakSectorIsNotCovered` and `theAnnulusRefusalReadsTheReportsRadii` added, both mutants then red) |
| `mutants-lane-a.py` | §5, §6: tetragonal and the error texts | 12 of 12 red |
| `mutants-lane-b.py` | §3: q and radius | 6 of 6 red (the Bragg-guard mutant was green until the inverse's last check became `isFinite && != 0`) |

Run, for example: `python3 tools/kit-0.0.4/mutants-ringfit.py <copy-dir> <build-dir> [mutant-id ...]` (release builds).
