# DiffractoKit benchmarks against the Python counterparts

Measured 2026-10-21. Small and typical sizes only; the **large** rows are not measured yet (they run with
`DIFFRACTO_BENCH_LARGE=1`). Best of five after one warm-up: the median is
the figure, the spread (min-max) is in brackets. Same inputs on both sides (`tools/benchmarks/gen_inputs.py` writes them),
Float64 unless the function is Float32 by design (ACOM, ellipse pattern pixels).

| | |
|---|---|
| Machine | Apple M5 Pro, 64 GB, macOS 27.0.1 (Darwin 27.0.0) |
| Toolchain | Apple Swift 6.4, `swift build -c release` (separate package `tools/benchmarks`) |
| Python | 3.12.14, NumPy 2.5.4, SciPy 1.18.1, in two virtual environments: py4DSTEM 0.14.19 in one; pymatgen-core 2026.10.16 (the code; its `pymatgen` metapackage dist-info reads 2026.9.24), orix 0.15.0 and spglib 2.8.0 in the other |
| Load | **10.5 on a 12-core machine while other builds ran** (uptime 17:39): every figure carries that noise, the single-core rows least. Re-run when the machine is quiet before quoting a ratio. |
| Pin note | The pymatgen environment runs pymatgen-core **2026.10.16** under a `pymatgen` 2026.9.24 metapackage, and every fixture that records pymatgen-core records 2026.10.16 (`tools/pins/requirements/pymatgen.txt`). Nothing here compares numbers, only time |
| 0.9.1 | The `symmetry_detection` and `cif_import` rows above are the 2026-10-21 measurement of 0.9.0. 0.9.1 replaced the closure check that took 95 % of them (see `docs/open-items.md`); measured in-process on the same `al_super2.cif` (machine load 3.6): `cif_import` typical median 67 ms (min 67, max 77), against 297 ms (one run at load 9) before the change and the table's 351 ms. Not yet re-run through `run.sh` |

Re-run from a fresh clone of the kit, on a quiet Mac: `gen_inputs.py` writes the shared inputs, `run.sh` builds the benchmark
package (`tools/benchmarks`, release) and prints the kit's JSON lines, `bench_python.py` prints the Python lines, and `merge.py`
turns them into the table below. Add `DIFFRACTO_BENCH_LARGE=1` to `run.sh` and `bench_python.py` for the large rows.

```sh
cd tools/benchmarks
python3.12 -m venv .venv-orix && .venv-orix/bin/pip install pymatgen==2026.9.24 orix==0.15.0 spglib==2.8.0
python3.12 -m venv .venv-py4dstem && .venv-py4dstem/bin/pip install py4DSTEM==0.14.19
.venv-orix/bin/python -I -B gen_inputs.py                                # writes inputs/
./run.sh > results/rerun.jsonl                                           # kit side
.venv-orix/bin/python -I -B bench_python.py pymatgen >> results/rerun.jsonl   # pymatgen, spglib, orix
.venv-py4dstem/bin/python -I -B bench_python.py py4dstem >> results/rerun.jsonl
python3 merge.py results/rerun.jsonl                                     # prints the per-function table
```

The raw lines of the run below are `tools/benchmarks/results/2026-10-21.jsonl`.

## Per-function table (ms per call; `acom_match_scan` is ms per pattern)

| Function | Size | n | Kit | Python | Kit / Python (first Python row) |
|---|---|---|---|---|---|
| `niggli_reduce` | small | 100 | 0.0679 ms (0.064-0.0839) | 11.3 ms (9.96-11.7) [cpu]<br>0.293 ms (0.29-0.296) [spglib] | 0.006 |
| `niggli_reduce` | typical | 1000 | 0.51 ms (0.508-0.512) | 163 ms (149-267) [cpu]<br>3.13 ms (2.88-3.5) [spglib] | 0.0031 |
| `niggli_reduce` | large | 10000 | not measured yet | not measured yet [cpu]<br>not measured yet [spglib] | - |
| `reflections` | small | 10 | 0.0606 ms (0.0581-0.065) | 0.0261 ms (0.024-0.0353) | 2.3 |
| `reflections` | typical | 20 | 0.753 ms (0.56-0.916) | 0.115 ms (0.0952-0.191) | 6.6 |
| `reflections` | large | 40 | not measured yet | not measured yet | - |
| `kinematic_reflections` | small | 10 | 0.98 ms (0.849-1.14) | 4.45 ms (3.96-55.3) | 0.22 |
| `kinematic_reflections` | typical | 20 | 7 ms (6.64-7.55) | 12.7 ms (11.4-34.7) | 0.55 |
| `kinematic_reflections` | large | 30 | not measured yet | not measured yet | - |
| `symmetry_detection` | small | 4 | 6.87 ms (6.45-7.47) | 5.38 ms (4.78-5.69) | 1.3 |
| `symmetry_detection` | typical | 32 | 361 ms (350-438) | 13.1 ms (11.9-16.5) | 28 |
| `symmetry_detection` | large | 108 | not measured yet | not measured yet | - |
| `cif_import` | small | 962 | 0.652 ms (0.632-0.713) | 1.28 ms (1.27-1.31) | 0.51 |
| `cif_import` | typical | 2364 | 351 ms (329-506) | 2.52 ms (2.32-2.64) | 1.4e+02 |
| `cif_import` | large | 6483 | not measured yet | not measured yet | - |
| `orientation_reduce` | small | 1000 | 1.88 ms (1.87-5.39) | 22.1 ms (17.8-26) | 0.085 |
| `orientation_reduce` | typical | 10000 | 19.3 ms (17.8-29.4) | 48.2 ms (42.2-50.7) | 0.4 |
| `orientation_reduce` | large | 100000 | not measured yet | not measured yet | - |
| `ellipse_fit_1d` | small | 128 | 1.25 ms (1.2-1.41) | 1.34 ms (1.22-1.8) | 0.94 |
| `ellipse_fit_1d` | typical | 256 | 5.86 ms (4.35-17.7) | 5.26 ms (4.87-5.62) | 1.1 |
| `ellipse_fit_1d` | large | 512 | not measured yet | not measured yet | - |
| `acom_library` | small | 100 | 7.81 ms (6.88-16.9) | 13.7 ms (12.4-40.2) | 0.57 |
| `acom_library` | typical | 600 | 60.6 ms (49.4-80.7) | 62.5 ms (39-132) | 0.97 |
| `acom_library` | large | 2400 | not measured yet | not measured yet | - |
| `acom_match_scan` | small | 64 | 1.85 ms (1.83-1.88) [cpu 1 worker]<br>0.245 ms (0.235-0.274) [cpu all cores]<br>0.294 ms (0.285-0.314) [metal] | 25.7 ms (14.4-37.9) | 0.072 / 0.0095 / 0.011 |
| `acom_match_scan` | typical | 1024 | 1.94 ms (1.91-1.98) [cpu 1 worker]<br>0.216 ms (0.211-0.236) [cpu all cores]<br>0.236 ms (0.215-0.285) [metal] | 17.5 ms (10.4-21.5) | 0.11 / 0.012 / 0.014 |
| `acom_match_scan` | large | 16384 | not measured yet [cpu 1 worker]<br>not measured yet [cpu all cores]<br>not measured yet [metal] | - | - |

## What the rows compare, and what they do not

- `niggli_reduce`: the kit's `UnitCell.niggliReduced` against pymatgen's `Lattice.get_niggli_reduced_lattice` (pure Python) and
  spglib's C `niggli_reduce`; the kit is faster than both (the first Python row is pymatgen). Cells are built outside the timer.
- `reflections`: **not like for like.** The kit's `reflections(kMax:)` also gives d, the absence flag, the rings and their
  orbits; pymatgen's `get_points_in_sphere` returns the points only. The kit is 2.3 to 6.6 times slower on the points.
- `kinematic_reflections`: the kit's `kinematicReflections(kMax:retainingAll:false)` against py4DSTEM's
  `Crystal.calculate_structure_factors` on the same C2/m cell; the kit is 1.8 to 4.5 times faster. Both use the Lobato-Van Dyck table.
- `symmetry_detection` and `cif_import`: the kit's `CIFImport.read` always detects the symmetry of the cell it reads;
  pymatgen's `CifParser` only parses and `SpacegroupAnalyzer` is spglib in C. On a 32-atom P1 cell (an fcc supercell with
  1 536 operations) the kit took 361 ms against 13 ms, **28 times slower**; on a four-atom cell they are level. The cause was the
  closure check, not the symmetry search (`docs/open-items.md`); 0.9.1 made the check about five times faster, the rest stays listed there.
- `orientation_reduce`: the kit reduces one matrix at a time (24 proper operators); orix reduces a vectorised batch
  (`Orientation.from_matrix` + `map_into_symmetry_reduced_zone`). The kit is 2.5 to 12 times faster at these sizes.
- `ellipse_fit_1d`: level with py4DSTEM's `fit_ellipse_1D` (0.94 and 1.1).
- `acom_library`: level with `orientation_plan` for about 600 zone axes (0.97); py4DSTEM's count comes from its angular step
  (630 at 1.6 degrees), the kit's from `zoneAxisCount`. 128 in-plane steps on both sides.
- `acom_match_scan`: the kit's `matchAll` over the whole scan against py4DSTEM's `match_single_pattern` per pattern (one
  core; Python's per-pattern time includes its own refinements, so this is a time, not a like-for-like algorithm). Per pattern
  the kit on one core is 14 times (8 x 8 scan) and 9 times (32 x 32) faster than py4DSTEM on one core, and on all cores 105 and 81 times.

## Backend decision per function

| Function | Backend | Why / break-even |
|---|---|---|
| `niggli_reduce`, `reflections`, `kinematic_reflections`, `symmetry_detection`, `cif_import` | CPU only | small arrays, branchy integer logic, readers; the kit is already faster than Python except where `docs/open-items.md` lists it |
| `orientation_reduce` | CPU only | 24 small matrix products per rotation; the batch (100 000 in 0.2 s) is bound by memory traffic, a GPU would pay the transfer for no gain at these sizes |
| `ellipse_fit_1d` | CPU only | one 1D least-squares fit per pattern; level with Python |
| `acom_library` | CPU only | built once per phase (60 ms for 600 axes), not per pattern |
| `acom_match_scan` | **CPU all cores; Metal not above any measured N** | Metal 0.294 ms per pattern at 64 patterns, 0.236 at 1 024, against 0.245 and 0.216 on all CPU cores: the break-even is **not reached up to 1 024 patterns** (the Metal path is 7 to 22 % slower); the 128 x 128 row (16 384 patterns) is not measured yet, and the Metal path stays an option (`ACOMBackend.metal`), not the default. A larger library (2 400 axes) is the other axis still to measure |

Numbers in `docs/parity.md` are untouched by this table (time only).
