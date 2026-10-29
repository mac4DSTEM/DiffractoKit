# tools/kit-0.0.4: pins of DiffractoKit 0.0.4

Each script works in a copy of the repository, never the working tree.

| Script | What | Result (2026-10-15) |
|---|---|---|
| `q-pins.py` | q and its inverse, 60-digit decimal truth computed two ways (algebraic and series, agreeing to 1e-50) beside NumPy 2.5.4; writes `Tests/DiffractoGeometryTests/Fixtures/q-pins.json` (byte-identical on rerun) | 100 forward, 120 inverse, 320 identity, 51 NaN plants |
| `tetragonal-pins.py` | pymatgen 2026.9.24 / spglib 2.8.0 (symprec 1e-4) on rutile, anatase, β-Sn in three settings and two controls; writes `Tests/DiffractoCoreTests/Fixtures/tetragonal-pins.json` | the kit equals spglib in every setting |

The 0.0.4 changes were also checked by breaking the code once per mutation (2026-10-15): the ring fit 18 of 18 red, the
tetragonal setting and the error texts 12 of 12, q and radius 6 of 6. Each test's `RED UNDER` line names the mutation that
turns it red.
