# CLAUDE.md — start here

DiffractoKit: the macEM suite's crystallography and diffraction library, a Swift port of what mac4DSTEM and macTEM compute
about crystals (cells, symmetry, CIF import, d-spacings, orbits; scattering and structure factors; pattern geometry;
orientation), checked against pymatgen, py4DSTEM and pyxem, imported by tagged version by both apps. GPL-3.0 (README.md
says which parts come from where). Modules: `DiffractoCore` (slice 1), `DiffractoScattering` (slice 2), `DiffractoGeometry`,
`DiffractoOrientation` (placeholders with their planned scope in a doc comment). Rules only here; state in `docs/status.md`.

## Read first

1. The suite's `../CLAUDE.md` (board, tiers, how we work) and `../ROADMAP.md`, "Kits" (the direction, the oracles, the licence
   rule, the Apple-infrastructure requirements).
2. `docs/registration-slice-1.md`: what slice 1 makes, the definitions, the pins, the planted truth, the gates.
3. `docs/status.md`: what is live, the hand-off, what is unverified.
4. `../SpectroKit/CLAUDE.md`: the model for a kit repo (this file follows it).

## Hard rules

- **mac4DSTEM and macTEM are read-only from here**: `git -C ../<repo> --no-optional-locks show <commit>:<path>`, never their
  working trees. Never build, test or write in them. The seeds of slice 1 were read at mac4DSTEM `7448401b` and macTEM `a58c84c`
  (their pushed heads on 2026-10-14 night).
- **Registration before code** for anything tier 1 (every number a scientist reads): `docs/registration-<slice>.md`, with the
  definitions, the pins, the planted truth, the gates and the planned mutants, committed before the first Swift line; a change
  after code starts is a dated amendment with its reason.
- **The gate is `tools/test.sh`** (`swift test -c release "$@"`). Never run `swift build --build-tests -c release` (fails on
  this toolchain). The manifest carries no unsafe flags, so a consumer builds the kit optimised itself.
- **Green stamp**: a full, passing `tools/test.sh` run on a clean tracked tree writes `.build/last-green` (commit, UTC time);
  `tools/last-green.sh` prints it and says whether HEAD equals it.
- **Hook**: `.claude/settings.json` runs `tools/hooks/no-sleep-poll.sh` before every Bash call; it refuses a sleep of 30 s or
  more and any `while`/`until` loop with a sleep in it. Wait on background-task notifications or Monitor. Proof:
  `tools/hooks/test-no-sleep-poll.sh`.
- **A reference for every number.** Pins: pymatgen 2026.9.24 (tag `v2026.9.24` = `8461938`, the umbrella package) with
  pymatgen-core 2026.9.23 (the package that holds `core`, `io.cif` and `symmetry`; its tag is cloned beside the venv for
  file:line), spglib 2.8.0, NumPy 2.5.4, Python 3.14.7, all in the gitignored `References/` (`pymatgen-venv`,
  `pymatgen-8461938`, `pymatgen-core-<sha7>`, `pymatgen-venv-freeze.txt`). Oracles run as `References/pymatgen-venv/bin/python
  -I -B`. An unpinned copy is never a reference. Agreement with pymatgen is agreement, not truth: truth is a closed form, a
  group-theory count or a cell content the file's formula states.
- **Pins are generated before the Swift code** (`tools/slice-1/pins.py` writes `Tests/DiffractoCoreTests/Fixtures/
  slice-1-pins.json`); the Swift tests read the pins and the planted truth; a difference is a failure of the port, never a
  tolerance to widen.
- **Public API only, no app types.** Nothing an app needs is `package`; nothing here names `AppState`, `CrystalModel`,
  `ACOMCrystalSymmetry` or any other app type. `tools/consumer-check` is a separate package that depends on this one by path
  and compiles only while every symbol it names is `public`. Error messages never say "mac4DSTEM", "macTEM" or "ACOM".
- **Rigour follows risk** (`../CLAUDE.md`): tier named in every commit; tier 1 is registration, every new test broken once, one
  Gate B refuter on the frozen commit, Gate D when a cause is unproven, a fixture. Every accelerated path, when one comes,
  keeps a CPU reference path with a registered tolerance.
- **Names** (the suite's standard): the library is domain + `Kit`, modules carry its prefix: `DiffractoCore`,
  `DiffractoScattering`, `DiffractoGeometry`, `DiffractoOrientation`.
- **Pushing, tagging and creating the remote are the owner's.** One session per repo; `main` only; commit by explicit paths at
  every green checkpoint. At each tag, `docs/consumers-x.y.z.md` says what each consumer changes.
- **The owner's data** lives in the suite's `data/` folder, read as `$MACEM_DATA` (default `../data`): never copied or
  committed here; a manifest line in `../data/README.md` per file used.
- **Suite board**: this repo's row is `lanes/diffractokit`; how a session opens and closes it is in `../CLAUDE.md`.
