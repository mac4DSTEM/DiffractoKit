# Contributing to DiffractoKit

Thank you for your interest. Bug reports, questions and pull requests are welcome.

## Build and test

You need a Swift 6.2 toolchain on macOS 15 or later.

```sh
swift test -c release
```

The package builds optimised code in release mode; a Debug build of the reflection loops and the ring fit is slow. The tests
read committed fixtures and need no Python.

## How numbers are checked

Every number a scientist reads is covered by a test against a reference: pymatgen and spglib for cells, symmetry and
d-spacings, py4DSTEM for scattering factors and diffraction geometry, SciPy and NumPy for numerical results, or a closed form
or planted truth where no package applies. The reference outputs are committed as JSON fixtures under `Tests/*/Fixtures/`,
with the generator scripts under `tools/`. `docs/parity.md` lists what is held against which reference.

A test's tolerance is the one the Python package's own tests use for that function; a result outside it is a bug. A change
that moves a fixture says why in the pull request and in the CHANGELOG; regenerate a fixture with its script and the pinned
versions of the reference.

## Regenerating fixtures

The generators run in pinned Python virtual environments outside the repository (Python 3.12, NumPy 2.5.4, SciPy 1.18.1;
pymatgen, spglib and orix in one, py4DSTEM in the other), placed at `../References/python/orix` and
`../References/python/py4dstem` beside the checkout; each script's header names the one it needs. The current `orix` venv
holds pymatgen-core 2026.10.16 under a `pymatgen` 2026.9.24 metapackage; the README and `docs/parity.md` say 2026.9.24, and the
older fixtures record pymatgen-core 2026.9.23 (a previous venv). `cif-symmetry-pins.json` was regenerated on the current
install on 2026-10-21 and came out byte-identical; the other fixtures were not regenerated. Regenerate and re-pin them
together at the next fixture change.

## Pull requests

- Keep a change focused, with a test that fails before the change and passes after it.
- Tests use Swift Testing (`import Testing`, `@Test`, `#expect`).
- Public API only: new types and functions are `public`, and the separate package in `tools/consumer-check` must still build.
- Add a line to `CHANGELOG.md` for anything a user can see.
- By contributing you agree that your contribution is licensed under GPL-3.0-only, like the rest of the library.

## Reporting issues

Open an issue with the DiffractoKit version, your platform and Swift version, a small CIF or code sample that reproduces the
problem, what you expected and what you got. For security problems see `SECURITY.md`.
