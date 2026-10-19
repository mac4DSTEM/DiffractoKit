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

If the Swift code and a fixture disagree, that is a bug in the code, not a tolerance to widen. If you believe a fixture is
wrong, say why in the pull request, and regenerate it with its script and the pinned versions of the reference.

## Pull requests

- Keep a change focused, with a test that fails before the change and passes after it.
- Tests use Swift Testing (`import Testing`, `@Test`, `#expect`).
- Public API only: new types and functions are `public`, and the separate package in `tools/consumer-check` must still build.
- Add a line to `CHANGELOG.md` for anything a user can see.
- By contributing you agree that your contribution is licensed under GPL-3.0-or-later, like the rest of the library.

## Reporting issues

Open an issue with the DiffractoKit version, your platform and Swift version, a small CIF or code sample that reproduces the
problem, what you expected and what you got. For security problems see `SECURITY.md`.
