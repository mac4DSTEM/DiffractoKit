# DiffractoKit

The macEM suite's crystallography and diffraction library for macOS: a Swift package of numbers and readers, never UI, that
mac4DSTEM and macTEM import by tagged version and that other people can build their own apps on. Mac-native (Swift 6 strict
concurrency, Accelerate where a path is accelerated, every accelerated path with a CPU reference and a registered tolerance),
headless, one SwiftPM product per module.

| Module | Scope | State |
|---|---|---|
| `DiffractoCore` | cells of all seven crystal systems, symmetry operations and Laue groups, CIF import, d-spacings, the reflection list with systematic absences, orbits and families (rings) | slice 1, `docs/registration-slice-1.md` |
| `DiffractoScattering` | electron scattering factors, structure factors, kinematic intensities | slice 2, to be registered |
| `DiffractoGeometry` | pattern geometry, ring and ellipse calibration, camera length, zone axes, indexing | later |
| `DiffractoOrientation` | orientation results, symmetry reduction, fundamental zones, IPF keys, Euler conventions | later |

## Licence

GPL-3.0 (`LICENSE`, `NOTICE`). The code of slice 1 is ported from mac4DSTEM and macTEM (both GPL-3.0; `github.com/mac4DSTEM`),
whose lattice and reflection conventions follow py4DSTEM's `Crystal` class (GPL-3.0): the real-lattice matrix with a along x
and b in the xy plane, reciprocal vectors without 2π, the per-axis index bound. The electron scattering table of slice 2 will be
py4DSTEM's (Lobato & Van Dyck 2014), also GPL-3.0. **pymatgen** (MIT) is the oracle of slice 1, cells, symmetry and
d-spacings are checked against it, and no pymatgen code is copied or ported; a part ported from pymatgen in the future would be
MIT on its own, but one package with any GPL part is GPL, so the whole is GPL-3.0 and an app built on it must be
GPL-compatible. The CIF fixtures under `Tests/` come from mac4DSTEM's tests (GPL-3.0) with their provenance in a comment line;
the β″ cell content (Andersen et al. 1998) is written by a generator into the gitignored `References/`, as the apps do.

## Build and test

```sh
tools/test.sh                      # the gate: swift test -c release, stamps .build/last-green on a clean tree
tools/last-green.sh                # the last green commit and whether HEAD equals it
tools/consumer-check/run.sh        # a separate package using the public API only (slice 1's skeleton)
```

The oracle venv is not in the repository: `docs/registration-slice-1.md` section 2 and `CLAUDE.md` say how it is pinned and
provisioned (`References/pymatgen-venv`, pymatgen 2026.9.24 with pymatgen-core 2026.9.23).

## Status

`docs/status.md`. No tag yet; the first tag is the owner's after slice 1's refuter.
