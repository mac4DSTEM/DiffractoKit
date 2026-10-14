# consumer-check

A separate SwiftPM package that depends on DiffractoKit **by path** (`.package(path: "../..")`) and uses, from outside the
package, the public API of slice 1: `UnitCell`, `MillerIndex`, `IntegerMatrix`, `SymmetryOperation`, `PointSymmetry`,
`LaueClass`, `CrystalSystem`, `Element`, `AtomSite`, `CrystalStructure`, `SymmetryDetection`, `CIFImport`, `CIFDocument`,
`ReflectionList`, `Reflection`, `ReflectionFamily`, the typed errors, and the three placeholder modules' scope statements.
`package` access does not cross a package boundary, so the target compiles only while each symbol it names is `public`;
the values it checks are ones a reader verifies by hand (β″'s d(200) = a sin β / 2, gold's eight families).

This is the slice-1 **skeleton**: it proves the surface is public and reachable. The negative controls (SpectroKit's
`controls.py` pattern: revert one declaration to `package` in a copy, require the build to fail naming the symbol) and the
per-consumer symbol table come with the first tag's `docs/consumers-x.y.z.md`.

```
tools/consumer-check/run.sh
```
