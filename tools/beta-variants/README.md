# beta'' variants: planted truth (DiffractoKit task 35)

`planted.py` builds the 12 orientation-relationship variants of beta'' (Mg5Si6, C2/m) in Al as rotation matrices, with their
pairwise disorientation angles in closed form and a set of seeded "measured" orientations whose nearest variant is known.
numpy only; no orix or diffsims is run and no upstream number is claimed. The Mac pin script (kit audit step 0) adds
`Misorientation.equivalent()` / `reduce()` at atol 1e-4 from the same matrices.

```sh
python3 -I tools/beta-variants/planted.py --out "$TMPDIR/beta-planted"    # under 1 s
```

Output: `variants.npy` (12, 3, 3), `measured.npy` (48, 3, 3), `truth.json`. Do not commit the output.

## Truth it carries
- 12 variants (24 cubic rotations / the 2-element precipitate group; the precipitate's C2 about b is the cubic 180 degree
  rotation about [001]Al), four needle variants for each of the three <001> axes.
- Variant-to-variant disorientation (2/m on both sides): only 0, 90, 120, 180 degrees; ordered-pair counts 60 / 48 / 24.
- The OR's own reduced disorientation (432 on the left, 2/m on the right): arctan(2/3) = 33.690068 degrees.
- 48 measured orientations (0.5 degree seeded noise, random 2/m setting) with the true variant; the nearest by disorientation
  is unique (second best at least 89 degrees away).

## How the Mac session uses it
1. `CrystalPresets.betaDoublePrime` plus `OrientationRelationship.variants(matrix:precipitate:)` must return the 12 classes of
   `variants.npy` (compare as classes under the right group, atol 1e-9) and the angle matrix in `truth.json`.
2. `nearestVariant(of:)` on `measured.npy` returns `truth.json.measured.true_variant` exactly.
3. The pin script runs orix 0.15.0 `Misorientation.equivalent()` on the same R0 and compares at atol 1e-4 (the orix test's own).
4. Frame: precipitate Cartesian e1 || a, e2 || b, e3 = e1 x e2; if the kit's `UnitCell` Cartesian setting differs, convert R0
   in the test wrapper and note it in `docs/parity.md`.

## Open point for Paul
The OR (a || [3 2 0]Al, b || [001]Al, c || [1 -3 0]Al) is derived in the script from the lattice parameters (105.25 degrees
against beta = 105.3, a misfit 3.7 %); it is not copied from the paper. Scope question 1 (which OR the study used) decides it;
replacing it means editing `R0` in `planted.py` only. Check against Andersen et al. 1998 before the figure quotes it.
