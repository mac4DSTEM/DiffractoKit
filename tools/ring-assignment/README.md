# Ring assignment: planted cases

`planted_cases.py` writes 64 planted sets of ring radii for `LinearScale.assign`: fcc Al, bcc Fe, diamond Si (one with the
forbidden {200} and {222} spots) and hcp Mg; single zone axes, three grains and polycrystals; tight and broad rings; with and
without a 9 % elliptical distortion. The planted scale of each case is in the file. numpy only; no upstream number is claimed
(py4DSTEM does not assign the first ring).

```sh
python3 -I tools/ring-assignment/planted_cases.py "$TMPDIR/ring-cases.json"    # a few seconds
```

Read the JSON, call `LinearScale.estimateAssigned` on each case and compare the scale with `dq`. Measured 2026-11-10: 20
right (within 5 %), 37 refused, 7 wrong; the 7 are single crystals with tight spots under the 9 % ellipse, the Al [111]
zone (two rings, {220} and {422}, whose ratio is also {111}, {220} within the tolerance) and one broad Si polycrystal
under the ellipse, and all 7 were wrong before the two readings had to agree. Do not commit the output.
