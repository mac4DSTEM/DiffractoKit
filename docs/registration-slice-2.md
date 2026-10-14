# DiffractoKit slice 2: electron scattering factors, structure factors and the kinematic ring list (registration)

Registered 2026-10-15, about 01:10 CEST, before any code, by the DiffractoKit session (the owner's instruction relayed by
the coordinator at 00:55: "build slice 2, then slice 3"). Tier 1: |F|, |F|² and which rings survive the cut are numbers a
scientist reads (and slice 3's indexing and vector matching consume). What is written here is fixed before the code; a
change after code starts is a dated amendment with its reason.

## 1. What slice 2 makes

The module `DiffractoScattering` (depends on `DiffractoCore` only). Public API only, no app types.

- `ElectronScatteringFactors`: the Lobato & Van Dyck (Acta Cryst. A70, 636–649, 2014) parameterisation as py4DSTEM
  ships it, Z = 1 ... 103; `f(z:gSquared:)` in Å (nil outside the table), `isSupported(z:)`, the raw parameters.
  f_e(g²) = Σ_{i<5} a_i (2 + b_i g²) / (1 + b_i g²)², g² = |g|² in Å⁻² (no 2π).
- `StructureFactor` (h, k, l, the Cartesian g, |g|, Re F, Im F, |F|, |F|²) and the function that computes it for a
  `CrystalStructure`: F(hkl) = (1/V) Σ_sites f(Z, |g|²) · occ · exp(2πi (h x + k y + l z)), in Å (py4DSTEM's "A"
  convention). The form factor is a parameter (`formFactor: (Int, Double) -> Double?`, default the electron one), so
  the phase sum can be checked with another set of atomic form factors (section 4.3) and slice 3 can reuse it.
- `KinematicReflectionList`: every reflection of slice 1's list (0 < |g| ≤ k_max, the same enumeration and order) with
  its F; `kept` = the reflections with |F| > `cut` (default 1e-4 Å, the apps'), whether or not the space-group
  operations allow them (an operation-allowed reflection whose F vanishes by the basis, Si {2 2 2}, is dropped by the
  cut; an operation-absent one has |F| ≤ 1e-12 and is dropped too, a self-check); `families` = slice 1's ring rule
  (1e-9 relative, the structure's Laue group for orbits and names) applied to `kept`, each ring carrying Σ|F|² over its
  members. Refusals: slice 1's k_max and size refusals, and `unsupportedElements([Int])` when a site's Z has no row
  (the apps' `unsupportedElements`, checked before any F is computed).

Not in slice 2: Debye–Waller factors, absorption, dynamical effects, relativistic or wavelength-dependent corrections
(py4DSTEM's "VA" unit, the interaction constant), excitation errors and the Ewald construction (slice 3), X-ray or
neutron scattering factors as public API.

## 2. The references

- **The table.** py4DSTEM `f050d207851f8accd3aa550365baaf1f38230428` (dev, 2026-03-26, 0.14.19; mac4DSTEM's lock and the pin
  macTEM's `tools/lib/fetch-py4dstem-scattering.sh` uses), `py4DSTEM/process/utils/scattering_factors.txt`, sha256
  `be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4`, 103 rows, row Z−1 = a0 b0 … a4 b4; fetched into the
  gitignored `References/py4dstem-scattering/` by `tools/slice-2/fetch-table.sh` (macTEM's script, adapted) and refused on
  any other sha256. The Swift table is mac4DSTEM `7d6fdbd1`'s `ScatteringFactors.swift` (identical to `7448401b`'s), GPL-3.0.
- **The seeds.** mac4DSTEM `7d6fdbd1` `mac4DSTEM/Core/Crystal/Crystal.swift` (`reflections(kMax:tolerance:)`, the
  structure-factor loop) and `ScatteringFactors.swift`; macTEM `a58c84c` `tools/crystal-parity/` (the NumPy reference
  written from the definitions, its cases, its controls), read only by `git show`.
- **Oracle 1, NumPy from the definitions** (`tools/slice-2/pins.py`, run in `References/pymatgen-venv`, NumPy 2.5.4):
  the table read from the fetched file, F from the definition on the sites pymatgen's `CifParser` expands (or the
  fixture's own list), every reflection with |g| ≤ k_max (1 + 1e-6) and |F| > 1e-6 listed so the cut's edge band is
  visible. Written from this text, not from the Swift.
- **Oracle 2, pymatgen's X-ray structure factor** (pymatgen-core 2026.9.23 `analysis/diffraction/xrd.py`,
  `XRDCalculator.get_pattern(structure, scaled=False)`): an independent implementation of the phase sum (its own site
  expansion, its own reciprocal-space enumeration, its own sum), with X-ray form factors. Each peak's unscaled intensity
  divided by the Lorentz–polarisation factor and by 2 (pymatgen's Friedel doubling) is Σ|F_X|² over the half-space
  members of the ring; the kit's phase sum with pymatgen's X-ray form factor (its `atomic_scattering_params.json`
  coefficients, f(s) = Z − 41.78214 s² Σ a_k exp(−b_k s²), s = |g|/2) over the same members must give the same sum.
  This checks the phase sums and the enumeration, not f_e.
- **Oracle 3, the apps' code** for bit identity: mac4DSTEM `7d6fdbd1`'s `Crystal.swift` and `ScatteringFactors.swift`,
  compiled by `git show` into a scratch tool (`tools/slice-2/apps-parity/`), on the apps' presets.

Agreement with py4DSTEM's table, NumPy, pymatgen or the apps is agreement, not truth; truth is section 5.

## 3. The definitions, pinned

### 3.1 The scattering factor
As section 1, py4DSTEM `single_atom_scatter.electron_scattering_factor` with units "A" (mac4DSTEM's DEVIATION note: the "VA"
branch is not ported). The parameters are the table's doubles: every Swift literal equals `np.loadtxt` of the fetched file
bit for bit (a test over all 103 × 10).

### 3.2 The structure factor
For each reflection of slice 1's list, in this order of operations (the apps' expression, so the result can be bit
identical): g = h b₁ + k b₂ + l b₃ (slice 1's `gVector`), |g| = `length(g)`, g² = |g|·|g|; for each site in the
structure's site order: f = f_e(Z, g²); phase = 2π (h x + k y + l z) with the products and sums in that order; Re +=
f · occ · cos(phase), Im += f · occ · sin(phase); then Re /= V, Im /= V with V = `UnitCell.volume`; |F|² = Re² + Im²;
|F| = √|F|². The cut keeps |F| > cut (strict).

### 3.3 Families after the cut
Slice 1's rule on the kept reflections, unchanged: increasing |g|, ties by descending (h, k, l); a reflection joins the
current ring while |g| ≤ g₀(1 + 1e-9); the ring's g is its smallest member's; members, multiplicity (member count),
orbits (canonical forms under the structure's Laue group, descending), name (Miller–Bravais in hexagonal axes only), and
Σ|F|² over the members. A ring may lose members to the cut when its orbits have different |F| (two coincident orbits,
Si {5 1 1} + {3 3 3}, keep both: both are strong); a single orbit has one |F| (the Laue group preserves |F| for a
centrosymmetric or Friedel-paired structure, kinematically), so the cut takes a whole orbit or none: a self-check.

## 4. Pins

### 4.1 NumPy (`slice-2-pins.json`, written before the Swift)
For the structures: the apps' presets Au (fcc 4.0782), Al (fcc 4.0495), Si (diamond 5.4309), Mg (hcp 3.2094 / 5.2108), WS₂
(2H, Schutte 1987), β″ (the generator's CIF), Cu3Au-occ (macTEM's L1₂-like disordered cell, a 3.75 Å, Au 0.6 / Cu 0.4 on
the origin, Cu on the faces), Weak-cut (macTEM's CsCl-like Cu cell, a 3 Å, the second site's occupancy solved so (1 0 0)
has |F| = 2.5e-4), and slice 1's CIF fixtures with elements in the table: every listed reflection's h, k, l, |g|, Re F,
Im F, |F|; the kept count; the families after the cut (count, g, members, multiplicity, orbits, name, Σ|F|²).

Agreement: the kept set equal except reflections whose NumPy |F| is within a factor 10 of the cut (listed, compared
separately: the same side of the cut unless within 1e-12 relative of it); Re F, Im F and |F| within 1e-9 relative of
max(|F|, 1e-4); |F|² within 2e-9 relative; families equal (count, g 1e-12, members, multiplicity, orbits, name; Σ|F|²
1e-9 relative). A difference is a failure of the port, never a tolerance to widen.

### 4.2 The table
All 103 rows × 10 parameters equal to the fetched file bit for bit; f_e at g² ∈ {0, 0.25, 1, 4} for every Z within
1e-15 relative of NumPy's value from the file.

### 4.3 pymatgen X-ray (independent phase sum)
For Au, Si, Mg, β″, Cu3Au-occ and WS₂ at λ = 0.5 Å (so 2θ reaches |g| = 1.2 Å⁻¹ well inside 180°): for every pymatgen
peak, its Σ|F_X|² (as section 2) equals the kit's Σ over the same half-space members of |F_X|² computed through the
kit's phase sum with pymatgen's X-ray form factor, within 1e-6 relative (pymatgen's 2θ binning and arcsin are the
limiting precision; recorded in the pins). The peak's member set is matched by the pins' (h, k, l) list.

### 4.4 The apps, bit for bit
On the presets Au, Al, Ni, Cu, Fe (bcc), Si, Mg, WS₂ and β″ at k_max 1.2 Å⁻¹ (and β″ at 2.0): the kit's kept set equals
the apps' `Crystal.reflections(kMax:tolerance: 1e-4)` set, and every kept reflection's |g|, Re F, Im F and |F|² are bit
identical (Double bit patterns) to the apps'. The apps' tool reads its own presets; the kit builds the same cells and
sites (β″ through `CIFImport` of the generator's text: its site ORDER differs from the apps' hand list, so for β″ the
gate is 1e-14 relative on Re F and Im F, a summation-order difference, recorded; bit identity is asserted for the presets
whose site order is the same).

## 5. Truth (independent of every oracle)

- **Extinctions, exact.** fcc (Au, Al): |F| ≤ 1e-12 for h, k, l of mixed parity, |F| > 1e-3 otherwise within 1.2 Å⁻¹.
  Diamond (Si): additionally h + k + l = 4n + 2 with all even ≤ 1e-12 ({2 0 0}, {2 2 2}, {4 2 0}, {4 4 2}); hcp (Mg): 000l
  with l odd and h − k = 3n with l odd ≤ 1e-12.
- **Rings after the cut, the apps' counts** (k_max 1.2 Å⁻¹): Au 8, Al 8, Si 10 (12 by symmetry minus {2 2 2} and {4 4 2}),
  Mg 24 (25 by symmetry minus {3 0 −3 1}), ideal hcp 23 (macTEM's T: 24 by symmetry minus one), β″ **417 rings, 1442
  reflections** (1448 by symmetry minus {3 3 3} (4 members, |F| 9.3e-5) and {14 0 1} (2 members, 7.9e-5), macTEM's S2c
  record, gemmi-checked); β″ at k_max 1.6 / 2.0 Å⁻¹: 953 / 1805 rings (macTEM's S2c predictions).
- **Closed forms.** fcc: F = 4 f_e(Z, g²) / V for h, k, l unmixed. Diamond: F = 4 f_e (1 + i^{h+k+l}) / V, so |F| = 4√2 f_e/V
  for h + k + l odd and 8 f_e / V for h + k + l = 4n. hcp Mg: |F| = 2 f_e |cos(π((h + 2k)/3 + l/2))| / V. Cu3Au-occ:
  unmixed h, k, l: (0.6 f_Au + 3.4 f_Cu)·…/V by the definition; mixed: (0.6 f_Au + 0.4 f_Cu − f_Cu)/V (the superlattice
  reflections exist only through the occupancy). Each to 1e-13 relative.
- **The weak cut.** In Weak-cut, the h + k + l odd reflections straddle 1e-4: at least 3 kept in (1e-4, 1e-3) and 3 dropped
  in (1e-5, 1e-4), each on the side the closed form (1 − q) f_e / V puts it.
- **f_e itself.** Positive and strictly decreasing in g on [0, 1.2] Å⁻¹ for Z ∈ {1, 6, 12, 13, 14, 29, 74, 79, 103};
  f_e(Z, 0) = 2 Σ a_i.

## 6. Gates and predicted outcomes

1. Pins (4.1–4.4) within their bounds; predicted: equal (the code is the apps', the table the file's).
2. Truth (section 5) exact or to its stated bound; predicted: the counts above.
3. Every new test broken once (mutants planned: the table row Z−1 → Z; f_e of |g| instead of |g|²; no 1/V; the cut ≥
   instead of >; the cut at 1e-3; occupancy ignored; the phase sign flipped (only Im F and the X-ray pin see it on
   non-centrosymmetric structures: WS₂, quartz); the diamond basis dropped; the sum over sites skipping the last site;
   families on the uncut list).
4. One Gate B refuter on the frozen slice-2 commit (a second only on a named doubt); Gate D if a cause is unproven.
5. `tools/test.sh` green with the stamp; consumer-check extended to `DiffractoScattering`.

## 7. What is not claimed

Kinematic, single-scattering, no thermal motion: a measured ring intensity can differ by large factors (dynamical
scattering, Debye–Waller, texture); the cut is a numerical threshold, not a visibility claim. Agreement with py4DSTEM's
table is agreement with Lobato & Van Dyck's fit, whose own accuracy (about 1e-3 of the Hartree–Fock values it fits) is
theirs. No real pattern is indexed here; that is slice 3's and the apps'.

## 8. Amendments

### 2026-10-15, after the pin generator ran, before the pins test (no registered number changes)

`tools/slice-2/pins.py` (7138 self-checks) and the apps' dump (`tools/slice-2/apps-parity`, mac4DSTEM 7d6fdbd1 compiled -O)
held every predicted count: rings after the cut Au 8, Al 8, Si 10, Mg 24, ideal hcp 23; β″ 417 rings, 1442 reflections,
953 / 1805 rings (3420 / 6596 reflections) at 1.6 / 2.0, dropping exactly {3 3 3} and {14 0 1}; the apps' kept sets equal
NumPy's; the X-ray phase sums equal pymatgen's to 1.8e-13 (β″) and ≤ 2.5e-15 elsewhere. Five readings are settled:

1. **Bit identity (4.4) excludes Im F below 1e-15 Å absolute.** The apps' own -O and -Onone builds differ only there (β″:
   6 reflections at 1.2, 16 at 2.0, Im F about 1e-18, rounding noise of a centrosymmetric sum); |g|, Re F and |F|² are bit
   identical between them. The gate: bit identity of |g|, Re F and |F|², and of Im F where |Im F| ≥ 1e-15; below that,
   |Δ Im F| ≤ 1e-15. The -O dump is the reference.
2. **The phase-sign mutant is seen by Im F on a non-centrosymmetric structure, not by WS₂ or the X-ray pin.** WS₂ (and
   the Mg and WS₂ CIF fixtures) are centrosymmetric about the origin (Im F/|F| ~ 1e-14), and Σ|F|² is invariant under
   conjugation. The mutant is aimed at quartz, the Si and Mg presets (atom at the origin, not the centre) and an added
   preset, zinc-blende GaAs (a = 5.6533 Å, F = 4(f_Ga + i^{h+k+l} f_As)/V, a closed form for the sign).
3. **CIF structures are pinned on their sites as written** (the kit keeps 0.333333; pymatgen snaps to ⅓): there the
   operation-absent |F| is up to 3.0e-8 (Mg) and 7.9e-7 (WS₂), far below the cut, and ≤ 1e-12 holds on the snapped sites
   and on the presets (1.6e-16). Section 5's "≤ 1e-12" is for exact sites; for written sites the bound is the cut's
   side (no absent reflection is kept).
4. **The cut's strictness gets its own case.** No pinned structure has |F| exactly 1e-4, so `>` versus `≥` was invisible.
   `cut_edge`: a = b = 100, c = 1 Å, one site at the origin, a constant unit form factor (the `formFactor` parameter):
   every reflection has |F| = 1/V = 1e-4 exactly (V by the cofactor determinant, as the kit and the apps); `>` keeps none.
6. **β″ against the apps (4.4).** Bit identity holds when the kit builds β″ from the apps' own hand site list (a test, at
   1.2 and 2.0 Å⁻¹). Built through `CIFImport` its 22 sites come in another order, and the sum's rounding differs by at
   most 2.7e-13 of max(|F|, 1e-4) on Re F and 5.5e-13 on |F|² (measured; |g|, V and the kept sets are bit equal). The
   registered 1e-14 for the CIF-order comparison was a guess below a summation-order effect; it is 1e-12 (Re F, Im F)
   and 2e-12 (|F|²). The `cut_edge` count may differ from NumPy's by the reflections exactly on |g| = 1.2 (2 of 72890),
   where rounding decides inclusion; the test allows exactly those.
5. **pymatgen's 2θ range** is limited to 0.1 % beyond |g| = 1.2 (1 + 1e-6): with no range the Lorentz–polarisation
   divergence near 180° sets the scale and pymatgen drops every low-angle peak. Recorded in the pins.
