# DiffractoKit slice 1: cells, symmetry, CIF import, d-spacings and reflection orbits (registration)

Registered 2026-10-14 night, before any code, by the DiffractoKit session (the coordinator's side session; the owner's
direction in the suite `ROADMAP.md`, "Kits": DiffractoKit after β″ lands in both apps, pymatgen for cells, symmetry and
d). Tier 1: every number here is one a scientist reads (a d-spacing, a multiplicity, an orbit name, the atoms of a cell).
What is written here is fixed before the code; a change after code starts is an amendment, dated, with its reason.

## 1. What slice 1 makes

One module and product, `DiffractoCore`, with no dependency beyond Foundation and simd. Public API only (the apps are
not clients yet; nothing is `package`). Cells of all seven crystal systems, **monoclinic first** in every list, test
and example, because β″ (Mg5Si6, C2/m) is what the owner's study needs.

- `UnitCell`: a, b, c (Å), α, β, γ (°); the real lattice matrix in the apps' (py4DSTEM / VESTA) convention, a along x,
  b in the xy plane; the reciprocal lattice (no 2π: |g| in Å⁻¹, d = 1/|g| in Å); the metric tensor and its inverse; the
  volume; `d(h, k, l)`, `g(h, k, l)` (the Cartesian vector), `gLength`; the crystal system of the metric at a stated
  tolerance; the angle between two plane normals. Constructors per system that take only the free parameters
  (`cubic(a:)`, `tetragonal(a:c:)`, `orthorhombic(a:b:c:)`, `hexagonal(a:c:)`, `rhombohedral(a:alpha:)`,
  `monoclinic(a:b:c:beta:)` (unique axis b), `triclinic(...)`).
- `CrystalSystem` (the seven) and `LaueClass` (the eleven: -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm, m-3, m-3m).
- `SymmetryOperation`: an affine operation on fractional coordinates, p' = R p + t, R an integer 3 × 3 matrix, t a
  translation; parsed from and printed as the CIF `x,y,z` form (`-x,y+1/2,-z`); apply, compose, identity test, the
  rotation part's order and determinant.
- `PointSymmetry`: a finite set of integer rotation matrices acting on fractional coordinates; the Laue group of a
  structure (the rotation parts of its operations, closed under inversion); its `LaueClass`; the orbit of (h, k, l)
  under it and the orbit's canonical form (section 3.6).
- `AtomSite` (Z, fractional position, occupancy, the CIF label when it has one), `Element` (symbols 1 ... 118).
- `CrystalStructure`: a cell, its sites, its symmetry operations (from the file, or detected from the sites, section
  3.4), its Laue group, the IT space-group number when the file states one, and the validation issues the apps'
  `CrystalModel` reports (cell lengths and angles, volume, occupancies, finite coordinates, the 0.5 Å close contact).
- `CIFImport`: the apps' parser, generalised: every crystal system is admitted; the refusals that protect the structure
  (section 3.7) stay; the two app-specific refusals go (no point-group refusal; no scattering-factor refusal, which
  moves to slice 2).
- `ReflectionList`: every (h, k, l) with 0 < |g| ≤ k_max, each with g, |g|, d and whether the space-group operations
  allow it (section 3.8); in increasing |g|. `ReflectionFamily`: the list grouped into rings by |g| (section 3.9),
  each with its members, their count (the multiplicity), its orbits by canonical form and its name.

Not in slice 1: scattering factors, structure factors, the |F| cut and kinematic intensities (slice 2); the Materials
Project importer; pattern geometry, calibration, zone axes, indexing (`DiffractoGeometry`); orientation
(`DiffractoOrientation`). The three modules exist as empty products with a doc comment naming their scope.

## 2. The reference

pymatgen 2026.9.24, the release tag `v2026.9.24` = commit `8461938451eedd7591b96998de4b48e4b1901ea0` (MIT), installed
from PyPI into the gitignored `References/pymatgen-venv` (Python 3.14.7) with its own dependencies as resolved on
2026-10-14 (spglib and NumPy versions recorded in `References/pymatgen-venv-install.log` and in the pins file), and the
same tag cloned into `References/pymatgen-8461938` so the pins below can cite file:line. The oracle runs as
`References/pymatgen-venv/bin/python -I -B`. An unpinned copy is never a reference. Agreement with pymatgen is agreement,
not truth; truth is section 5.

What of pymatgen is used, and for what:

- `pymatgen.io.cif.CifParser(...).parse_structures(primitive=False)`: the cell and the expanded sites of a CIF.
- `pymatgen.core.Lattice.from_parameters(..., vesta=True)` (a along x, b in xy: the apps' convention), `.matrix`,
  `.metric_tensor`, `.volume`, `.reciprocal_lattice_crystallographic` (no 2π), `.d_hkl`, `.get_points_in_sphere`.
- `pymatgen.symmetry.analyzer.SpacegroupAnalyzer(structure, symprec=1e-4)`: `get_space_group_number()`,
  `get_point_group_symbol()`, `get_symmetry_operations(cartesian=False)` (spglib behind it).
- `pymatgen.core.operations.SymmOp.from_xyz_str` / `.as_xyz_str()`: the `x,y,z` grammar.
- The Laue orbit and the systematic-absence test are computed in the pin generator from pymatgen's operations with the
  definitions of section 3, in NumPy, not by a pymatgen method (pymatgen has no reciprocal-space orbit function); they
  are written from this text, not from the Swift.

## 3. The definitions, pinned

### 3.1 The lattice

Real lattice rows (the apps' `Crystal.latReal`, py4DSTEM `Crystal.__init__`): a₁ = (a, 0, 0), a₂ = (b cos γ, b sin γ, 0),
a₃ = (c cos β, −c (cos β cos γ − cos α) / sin γ, V / (a b sin γ)) with V = abc √(1 − cos²α − cos²β − cos²γ + 2 cos α cos β
cos γ). Metric G = A Aᵀ; reciprocal rows B = G⁻¹ A so aᵢ · bⱼ = δᵢⱼ; g(hkl) = h b₁ + k b₂ + l b₃; d = 1/|g|. The volume
is |det A| and must equal the closed form to 1e-12 relative (a self-check in the tests). The reciprocal metric
G* = G⁻¹ gives |g|² = (h k l) G* (h k l)ᵀ: `d(h, k, l)` is computed from G*, and `gLength` from the vector; the two agree
to 1e-12 relative (another self-check).

### 3.2 The crystal system of a metric

At a tolerance ε (default 1e-9 relative on lengths, 1e-9° on angles; the apps' indexing uses the same; a looser one is
the caller's choice): cubic a = b = c, α = β = γ = 90°; hexagonal a = b, α = β = 90°, γ = 120°; tetragonal a = b ≠ c,
all 90°; rhombohedral (trigonal, rhombohedral axes) a = b = c, α = β = γ ≠ 90°; orthorhombic a ≠ b ≠ c, all 90°;
monoclinic two angles 90°, one not (the unique axis is the one whose angle is not 90°: b when β ≠ 90°); triclinic
otherwise. A trigonal structure in hexagonal axes has a hexagonal metric: the system is then decided by the symmetry
(section 3.4), never by the metric alone, which is the apps' lesson (quartz, calcite). Tested in this order, first that
applies.

### 3.3 Symmetry operations

Parsed from the CIF grammar as the apps do (`CIFImport.parseSymmetryOperation` at mac4DSTEM 7448401b): three
comma-separated components, each a sum of signed terms ± n/m·{x, y, z} or ± n/m; the rotation part is required to be
an integer matrix (a coefficient that is not an integer is a refusal); translations are kept as doubles, compared mod 1
within 1e-9. The identity is x,y,z with whole-number translations. The order of an operation is the smallest n with
Rⁿ = I; the determinant is ±1. Distinct operations: rotation exact, translation mod 1 within 1e-3 (the apps' rule for
the incomplete-list check, kept).

### 3.4 The symmetry of a structure

(a) A file with an operation loop: the structure's operations ARE that loop (after the apps' completeness check,
section 3.7). (b) A structure given as a site list (a P1 file, or a cell built in code): its operations are detected.
The candidate rotations are the holohedry of the metric's system (section 3.2) in its conventional setting: the closure
of the generators {−I} ∪ {C₂ ∥ b} (monoclinic, unique axis b; C₂ ∥ c or ∥ a when that is the unique axis), {C₂ ∥ a, C₂
∥ b, C₂ ∥ c} (orthorhombic), {C₄ ∥ c, C₂ ∥ a} (tetragonal), {C₆ ∥ c, C₂ ∥ a} (hexagonal, also for trigonal in hexagonal
axes), {C₃ ∥ [111], C₂ ∥ [1-10]} (rhombohedral axes), {C₃ ∥ [111], C₄ ∥ c} (cubic), {} (triclinic); each with −I. A
candidate R is a symmetry when some translation t makes p ↦ R p + t a permutation of the site set preserving Z and
occupancy (the apps' `isSymmetry`: t from mapping site 0 onto every compatible site, positions within 5e-3 fractional,
occupancies within 0.02, targets distinct). Every such (R, t) with t mod 1 is an operation, so a centred cell yields its
centring translations (R = I, t = (½, ½, 0)) and a screw axis its translation. The set found is closed (a self-check).
The Laue group is {±R} over the rotation parts. The position tolerance 5e-3 is fixed on purpose (the apps' reason: a
file must not widen its own gate); when a candidate fails at 5e-3 but passes at 10 × the coarsest written half-step, the
structure reports `symmetryUnresolvableAtWrittenPrecision` for that candidate and the operation is not taken.

### 3.5 Laue class from a group

By the order and the axis orders present: 2 → -1; 4 → 2/m; 8 with a 4-fold → 4/m, else mmm; 16 → 4/mmm; 6 → -3; 12
with a 6-fold → 6/m, else -3m; 24 with a 6-fold → 6/mmm, else m-3; 48 → m-3m. The order of a rotation is read from its
matrix, so the setting does not matter. The table is checked against pymatgen's point-group symbol for every pinned
structure, and against the IT-number ranges of section 3.7 for the declared groups.

### 3.6 Orbits and the canonical form

The orbit of (h, k, l) under a Laue group is {(h k l) R : R in the group} (hkl is covariant: the row vector times the
matrix that acts on fractional coordinates). The canonical form is the member that is greatest under the key (number of
non-negative indices, then lexicographic (h, k, l)). This one rule reproduces every name the apps use today, which is
why it is chosen over the plain lexicographic maximum: cubic (|h|, |k|, |l|) sorted descending ({3 1 1}); hexagonal the
Miller–Bravais form with two non-negative in-plane indices and H ≥ K ({1 1 -2 0}, {1 0 -1 1}; the plain maximum would
name it {2 -1 -1 0}); monoclinic macTEM's S2c form (|h|, |k|, L) with L = l · sign(h), |l| when h = 0 ({2 0 -1} and
{2 0 1} are different rings in β″, 5.6555 and 4.3246 Å). The name prints the canonical form in braces, negatives as
"-1", the hexagonal-axes classes (-3, -3m, 6/m, 6/mmm) in four Miller–Bravais indices with i = −(h + k). Orbits under
4/m, 6/m, -3 and m-3 are smaller than their holohedry's and their names keep the distinction ({2 1 0} and {1 2 0} are
two orbits in 4/m): the rule is the same, the group is smaller.

### 3.7 CIF import, what is carried from the apps and what changes

Carried verbatim in behaviour (mac4DSTEM `CIFImport.swift` at 7448401b, 1219 lines; the five-line tokenizer rules, the
`(3)` uncertainty suffix, the element symbol from `_atom_site_type_symbol` else the label with trailing digits and
charges stripped, `_atom_site_occupancy` default 1, the written-precision half-step 10⁻ᵈ capped at 1e-2 with dyadic
rationals exact, the expansion that merges images within the propagated half-step AND 0.10 Å, the refusals:
`missingTag`, `unparsableNumber`, `noAtomSites`, `unknownElementSymbol`, `unsupportedSymmetryOperation`,
`missingSymmetryOperations` (a non-P1 declaration with no expanding operations), `truncatedLoop`, `multipleStructures`,
`symmetryOperationsOutsideStructure`, `incompleteSymmetryOperations` (the point-group-order table for IT 1 ... 230 ×
the centring, with the primitive-cell exception by angles), `invalidModel`), and the display-name rule (`VESTA_phase_n`
falls back to the file name). The error messages lose "mac4DSTEM" and "ACOM" (the apps' recorded wart).

Changed: `unsupportedPointGroup` and `symmetryNotSupportedByStructure` are gone. Every system is admitted; a trigonal
structure in hexagonal axes imports with Laue class -3 or -3m, not as a refusal and not as 6/mmm. The symmetry comes
from section 3.4; `symmetryUnresolvableAtWrittenPrecision` survives as the report on a candidate operation, not a
refusal of the file. `unsupportedElement` is gone from the importer (every symbol 1 ... 118 resolves); slice 2's
scattering layer reports what its table lacks.

### 3.8 Reflections and systematic absences

Every (h, k, l) ≠ 0 with |g| ≤ k_max, enumerated with the apps' per-axis bound |h| ≤ ⌈k_max |a₁|⌉ etc. (the DEVIATION
of 2026-09-14, complete for every cell; the old ten-direction bound lost reflections in oblique cells). k_max must be
finite, positive and at most 4 Å⁻¹ (macTEM's cap). Order: increasing |g|; equal |g| (within 1e-12 relative) by
descending lexicographic (h, k, l), so the list is deterministic. A reflection is **allowed** unless some operation
(R, t) of the structure has (h k l) R = (h k l) and (h k l) · t not an integer (within 1e-9): the general reflection
condition. For a structure whose operations were detected (3.4b) this reproduces the centring and glide/screw rules of
the operations found; a reflection that the operations allow but whose structure factor is zero by the basis (diamond
{2 2 2}, {4 4 2}) is allowed here and is slice 2's to score.

### 3.9 Families (rings)

macTEM's S2 rule (ADR 022, `DSpacingIndexing.families` at a58c84c), on the allowed reflections: swept in increasing
|g|, a reflection joins the current family while |g| ≤ g₀ (1 + 1e-9), g₀ the family's first; the family's g is its
smallest member's |g|; its members (Friedel mates included), their count (the multiplicity), its orbits (the distinct
canonical forms, descending) and its name (the orbits' names joined by " + "; Si {3 3 3} + {5 1 1} and ideal-hcp
{3 0 -3 2} + {0 0 0 6} are one ring each). A family built on the full list instead (allowed and absent alike) is
offered for a caller that wants the lattice rings.

## 4. Pins (agreement with pymatgen)

A generator, `tools/slice-1/pins.py` (run in `References/pymatgen-venv`), writes `Tests/DiffractoCoreTests/Fixtures/
slice-1-pins.json` before the Swift code exists, for every structure of section 5.3 and the parameter cells of 5.1:

- the cell as pymatgen reads it (a, b, c, α, β, γ), the lattice matrix (vesta=True), the metric tensor, the volume,
  the reciprocal rows;
- the expanded sites (Z, fractional, occupancy) from `CifParser`;
- `get_space_group_number`, `get_point_group_symbol`, and every symmetry operation (rotation, translation) from
  `SpacegroupAnalyzer` at symprec 1e-4;
- every (h, k, l) with |g| ≤ k_max (1.2 Å⁻¹; β″ also at 1.6 and 2.0) from `get_points_in_sphere` of the
  crystallographic reciprocal lattice, with |g| and `d_hkl`;
- for each, its Laue orbit (the set), its canonical form by section 3.6's rule written in Python, and whether it is
  allowed by section 3.8's rule applied to pymatgen's operations;
- the families of section 3.9 on the allowed reflections: g, members, multiplicity, orbits, name.

Agreement required: cell parameters equal to 1e-12 relative; matrix elements, metric, volume and reciprocal rows within
1e-12 · max(a, b, c) (or relative 1e-12 where non-zero); the expanded site multiset equal (same count; each site matched
on Z, occupancy within 1e-9, position equal mod 1 within 1e-9); space-group number equal; the Laue class equal to the
class of pymatgen's point-group symbol; the operation set equal (same count; rotations equal as integers; translations
equal mod 1 within 1e-6, spglib's own rounding); the (h, k, l) set equal; |g| and d relative 1e-12; orbit sets and
canonical forms equal; allowed flags equal; families equal (count, g to 1e-12, members, multiplicity, orbits, name). A
difference is a failure of the port, never a tolerance to widen. The pins carry the versions (pymatgen, spglib, NumPy,
Python) and the generator's own sha256 of each CIF text.

Where pymatgen and the registered rule are known to differ, the pin is not taken and the case is planted truth instead:
`CifParser` merges sites within a fixed 1e-4 fractional (`site_tolerance`), ours within the written precision, so
coarse-coordinate files (section 5.4) are planted cases; `SpacegroupAnalyzer` works in Cartesian Å at symprec, ours in
fractional at 5e-3, so the detection pins are taken on exactly written structures only.

## 5. Truth (independent of pymatgen)

### 5.1 Closed-form d-spacings, one cell per system (hand formulas, rational parameters)

- cubic a = 4: d = a / √(h² + k² + l²);
- tetragonal a = 3, c = 5: 1/d² = (h² + k²)/a² + l²/c²;
- orthorhombic 3, 4, 5: 1/d² = h²/a² + k²/b² + l²/c²;
- hexagonal a = 3, c = 5: 1/d² = (4/3)(h² + hk + k²)/a² + l²/c²;
- rhombohedral a = 4, α = 70°: 1/d² = [(h² + k² + l²) sin²α + 2(hk + kl + hl)(cos²α − cos α)] / [a² (1 − 3cos²α +
  2cos³α)];
- monoclinic a = 15.16, b = 4.05, c = 6.74, β = 105.3° (β″'s cell): 1/d² = (h²/a² + l²/c² − 2hl cos β/(ac))/sin²β +
  k²/b²; the hand anchors from macTEM's S2c: d(200) = a sin β/2 = 7.31135, d(001) = c sin β = 6.50112, d(020) = 2.025,
  d(110) = 3.90306, d(201) = 4.32457, d(20-1) = 5.65554 Å;
- triclinic a = 4, b = 5, c = 6, α = 80°, β = 85°, γ = 95°: the general formula 1/d² = (1/V²)[S₁₁h² + S₂₂k² + S₃₃l² +
  2S₁₂hk + 2S₂₃kl + 2S₁₃hl] with S₁₁ = b²c² sin²α, S₂₂ = a²c² sin²β, S₃₃ = a²b² sin²γ, S₁₂ = abc²(cos α cos β − cos γ),
  S₂₃ = a²bc(cos β cos γ − cos α), S₁₃ = ab²c(cos α cos γ − cos β).

Every (h, k, l) with indices in −3 ... 3: the metric-tensor d of section 3.1 equals the formula to 1e-12 relative.
Every other system's formula is also the triclinic one specialised, a second self-check.

### 5.2 Orbit sizes and names (group theory, no reference)

General position multiplicities: -1 2, 2/m 4, mmm 8, 4/m 8, 4/mmm 16, -3 6, -3m 12, 6/m 12, 6/mmm 24, m-3 24, m-3m 48.
Special: m-3m {100} 6, {110} 12, {111} 8, {210} 24; 6/mmm {10-10} 6, {0002} 2, {10-11} 12, {11-20} 6; 2/m (b unique)
{h0l} 2, {0k0} 2, {hkl} 4, {0kl} 4; 4/m {210} 8 and {120} 8 as two orbits, one orbit of 16 under 4/mmm. The holohedry
of each system has the registered order (2, 4, 8, 16, 12 rhombohedral / 24 hexagonal, 48) and is closed under
multiplication, every element of determinant ±1 and finite order. The canonical forms named in section 3.6 for the apps'
cases: {3 1 1}, {1 1 -2 0}, {1 0 -1 1}, {2 0 -1}, {2 0 1}, {0 2 0}.

### 5.3 Structures with known cell contents (the CIF texts are committed fixtures)

Monoclinic first: **β″** from the pinned generator's values (mac4DSTEM `tools/crystal-structures/
make_beta_double_prime.py` at 7448401b; Andersen et al. 1998, Acta Mater. 46, 3283; C2/m, 8 operations, six sites on
y = 0): 22 atoms, Mg10Si12; Laue 2/m; h + k odd absent; at k_max 1.2 Å⁻¹ **1448** allowed reflections (gemmi's count in
macTEM's S2c record; the apps' 1442 is after the |F| cut, which is slice 2's) and **419** families (S2c: "with |F| >
1e-9 419"; 417 + {3 3 3} + {14 0 1}); the first families {2 0 0} 7.3113 (2), {0 0 1} 6.5011 (2), {2 0 -1} 5.6555 (2),
{2 0 1} 4.3246 (2), {1 1 0} 3.9031 (4), {4 0 0} 3.6557 (2), {4 0 -1} 3.6207 (2), {1 1 -1} 3.4554 (4), {2 0 -2}
3.3124 (2), {0 0 2} 3.2506 (2) (d in Å, multiplicity). The same cell as a P1 list of the 22 atoms detects the 8
operations and the same Laue class, absences and families.

Then, from the apps' test fixtures (their CIF texts at 7448401b, expectations re-derived here for the general case):
gold P1 4 sites (detected: Fm-3m's 192 operations, m-3m, mixed parity absent, 8 families at 1.2 Å⁻¹: {1 1 1}, {2 0 0},
{2 2 0}, {3 1 1}, {2 2 2}, {4 0 0}, {3 3 1}, {4 2 0}); gold with one site and the four centring operations (4 sites);
silicon Fd-3m (8 sites, m-3m; {2 0 0}, {4 2 0}, {6 0 0} absent by the d-glide and F-centring; {2 2 2} and {4 4 2}
allowed by the operations, so the allowed families at 1.2 are the apps' 10 plus those two: 12); magnesium P6₃/mmc
(2 sites, 6/mmm, {0 0 0 1} absent, 24 families at 1.2 as macTEM counts, {1 0 -1 0} 6, {0 0 0 2} 2, {1 0 -1 1} 12,
{1 1 -2 0} 6); rock salt Fm-3m with its full 192 operations (8 sites, 225); α-quartz P3₂21 (9 sites; Laue -3m, NOT
6/mmm; 00l with l ≠ 3n absent); calcite R-3c in hexagonal axes (30 sites; -3m; −h + k + l ≠ 3n absent); fluorapatite
P6₃/m (42 sites; 6/m, not 6/mmm; {2 1 0}-type orbits of 12, not 24); pyrite Pa-3 (12 sites; m-3, not m-3m; {2 1 0} and
{1 2 0} two orbits of 12); 2H-WS₂ (6 sites; 6/mmm); CdI₂ P-3m1 (3 sites; -3m); the apps' "tetragonal a = c" cell (its
Laue class is read from its atoms by section 3.4, not refused); ideal-c/a hcp (c = a √(8/3)): {3 0 -3 2} + {0 0 0 6}
one ring of 14. Every site count above is the cell content the file's formula states.

### 5.4 Planted refusals and precision cases (the apps' tests, same expectations where the rule is unchanged)

Missing cell tag; unparsable number; no atom loop; unknown element; unparsable operation; non-P1 declaration without
operations (named); identity-only loop with a non-P1 number; truncated loop; two structure blocks; a symmetry-only block
before or after the structure; Fm-3m listing only its 48 point operations (refused, 192 required); a primitive cell of
a centred group not refused; duplicated rows not counted; a garbled 1e21 coefficient refused, not trapped; gold with a
split site (the symmetry map must be a permutation: Laue class not m-3m); magnesium at 3 decimals admitted (6/mmm);
gold at 2 decimals (`0.50` is exact, admitted); a coarse ⅓ (`0.33`) merged by the written precision; two distinct sites
0.3 Å apart kept and then refused by the close-contact rule; the standard-uncertainty suffix stripped; VESTA block names.

## 6. Gates and predicted outcomes (registered now)

1. Pins (section 4): all equal within the stated bounds, on every structure of 5.3 and cell of 5.1.
2. Truth (section 5): every closed form to 1e-12 relative; every count, class, multiplicity, name and absence exact;
   predicted β″ 1448 / 419, Si 12 families, Au 8, Mg 24.
3. Every new test broken once by a mutant aimed at its rule (planned: the reciprocal lattice as G A instead of G⁻¹ A;
   β read from α; the per-axis tile from the shortest reciprocal direction; the orbit as hkl R⁻¹ for a non-symmetric R
   (hexagonal 6-fold); the canonical key without the non-negative count (hexagonal names change); the mirror dropped
   from 2/m (orbits halve); the C₂ about a instead of b; the absence test without the "R fixes h" guard; the family
   tolerance at 1e-4 and at 0; the multiplicity as the orbit count; the dedup tolerance 0.5; the centring translation
   skipped in the expansion (β″ 11 sites); "1/2+x" read as "x"; the Laue table's 8-with-4-fold branch inverted; the
   detection's target-distinct guard dropped (the split gold passes as m-3m)).
4. Self-checks in the tests: V by both formulas, d by G* and by the vector, every holohedry closed, every detected set
   closed.
5. One Gate B refuter on the frozen slice-1 commit (a second only on a named doubt); Gate D only if a discrepancy's
   cause is unproven. The refuter is asked to compare the port with pymatgen on structures the fixtures do not hold
   (other space groups, other settings: unique axis c, rhombohedral axes, origin choice 2) and to try to break the
   detection with near-symmetric cells.
6. `tools/test.sh` green with the stamp; a cold build with zero warnings; the consumer-check skeleton compiles against
   the public API from a separate package.

## 7. What is not claimed

Agreement with pymatgen and recovery of planted cells do not make a d-spacing right on a real pattern: that needs a
calibration (`DiffractoGeometry`, later) and is the apps' validation to register. No intensity is computed here. The
symmetry detection of section 3.4 is a holohedry test at a fixed fractional tolerance, not spglib: it finds the
operations of a structure written to ordinary precision in its conventional cell; a structure in a non-conventional
cell (a monoclinic cell with its unique axis along c is conventional too and is handled; a rhombohedral structure given
in a primitive triclinic-looking cell is not) detects only what the metric's holohedry contains, and says so by its
Laue class. The family rule's 1e-9 is macTEM's; near-coincident rings are the caller's to judge (ADR 022). Not
claimed (Gate B, 2026-10-15): the symmetry of a P1 cell in a non-conventional setting (no lattice reduction yet);
symmetry finer than the 5e-3 fractional tolerance in a large cell; the speed of a large symmetric P1 cell.

## 8. Amendments

### 2026-10-14, before any Swift code: the oracle's package is pymatgen-core

Found while provisioning `References/`. pymatgen split in 2026-09: the PyPI package `pymatgen` 2026.9.24 (tag `v2026.9.24` =
`8461938`, cloned to `References/pymatgen-8461938`) is an umbrella that depends on `pymatgen-core>=2026.9.23`, and the code
section 2 uses (`pymatgen.core.Lattice`, `pymatgen.io.cif.CifParser`, `pymatgen.symmetry.analyzer.SpacegroupAnalyzer`,
`pymatgen.core.operations.SymmOp`) lives in pymatgen-core (the monorepo's `pymatgen-core` submodule,
`github.com/materialsproject/pymatgen-core`). pip first resolved pymatgen-core 2026.10.2 (12 days old); the venv is pinned to
**pymatgen-core 2026.9.23** (tag `v2026.9.23` = `ae9a8348177b92fe33b214e54e61a3b8f4244b6e`, cloned to
`References/pymatgen-core-ae9a834`), the version the umbrella tag was released against; the installed `core/lattice.py`,
`io/cif.py`, `symmetry/analyzer.py` and `core/operations.py` are byte-identical to that clone (sha256 compared 2026-10-14).
File:line pins cite `pymatgen-core-ae9a834/src/pymatgen/...`. spglib 2.8.0, NumPy 2.5.4 (`References/pymatgen-venv-freeze.txt`).
No registered number changes.

### 2026-10-14, at the first test run: a hand error in section 5.2's tetragonal multiplicities

Section 5.2 wrote "4/m {210} 8 and {120} 8 as two orbits, one orbit of 16 under 4/mmm". The group code, written from
section 3.6's definition, gives 4, 4 and 8, and that is the International Tables' value: hk0 is a special form, fixed
by the mirror perpendicular to c (in 4/m the inversion composed with the 2-fold about c), so its multiplicity is half
the general one: 4 under 4/m and 8 under 4/mmm; the general hkl are 8 and 16 as registered. The two-orbit statement
stands ({2 1 0} and {1 2 0} are two orbits of 4/m, one of 4/mmm). The test asserts the corrected numbers and adds the
general (2 1 3) at 8 and 16. The rule did not change; a planted number was wrong by hand.

### 2026-10-14, after the pin generator ran, before the pins test: what pymatgen and the fixtures found

Found by `tools/slice-1/pins.py` (6960 self-checks passed) on 19 structures. The β″ predictions held exactly (1448
allowed reflections, 419 families, the first ten rings' names, d and multiplicities); Au 8, Si 12 and the ideal-hcp
ring of 14 held. Seven readings are settled here; the rules of sections 3.1, 3.6, 3.8 and 3.9 do not change.

1. **Section 3.5 counts the orders of the proper rotations** (determinant +1). Read over every element, −3 has order 6
   and -3m reads as 6/m, m-3 as 6/mmm; the Swift was written over the proper rotations from the start.
2. **Section 3.4 is amended: the operations of a structure are always detected from its sites.** A file's operation
   loop expands the asymmetric unit (3.7's checks unchanged) and is kept as `fileOperations`; the structure's
   `operations`, Laue group, absences and orbits come from the detection of 3.4b on the expanded cell, with the loop's
   rotations added to the holohedry's candidates (so a non-conventional setting is not under-detected). Reason: a loop
   may legitimately declare less symmetry than the atoms have (the apps' gold fixture lists one atom and the four
   centring translations, a legal CIF whose loop alone reads as Laue -1, while its atoms are m-3m; spglib finds 192
   operations), and the apps' importer already confirmed the family from the atoms, never from the loop. The pins
   record whether the loop equals or is a subset of the detected set (`file_operations`); the port reports
   `fileOperationsAreSymmetries`. The β″ loop equals its detected set.
3. **Detected translations are snapped to the nearest k/24** when within the position tolerance (5e-3): a translation
   found from coordinates written to six or seven digits carries their error (⅔ written as 0.666667), and the absence
   test of 3.8 needs h · t exact. Every translation of a space-group operation in a conventional or centred setting is
   a multiple of 1/24 (halves, thirds, quarters, sixths, eighths). pymatgen does the same (`limit_denominator(1000)`);
   the pins carry k/24.
4. **pymatgen's `CifParser` snaps a coordinate within 1e-4 of ⅓ or ⅔ to the exact fraction** before expanding
   (`pymatgen-core-ae9a834/src/pymatgen/io/cif.py`, about line 610). The port keeps the written value, as the apps do.
   For magnesium, fluorapatite, WS₂ and CdI₂ the pinned sites therefore differ from the written ones by up to 3.3e-7;
   the site pin compares the port against `sites_as_written`, and `site_snap_max` records pymatgen's displacement.
   Symmetry, absences and families are identical either way.
5. **Magnesium has 25 symmetry-allowed families at 1.2 Å⁻¹, not 24**, and ideal hcp 24, not 23: {3 0 -3 1} is allowed
   by P6₃/mmc's operations and has F = 0 by the basis (unit scatterers: 24 survive), exactly as Si's {2 2 2} and
   {4 4 2}. macTEM's 24 is after the |F| cut, slice 2's. Section 5.3's "24 families as macTEM counts" reads 25.
6. **Names of a multi-orbit ring follow the descending lexicographic order of the orbits**, as 3.9 says: Si's ring is
   `{5 1 1} + {3 3 3}` (macTEM's `orbits` are `[(5 1 1), (3 3 3)]` too); 5.3's prose had them the other way round.
7. **The apps' "tetragonal a = c" fixture** (the P1 two-site cell) is Amm2 by its atoms, Laue mmm, not tetragonal:
   section 5.3 said only that its class is read from its atoms; it is. The seeded quartz, calcite and pyrite texts of
   `cif-symmetry-test` carry incomplete operation loops (the apps used them to show a refusal, not a cell content): the
   fixtures here are completed to P3₂21 (6 operations, 9 atoms), R-3c (36, 30) and Pa-3 (24, 12), marked `# MODIFIED`
   with the seeded texts kept beside them (`*_as_seeded.cif`, pinned as what they are: 12, 5 and 8 atoms).

### 2026-10-15, after Gate B: the refuter's findings and the fixes they require (tier 1)

One fresh Sonnet refuter on the frozen `b4cfe17` (its scripts and CIFs in the session's scratchpad, `refuter/`). It
confirmed agreement with pymatgen and spglib for all 230 space groups in their standard settings (two atoms at generic
coordinates, as a full loop and as P1: 458 of 460 comparisons; the two others are F6 below), for β″ in the unique-axis
b, c and a settings and in I2/m (8 operations, 1448, 419 each), and for the edge inputs of section 5.4. It found, with
repros (now fixtures `refuter_*.cif` and `RefuterRegressionTests`; truth from pymatgen on the UNSHIFTED structures, an
origin shift changes no absence):

- **F1, port defect (silent wrong numbers): translations off the standard origin.** (a) Amendment 3's snap to k/24
  within 5e-3 moved genuine off-grid translations (a shifted origin, a lone atom at a general position): operations lost
  or merged, sets that were not groups, wrong absences (rock salt shifted by (0.9172, 0.6292, 0.5141): 160 operations,
  no Laue class, 24 flags wrong; one Fe atom in a 5 Å cube: 43 operations, ring {1 1 0} of 10). (b) The absence test's
  1e-9 failed on translations carrying the coordinates' rounding (Mg shifted, 8 decimals: 262 allowed of 282; a 7×
  supercell's k/7 translations: 28 allowed of 122). **Fix:** each translation is refined by least squares over the
  matched sites; it is replaced by the file's exact translation when the loop holds that operation, else by the k/24
  grid point only when within three times its own uncertainty (residual spread + 2 × written half-step); the final
  translation is re-verified as a permutation (falling back to the refined, then the trial one); the absence test
  allows (|h| + |k| + |l|) × the structure's largest translation uncertainty + 1e-9, capped at 0.04 (below 1/24): for
  an operation that fixes h, h · t does not depend on the origin. **Amendment 3 above is replaced by this.**
- **Found while fixing F1:** detection must return a group. Pseudo-symmetry within the fractional tolerance passed
  single operations whose products do not pass (the apps' split-atom gold: 40 operations with m-3m's rotations).
  **Fix:** operations that compose out of the set (within 5e-3) are dropped until it is closed (sets of at most 2048
  operations; a larger one, a big symmetric supercell, is returned as found). The split gold is mmm.
- **F2, port defect: co-sited rows merged.** Two rows on one position (Fe²⁺ 0.5 and Fe³⁺ 0.5) under a loop merged into
  one (cell content Fe 2 instead of 4). **Fix:** only images of the same row merge; co-sited rows stay, and at full
  occupancy are refused by the close-contact rule.
- **F3, registration gap, recorded in section 7 (not fixed in slice 1):** a P1 cell in a non-conventional setting (a
  primitive cell of a centred lattice, unique axis a for a tetragonal structure, γ = 60° for a hexagonal one) is tested
  only against its own metric's holohedry and detects a lower Laue class, so its absences and ring names follow the lower
  class. The Laue class says so; a lattice reduction (Niggli or Delaunay, then the conventional cell) is slice 1b's.
  With the file's loop the rotations join the candidates and these cells agree.
- **F4, defect: Miller–Bravais names in rhombohedral axes.** **Fix:** four indices only for a hexagonal-axes Laue class
  in a hexagonal metric (`ReflectionList.usesMillerBravais`); corundum in rhombohedral axes keeps three.
- **F5, traps and resources.** A finite cell of 1e30 Å trapped (`Int(Double)`), 1e10 Å hung, 100 Å at 4 Å⁻¹ would need
  about 68 GB. **Fix:** a list is refused (`cellTooLargeOrInvalid`) when the cell is invalid or the expected count
  (4/3)π k_max³ V exceeds 10 million; a rotation with determinant other than ±1 (`x,x,z`) is refused
  (`notInvertible`). Not fixed, recorded: importing a large SYMMETRIC P1 cell scales about as N³ (500 atoms, 30 s).
- **F6, by design, recorded:** the 5e-3 fractional tolerance is 0.25 Å in a 50 Å cell (spglib's is 1e-4 Å), so a
  pseudo-centring of 0.013 Å reads as symmetry (sg 54 with x = 0.00247); a Cartesian tolerance is a candidate for
  slice 1b. **Changed:** the metric family for detection is proposed at 1e-9 (section 3.2's default) instead of the
  apps' 1e-3 / 0.05°, because a rotation must preserve every |g| to the ring rule's 1e-9 or an orbit splits across two
  rings (a cell with c/a = 1.0001 read m-3m with every ring named twice; now 4/mmm, 16 operations, as spglib).

### 2026-10-14, 23:30 CEST, from the coordinator: parsing comparisons named, no change of scope

No Swift crystallography package exists (d-spacings, structure factors, indexing): DiffractoKit is the first. A Swift
CIF parser does: SKCIFParser inside iRASPA (MIT), vendored by the macOS Quick Look extension CIFPreview
(`github.com/grebenyyk/CIFPreview`), with symmetry expansion. It and Gemmi (C++/Python, reads COD and small-molecule
CIFs) are **parsing comparisons** for section 3.7 (which files parse, which symmetry lists expand), never a numbers
oracle; pymatgen stays the oracle for cells, symmetry, d and F. Neither is run in this session; a later session may
add them as a tier-2 check of the parser on database CIFs. Dependencies: none beyond Foundation and simd in slice 1;
when a linear-algebra path comes, Apple's swift-numerics and Accelerate directly, never the small single-maintainer
wrappers (LASwift, LANumerics, Surge). No registered number changes.

### 2026-10-14, at the first test run: two parser readings tightened (section 3.3)

Found by the planted refusals. (1) The apps' fraction parser split "/2" into ["2"] and read `x/2` as x + 2; the kit's
refuses an empty numerator or denominator as unparsable. (The coordinator relayed the same reading at 23:40 CEST as
macTEM's refuter finding F5 on its own CIF lane, upstream in mac4DSTEM; the planted cases here hold `x/2`, `0.5x`,
`1/3x` refused and `-y+1/2`, `x-y`, `2x`, `-z+1/3` parsed and round-tripped; a non-integer coefficient is refused
whatever pymatgen would read into such a string, since no space-group operation has one.) (2) A coefficient or constant that is non-finite or beyond
1000 in magnitude (`1000000000000000000000x`, which the apps' incomplete-list check guards against trapping on) is
refused as unparsable at the parser, never converted to Int. No registered number changes.
