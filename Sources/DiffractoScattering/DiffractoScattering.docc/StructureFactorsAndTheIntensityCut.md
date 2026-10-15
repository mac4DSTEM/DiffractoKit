# Structure factors and the intensity cut

How a reflection's F(hkl) is computed, which reflections are kept, and what a ring's intensity means.

## The structure factor

For a reflection hkl, ``StructureFactors/compute(_:in:formFactor:)`` returns a ``StructureFactor``:

F(hkl) = (1/V) sum over sites of f(Z, |g|^2) * occupancy * exp(2 pi i (h x + k y + l z))

with f the electron scattering factor in angstroms and V the cell volume in cubic angstroms, so F is in A^-2 and
``StructureFactor/intensity`` = |F|^2 is in A^-4. The default form factor is `ElectronScatteringFactors.f`; any
``FormFactor`` can be passed in its place. The function returns nil when a site's element has no scattering factor.

## The cut

``DiffractoCore/CrystalStructure/kinematicReflections(kMax:cut:formFactor:)`` lists every reflection with
0 < |g| <= k_max, computes F for each, and keeps those with |F| above `cut`. The default,
``StructureFactors/defaultCut`` = 1e-4 A^-2, is a fixed number against a quantity that scales as 1/V: a large cell
keeps fewer reflections at the same cut.

```swift
import DiffractoCore
import DiffractoScattering

let list: KinematicReflectionList = try gold.kinematicReflections(kMax: 1.2)
print(list.all.count, list.kept.count)          // every reflection, and those above the cut
for ring in list.families.prefix(3) {
    print(ring.family.name, ring.family.d, ring.family.multiplicity, ring.intensity)
}
```

``KinematicReflectionList/all`` holds every reflection of the Core list with its F, ``KinematicReflectionList/kept`` those
above the cut, and ``KinematicReflectionList/families`` the rings of the kept ones. A ``KinematicFamily`` carries the
Core `ReflectionFamily` and ``KinematicFamily/intensity``, the sum of |F|^2 over its members: the ring's
kinematic weight, comparable between rings of one structure, not an absolute count.

A structure that fails validation (a non-finite coordinate, an occupancy out of range, a close contact), an element the
table lacks, an invalid cut or a non-finite F throws a ``KinematicError`` before any partial result is returned.
Reflections kept by |F| alone although the detected operations forbid them (written coordinates that break the symmetry
slightly) are named by ``KinematicReflectionList/keptButOperationAbsent``.

## Memory: `retainingAll`

By default ``KinematicReflectionList/all`` holds every reflection before the cut, which for a large cell is most of the
memory the call uses (a 30 angstrom cell at k_max 2.5 holds about 1.7 million). Pass `retainingAll: false` to
`kinematicReflections` to compute F as each reflection is enumerated and keep only those above the cut: `all` is then empty,
the list before the cut is never held, and ``KinematicReflectionList/kept`` and ``KinematicReflectionList/families`` are
exactly the values the default gives.

## Polarity cannot be read from a kinematic pattern

A kinematic pattern cannot tell +c from -c on a polar structure. Without anomalous scattering, |F(hkl)|^2 equals
|F(-h-k-l)|^2 (Friedel's law), so a structure and its inversion give the same intensities, and so does every
kinematic simulation of this kit and of py4DSTEM. No consumer may claim a polarity (which way a polar axis points, which
of two enantiomorphs a crystal is) from a kinematic match, however good the match is; that needs dynamical diffraction or
a measurement that is not centrosymmetric in intensity.
