# ``DiffractoCore``

Crystals as numbers: cells, symmetry, CIF import, and the list of reflections with their rings.

## Overview

DiffractoCore holds the geometry of a crystal and nothing that depends on an electron. A ``UnitCell`` gives the real and
reciprocal lattice (reciprocal vectors without 2 pi, so d = 1/|g|), a ``CrystalStructure`` adds the atom sites and the
symmetry found from them, ``CIFImport`` reads a CIF into a structure with errors that name the tag or value at fault, and
``ReflectionList`` lists every hkl up to a reciprocal radius k_max with its systematic-absence flag and groups them into
rings (``ReflectionFamily``).

It uses `Foundation` and `simd` only. Start with <doc:FromACIFToRings>.

## Topics

### Essentials

- <doc:FromACIFToRings>
- ``CIFImport``
- ``CrystalStructure``
- ``UnitCell``

### Reading a CIF

- ``CIFImport``
- ``CIFDocument``
- ``CIFImportError``
- ``CIFImportWarning``

### Cells, sites and elements

- ``UnitCell``
- ``AtomSite``
- ``Element``
- ``MillerIndex``
- ``IntegerMatrix``

### Symmetry

- ``SymmetryOperation``
- ``PointSymmetry``
- ``SymmetryDetection``
- ``LaueClass``
- ``CrystalSystem``

### Reflections and rings

- ``ReflectionList``
- ``Reflection``
- ``ReflectionFamily``
- ``ReflectionListError``

### Reduced and conventional cells

- ``ReducedCell``
- ``ConventionalCell``
- ``ConventionalSetting``
- ``BravaisLattice``
