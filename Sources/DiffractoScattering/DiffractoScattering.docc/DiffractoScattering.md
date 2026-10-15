# ``DiffractoScattering``

Electron scattering factors, kinematic structure factors, and the rings that survive an intensity cut.

## Overview

DiffractoScattering turns the geometry of DiffractoCore into intensities. ``ElectronScatteringFactors`` is the
Lobato and Van Dyck parameterisation (Z = 1 to 103, the table py4DSTEM ships), ``StructureFactors`` sums it over the atom
sites to give F(hkl) in inverse square angstroms, and ``KinematicReflectionList`` is the result an application reads: the
reflections up to a reciprocal radius, those with |F| above a cut, and their rings with the summed |F|^2.

The model is kinematic: no dynamical diffraction, no thermal factor beyond what the sites carry. Start with
<doc:StructureFactorsAndTheIntensityCut>.

## Topics

### Essentials

- <doc:StructureFactorsAndTheIntensityCut>
- ``KinematicReflectionList``
- ``StructureFactors``

### Scattering factors

- ``ElectronScatteringFactors``
- ``FormFactor``

### Structure factors

- ``StructureFactors``
- ``StructureFactor``

### The reflection list after the cut

- ``KinematicReflectionList``
- ``KinematicReflection``
- ``KinematicFamily``
- ``KinematicError``
