//
//  DiffractoOrientation.swift
//  Role: the module overview. DiffractoOrientation holds (registration docs/registration-slice-4.md):
//
//    - `PhaseReferenceLibrary` and its value types (`PhaseDefinition`, `PhaseReferenceSettings`, `PhaseOrientationReference`,
//      `ReferenceVector`, `LatticeVector`, `OrientationRelationship`): for each phase, zone axis and in-plane rotation, the
//      reciprocal-lattice vectors a detector would show, from the kit's kinematic reflections;
//    - `PhaseVectorMatcher` and its value types (`PhaseVectorSettings`, `PhaseVectorResolution`, `PhaseVerdict`,
//      `PhaseVectorResult`, `PhaseMap`, `PhaseDiskClaim`, `PeakCoordinates`): the vector-matched phase map of a scan, after
//      Thronsen et al., Ultramicroscopy 255 (2024) 113861;
//    - `EulerAngles` (py4DSTEM's orientation matrix to and from SciPy's Bunge angles), `ZoneAxisFrame` (the detector frame of
//      a zone axis) and `OrientationResult` (one position's orientation result).
//
//  Template matching remains planned (slice 5). CrystalOrientationSymmetry implements eleven-class minimum-angle reduction,
//  native direction sampling and explicit IPFColorKey conventions (docs/orientation-symmetry.md).
//  Depends on DiffractoCore and DiffractoScattering.
//
