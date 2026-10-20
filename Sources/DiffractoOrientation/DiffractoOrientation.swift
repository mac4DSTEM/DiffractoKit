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
//    - `ACOMTemplateLibrary` and `ACOMMatcher` (+ `ACOMTemplateSettings`, `ACOMExcitationModel`, `ACOMReflection`, `ACOMPeak`,
//      `ACOMError`): the apps' simplified polar-correlation ACOM route (docs/parity.md P25); not py4DSTEM's SO(3) correlogram.
//
//  py4DSTEM's SO(3) correlogram matcher and pyxem's template matching are not ported (slice 5); the ACOM route above is the apps' simplified one. CrystalOrientationSymmetry implements eleven-class minimum-angle reduction,
//  native direction sampling, disorientation, the py4DSTEM direction families and explicit IPFColorKey conventions (native, py4DSTEM, orix); PhaseModel (validation, orientation class, Materials Project models) (docs/orientation-symmetry.md).
//  Depends on DiffractoCore, DiffractoScattering and ImagoKit's ImagoFFT (FFT1D, for the ring FFTs).
//
