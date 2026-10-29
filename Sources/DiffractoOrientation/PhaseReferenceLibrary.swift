//
//  PhaseReferenceLibrary.swift
//  Role: the LIBRARY half of vector-matched phase mapping: for each phase, each sampled beam direction and each in-plane
//        rotation, the set of reciprocal-lattice vectors that would appear in the detector plane.
//
//  Method: Thronsen, Hjelen, Holmestad et al., "Studying GPB zones and precipitates in Al-Cu-Li with scanning precession
//  electron diffraction", Ultramicroscopy 255 (2024) 113861, CC BY 4.0 (attribution in NOTICE). The method is
//  reusable under that licence. This implementation is written from the published description; no code
//  was copied from the authors' notebook (their repository carries no licence), and every parameter default below is derived
//  here.
//
//  Why this is composition and not new mathematics: `CrystalStructure.kinematicReflections` already builds the reciprocal
//  lattice of an arbitrary cell with kinematic structure factors. This file adds three things it had not: a zone axis named
//  in LATTICE indices rather than Cartesian, an explicit in-plane rotation axis of the search, and a visibility cut.
//
//  The visibility cut is the one that carries scientific weight, so it is stated here rather than buried in a default. A
//  real pattern shows the few strongest reflections. A library holding every kinematically allowed vector gives EVERY
//  experimental vector a near neighbour for EVERY phase, and a nearest-neighbour score saturates: the same saturation
//  that refuted a template-matched route on 2026-09-11, when a bare best-score argmax over phases returned gold for a
//  pattern containing only aluminium, because one radial bin of 0.05 Å⁻¹ is wider than the Al-Au (111) difference.
//  `minimumIntensityFraction` is the guard, and a deliberately over-dense library must fail a negative control.
//
//  Units: |g| in Å⁻¹ with no 2π factor, as in `UnitCell` and py4DSTEM.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

// MARK: - One reference vector

/// A reflection projected into the detector plane, in Å⁻¹.
///
/// `hkl` is carried for reporting only — the matcher never reads it. It is what
/// lets a caller say *which* reflection a matched peak is, which is the
/// difference between a phase map a user can check and one they must trust.
public struct ReferenceVector: Sendable, Equatable {
    public let h: Int, k: Int, l: Int
    /// In-plane reciprocal vector (Å⁻¹), after the in-plane rotation.
    public let q: SIMD2<Double>
    /// |q| (Å⁻¹). Precomputed because the matcher prunes on it.
    public let length: Double
    /// Excited intensity |F|² × exp(−(s_g/σ)²), relative to the strongest
    /// reflection in the same entry. In (0, 1].
    public let relativeIntensity: Double

    public init(h: Int, k: Int, l: Int, q: SIMD2<Double>,
                             length: Double, relativeIntensity: Double) {
        self.h = h; self.k = k; self.l = l
        self.q = q
        self.length = length
        self.relativeIntensity = relativeIntensity
    }
}

// MARK: - One candidate orientation of one phase

/// The reference vectors of one phase at one beam direction and one in-plane
/// rotation. `vectors` is sorted by increasing `length`, which the matcher's
/// pruning depends on.
public struct PhaseOrientationReference: Sendable {
    public let phaseIndex: Int
    /// Beam direction as lattice indices [u v w], the way a crystallographer
    /// names a zone axis. Cartesian is derivable but not what anyone reads.
    public let zoneAxis: SIMD3<Int>
    public let inPlaneRotationRad: Double
    /// Sorted by increasing `length`.
    public let vectors: [ReferenceVector]

    /// The fraction of the accessible plane this entry's reference vectors
    /// cover at a given pair radius — an upper bound on how often a vector
    /// pointing NOWHERE IN PARTICULAR matches this entry by chance.
    ///
    /// `V·ρ²/R²` by the union bound over EVERY reference — an upper bound,
    /// not an estimate. See the body for why restricting it to the references
    /// the data can reach was tried, measured, and reverted.
    /// It is the number that says whether a match carries information: an
    /// entry covering a third of the plane explains a third of any random set,
    /// and its "match" means nothing. Reported in provenance rather than
    /// enforced, because the threshold at which it becomes unacceptable
    /// depends on how many vectors a pattern has — and a limit invented here
    /// would be a guess where a printed number is a measurement.
    public func chanceMatchFraction(pairRadius: Double, accessibleRadius: Double) -> Double {
        guard accessibleRadius > 0, pairRadius > 0 else { return 0 }
        // EVERY reference counts, including the ones outside the data's reach,
        // and that is deliberate — it is what makes this an UPPER BOUND rather
        // than an estimate.
        //
        // A review (2026-09-12) called the whole count against a restricted area
        // a 12x overestimate and proposed restricting the numerator to match.
        // That remedy was implemented and MEASURED, and it made the guard
        // weaker in exactly the case the guard exists for: with the numerator
        // restricted, the one-distant-spurious-peak test
        // went red on its FIRST assertion — random vectors in a narrow annulus
        // were indexed with no outlier at all — because a lower chance
        // estimate is a lower bar for eligibility. Reverted.
        //
        // The reason the whole count is the right conservative choice: the
        // uniform-disc model puts most of its area at large r, while real and
        // spurious peaks cluster at small r where a reference set is densest,
        // so the model UNDERSTATES the true chance for the distributions that
        // matter. Counting every reference pushes the other way. An
        // overestimate of chance demands more evidence, which is the safe
        // direction for a guard whose whole job is refusing uninformative
        // matches. It is not a probability and this comment does not call it
        // one.
        let p = Double(vectors.count) * pairRadius * pairRadius
            / (accessibleRadius * accessibleRadius)
        return min(1, p)
    }

    public init(phaseIndex: Int, zoneAxis: SIMD3<Int>,
                             inPlaneRotationRad: Double, vectors: [ReferenceVector]) {
        self.phaseIndex = phaseIndex
        self.zoneAxis = zoneAxis
        self.inPlaneRotationRad = inPlaneRotationRad
        self.vectors = vectors
    }
}

// MARK: - The orientation relationship, stated the way a crystallographer writes it

/// A lattice vector named the way a crystallographer names it: a plane
/// "(hkl)" — a reciprocal vector, `h·a* + k·b* + l·c*` — or a direction
/// "[uvw]" — a real one, `u·a + v·b + w·c`. Both come out as Cartesian
/// vectors in the crystal's own frame (`UnitCell.reciprocalLattice` rows for a
/// plane, `UnitCell.realLattice` rows for a direction), which is all the
/// projection (`PhaseVectorMatcher.projectedAzimuth`) needs — it never sees
/// `hkl` or `uvw` again.
public enum LatticeVector: Sendable, Equatable {
    case plane(SIMD3<Int>)
    case direction(SIMD3<Int>)

    /// This vector's Cartesian form in `crystal`'s own frame. Not normalised
    /// — `projectedAzimuth` only needs its direction, and the caller never
    /// reads its length.
    public func cartesian(in crystal: CrystalStructure) -> SIMD3<Double> {
        switch self {
        case .plane(let hkl):
            return Double(hkl.x) * crystal.cell.reciprocalLattice[0]
                + Double(hkl.y) * crystal.cell.reciprocalLattice[1]
                + Double(hkl.z) * crystal.cell.reciprocalLattice[2]
        case .direction(let uvw):
            return Double(uvw.x) * crystal.cell.realLattice[0]
                + Double(uvw.y) * crystal.cell.realLattice[1]
                + Double(uvw.z) * crystal.cell.realLattice[2]
        }
    }
}

/// One statement of an orientation relationship: this phase's vector is
/// parallel to the matrix's — the way a paper states it, e.g.
/// "(002)θ′ ∥ (200)Al". List every symmetry-equivalent variant explicitly;
/// the matcher knows no symmetry.
public struct OrientationRelationship: Sendable, Equatable {
    public let candidate: LatticeVector
    public let matrix: LatticeVector

    public init(candidate: LatticeVector, matrix: LatticeVector) {
        self.candidate = candidate
        self.matrix = matrix
    }
}

// MARK: - A phase as the library sees it

/// One phase in the library: a structure, a name, and whether it is the matrix.
///
/// The matrix is not a label the method infers — it is stated by the user,
/// because "which phase is the bulk" is knowledge about the specimen and not
/// about the data. Their method leaves the matrix out of the candidate set and
/// assigns it by exclusion; `role` is what makes that possible.
public struct PhaseDefinition: Sendable {
    public enum Role: String, Sendable {
        /// The bulk. Its reflections are removed from every pattern, and it is
        /// the label assigned when too little survives removal to index.
        case matrix
        /// A candidate phase. Scored.
        case candidate
    }

    public let id: String
    public let displayName: String
    public let crystal: CrystalStructure
    public let role: Role
    /// Beam directions to sample, as lattice indices. Empty means the caller
    /// wants `PhaseReferenceLibrary.lowIndexZoneAxes` for this phase.
    public let zoneAxes: [SIMD3<Int>]
    /// The orientation relationship this phase is stated to have with the
    /// matrix, as pairs of parallel lattice vectors — "(002) ∥ (200)", not a
    /// pre-computed angle. Empty is free, and the caller
    /// lists every symmetry-equivalent variant explicitly; the matcher knows
    /// no symmetry.
    ///
    /// DERIVED, PER ENTRY (2026-09-15, replacing a library-frame degree
    /// list): each pair's candidate and matrix vectors are projected into
    /// their OWN entry's detector frame (`ZoneAxisFrame.detectorBasis` of
    /// that entry's zone axis, then that entry's own in-plane rotation) and
    /// compared as azimuths — so the same statement is checked against
    /// every zone axis and rotation a phase samples, rather than one library
    /// frame the caller had to work out by hand. The comparison is modulo
    /// 180°: a ZOLZ under the flat Ewald sphere this file already uses (see
    /// `projectedVectors`) is centrosymmetric, g and −g both excited, so a
    /// plane and its negative are indistinguishable and "parallel" only
    /// means so up to sign.
    public let orientationRelationships: [OrientationRelationship]
    /// Per-phase override of `PhaseReferenceSettings.excitationSlabInvAngstrom`.
    /// `nil` (the default) means "use the library's global setting".
    ///
    /// A thin plate streaks its reciprocal-lattice points along the plate
    /// normal, so the slab a beam direction admits is a property of THAT
    /// PHASE'S shape, not a single number every phase shares. Thronsen et al.
    /// use 0.030 Å⁻¹ for T1 and 0.300 Å⁻¹ for θ′ (their Table 2); recording
    /// which value a given library was built with is provenance, and is the
    /// job of whoever calls `PhaseReferenceLibrary.build`, not of this type.
    public let excitationSlabInvAngstrom: Double?

    public init(id: String, displayName: String, crystal: CrystalStructure,
                             role: Role, zoneAxes: [SIMD3<Int>],
                             orientationRelationships: [OrientationRelationship] = [],
                             excitationSlabInvAngstrom: Double? = nil) {
        self.id = id
        self.displayName = displayName
        self.crystal = crystal
        self.role = role
        self.zoneAxes = zoneAxes
        self.orientationRelationships = orientationRelationships
        self.excitationSlabInvAngstrom = excitationSlabInvAngstrom
    }
}

// MARK: - Generation parameters

/// The settings of a reference library: seven parameters, each with the tuned value as its documented default.
public struct PhaseReferenceSettings: Sendable, Equatable {
    /// The documented defaults, `PhaseReferenceSettings()`: `kMaxInvAngstrom` 1.6, `excitationSlabInvAngstrom` 0.05,
    /// `excitationWidthInvAngstrom` 0.03, `minimumIntensityFraction` 0.05, `maximumVectorsPerEntry` 48, `inPlaneStepDeg` 2,
    /// `maximumEntries` 4000.
    public static let `default` = PhaseReferenceSettings()

    /// Largest |q| kept, Å⁻¹. Should be the detector's calibrated reach: a
    /// reference vector beyond it can never be matched, only counted.
    public var kMaxInvAngstrom: Double = 1.6
    /// Half-thickness of the slab about the zero-order Laue zone, Å⁻¹
    /// (`sgMax`). Larger admits reflections the
    /// beam excites only weakly; smaller risks dropping ones precession
    /// excites. Thronsen et al. use phase-dependent values (0.030 and 0.300
    /// Å⁻¹); this default is one number and is derived, not copied.
    public var excitationSlabInvAngstrom: Double = 0.05
    /// Gaussian width of the excitation weight, Å⁻¹ (`sgWidth`): the weight is exp(−(s/w)²) with w = max(width, 1e-9).
    public var excitationWidthInvAngstrom: Double = 0.03
    /// Reflections weaker than this fraction of the entry's strongest are
    /// dropped. One half of the guard against a saturating library — see the
    /// file header, and `maximumVectorsPerEntry`, which is the half that
    /// actually bounds the density.
    public var minimumIntensityFraction: Double = 0.05
    /// At most this many reference vectors per entry, strongest kept.
    ///
    /// An intensity fraction alone does NOT control density, which is what
    /// matters: measured 2026-09-12, a β″ [010] entry at
    /// `minimumIntensityFraction` 0.02 holds 130 vectors out to 1.2 Å⁻¹, and
    /// with a 0.06 Å⁻¹ pair radius those cover a THIRD of the plane — so a
    /// third of any random vector set matched by chance and the test scan's
    /// random-vector positions were labelled β″ 100 % of the time. A cap is
    /// the direct control, and it is also the physical statement: a real
    /// pattern from a thin precipitate shows a few dozen clear reflections,
    /// not every kinematically allowed one. See
    /// `PhaseOrientationReference.chanceMatchFraction`, which is what makes
    /// the trade visible rather than a matter of taste. A value of 0 or less
    /// disables the cap.
    public var maximumVectorsPerEntry: Int = 48
    /// In-plane rotation step, degrees. The search is over [0, 360).
    public var inPlaneStepDeg: Double = 2
    /// Refuse to build a library with more entries than this. A phase map is
    /// entries × scan positions × vectors of work; a library that takes an hour
    /// is a mistake in the settings, and saying so beats discovering it.
    public var maximumEntries: Int = 4000

    public init() {}
}

// MARK: - The library

public struct PhaseReferenceLibrary: Sendable {

    public enum Failure: Error, Equatable {
        /// A phase's structure has elements with no scattering factors.
        case unsupportedElements(phase: String, z: [Int])
        /// A zone axis of [0 0 0], which names no direction.
        case degenerateZoneAxis(phase: String)
        /// A phase produced no visible reflection at any sampled orientation.
        case noVisibleReflections(phase: String)
        /// `maximumEntries` exceeded. Carries what it would have been.
        case libraryTooLarge(requested: Int, limit: Int)
        /// No phase was marked `.matrix`, or more than one was.
        case matrixPhaseNotUnique(count: Int)
        /// Nothing to score.
        case noCandidatePhases
        /// The kit refused a phase's structure before any reflection was made (non-finite or out-of-range occupancy or
        /// coordinate, a degenerate cell, a close contact, a k-max the reflection list is not built to): the
        /// `KinematicError` text. A consequence of taking a `CrystalStructure`, whose validation a bare lattice and site
        /// list in an application did not have.
        case structureRefused(phase: String, reason: String)
    }

    public let phases: [PhaseDefinition]
    public let settings: PhaseReferenceSettings
    /// Every candidate entry, matrix included. `matrixEntryIndices` selects the
    /// matrix's; `candidateEntryIndices` the rest.
    public let entries: [PhaseOrientationReference]
    public let matrixPhaseIndex: Int

    /// Explicit because the implicit memberwise initializer of a public
    /// struct is `internal` (SE-0386), which `build()` can call from this
    /// same file but a consumer cannot. Every sibling struct in
    /// `PhaseVectorTypes.swift` this library feeds (`PhaseVectorSettings`,
    /// `PhaseVectorResult`, `PhaseOrientationReference`) states its own
    /// public init for the same reason: tests build small fixtures
    /// directly, bypassing `build()`'s crystallography so a `.knownVariants`
    /// test's distances are exact rather than derived from a projected cell.
    /// `build()` below is unaffected: same four fields, same order.
    public init(phases: [PhaseDefinition], settings: PhaseReferenceSettings,
                             entries: [PhaseOrientationReference], matrixPhaseIndex: Int) {
        self.phases = phases
        self.settings = settings
        self.entries = entries
        self.matrixPhaseIndex = matrixPhaseIndex
    }

    /// Does every entry name a phase of this library, and is `matrixPhaseIndex` one of them? A library from `build` always
    /// is; a hand-built one may not be, and the matcher refuses it rather than index out of range.
    public var isConsistent: Bool {
        phases.indices.contains(matrixPhaseIndex) && entries.allSatisfy { phases.indices.contains($0.phaseIndex) }
    }

    public var matrixEntryIndices: [Int] {
        entries.indices.filter { entries[$0].phaseIndex == matrixPhaseIndex }
    }
    public var candidateEntryIndices: [Int] {
        entries.indices.filter { entries[$0].phaseIndex != matrixPhaseIndex }
    }

    /// Low-index beam directions, one per ±pair, |u|,|v|,|w| ≤ 2 and coprime.
    /// Deliberately NOT reduced by crystal symmetry: a monoclinic precipitate
    /// read from a file is reported unreduced (2026-09-11), so a symmetry
    /// reduction would be a no-op where it matters most and an inconsistency
    /// everywhere else.
    public static let lowIndexZoneAxes: [SIMD3<Int>] = {
        var seen = Set<SIMD3<Int>>()
        var out: [SIMD3<Int>] = []
        for u in -2...2 {
            for v in -2...2 {
                for w in -2...2 {
                    if u == 0 && v == 0 && w == 0 { continue }
                    // max(1, …) at the END, not inside `gcd`: a `gcd`
                    // that floors itself at 1 returns gcd(0, 0) = 1, so
                    // [0 0 2] never reduced and the list held it beside
                    // [0 0 1] — 50 entries for 49 directions, and 180
                    // redundant library entries per phase at the default
                    // step (a review finding, 2026-09-12).
                    let g = max(1, gcd(gcd(abs(u), abs(v)), abs(w)))
                    let r = SIMD3(u / g, v / g, w / g)
                    // One of each ±pair: [uvw] and [-u-v-w] are the same axis
                    // viewed from opposite sides, and the projected vector set
                    // of one is the other's reflected through the origin —
                    // which a set containing both ±g cannot tell apart.
                    let neg = SIMD3(-r.x, -r.y, -r.z)
                    if seen.contains(r) || seen.contains(neg) { continue }
                    seen.insert(r)
                    out.append(r)
                }
            }
        }
        return out
    }()

    /// A true gcd: `gcd(0, 0)` is 0, not 1. Callers floor the RESULT.
    private static func gcd(_ a: Int, _ b: Int) -> Int {
        var a = a, b = b
        while b != 0 { (a, b) = (b, a % b) }
        return a
    }

    // MARK: Build

    /// Builds the library. A cancelled build (`isCancelled` returns true) throws `CancellationError`, never a library
    /// that lacks phases; a phase the kit's scattering table does not cover, or a structure the kit refuses, throws.
    public static func build(phases: [PhaseDefinition],
                              settings: PhaseReferenceSettings = .init(),
                              isCancelled: (@Sendable () -> Bool)? = nil)
        throws -> PhaseReferenceLibrary {

        let matrixIndices = phases.indices.filter { phases[$0].role == .matrix }
        guard matrixIndices.count == 1 else {
            throw Failure.matrixPhaseNotUnique(count: matrixIndices.count)
        }
        guard phases.contains(where: { $0.role == .candidate }) else {
            throw Failure.noCandidatePhases
        }
        for phase in phases {
            let unsupported = unsupportedElements(of: phase.crystal)
            guard unsupported.isEmpty else {
                throw Failure.unsupportedElements(phase: phase.displayName, z: unsupported)
            }
        }

        // Size the library before building it, so a refusal costs nothing: the count is a Double (a step of 1e-300 degrees
        // is 3.6e302 rotations) and no array of steps exists until the size is known to be within `maximumEntries`.
        let stepCount = inPlaneStepCount(settings)
        var requested = 0.0
        for phase in phases {
            let axes = phase.zoneAxes.isEmpty ? lowIndexZoneAxes : phase.zoneAxes
            guard !axes.contains(where: { $0 == SIMD3(0, 0, 0) }) else {
                throw Failure.degenerateZoneAxis(phase: phase.displayName)
            }
            requested += Double(axes.count) * stepCount
        }
        guard requested <= Double(settings.maximumEntries), let steps = Int(exactly: stepCount) else {
            // Saturated at Int.max when the count does not fit an Int.
            let asked = requested < 9.0e18 ? Int(requested) : Int.max
            throw Failure.libraryTooLarge(requested: asked, limit: settings.maximumEntries)
        }

        var entries: [PhaseOrientationReference] = []
        entries.reserveCapacity(Int(requested))

        for (phaseIndex, phase) in phases.enumerated() {
            if isCancelled?() == true { throw CancellationError() }
            let reflections = try Self.reflections(of: phase.crystal, kMaxInvAngstrom: settings.kMaxInvAngstrom,
                                                   phase: phase.displayName)
            let axes = phase.zoneAxes.isEmpty ? lowIndexZoneAxes : phase.zoneAxes
            var phaseHasAny = false
            // This phase's own slab, or the library's global one — applied to
            // THIS phase's entries only, never the others (see
            // `PhaseDefinition.excitationSlabInvAngstrom`).
            var phaseSettings = settings
            phaseSettings.excitationSlabInvAngstrom =
                phase.excitationSlabInvAngstrom ?? settings.excitationSlabInvAngstrom

            for axis in axes {
                if isCancelled?() == true { throw CancellationError() }
                let base = projectedVectors(reflections: reflections,
                                            crystal: phase.crystal,
                                            zoneAxis: axis,
                                            settings: phaseSettings)
                if base.isEmpty { continue }
                phaseHasAny = true
                for stepIndex in 0..<steps {
                    let theta = inPlaneAngle(stepIndex, of: steps)
                    entries.append(PhaseOrientationReference(
                        phaseIndex: phaseIndex, zoneAxis: axis, inPlaneRotationRad: theta,
                        vectors: rotate(base, by: theta)
                    ))
                }
            }
            guard phaseHasAny else {
                throw Failure.noVisibleReflections(phase: phase.displayName)
            }
        }

        return PhaseReferenceLibrary(phases: phases, settings: settings,
                                     entries: entries,
                                     matrixPhaseIndex: matrixIndices[0])
    }

    /// Atomic numbers of `structure` the kit's scattering table has no form factor for, ascending: the elements for
    /// which structure factors would be wrong (`build` refuses them).
    public static func unsupportedElements(of structure: CrystalStructure) -> [Int] {
        Array(Set(structure.sites.map(\.z).filter { !ElectronScatteringFactors.isSupported(z: $0) })).sorted()
    }

    /// The reflections a library is made from, in the order the matcher's sorts rely on: the kit's kinematic reflections
    /// of `structure` with 0 < |g| ≤ `kMaxInvAngstrom` and |F| above the default cut (`StructureFactors.defaultCut`),
    /// `g` Cartesian in Å⁻¹ without 2π and `structureFactor.intensity` = |F|².
    ///
    /// THE ORDER. The density cap's sort and the final length sort break ties on input order,
    /// so the input order is part of the definition. The original reflection list was enumerated h, k, l ascending
    /// lexicographically, filtered, then sorted by increasing |g| with a stable sort: equal |g| (to the last bit) kept
    /// the enumeration order, so ties are ordered by (h, k, l) ASCENDING. The kit's own list orders ties within 1e-12
    /// relative by DESCENDING (h, k, l), so it is not used as it is: the kept reflections are put in the order
    /// (|g| ascending, then h, k, l ascending), a total order that equals the stable sort of the enumeration.
    public static func reflections(of structure: CrystalStructure, kMaxInvAngstrom: Double, phase: String = "")
        throws(Failure) -> [KinematicReflection] {
        let list: KinematicReflectionList
        do {
            list = try structure.kinematicReflections(kMax: kMaxInvAngstrom)
        } catch {
            switch error {
            case .unsupportedElements(let z): throw Failure.unsupportedElements(phase: phase, z: z)
            default: throw Failure.structureRefused(phase: phase, reason: error.description)
            }
        }
        return list.kept.sorted {
            if $0.reflection.gLength != $1.reflection.gLength { return $0.reflection.gLength < $1.reflection.gLength }
            return $0.reflection.hkl < $1.reflection.hkl
        }
    }

    /// The most in-plane rotations `inPlaneSteps` returns and `fitZoneAxis` accepts: 360 000, a step of 0.001 degrees. A finer
    /// step is not a sampling choice but a mistake (a step of 1e-9 degrees asks for 3.6e11 rotations); `inPlaneSteps` holds
    /// it at this count, `fitZoneAxis` refuses it, and `build` refuses it through `maximumEntries`.
    public static let maximumInPlaneStepCount = 360_000

    /// How many in-plane rotations `settings.inPlaneStepDeg` asks for, as a `Double` so that a tiny step is a large number
    /// and not a conversion that traps: 1 for a step that is not finite, not positive or above 360, else 360 / step rounded.
    static func inPlaneStepCount(_ settings: PhaseReferenceSettings) -> Double {
        let step = settings.inPlaneStepDeg
        guard step.isFinite, step > 0, step <= 360 else { return 1 }
        return max(1, (360.0 / step).rounded())
    }

    /// The `index`-th of `count` evenly spaced rotations, radians (the one expression `inPlaneSteps` and `build` share).
    @inline(__always)
    static func inPlaneAngle(_ index: Int, of count: Int) -> Double {
        Double(index) * (2 * Double.pi / Double(count))
    }

    /// In-plane rotations sampled, radians, k 2π/N for N = 360 / step rounded. Always includes 0. A step that is not finite,
    /// not positive or above 360 gives `[0]`; a step finer than 360 / `maximumInPlaneStepCount` degrees gives
    /// `maximumInPlaneStepCount` rotations (the cap is on this array, which a caller holds in memory).
    public static func inPlaneSteps(_ settings: PhaseReferenceSettings) -> [Double] {
        let n = Int(min(inPlaneStepCount(settings), Double(maximumInPlaneStepCount)))
        return (0..<n).map { inPlaneAngle($0, of: n) }
    }

    // MARK: The projection itself

    /// The beam direction [u v w] as a Cartesian unit vector.
    ///
    /// A zone axis is a REAL-space direction, so it is `u·a₁ + v·a₂ + w·a₃`
    /// over `latReal`, not the Cartesian triple (u, v, w). For a cubic cell the
    /// two agree, which is exactly why getting it wrong stays invisible until
    /// the first monoclinic phase — the case this whole feature exists for.
    public static func cartesianZoneAxis(_ uvw: SIMD3<Int>, crystal: CrystalStructure)
        -> SIMD3<Double>? {
        let v = Double(uvw.x) * crystal.cell.realLattice[0]
            + Double(uvw.y) * crystal.cell.realLattice[1]
            + Double(uvw.z) * crystal.cell.realLattice[2]
        let n = simd_length(v)
        guard n.isFinite, n > 1e-12 else { return nil }
        return v / n
    }

    /// Reflections excited at `zoneAxis`, projected into the detector plane at
    /// in-plane rotation 0 and cut to the visible ones. Sorted by |q|.
    public static func projectedVectors(reflections: [KinematicReflection],
                                         crystal: CrystalStructure,
                                         zoneAxis uvw: SIMD3<Int>,
                                         settings: PhaseReferenceSettings)
        -> [ReferenceVector] {
        guard let n = cartesianZoneAxis(uvw, crystal: crystal) else { return [] }
        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: n)
        let e1 = basis.columns.0, e2 = basis.columns.1
        let sgMax = settings.excitationSlabInvAngstrom
        let sgWidth = max(settings.excitationWidthInvAngstrom, 1e-9)

        // MEASURED 2026-09-12 (a review finding), and it changes how the two
        // excitation settings should be read: with a flat sphere
        // s_g = g·n = (hu + kv + lw)/|r_uvw|, which is EXACTLY 0 for every
        // zone-law reflection and at least 1/|r_uvw| otherwise. For the axes
        // measured that day that spacing was 0.066 Å⁻¹ or more — above
        // `excitationSlabInvAngstrom` — so no non-ZOLZ reflection was admitted
        // and the Gaussian was identically 1 on everything that was. The slab
        // is then a ZONE SELECTOR, not a weighting, and
        // `excitationWidthInvAngstrom` cannot change any output until the
        // Ewald curvature is restored or the slab is widened past a Laue-zone
        // spacing. Both settings are kept because that is exactly what changes
        // for a long-axis cell, where 1/|r_uvw| falls below the slab.
        //
        // NOT true for every axis in use (measured 2026-09-30): T1's [0 -4 1] entry
        // (c = 14.145 Å cell) holds 48 vectors, 30 of them first-order-Laue-
        // zone reflections at |s_g| = 0.0411 Å⁻¹ — below the 0.05 slab, so
        // admitted, at Gaussian weight e^-1.9 at the default 0.03 width
        // (e.g. (0 0 ±1), |g| 0.0575 Å⁻¹, relative intensity 0.173). For that
        // entry the slab is a weighting too, and the width does move output.
        //
        // Flat Ewald sphere: s_g = g·n. Deliberately flat, and NOT a
        // curvature-corrected form. The
        // curvature term exists in template matching to break a 180° ambiguity in an
        // AZIMUTHAL CORRELATION by making g and −g excite differently. Here the
        // score is a distance between vectors, the in-plane rotation is an
        // explicit axis of the search rather than something recovered from an
        // FFT, and there is no ambiguity to break — so the curvature would only
        // add a wavelength the caller may not have and an asymmetry nothing
        // reads. Under precession, which is how SPED data is taken, the
        // effective excitation is symmetrised anyway.
        var raw: [(ReferenceVector, Double)] = []
        raw.reserveCapacity(reflections.count)
        var strongest = 0.0
        for refl in reflections {
            let sg = simd_dot(refl.reflection.g, n)
            if abs(sg) > sgMax { continue }
            let excited = refl.structureFactor.intensity * exp(-(sg / sgWidth) * (sg / sgWidth))
            guard excited.isFinite, excited > 0 else { continue }
            let x = simd_dot(refl.reflection.g, e1)
            let y = simd_dot(refl.reflection.g, e2)
            let q = SIMD2(x, y)
            let len = simd_length(q)
            guard len.isFinite, len > 1e-9, len <= settings.kMaxInvAngstrom else { continue }
            strongest = max(strongest, excited)
            raw.append((ReferenceVector(h: refl.reflection.hkl.h, k: refl.reflection.hkl.k, l: refl.reflection.hkl.l,
                                        q: q, length: len, relativeIntensity: excited), excited))
        }
        guard strongest > 0 else { return [] }

        let floor = settings.minimumIntensityFraction * strongest
        var kept = raw.filter { $0.1 >= floor }
        // Strongest first, then cap, then sort by length for the matcher's
        // prune. Ties broken on |q| so the cap is deterministic: a set that
        // depends on hash or enumeration order is a gate that drifts.
        if settings.maximumVectorsPerEntry > 0, kept.count > settings.maximumVectorsPerEntry {
            kept.sort {
                $0.1 != $1.1 ? $0.1 > $1.1 : $0.0.length < $1.0.length
            }
            kept.removeLast(kept.count - settings.maximumVectorsPerEntry)
        }
        var out = kept.map { pair in
            ReferenceVector(h: pair.0.h, k: pair.0.k, l: pair.0.l, q: pair.0.q,
                            length: pair.0.length, relativeIntensity: pair.1 / strongest)
        }
        out.sort { $0.length < $1.length }
        return out
    }

    /// Rotate a projected set in the detector plane. Keeps the sort, because
    /// a rotation does not change any |q|.
    public static func rotate(_ vectors: [ReferenceVector], by theta: Double)
        -> [ReferenceVector] {
        guard theta != 0 else { return vectors }
        let c = cos(theta), s = sin(theta)
        return vectors.map {
            ReferenceVector(h: $0.h, k: $0.k, l: $0.l,
                            q: SIMD2(c * $0.q.x - s * $0.q.y, s * $0.q.x + c * $0.q.y),
                            length: $0.length, relativeIntensity: $0.relativeIntensity)
        }
    }
}
