//
//  PhaseVectorMatching.swift
//  Role: the MATCHER half of vector-matched phase mapping: for each scan position, remove the matrix's reflections, score
//        what survives against every candidate entry of a `PhaseReferenceLibrary`, and return a phase label, a distance, and
//        an explicit "not indexed" (registration docs/registration-slice-4.md, section 1).
//
//  Method: Thronsen et al., Ultramicroscopy 255 (2024) 113861, CC BY 4.0 (attribution in NOTICE). The decision rule of
//  `.knownVariants` and its two constants (direct-matrix count 1, residual cutoff 0.07 Å⁻¹) are as the paper describes them;
//  no code was copied from the authors' notebook (their repository carries no licence). Three deviations from what the paper describes are deliberate and are named where
//  they occur: the matrix mask is applied in VECTOR space rather than to the image, the matrix orientation is fitted ONCE
//  for the scan rather than per pattern, and the verdict carries a completeness requirement the paper's mean-distance
//  score does not have.
//
//  WHY THIS EXISTS AT ALL. On 2026-09-11 a template-matched route was refuted on a measurement: a bare best-score argmax
//  over phases returned GOLD for a pattern containing only ALUMINIUM, because one radial bin is 0.05 Å⁻¹ and Al-Au (111)
//  differ by 6 % of one bin. Vector matching works on the peak positions themselves, at the precision the sub-pixel
//  refinement already delivers, so the discriminating quantity is a distance in Å⁻¹ that a reader can compare against the
//  spacing they care about, not a correlation score that hides its own resolution.
//
//  WHAT IT CANNOT DO, stated here rather than discovered later: it is zone-axis only; it cannot separate phases
//  overlapping along the beam; and strain shifts Bragg positions, which Thronsen et al. name as the main cause of their
//  interface errors. One matrix orientation is fitted for the whole scan, so a second matrix grain would be mis-removed:
//  that is a limit of this implementation, not of the method.
//
//  NOT SCORED AGAINST THE PUBLISHED GROUND TRUTH. Scoring this against the published ground truth of Thronsen et al. needs
//  a ~7.4 GB download and has not been run; the verdict rules have no reference implementation (py4DSTEM's phase tool is
//  a different method), so every number this file produces is an output of the definitions in the registration, not a
//  measurement against an independent oracle.
//
//  THE INPUT. A scan is `[[SIMD2<Double>]]`, one array of peak positions per scan position, row-major (`scanWidth` ×
//  `scanHeight`), in one of two forms (`PeakCoordinates`): detector pixels with an origin and a scale (converted with a
//  `Float` subtraction first, so peaks held as `Float` positions give bit-identical vectors), or calibrated Å⁻¹ vectors
//  about the origin.
//

import Foundation
import simd
import DiffractoCore
import DiffractoScattering

// MARK: - Settings

/// The 18 settings of the matcher. Each is a tuned policy, not physics: the default is the value measured on the data
/// named in its comment, and a different dataset may want a different value.
public struct PhaseVectorSettings: Sendable, Equatable {
    /// The documented defaults, `PhaseVectorSettings()`: `directBeamRadiusInvAngstrom` 0.15, `maximumVectorInvAngstrom` 0
    /// (no outer reach), `matrixToleranceInvAngstrom` 0.02, `pairRadiusInvAngstrom` 0.02, `minimumVectors` 2,
    /// `matrixFallbackExplainedFraction` 0 (off), `minimumMatchedVectors` 3, `friedelPairMinimumMatchedVectors` 2,
    /// `chanceMatchMultiple` 5, `notIndexedAboveInvAngstrom` 0.015, `minimumPhaseContrastInvAngstrom` 0 (off),
    /// `orientationRelationshipToleranceDeg` 10, `completenessAwareCrossPhaseRanking` false, `classificationRule` `.search`,
    /// `directMatrixMaximumVectors` 1, `residualCutoffInvAngstrom` 0.07, `knownVariantsMinimumSpecificReflections` 1,
    /// `minimumMatrixEvidence` 0 (off).
    public static let `default` = PhaseVectorSettings()

    /// Peaks within this of the pattern origin are the direct beam. Å⁻¹.
    public var directBeamRadiusInvAngstrom: Double = 0.15
    /// Peaks beyond this are not reflections: the data's own reach — a
    /// detector edge, or a mask. 0 means the detector is the limit. Å⁻¹.
    /// MEASURED (Thronsen step 3, 2026-09-15): their patterns are masked
    /// beyond 0.70 Å⁻¹, the mask edge is a ring of maxima no phase explains,
    /// and without this every Al position came back "not indexed" (86 %
    /// mislabelled against 13 % with it at the same settings).
    public var maximumVectorInvAngstrom: Double = 0
    /// An experimental vector this close to a matrix reference vector is the
    /// matrix's, and is removed. Their image-space DoG masking, done in vector
    /// space. Å⁻¹.
    public var matrixToleranceInvAngstrom: Double = 0.02
    /// Largest distance at which an experimental and a reference vector count
    /// as the same reflection. Å⁻¹.
    ///
    /// The default is a budget, not a preference: sub-pixel peak refinement
    /// lands well inside 0.5 detector pixels, and on a typical calibration a
    /// pixel is ~0.008 Å⁻¹, so measurement contributes ~0.004. What dominates
    /// is Q-calibration error — 1 % at |q| = 0.7 Å⁻¹ is 0.007 Å⁻¹ — and
    /// strain, which shifts real Bragg positions and which Thronsen et al.
    /// name as their main source of interface error. 0.02 covers those with
    /// room. It is deliberately NOT generous: every increase widens the
    /// reference discs quadratically, and `chanceMatchFraction` is what that
    /// costs. Measured 2026-09-12: at 0.06 a large-cell precipitate entry
    /// matched a THIRD of any random vector set.
    public var pairRadiusInvAngstrom: Double = 0.02
    /// Fewer surviving vectors than this → the matrix, by exclusion.
    public var minimumVectors: Int = 2
    /// "The matrix, by exclusion" (`minimumVectors` in `.search`, `directMatrixMaximumVectors` in `.knownVariants`) needs
    /// at least this many vectors REMOVED as the matrix's; a position with fewer is "not indexed" instead. 0 disables it,
    /// the shipped default: with it off a position showing no vector at all of the matrix (one stray peak, a vacuum
    /// edge, an amorphous patch) is called matrix because too little survived, and the phase fraction counts it.
    /// Diagnosed 2026-10-21: the exclusion rule is the paper's (one or no survivor is the matrix) and the survivor
    /// count cannot tell "the matrix explained everything" from "nothing was there"; `removedCount` can. Not measured
    /// on a labelled scan, so it ships off; an empty pattern is already neither (`PhaseVectorResult()`).
    public var minimumMatrixEvidence: Int = 0

    /// When NO candidate phase clears its guards, fall back to the matrix
    /// verdict if the matrix already explained at least this fraction of the
    /// position's detected vectors. 0 disables it — the shipped default, and
    /// nothing moves until the objection recorded with it is ruled on.
    /// Measured 2026-09-16: improving detection turns matrix positions into
    /// refusals, because a position is otherwise called matrix only by exclusion.
    public var matrixFallbackExplainedFraction: Double = 0
    /// An entry must account for at least this many surviving vectors. Three,
    /// not two: two points fix a lattice only if you already know which two.
    public var minimumMatchedVectors: Int = 3
    /// … except when the survivors hold a Friedel pair — u and −u within the
    /// pair radius — which IS knowing which two: one lattice row, and the
    /// matrix removal has already said it is not the matrix's. Two, then.
    /// MEASURED on Thronsen et al.'s dataset A (2026-09-15): along [001]Al a T1 variant leaves exactly
    /// that pair and nothing else inside their mask; with a floor of three
    /// T1 recall was 6 %, with two 59 %, and the Al class paid 28 positions
    /// of 21 494 (0.13 %). Two limits, both from a review (2026-09-15): the
    /// chance guard still applies, and for a 48-vector entry it admits a pair
    /// only when the accessible radius is above ≈ 0.3 Å⁻¹ (below that the
    /// floor is silently a no-op — recall lost, not safety); and the matrix
    /// challenge cannot reach a two-of-two winner (it must explain strictly
    /// more than two from two), so a second matrix grain whose residual is
    /// exactly a candidate's pair is labelled that candidate. Not seen on
    /// the dataset; not tested.
    public var friedelPairMinimumMatchedVectors: Int = 2
    /// … and must beat CHANCE by this multiple.
    ///
    /// THIS IS THE DEVIATION from the published mean-distance score, and it is
    /// the second thing tried, not the first. A pure mean distance lets a
    /// phase explaining one vector out of ten beat one explaining nine,
    /// because the mean of a single small number is a small number — so the
    /// first attempt required a minimum matched FRACTION. Measured 2026-09-12,
    /// that fails on real geometry: the reference library is capped at
    /// `maximumVectorsPerEntry`, so a pattern showing more spots than the
    /// library holds can never reach any fraction, and the test scan's β″
    /// positions went to 0 % indexed while being perfectly matched.
    ///
    /// What actually needs excluding is a match that carries no information,
    /// and that has a number:
    /// `PhaseOrientationReference.chanceMatchFraction` is how often a vector
    /// pointing nowhere in particular lands within the pair radius of some
    /// reference. An entry is eligible only when its matched count exceeds
    /// that expectation by this multiple, so a dense library has to clear a
    /// proportionally higher bar — which is the honest way to let one phase
    /// carry forty reflections and another two.
    ///
    /// **WHEN IT ACTUALLY BINDS — corrected 2026-09-12, after a review refuted
    /// the first claim made for it.** The bar is `multiple · V · ρ² · n / R²`,
    /// which at the shipped settings crosses `minimumMatchedVectors` only
    /// above about **45 surviving vectors per pattern** — and the repo's own
    /// real cube has a median of **7**. So on typical SPED data the
    /// matched-vector floor is the binding constraint and this is inert.
    /// Measured: removing this guard ALONE
    /// leaves random-vector accuracy at 99.2 %, unchanged; removing the FLOOR
    /// alone takes it to 48.8 %; removing both takes it to 5.5 %. The note
    /// that stood here read the 99.2 → 5.5 collapse as this guard's doing,
    /// and it is the floor's.
    ///
    /// It is kept rather than deleted as dead weight, because it is what
    /// bounds a dense library on a pattern rich enough to need bounding —
    /// exactly the case a fixed `maximumVectorsPerEntry` cannot see.
    public var chanceMatchMultiple: Double = 5
    /// Best mean distance above this → "not indexed" rather than a forced
    /// label. Å⁻¹. Under `pairRadiusInvAngstrom`: the pair radius says what
    /// could be the same reflection, this says what is close enough to call
    /// it. THREE QUARTERS of the pair radius, not half (MEASURED 2026-09-15,
    /// Thronsen step 3): a correct 6–11-vector fit under strain and
    /// sub-pixel jitter has a mean residual of 0.5–0.8 pair radii, so at
    /// half it fell over the cliff while a lucky two-vector pair passed;
    /// at 0.75 and at 1.0 both datasets read the same (Thronsen 7.96 %
    /// free / 6.64 % with the OR at 0.75, 7.94 / 6.61 at 1.0, seven
    /// positions apart; the demonstration scan identical at all three), Al and the
    /// demonstration scan's false labels unchanged. Not 1.0:
    /// a mean of distances each within the pair radius is within it, so
    /// the refusal would never fire and the setting would be dead.
    public var notIndexedAboveInvAngstrom: Double = 0.015
    /// The winner must beat the best entry of every OTHER candidate phase by
    /// this margin, else "not indexed". 0 disables the requirement. This is
    /// the contrast margin the 2026-09-11 refutation asked for; it is off by
    /// default because a margin that is not measured on the data at hand is a
    /// guess, and it is to be measured on the data at hand.
    public var minimumPhaseContrastInvAngstrom: Double = 0
    /// A candidate entry whose phase lists an orientation relationship
    /// (`PhaseDefinition.orientationRelationships`) is scored only when its
    /// derived azimuth agrees with a listed pair's, relative to the fitted
    /// matrix entry, within this many degrees. MEASURED 2026-09-15 (Thronsen step 3):
    /// correct face-on matches land within ±15° of the orientation
    /// relationship and the confusing ones — face-on taken for edge-on —
    /// at 45°; but half of the correct edge-on matches sit near 22° and
    /// 67°, which a {0, 90} list at this tolerance cuts (edge-on recall
    /// 50 → 31 %). Open; no phase lists angles by default.
    public var orientationRelationshipToleranceDeg: Double = 10

    /// OFF BY DEFAULT (the best entry per phase is still chosen by
    /// mean distance alone). When true, the CROSS-phase winner in
    /// `classify` step 3 is ranked the same way `fitMatrixOrientation`
    /// already ranks candidate zone axes in this same file — by matched
    /// count first, mean distance only as the tiebreak — instead of by mean
    /// distance alone. Measured, not shipped: a two-vector Friedel pair at
    /// 0.004 Å⁻¹ beating a ten-vector match of another phase at 0.012 Å⁻¹
    /// is the observed case (a review, 2026-09-15) this exists to close, but
    /// it moves a per-position verdict and has not been cleared on real
    /// data. Measure it against a labelled scan before flipping it (the
    /// Thronsen scan is confounded by its unfixed T1 reference).
    public var completenessAwareCrossPhaseRanking = false

    // MARK: The known-variants rule (Thronsen et al.'s own, beside the search above)

    /// Which per-position rule decides a verdict. `.search` (shipped) is
    /// everything above: a library search ranked by mean distance and gated
    /// by the floors and guards this file measures. `.knownVariants` is
    /// the rule Thronsen et al. describe — score every non-matrix library entry
    /// against the survivors, no floors, argmin wins
    /// (`PhaseVectorMatcher.classifyKnownVariants`, written from the paper).
    public enum ClassificationRule: String, Sendable, Equatable, CaseIterable {
        case search
        case knownVariants
    }
    /// Shipped. Adding `.knownVariants` changes nothing about this rule.
    public var classificationRule: ClassificationRule = .search
    /// `.knownVariants` only: survivors at or below this many → the matrix,
    /// by exclusion. The paper's value (1).
    public var directMatrixMaximumVectors: Int = 1
    /// `.knownVariants` only: the winning entry's residual mean above this →
    /// "not indexed". Å⁻¹. The paper's value (0.07) on
    /// the paper's dataset — shipped as a setting, not retuned.
    public var residualCutoffInvAngstrom: Double = 0.07
    /// `.knownVariants` only: an `.indexed` winner must match at least this
    /// many PHASE-SPECIFIC reference reflections, or the position falls back
    /// to the matrix. 0 disables it.
    ///
    /// "Specific" = a reference vector of the winning entry that some
    /// survivor hit within `pairRadiusInvAngstrom`, and whose own `q` lies
    /// outside `matrixToleranceInvAngstrom` of every matrix reference, i.e. a
    /// reflection the matrix cannot also explain. The fallback is the matrix,
    /// not "not indexed": this rule has no matrix challenge, and too little
    /// evidence for a precipitate is what the matrix means here.
    ///
    /// DEVIATION from the paper's rule: it has no such guard (the argmin wins however little the
    /// winner matched). MEASURED on Thronsen et al.'s dataset A, stride 3,
    /// at this rule's defaults (pair radius 0.020 Å⁻¹): 529 → 423 of 29 241
    /// mislabelled at the 0.1 % detection floor, 424 → 383 at 0.15 %,
    /// 1189 → 1168 at 0.2 %. It fixes Al → precipitate false calls at a cost
    /// of 6 of 7 625 correct precipitate calls, and it relabels already-wrong
    /// precipitate calls as Al (precipitate → Al errors 39 → 85).
    /// Held out on spatial halves, k = 1 was chosen 170 of 200 times. It is
    /// NOT parameter-free: with k fixed, the result moves with the pair
    /// radius (2.18 % at 0.010 Å⁻¹, 1.35 % at 0.015). One dataset.
    /// Shipped on at 1 (2026-09-23).
    public var knownVariantsMinimumSpecificReflections: Int = 1

    public init() {}
}

/// What the tolerances mean on a PARTICULAR detector, and whether they can be
/// met on it.
///
/// WHY THIS EXISTS. The settings are in Å⁻¹ because that is the physics; the
/// measurement is on a detector grid. The two only meet through the Q
/// calibration, and nothing showed the user the conversion. Measured on a
/// run of 2026-09-12: on a 4x-binned cube at 0.045741 Å⁻¹ per
/// detector pixel, the shipped 0.020 Å⁻¹ tolerances are **0.44 of one pixel**.
/// Matrix removal then removed NOTHING — 802 752 detected peaks across
/// 108 900 patterns, matrix verdict count **zero** — and 108 899 of 108 900
/// positions came back "not indexed". Every number needed to say so was known
/// before the run.
///
/// It is arithmetic, not a threshold anyone tuned, and it is a value type so it
/// is testable and so every caller reads the same numbers.
public struct PhaseVectorResolution: Sendable, Equatable {
    public let invAngstromPerPixel: Double
    public let pairRadiusPixels: Double
    public let matrixRemovalPixels: Double
    public let notIndexedAbovePixels: Double

    /// A sentence for the user when a tolerance is smaller than the detector
    /// grid it will be measured on, else nil.
    ///
    /// Deliberately worded as DEMANDING and not as impossible: sub-pixel
    /// refinement genuinely beats one pixel on a well-sampled disk, as
    /// measured. What it cannot do is beat one pixel on a
    /// heavily binned detector whose disks are two or three pixels across —
    /// which is the case that produced the run above. The number is given so
    /// the user can judge; nothing refuses on it.
    public var advice: String? {
        guard invAngstromPerPixel > 0, pairRadiusPixels.isFinite else { return nil }
        guard pairRadiusPixels < 1 else { return nil }
        return String(
            format: "Pair radius is %.2f of one detector pixel (%.4f Å⁻¹ each). "
            + "Matching that tightly needs sub-pixel peak positions AND an "
            + "accurate Q scale; on a binned detector it will match nothing.",
            pairRadiusPixels, invAngstromPerPixel)
    }

    /// The same settings expressed on this detector: one pixel for the two
    /// tolerances, three quarters of a pixel for the verdict distance — the
    /// shipped ratio (0.02 / 0.02 / 0.015, measured 2026-09-15), carried
    /// onto the grid the data is on.
    public func scaledToDetector(_ settings: PhaseVectorSettings) -> PhaseVectorSettings {
        guard invAngstromPerPixel > 0 else { return settings }
        var out = settings
        out.pairRadiusInvAngstrom = invAngstromPerPixel
        out.matrixToleranceInvAngstrom = invAngstromPerPixel
        out.notIndexedAboveInvAngstrom = invAngstromPerPixel * 0.75
        out.directBeamRadiusInvAngstrom = max(settings.directBeamRadiusInvAngstrom,
                                              3 * invAngstromPerPixel)
        return out
    }

    public init(settings: PhaseVectorSettings, invAngstromPerPixel: Double) {
        self.invAngstromPerPixel = invAngstromPerPixel
        let scale = invAngstromPerPixel > 0 ? invAngstromPerPixel : .nan
        pairRadiusPixels = settings.pairRadiusInvAngstrom / scale
        matrixRemovalPixels = settings.matrixToleranceInvAngstrom / scale
        notIndexedAbovePixels = settings.notIndexedAboveInvAngstrom / scale
    }
}

// MARK: - Result

public enum PhaseVerdict: UInt8, Sendable, CaseIterable {
    /// Too little survived matrix removal to index — the matrix, by exclusion.
    case matrix = 0
    /// A candidate phase won, within every threshold.
    case indexed = 1
    /// Vectors survived, and nothing explained them well enough.
    case notIndexed = 2
    /// No usable peaks at this position at all.
    case noData = 3

    public var displayName: String {
        switch self {
        case .matrix: "Matrix"
        case .indexed: "Indexed"
        case .notIndexed: "Not indexed"
        case .noData: "No peaks"
        }
    }
}

/// One scan position's verdict and the evidence behind it. Every field a
/// caller needs to draw the evidence is here, so a label can always be taken apart.
public struct PhaseVectorResult: Sendable, Equatable {
    public var verdict: PhaseVerdict = .noData
    /// Winning phase, or the matrix's index for `.matrix`. −1 otherwise.
    public var phaseIndex: Int32 = -1
    /// Winning library entry (phase × zone axis × in-plane rotation). −1 when
    /// the verdict did not come from an entry.
    public var entryIndex: Int32 = -1
    /// Mean |u − v| over the unique reference vectors matched, Å⁻¹. NaN when
    /// no entry was scored.
    public var score: Float = .nan
    /// Surviving experimental vectors that found a reference within
    /// `pairRadius` in the winning entry.
    public var matchedCount: Int32 = 0
    /// Experimental vectors left after the direct beam and the matrix.
    public var survivingCount: Int32 = 0
    /// Experimental vectors removed as the matrix's.
    public var removedCount: Int32 = 0
    /// Best entry of the best OTHER candidate phase — the contrast the
    /// 2026-09-11 refutation showed a bare argmax hides. −1 / NaN if there is
    /// only one candidate phase.
    public var runnerUpPhaseIndex: Int32 = -1
    public var runnerUpScore: Float = .nan

    public init() {}

    /// Margin over the best other phase, Å⁻¹. NaN when there is no other.
    public var phaseContrast: Float {
        guard runnerUpScore.isFinite, score.isFinite else { return .nan }
        return runnerUpScore - score
    }
}

/// A verdict per scan position, row-major.
public struct PhaseMap: Sendable {
    public let width: Int, height: Int
    public var results: [PhaseVectorResult]
    /// The one matrix entry fitted for the whole scan. −1 if none was fitted.
    public let matrixEntryIndex: Int
    /// Phase display names, indexed by `phaseIndex`.
    public let phaseNames: [String]
    public let matrixPhaseIndex: Int

    public init(width: Int, height: Int, matrixEntryIndex: Int,
                             phaseNames: [String], matrixPhaseIndex: Int) {
        self.width = width; self.height = height
        // A size whose product overflows holds no results (`map` refuses such a scan; this initialiser must not trap).
        let (positions, overflow) = width.multipliedReportingOverflow(by: height)
        self.results = Array(repeating: PhaseVectorResult(), count: overflow ? 0 : max(0, positions))
        self.matrixEntryIndex = matrixEntryIndex
        self.phaseNames = phaseNames
        self.matrixPhaseIndex = matrixPhaseIndex
    }

    public func count(of verdict: PhaseVerdict) -> Int {
        results.reduce(0) { $0 + ($1.verdict == verdict ? 1 : 0) }
    }

    /// Positions labelled with each phase, by phase index.
    public var phaseCounts: [Int] {
        var out = [Int](repeating: 0, count: phaseNames.count)
        for r in results where r.phaseIndex >= 0 && Int(r.phaseIndex) < out.count {
            if r.verdict == .indexed || r.verdict == .matrix { out[Int(r.phaseIndex)] += 1 }
        }
        return out
    }
}

// MARK: - The two forms of input

/// How the positions of a scan's peaks become calibrated in-plane vectors (Å⁻¹).
public enum PeakCoordinates: Sendable, Equatable {
    /// Detector pixels: u = (Double(Float(x) − originX) · scale, Double(Float(y) − originY) · scale). The subtraction is
    /// done in `Float`, exactly as the apps' original did on its `Float` peak positions, so an application holding `Float`
    /// positions passes `SIMD2<Double>(Double(x), Double(y))` and gets bit-identical vectors.
    case detectorPixels(originX: Float, originY: Float, invAngstromPerPixel: Double)
    /// Calibrated vectors in Å⁻¹ about the origin, used as given (no `Float` rounding, no scale).
    case calibrated

    /// The calibrated vector of one position.
    @inline(__always)
    func vector(of p: SIMD2<Double>) -> SIMD2<Double> {
        switch self {
        case .detectorPixels(let originX, let originY, let invAngstromPerPixel):
            return SIMD2(Double(Float(p.x) - originX) * invAngstromPerPixel,
                         Double(Float(p.y) - originY) * invAngstromPerPixel)
        case .calibrated:
            return p
        }
    }
}

// MARK: - The matcher

public enum PhaseVectorMatcher {

    /// One pattern's peaks as calibrated in-plane vectors (Å⁻¹), with the
    /// direct beam removed.
    ///
    /// The frame is `dx = x − originX`, `dy = y − originY`, times
    /// `invAngstromPerPixel`, the frame of template matching. Ellipse
    /// correction and per-position origin collapse are the caller's job: they
    /// belong to calibration, not to matching. The R-Q rotation is deliberately
    /// NOT applied: the in-plane rotation is an axis of this search, fitted from
    /// the data, so applying a stored rotation as well would double-count it.
    ///
    /// Detector-pixel form: `peaks` are pixel positions (an application holding
    /// `Float` positions passes `SIMD2<Double>(Double(x), Double(y))` and gets
    /// bit-identical vectors); the subtraction is done in `Float`, then `Double`,
    /// then the scale, as in the apps' original.
    public static func experimentalVectors(peaks: [SIMD2<Double>],
                                            originX: Float, originY: Float,
                                            invAngstromPerPixel: Double,
                                            directBeamRadiusInvAngstrom: Double,
                                            maximumVectorInvAngstrom: Double = 0)
        -> [SIMD2<Double>] {
        experimentalVectors(
            peaks: peaks,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
            maximumVectorInvAngstrom: maximumVectorInvAngstrom)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func experimentalVectors(peaks: [SIMD2<Double>],
                                            coordinates: PeakCoordinates,
                                            directBeamRadiusInvAngstrom: Double,
                                            maximumVectorInvAngstrom: Double = 0)
        -> [SIMD2<Double>] {
        var out: [SIMD2<Double>] = []
        out.reserveCapacity(peaks.count)
        for p in peaks {
            if let q = experimentalVector(
                peak: p, coordinates: coordinates,
                directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: maximumVectorInvAngstrom) {
                out.append(q)
            }
        }
        return out
    }

    /// One peak as a calibrated vector, or nil when the matcher leaves it
    /// out (the direct beam, past the maximum vector, not finite). The single
    /// definition of "usable", shared by `experimentalVectors` and
    /// `claims(...)`, so a caller's idea of which disks were matched cannot
    /// drift from the map's.
    public static func experimentalVector(peak p: SIMD2<Double>,
                                           originX: Float, originY: Float,
                                           invAngstromPerPixel: Double,
                                           directBeamRadiusInvAngstrom: Double,
                                           maximumVectorInvAngstrom: Double = 0)
        -> SIMD2<Double>? {
        experimentalVector(
            peak: p,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
            maximumVectorInvAngstrom: maximumVectorInvAngstrom)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func experimentalVector(peak p: SIMD2<Double>,
                                           coordinates: PeakCoordinates,
                                           directBeamRadiusInvAngstrom: Double,
                                           maximumVectorInvAngstrom: Double = 0)
        -> SIMD2<Double>? {
        let q = coordinates.vector(of: p)
        let len = simd_length(q)
        guard len.isFinite, len > directBeamRadiusInvAngstrom,
              maximumVectorInvAngstrom <= 0 || len < maximumVectorInvAngstrom else { return nil }
        return q
    }

    // MARK: Nearest reference, with a length-band prune

    /// Index of and distance to the nearest vector of `refs` within `radius`.
    ///
    /// `refs` must be sorted by `length`, which `PhaseReferenceLibrary`
    /// guarantees: the reverse triangle inequality gives
    /// `| |u| − |v| | ≤ |u − v|`, so only the band `|u| ± radius` can hold a
    /// match and the rest is skipped without a distance computation.
    public static func nearest(_ u: SIMD2<Double>, in refs: [ReferenceVector],
                                radius: Double) -> (index: Int, distance: Double)? {
        guard !refs.isEmpty else { return nil }
        let uLen = simd_length(u)
        var lo = 0, hi = refs.count
        let lowerLength = uLen - radius
        while lo < hi {
            let mid = (lo + hi) / 2
            if refs[mid].length < lowerLength { lo = mid + 1 } else { hi = mid }
        }
        let upperLength = uLen + radius
        var best = -1
        var bestD = Double.infinity
        var i = lo
        while i < refs.count, refs[i].length <= upperLength {
            let d = simd_distance(u, refs[i].q)
            if d < bestD { bestD = d; best = i }
            i += 1
        }
        return (best >= 0 && bestD <= radius) ? (best, bestD) : nil
    }

    /// Per-worker scratch: one slot per reference vector, so scoring an entry
    /// allocates nothing.
    public final class Scratch: @unchecked Sendable {
        var bestPerRef: [Double]
        var touched: [Int]
        /// Room for `capacity` reference vectors (it grows on demand).
        public init(capacity: Int) {
            bestPerRef = [Double](repeating: .infinity, count: max(1, capacity))
            touched = []
            touched.reserveCapacity(max(1, capacity))
        }
    }

    /// Score `vectors` against one entry.
    ///
    /// For each experimental vector keep only its closest reference vector
    /// (their per-pattern reference subset, so a thin precipitate showing two
    /// spots is not penalised against a phase with forty), then average over
    /// the UNIQUE references hit — a reference claimed twice contributes its
    /// better distance once. Returns nil when nothing matched.
    public static func score(vectors: [SIMD2<Double>],
                              against entry: PhaseOrientationReference,
                              pairRadius: Double,
                              scratch: Scratch)
        -> (score: Double, matched: Int, uniqueReferences: Int)? {
        if scratch.bestPerRef.count < entry.vectors.count {
            scratch.bestPerRef = [Double](repeating: .infinity, count: entry.vectors.count)
        }
        scratch.touched.removeAll(keepingCapacity: true)
        var matched = 0
        for u in vectors {
            guard let hit = nearest(u, in: entry.vectors, radius: pairRadius) else { continue }
            matched += 1
            if scratch.bestPerRef[hit.index].isInfinite { scratch.touched.append(hit.index) }
            if hit.distance < scratch.bestPerRef[hit.index] {
                scratch.bestPerRef[hit.index] = hit.distance
            }
        }
        defer { for i in scratch.touched { scratch.bestPerRef[i] = .infinity } }
        guard !scratch.touched.isEmpty else { return nil }
        var sum = 0.0
        for i in scratch.touched { sum += scratch.bestPerRef[i] }
        return (sum / Double(scratch.touched.count), matched, scratch.touched.count)
    }

    /// Do `vectors` hold a Friedel pair — some u and v with |u + v| within
    /// `radius`? A ±g pair is one reciprocal-lattice row, the smallest thing
    /// that is evidence of a lattice rather than of two coincidences.
    public static func containsFriedelPair(_ vectors: [SIMD2<Double>], radius: Double) -> Bool {
        guard vectors.count >= 2 else { return false }
        for i in 0..<(vectors.count - 1) {
            for j in (i + 1)..<vectors.count
            where simd_length(vectors[i] + vectors[j]) <= radius && simd_length(vectors[i]) > radius {
                return true
            }
        }
        return false
    }

    // MARK: Which zone axis is the specimen on?

    /// One candidate beam direction, scored against the scan's own peaks.
    public struct ZoneAxisFit: Sendable, Equatable {
        public init(zoneAxis: SIMD3<Int>, inPlaneRotationRad: Double,
                     matchedVectors: Int, totalVectors: Int, meanDistance: Double,
                     chanceMatchedVectors: Double, sweepMedianFraction: Double) {
            self.zoneAxis = zoneAxis
            self.inPlaneRotationRad = inPlaneRotationRad
            self.matchedVectors = matchedVectors
            self.totalVectors = totalVectors
            self.meanDistance = meanDistance
            self.chanceMatchedVectors = chanceMatchedVectors
            self.sweepMedianFraction = sweepMedianFraction
        }

        public let zoneAxis: SIMD3<Int>
        public let inPlaneRotationRad: Double
        /// Experimental vectors this axis accounts for, across the sample.
        public let matchedVectors: Int
        /// Experimental vectors in the sample, total.
        public let totalVectors: Int
        /// Mean |u − v| over the pairs it made, Å⁻¹.
        public let meanDistance: Double
        /// How many of `totalVectors` an entry this dense would match by
        /// accident, summed over the sampled patterns at each pattern's own
        /// reach. The number that says whether `matchedVectors` is evidence.
        ///
        /// WHY IT IS HERE (2026-09-14). Without it the sweep presented a
        /// ⟨112⟩ family at 8 % exactly the way it presents a ⟨110⟩ family at
        /// 38 %, and a reader took a coin toss for a fit. It is computed from
        /// the same `chanceMatchFraction` the matcher's own guard uses, so
        /// there is one definition of chance in this file rather than two.
        ///
        /// **IT DOES NOT MARK THAT MOTIVATING CASE, and saying so is the
        /// point** (a review, 2026-09-14). At the default reference
        /// settings aluminium's ⟨112⟩ entries carry 12 to 16 vectors, not the
        /// 48 the cap allows, so the expectation is 0.37–1.3 % and 8 % clears
        /// five times it. The mark fires only while each pattern's
        /// second-largest |u| stays under about 0.55–0.63 Å⁻¹, and Al {220}
        /// alone is at 0.699. What this buys is the number itself, reported
        /// instead of absent. That case is caught by the SECOND null
        /// below, `sweepMedianFraction`, not by a threshold on this one.
        public let chanceMatchedVectors: Double
        /// The sweep's own null (a diagnostic run, 2026-09-15): the explained
        /// fraction of the MEDIAN axis in the sweep this fit came from.
        ///
        /// WHY A SECOND NULL. On a real crystal every wrong axis explains a
        /// large share of the vectors by SHARED reflections, not by chance:
        /// measured on planted aluminium at 2° steps with 30 % of spots
        /// missing, 0.004 Å⁻¹ of jitter and three spurious peaks per
        /// pattern, the median wrong axis explains 11–16 % (⟨110⟩, ⟨112⟩
        /// and ⟨001⟩ plants) and the worst wrong family up to 25 %, while
        /// every one of the 49 axes clears five times the disc-chance
        /// expectation. That is why a ⟨112⟩ at 8 % could never be
        /// marked by `chanceMatchedVectors`: it was not chance, it was the
        /// true axis's reflections seen through a wrong projection. Against
        /// the sweep's median those true families sat at 4.8–6.8× and every
        /// wrong family at 0.7–1.6×; the bar is placed at 2×. On vectors
        /// pointing nowhere the sweep median is ~0.5 % and the disc rule is
        /// what refuses; the two rules cover different failures and both
        /// must pass.
        ///
        /// WHAT THE BAR IS NOT (a review, 2026-09-15): a clean separator
        /// for every truth. ⟨111⟩, ⟨012⟩ and ⟨210⟩ plants put the true family
        /// at 2.7–2.8× and the worst wrong at 1.0–2.0×; a ⟨122⟩ plant puts
        /// the true family at 2.8× and the ⟨100⟩ family at 2.7–2.9× — a
        /// wrong family clearing the bar, sometimes above the truth. A
        /// high-index truth with few reflections shares too much with a
        /// low-index axis for a ratio to tell them apart; both rows then
        /// read as informative and the tie caption is the honest thing on
        /// screen. Heavier degradation (to 92 % missing, 15 spurious) never
        /// made the sweep rule bind before the disc rule; a two-grain scan
        /// left both true families at 5×. DEVIATION: py4DSTEM's zone-axis
        /// tools have no sweep null and no informativeness verdict.
        ///
        /// Fewer than five axes in the sweep is no null at all — with two,
        /// the "median" is the top axis's own fraction and the ratio
        /// saturates at 1 — so the median is reported as 0 (ratio infinite,
        /// the sweep rule inert) below that count.
        public let sweepMedianFraction: Double
        /// Explained fraction relative to the sweep's median axis.
        public var sweepRatio: Double {
            sweepMedianFraction > 0 ? explainedFraction / sweepMedianFraction : .infinity
        }
        /// Below this ratio an axis explains no more than a wrong axis does
        /// on the same data (measured: wrong ≤ 1.6×, true ≥ 4.8×).
        public static let informativeSweepRatio = 2.0
        public var isAboveSweep: Bool { sweepRatio >= Self.informativeSweepRatio }
        /// Both nulls: chance vectors (the disc model) and wrong axes on a
        /// real crystal (the sweep median).
        public func isInformative(multiple: Double) -> Bool {
            isAboveChance(multiple: multiple) && isAboveSweep
        }

        public var explainedFraction: Double {
            totalVectors > 0 ? Double(matchedVectors) / Double(totalVectors) : 0
        }
        /// The same fraction chance alone would explain.
        public var chanceFraction: Double {
            totalVectors > 0 ? chanceMatchedVectors / Double(totalVectors) : 0
        }
        /// Does this fit carry information, by the matcher's own standard?
        /// `multiple` is `PhaseVectorSettings.chanceMatchMultiple`.
        public func isAboveChance(multiple: Double) -> Bool {
            Double(matchedVectors) >= multiple * chanceMatchedVectors
        }
        public var inPlaneDegrees: Double { inPlaneRotationRad * 180 / .pi }

        /// Mean distances closer than this, in Å⁻¹, rank as equal.
        ///
        /// WHY (2026-09-30): symmetry-equivalent axes tie on matched count
        /// by construction, and their mean distances then differ only in the
        /// last bits of the arithmetic that rotated their references — so a
        /// tied family's order was set by rounding: deterministic for one build
        /// and peak set, but not a property of the data. Hardening, not the
        /// fix of a recorded "run-dependent tie order": that came from a
        /// caller ranking through a Dictionary (hash order per process).
        /// The bar is far under
        /// anything a detector resolves: float32 peaks carry ~1e-7 Å⁻¹.
        public static let meanDistanceTieResolution = 1e-9

        /// The ranking `fitZoneAxis` returns, as a TOTAL order: more matched
        /// vectors first; then the lower mean distance (compared on the
        /// `meanDistanceTieResolution` grid); then the axis indices,
        /// lexicographically ascending — the order `lowIndexZoneAxes` already
        /// lists one symmetry family in, so a tied family reads as it did
        /// when it tied exactly (an exact tie ACROSS families may not). The
        /// last key only ever decides between fits the first two cannot
        /// tell apart; it never outranks a better count or a closer fit.
        public static func ranksBefore(_ a: ZoneAxisFit, _ b: ZoneAxisFit) -> Bool {
            if a.matchedVectors != b.matchedVectors { return a.matchedVectors > b.matchedVectors }
            let da = rankedDistance(a.meanDistance), db = rankedDistance(b.meanDistance)
            if da != db { return da < db }
            if a.zoneAxis.x != b.zoneAxis.x { return a.zoneAxis.x < b.zoneAxis.x }
            if a.zoneAxis.y != b.zoneAxis.y { return a.zoneAxis.y < b.zoneAxis.y }
            return a.zoneAxis.z < b.zoneAxis.z
        }

        /// A grid index, so equality is transitive (a tolerance comparison is
        /// not, and `sort` needs a strict weak order). Non-finite ranks last.
        private static func rankedDistance(_ d: Double) -> Double {
            d.isFinite ? (d / meanDistanceTieResolution).rounded() : .infinity
        }
    }

    /// Rank beam directions by how much of the scan's measured peaks each one
    /// explains — the question "what orientation is this specimen on?", asked
    /// of the data instead of assumed.
    ///
    /// WHY IT EXISTS. A matrix viewed down an axis it is not on presents no
    /// reflections to remove, nothing is removed, and every position comes
    /// back "not indexed" — which on screen is indistinguishable from a method
    /// that does not work. That happened on 2026-09-12: Al set to
    /// [001], matrix verdict count **zero** across 108 900 positions. Measured
    /// afterwards, the specimen was on **⟨110⟩**, and
    /// all five sampled ⟨110⟩ equivalents tied at 39.0 % — an exact tie across
    /// the family being what cubic symmetry requires, and therefore a check
    /// that the sweep is behaving rather than merely returning something.
    ///
    /// Ranked by explained fraction, tie-broken on the lower mean distance,
    /// then on the axis indices (`ZoneAxisFit.ranksBefore`, a total order).
    /// NOT by mean distance alone, for the same reason `fitMatrixOrientation`
    /// is not: an axis explaining one vector at 0.001 Å⁻¹ would beat one
    /// explaining nine at 0.01.
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`); a count that is not
    /// width · height gives no fits, as no peaks do. Throws when the kit refuses
    /// the structure (`PhaseReferenceLibrary.Failure.structureRefused`,
    /// `.unsupportedElements`), and `.libraryTooLarge(requested:limit:)` when `inPlaneStepDeg` asks for more than
    /// `PhaseReferenceLibrary.maximumInPlaneStepCount` rotations (a step finer than 0.001 degrees), with the count asked
    /// for and the cap.
    public static func fitZoneAxis(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                    crystal: CrystalStructure,
                                    referenceSettings: PhaseReferenceSettings,
                                    settings: PhaseVectorSettings,
                                    originX: Float, originY: Float,
                                    invAngstromPerPixel: Double,
                                    candidateAxes: [SIMD3<Int>]
                                        = PhaseReferenceLibrary.lowIndexZoneAxes,
                                    inPlaneStepDeg: Double = 5,
                                    sampleLimit: Int = 400,
                                    isCancelled: (@Sendable () -> Bool)? = nil)
        throws(PhaseReferenceLibrary.Failure) -> [ZoneAxisFit] {
        try fitZoneAxis(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, crystal: crystal,
            referenceSettings: referenceSettings, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            candidateAxes: candidateAxes, inPlaneStepDeg: inPlaneStepDeg,
            sampleLimit: sampleLimit, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func fitZoneAxis(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                    crystal: CrystalStructure,
                                    referenceSettings: PhaseReferenceSettings,
                                    settings: PhaseVectorSettings,
                                    coordinates: PeakCoordinates,
                                    candidateAxes: [SIMD3<Int>]
                                        = PhaseReferenceLibrary.lowIndexZoneAxes,
                                    inPlaneStepDeg: Double = 5,
                                    sampleLimit: Int = 400,
                                    isCancelled: (@Sendable () -> Bool)? = nil)
        throws(PhaseReferenceLibrary.Failure) -> [ZoneAxisFit] {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight),
              !scan.isEmpty, !candidateAxes.isEmpty else { return [] }

        // A subsample, spread across the scan: this answers a question about
        // the GRAIN, and a grain does not change between neighbouring probe
        // positions. Sampling is what keeps a 62-axis sweep interactive.
        let step = max(1, scan.count / max(1, sampleLimit))
        var sample: [[SIMD2<Double>]] = []
        var index = 0
        while index < scan.count {
            let v = experimentalVectors(
                peaks: scan[index], coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            if !v.isEmpty { sample.append(v) }
            index += step
        }
        guard !sample.isEmpty else { return [] }
        let totalVectors = sample.reduce(0) { $0 + $1.count }

        // Each pattern's own reach, for the chance expectation: the
        // SECOND-largest |u|, the same rule `classify` uses and for the same
        // measured reason — with `max`, one spurious maximum far out inflates
        // the area and weakens the estimate for every other vector.
        let reach: [Double] = sample.map { vectors in
            var largest = 0.0, second = 0.0
            for u in vectors {
                let r = simd_length(u)
                if r > largest { second = largest; largest = r }
                else if r > second { second = r }
            }
            return second > 0 ? second : largest
        }

        var reference = referenceSettings
        reference.inPlaneStepDeg = inPlaneStepDeg
        // A step finer than 360 / maximumInPlaneStepCount degrees is refused, with the count asked for (saturated at
        // Int.max) and the cap, not held at the cap: a sweep at a different step than the one asked for would be a
        // different answer.
        let stepCount = PhaseReferenceLibrary.inPlaneStepCount(reference)
        guard stepCount <= Double(PhaseReferenceLibrary.maximumInPlaneStepCount) else {
            throw PhaseReferenceLibrary.Failure.libraryTooLarge(
                requested: stepCount < 9.0e18 ? Int(stepCount) : Int.max,
                limit: PhaseReferenceLibrary.maximumInPlaneStepCount)
        }
        let reflections = try PhaseReferenceLibrary.reflections(of: crystal, kMaxInvAngstrom: reference.kMaxInvAngstrom)
        let rotations = PhaseReferenceLibrary.inPlaneSteps(reference)
        let scratch = Scratch(capacity: max(1, reference.maximumVectorsPerEntry))

        // Every axis's best rotation first; the sweep median is known only
        // once all of them are, and each fit carries it.
        var winners: [(axis: SIMD3<Int>, theta: Double, matched: Int, meanDistance: Double, chance: Double)] = []
        for axis in candidateAxes {
            if isCancelled?() == true { break }
            let base = PhaseReferenceLibrary.projectedVectors(
                reflections: reflections, crystal: crystal,
                zoneAxis: axis, settings: reference)
            guard !base.isEmpty else { continue }
            // Constant across this axis's rotations: rotating an entry moves
            // its references, it does not change how many there are.
            let unrotated = PhaseOrientationReference(
                phaseIndex: 0, zoneAxis: axis, inPlaneRotationRad: 0, vectors: base)
            var chanceMatched = 0.0
            for (index, vectors) in sample.enumerated() {
                chanceMatched += unrotated.chanceMatchFraction(
                    pairRadius: settings.matrixToleranceInvAngstrom,
                    accessibleRadius: reach[index]) * Double(vectors.count)
            }
            var best: (theta: Double, matched: Int, meanDistance: Double)?
            for theta in rotations {
                let entry = PhaseOrientationReference(
                    phaseIndex: 0, zoneAxis: axis, inPlaneRotationRad: theta,
                    vectors: PhaseReferenceLibrary.rotate(base, by: theta))
                var matched = 0, pairs = 0
                var total = 0.0
                for vectors in sample {
                    guard let sc = score(vectors: vectors, against: entry,
                                         pairRadius: settings.matrixToleranceInvAngstrom,
                                         scratch: scratch) else { continue }
                    matched += sc.matched
                    total += sc.score * Double(sc.uniqueReferences)
                    pairs += sc.uniqueReferences
                }
                guard pairs > 0 else { continue }
                let meanDistance = total / Double(pairs)
                if best == nil || matched > best!.matched
                    || (matched == best!.matched && meanDistance < best!.meanDistance) {
                    best = (theta, matched, meanDistance)
                }
            }
            if let best {
                winners.append((axis, best.theta, best.matched, best.meanDistance, chanceMatched))
            }
        }
        let fractions = winners.map { Double($0.matched) / Double(totalVectors) }.sorted()
        let sweepMedian = fractions.count >= 5 ? fractions[fractions.count / 2] : 0
        var out = winners.map {
            ZoneAxisFit(zoneAxis: $0.axis, inPlaneRotationRad: $0.theta,
                        matchedVectors: $0.matched, totalVectors: totalVectors,
                        meanDistance: $0.meanDistance, chanceMatchedVectors: $0.chance,
                        sweepMedianFraction: sweepMedian)
        }
        out.sort(by: ZoneAxisFit.ranksBefore)
        return out
    }

    // MARK: The matrix orientation, fitted once

    /// The single matrix entry that best explains the whole scan.
    ///
    /// DEVIATION, deliberate: the paper masks the matrix per pattern. Here the
    /// matrix is one grain across the scan, so its in-plane rotation is one
    /// number — fitting it once is both faster and steadier than re-deciding
    /// it at every position from a handful of peaks. The cost is stated in the
    /// file header: a second matrix grain would be mis-removed.
    ///
    /// Chosen by total matched vectors, tie-broken by lower mean distance. Not
    /// by mean distance alone: an entry matching one vector at 0.001 Å⁻¹ would
    /// win over one matching nine at 0.01.
    ///
    /// The candidate side is guarded by `minimumMatchedVectors`
    /// and `chanceMatchMultiple`, which are floors on COUNT, not on the
    /// fraction of an entry's accessible vectors that matched. So the trap IS
    /// open on the candidate side: the cross-phase winner below is still
    /// chosen by mean distance alone, and a 2-vector Friedel pair at 0.004
    /// beats a 10-vector match at 0.012. That is measured, not hypothetical
    /// (measured, and recorded as an open item).
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`); a count that is not
    /// width · height is refused (nil), as is a scan with no peaks.
    public static func fitMatrixOrientation(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                             library: PhaseReferenceLibrary,
                                             settings: PhaseVectorSettings,
                                             originX: Float, originY: Float,
                                             invAngstromPerPixel: Double,
                                             sampleLimit: Int = 2000,
                                             isCancelled: (@Sendable () -> Bool)? = nil)
        -> (entryIndex: Int, matched: Int, meanDistance: Double)? {
        fitMatrixOrientation(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            sampleLimit: sampleLimit, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func fitMatrixOrientation(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                             library: PhaseReferenceLibrary,
                                             settings: PhaseVectorSettings,
                                             coordinates: PeakCoordinates,
                                             sampleLimit: Int = 2000,
                                             isCancelled: (@Sendable () -> Bool)? = nil)
        -> (entryIndex: Int, matched: Int, meanDistance: Double)? {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight) else { return nil }
        let matrixEntries = library.matrixEntryIndices
        guard !matrixEntries.isEmpty, !scan.isEmpty else { return nil }

        let stride = max(1, scan.count / max(1, sampleLimit))
        var sample: [[SIMD2<Double>]] = []
        sample.reserveCapacity(min(scan.count, sampleLimit) + 1)
        var index = 0
        while index < scan.count {
            let v = experimentalVectors(
                peaks: scan[index], coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom
            )
            if !v.isEmpty { sample.append(v) }
            index += stride
        }
        guard !sample.isEmpty else { return nil }

        let capacity = matrixEntries.map { library.entries[$0].vectors.count }.max() ?? 1
        let scratch = Scratch(capacity: capacity)
        var best: (entryIndex: Int, matched: Int, meanDistance: Double)?
        for entryIndex in matrixEntries {
            if isCancelled?() == true { return best }
            let entry = library.entries[entryIndex]
            var matched = 0
            var sum = 0.0
            var pairs = 0
            for vectors in sample {
                guard let s = score(vectors: vectors, against: entry,
                                    pairRadius: settings.matrixToleranceInvAngstrom,
                                    scratch: scratch) else { continue }
                matched += s.matched
                sum += s.score * Double(s.uniqueReferences)
                pairs += s.uniqueReferences
            }
            guard matched > 0, pairs > 0 else { continue }
            let mean = sum / Double(pairs)
            if let b = best {
                if matched > b.matched || (matched == b.matched && mean < b.meanDistance) {
                    best = (entryIndex, matched, mean)
                }
            } else {
                best = (entryIndex, matched, mean)
            }
        }
        return best
    }

    /// Is `scan` a scan: width and height positive, and exactly width · height positions? A malformed scan is refused by
    /// `map` and `fitMatrixOrientation` (nil) and gives no zone-axis fits, as the apps' original refused one. A width · height that
    /// overflows is not a scan.
    static func scanIsWellFormed(_ scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int) -> Bool {
        guard scanWidth > 0, scanHeight > 0 else { return false }
        let (positions, overflow) = scanWidth.multipliedReportingOverflow(by: scanHeight)
        return !overflow && scan.count == positions
    }

    /// True when `a` outranks `b` as the cross-phase winner. Shipped
    /// (`completenessAware: false`): mean distance alone, ascending —
    /// `a.score < b.score`. Candidate (`completenessAware: true`): matched
    /// count first, descending, mean distance only as the tiebreak — the
    /// same rule `fitMatrixOrientation` already uses above for candidate
    /// zone axes. Extracted from `classify`'s sort so it is testable without
    /// a crystal library: `PhaseVectorSettings.completenessAwareCrossPhaseRanking`'s
    /// own doc comment has the measurement and the open item.
    public static func crossPhaseWinsOver(
        _ a: (matched: Int, score: Double), _ b: (matched: Int, score: Double),
        completenessAware: Bool
    ) -> Bool {
        guard completenessAware else { return a.score < b.score }
        return a.matched != b.matched ? a.matched > b.matched : a.score < b.score
    }

    // MARK: One pattern

    /// Classify one pattern's already-calibrated vectors.
    public static func classify(vectors: [SIMD2<Double>],
                                 library: PhaseReferenceLibrary,
                                 settings: PhaseVectorSettings,
                                 matrixEntry: PhaseOrientationReference?,
                                 candidateEntryIndices: [Int],
                                 scratch: Scratch,
                                 matrixChallenge: [PhaseOrientationReference] = []) -> PhaseVectorResult {
        // Dispatch only, and only here: everything below this line is
        // `.search`, untouched by `.knownVariants` existing. The second rule
        // is a second function, not a branch threaded through this one, so
        // `.search`'s numbers cannot drift by way of a shared code path —
        // pinned by `testKnownVariantsRuleDoesNotChangeSearchsOwnResult`.
        if settings.classificationRule == .knownVariants {
            return classifyKnownVariants(vectors: vectors, library: library, settings: settings,
                                         matrixEntry: matrixEntry,
                                         candidateEntryIndices: candidateEntryIndices)
        }
        var result = PhaseVectorResult()
        guard !vectors.isEmpty else { return result }
        // Caller-supplied indices outside the library give no verdict (`map` checks the library's own consistency once).
        guard candidateEntryIndices.allSatisfy({ library.entries.indices.contains($0) }) else { return result }

        // 1 — matrix removal in vector space.
        var surviving: [SIMD2<Double>] = []
        surviving.reserveCapacity(vectors.count)
        if let matrixEntry {
            for u in vectors {
                if nearest(u, in: matrixEntry.vectors,
                           radius: settings.matrixToleranceInvAngstrom) == nil {
                    surviving.append(u)
                }
            }
        } else {
            surviving = vectors
        }
        result.survivingCount = Int32(surviving.count)
        result.removedCount = Int32(vectors.count - surviving.count)

        // 2 — the matrix, by exclusion.
        if surviving.count < settings.minimumVectors {
            guard Int(result.removedCount) >= settings.minimumMatrixEvidence else {
                result.verdict = .notIndexed
                return result
            }
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            return result
        }

        // 3 — score every candidate entry; keep the best per phase.
        //
        // The area chance is measured against is the one the DATA occupies,
        // not the library's kMax: if every surviving vector sits inside
        // 0.7 Å⁻¹, references beyond it cannot be hit by accident and must not
        // dilute the estimate.
        // The radius is the SECOND largest |u|, not the largest. Measured
        // 2026-09-12 (a review finding): with `max`, one spurious maximum far
        // out inflates the area the chance expectation is computed over,
        // weakens the guard for every other vector in the pattern, and took a
        // 1024-pattern false-positive rate from 0 to 10. Second-largest costs
        // nothing, is exact whenever the outermost peak is real, and cannot be
        // moved by a single outlier.
        var largest = 0.0, accessibleRadius = 0.0
        for u in surviving {
            let r = simd_length(u)
            if r > largest { accessibleRadius = largest; largest = r }
            else if r > accessibleRadius { accessibleRadius = r }
        }
        if accessibleRadius <= 0 { accessibleRadius = largest }
        var bestPerPhase: [Int: (entryIndex: Int, score: Double, matched: Int)] = [:]
        // The floor a candidate must clear: `minimumMatchedVectors`, or the
        // Friedel-pair floor when the survivors hold u and −u (see the
        // setting's note). Computed once per pattern, not per entry.
        let matchedFloor = containsFriedelPair(surviving, radius: settings.pairRadiusInvAngstrom)
            ? min(settings.minimumMatchedVectors, settings.friedelPairMinimumMatchedVectors)
            : settings.minimumMatchedVectors
        for entryIndex in candidateEntryIndices {
            let entry = library.entries[entryIndex]
            guard let s = score(vectors: surviving, against: entry,
                                pairRadius: settings.pairRadiusInvAngstrom,
                                scratch: scratch) else { continue }
            let chance = entry.chanceMatchFraction(
                pairRadius: settings.pairRadiusInvAngstrom,
                accessibleRadius: accessibleRadius) * Double(surviving.count)
            guard s.matched >= matchedFloor,
                  Double(s.matched) >= settings.chanceMatchMultiple * chance
            else { continue }
            let current = bestPerPhase[entry.phaseIndex]
            if current == nil || s.score < current!.score {
                bestPerPhase[entry.phaseIndex] = (entryIndex, s.score, s.matched)
            }
        }
        guard !bestPerPhase.isEmpty else {
            if settings.matrixFallbackExplainedFraction > 0, !vectors.isEmpty,
               Double(result.removedCount) / Double(vectors.count)
                   >= settings.matrixFallbackExplainedFraction {
                result.verdict = .matrix
                result.phaseIndex = Int32(library.matrixPhaseIndex)
                result.entryIndex = -1
                return result
            }
            result.verdict = .notIndexed
            return result
        }

        // An exact tie goes to the lower phase index (the order the library lists the phases), so the label never depends on
        // the Dictionary's per-process hash order; identical to the apps' original
        // wherever it is deterministic.
        let ranked = bestPerPhase.sorted {
            let a = (matched: $0.value.matched, score: $0.value.score), b = (matched: $1.value.matched, score: $1.value.score)
            let aWins = crossPhaseWinsOver(a, b, completenessAware: settings.completenessAwareCrossPhaseRanking)
            let bWins = crossPhaseWinsOver(b, a, completenessAware: settings.completenessAwareCrossPhaseRanking)
            return aWins != bWins ? aWins : $0.key < $1.key
        }
        let winner = ranked[0]
        result.score = Float(winner.value.score)
        result.matchedCount = Int32(winner.value.matched)
        if ranked.count > 1 {
            result.runnerUpPhaseIndex = Int32(ranked[1].key)
            result.runnerUpScore = Float(ranked[1].value.score)
        }

        // 4 — the two refusals. Both are recorded with the winner's numbers
        // still attached, so "not indexed" can be inspected rather than
        // guessed at.
        if winner.value.score > settings.notIndexedAboveInvAngstrom {
            result.verdict = .notIndexed
            return result
        }
        if settings.minimumPhaseContrastInvAngstrom > 0, ranked.count > 1 {
            let margin = ranked[1].value.score - winner.value.score
            if margin < settings.minimumPhaseContrastInvAngstrom {
                result.verdict = .notIndexed
                return result
            }
        }

        // 5 — the matrix gets the last word.
        //
        // WHY, and what it fixes (2026-09-14, a diagnostic run). Steps 1-4 ask only
        // which CANDIDATE explains the surviving vectors best. They never ask
        // whether those vectors are unexplained at all. At a position where
        // the matrix crystal is present on a zone axis other than the fitted
        // one, step 1 removes almost nothing and the whole pattern is offered
        // to the candidates — and a candidate coherent with the matrix wins by
        // default, because nothing else is in the competition. Measured on the
        // demonstration scan: 2 250 positions of pure aluminium on [011] came back
        // 100 % "β″ [001]". β″ [001] covers 10 of the 16 [011]Al reflections
        // at 0.0059 Å⁻¹; Al [011] covers all 16 at 0.0000 Å⁻¹, and was never
        // asked.
        //
        // So before a position may be called a candidate phase, the matrix
        // crystal is offered every low-index zone axis
        // (`matrixChallengeBases`, the fixed 49-direction list — NOT a list
        // the scan selects, and it ignores any zone axis the user restricted
        // the matrix phase to), scored by the same rule and held to the same
        // eligibility as any candidate.
        //
        // THE RULE: the challenger must explain STRICTLY MORE vectors, at a
        // smaller mean distance. See `challengeByMatrix` for why "strictly
        // more" rather than "at least as many" — it is what makes a fully
        // explained precipitate impossible to erase.
        //
        // THE IN-PLANE ROTATION IS DERIVED HERE, NOT FITTED FOR THE SCAN. An
        // earlier implementation took one rotation per axis from
        // `fitZoneAxis` — a whole-scan fit — and the fixture refuted it: at a
        // grain-B position the pool carried [-1 1 0] at 100°, matching NOTHING,
        // while the same axis at 130° matches 8 of 8 at 0.0001 Å⁻¹. A
        // whole-scan fit ranks a rotation by matches summed over every grain,
        // so a rotation picking up scattered matches everywhere outranks the
        // one that fits the grain exactly. A grain's rotation is a property of
        // the GRAIN, and must be asked at the position.
        //
        // Searching 180 rotations × 49 axes per position would cost more than
        // the whole map. Instead the rotations that can possibly win are
        // enumerated: if an orientation explains this pattern, it puts each
        // reference under an observed vector of the same length, so taking the
        // longest few observed vectors and every reference of matching length
        // gives every rotation worth testing — exactly, not on a 2° grid.
        if !matrixChallenge.isEmpty,
           let challenged = challengeByMatrix(
               surviving: surviving, bases: matrixChallenge, settings: settings,
               accessibleRadius: accessibleRadius,
               beating: (matched: winner.value.matched, score: winner.value.score),
               scratch: scratch) {
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            result.entryIndex = -1
            // The CHALLENGER's numbers, not the rejected candidate's: they are
            // the ones that describe the verdict given. A review caught the first
            // version discarding them while a comment claimed the verdict could
            // be inspected.
            result.score = Float(challenged.score)
            result.matchedCount = Int32(challenged.matched)
            return result
        }

        result.verdict = .indexed
        result.phaseIndex = Int32(winner.key)
        result.entryIndex = Int32(winner.value.entryIndex)
        return result
    }

    // MARK: The known-variants rule (Thronsen et al.)

    /// The per-position rule Thronsen et al. describe
    /// (written from the paper, no code copied): score every non-matrix library entry against the
    /// survivors and take the argmin, no floors and no guards. `classify`
    /// dispatches here for `settings.classificationRule == .knownVariants`;
    /// this function never runs for `.search`.
    ///
    /// DEVIATION from the paper's rule (which assigns Al directly when few vectors survive): the paper's
    /// five references are θ′ and T1 only — Al is never among them. Here Al
    /// (the matrix) is decided ENTIRELY by step b, by survivor count, and a
    /// matrix-phase entry is excluded from scoring in step c for that same
    /// reason: it is not that it would score badly, it is that the rule
    /// never asks.
    public static func classifyKnownVariants(vectors: [SIMD2<Double>],
                                              library: PhaseReferenceLibrary,
                                              settings: PhaseVectorSettings,
                                              matrixEntry: PhaseOrientationReference?,
                                              candidateEntryIndices: [Int]) -> PhaseVectorResult {
        var result = PhaseVectorResult()
        guard !vectors.isEmpty else { return result }
        guard candidateEntryIndices.allSatisfy({ library.entries.indices.contains($0) }) else { return result }

        // a — matrix removal, identical to `classify` step 1: same
        // tolerance, same survivors, same removedCount. DEVIATION: the
        // paper masks in IMAGE space at a fixed 8 detector pixels; this
        // removes in Å⁻¹ VECTOR space at `matrixToleranceInvAngstrom` — the
        // same vector-space departure the file header already states for
        // `.search`, inherited here because this step IS `.search`'s step
        // 1, not a reimplementation of it.
        var surviving: [SIMD2<Double>] = []
        surviving.reserveCapacity(vectors.count)
        if let matrixEntry {
            for u in vectors {
                if nearest(u, in: matrixEntry.vectors,
                           radius: settings.matrixToleranceInvAngstrom) == nil {
                    surviving.append(u)
                }
            }
        } else {
            surviving = vectors
        }
        result.survivingCount = Int32(surviving.count)
        result.removedCount = Int32(vectors.count - surviving.count)

        // b — the matrix, by exclusion (survivors at or below `directMatrixMaximumVectors`).
        if surviving.count <= settings.directMatrixMaximumVectors {
            guard Int(result.removedCount) >= settings.minimumMatrixEvidence else {
                result.verdict = .notIndexed
                return result
            }
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            return result
        }

        // c — score every entry whose phase is NOT the matrix's. The paper's
        // per-entry scoring: for every survivor, the TRUE nearest
        // reference vector of the entry (the argmin over the WHOLE set, no
        // radius — see `nearestReferenceVector` below for why `nearest(_:in:
        // radius:)` above is not reused), summed over every survivor, and
        // divided by the count of DISTINCT reference vectors hit
        // (the count of unique reference vectors) — a reference claimed
        // twice costs the sum without buying a second denominator credit,
        // which is what makes two survivors splitting one spot score worse
        // than two survivors hitting two spots at the same distances.
        //
        // DEVIATION (the paper's loop carries none of these — each omission
        // named against the paper's rule, which has no such term):
        //  - no `minimumMatchedVectors` / Friedel-pair floor (the paper: every
        //    survivor is scored, however few clear any particular entry);
        //  - no `chanceMatchMultiple` guard (the paper: no chance expectation
        //    is computed at all);
        //  - no matrix challenge, `classify` step 5 (the paper: the matrix is
        //    never re-offered once step b has passed it by);
        //  - no `minimumPhaseContrastInvAngstrom` margin (the paper:
        //    the phase with the lowest score is taken
        //    outright, no runner-up margin considered);
        //  - no completeness-aware cross-phase ranking (the paper: the argmin
        //    is by score alone, matched count never enters the ordering).
        //
        // DEVIATION: the paper's five entries are hand-built at fixed
        // orientations (T1 at two in-plane rotations under the stated
        // orientation relationship, θ′ edge-on at two, face-on at one);
        // every non-matrix entry THIS library holds is scored here instead —
        // for this app, whatever zone axes and in-plane steps the phase list
        // and `PhaseReferenceSettings` produced. T1 is a free in-plane sweep
        // until a later session derives its stated relationship
        // ((0001)T1 ∥ (111)Al, [1‑10]Al ∥ [10‑10]T1) into the constrained
        // form θ′ already has.
        var bestPerPhase: [Int: (entryIndex: Int, score: Double, matched: Int)] = [:]
        for entryIndex in candidateEntryIndices {
            let entry = library.entries[entryIndex]
            guard entry.phaseIndex != library.matrixPhaseIndex, !entry.vectors.isEmpty
            else { continue }
            var sum = 0.0
            var uniqueHits = Set<Int>()
            var matchedWithinRadius = 0
            for u in surviving {
                guard let hit = nearestReferenceVector(u, in: entry.vectors) else { continue }
                sum += hit.distance
                uniqueHits.insert(hit.index)
                if hit.distance <= settings.pairRadiusInvAngstrom { matchedWithinRadius += 1 }
            }
            guard !uniqueHits.isEmpty else { continue }
            let score = sum / Double(uniqueHits.count)
            if bestPerPhase[entry.phaseIndex] == nil || score < bestPerPhase[entry.phaseIndex]!.score {
                bestPerPhase[entry.phaseIndex] = (entryIndex, score, matchedWithinRadius)
            }
        }
        guard !bestPerPhase.isEmpty else {
            result.verdict = .notIndexed
            return result
        }

        // d — argmin over ALL entries, ties to the first in library order.
        // Reduced through `bestPerPhase` above (the global minimum is always
        // some phase's own minimum), then picked by entry index ascending so
        // a tie keeps the earliest entry, the way an argmin over a flat
        // array would.
        let ranked = bestPerPhase.map {
            (phaseIndex: $0.key, entryIndex: $0.value.entryIndex,
             score: $0.value.score, matched: $0.value.matched)
        }.sorted { $0.entryIndex < $1.entryIndex }
        var winner = ranked[0]
        for candidate in ranked.dropFirst() where candidate.score < winner.score {
            winner = candidate
        }
        var runnerUp: (phaseIndex: Int, score: Double)?
        for candidate in ranked where candidate.phaseIndex != winner.phaseIndex {
            if runnerUp == nil || candidate.score < runnerUp!.score {
                runnerUp = (candidate.phaseIndex, candidate.score)
            }
        }

        // Filled before the verdict is decided: a refused position must
        // still carry the would-be winner's numbers, so a caller can print
        // the distribution of scores that missed the cutoff.
        result.phaseIndex = Int32(winner.phaseIndex)
        result.entryIndex = Int32(winner.entryIndex)
        result.score = Float(winner.score)
        result.matchedCount = Int32(winner.matched)
        if let runnerUp {
            result.runnerUpPhaseIndex = Int32(runnerUp.phaseIndex)
            result.runnerUpScore = Float(runnerUp.score)
        }

        // The paper's residual cutoff (0.07 Å⁻¹): a winner above it is not indexed.
        result.verdict = winner.score > settings.residualCutoffInvAngstrom ? .notIndexed : .indexed

        // e — the evidence guard (DEVIATION, `knownVariantsMinimumSpecificReflections`):
        // an accepted winner with too few phase-specific reflections is the
        // matrix. Only `.indexed` is judged; a refusal stays a refusal. As in
        // the other matrix fallbacks, `entryIndex` is cleared; the rejected
        // winner's score, matched count and runner-up stay for diagnosis.
        let minimumSpecific = settings.knownVariantsMinimumSpecificReflections
        if result.verdict == .indexed, minimumSpecific > 0,
           specificReflectionCount(surviving: surviving, entry: library.entries[winner.entryIndex],
                                   matrixEntry: matrixEntry, settings: settings) < minimumSpecific {
            result.verdict = .matrix
            result.phaseIndex = Int32(library.matrixPhaseIndex)
            result.entryIndex = -1
        }
        return result
    }

    /// Distinct reference vectors of `entry` hit by some survivor within
    /// `pairRadiusInvAngstrom`, whose own `q` is not within
    /// `matrixToleranceInvAngstrom` of any matrix reference. With no matrix
    /// entry every hit counts. Transcribes the inline guard of the measuring
    /// tool the rule was developed with, so this rule reproduces its measured
    /// numbers.
    public static func specificReflectionCount(surviving: [SIMD2<Double>],
                                                entry: PhaseOrientationReference,
                                                matrixEntry: PhaseOrientationReference?,
                                                settings: PhaseVectorSettings) -> Int {
        var hits = Set<Int>()
        for u in surviving {
            if let hit = nearest(u, in: entry.vectors, radius: settings.pairRadiusInvAngstrom) {
                hits.insert(hit.index)
            }
        }
        guard let matrixEntry else { return hits.count }
        return hits.filter { index in
            let q = entry.vectors[index].q
            return !matrixEntry.vectors.contains {
                simd_distance($0.q, q) <= settings.matrixToleranceInvAngstrom
            }
        }.count
    }

    /// The TRUE nearest reference vector of `refs` — brute force, no radius
    /// and no length-band prune.
    ///
    /// `nearest(_:in:radius:)` above is deliberately NOT reused: its
    /// length-band prune assumes a candidate outside `|u| ± radius` can be
    /// discarded without being computed, which is exactly the restriction
    /// this rule does not have — a survivor here must find its true nearest
    /// reference however far away that reference sits, because the residual
    /// itself is what the cutoff judges.
    private static func nearestReferenceVector(_ u: SIMD2<Double>, in refs: [ReferenceVector])
        -> (index: Int, distance: Double)? {
        guard !refs.isEmpty else { return nil }
        var best = 0
        var bestDistance = simd_distance(u, refs[0].q)
        for i in 1..<refs.count {
            let d = simd_distance(u, refs[i].q)
            if d < bestDistance { bestDistance = d; best = i }
        }
        return (best, bestDistance)
    }

    /// The MATRIX crystal projected down every low-index zone axis, unrotated —
    /// the pool `classify` challenges a candidate label with. One entry per
    /// axis (49 at the shipped list), built once per scan; the in-plane
    /// rotation is derived at each position by `challengeByMatrix`, because a
    /// grain's rotation is a property of the grain and not of the scan.
    public static func matrixChallengeBases(
        library: PhaseReferenceLibrary
    ) -> [PhaseOrientationReference] {
        guard library.phases.indices.contains(library.matrixPhaseIndex) else { return [] }
        let crystal = library.phases[library.matrixPhaseIndex].crystal
        // `build` made this very list for this phase and structure, so it cannot be refused here; a library assembled by
        // hand around a structure the kit refuses has no challenge pool.
        guard let reflections = try? PhaseReferenceLibrary.reflections(
            of: crystal, kMaxInvAngstrom: library.settings.kMaxInvAngstrom),
              !reflections.isEmpty else { return [] }
        var out: [PhaseOrientationReference] = []
        out.reserveCapacity(PhaseReferenceLibrary.lowIndexZoneAxes.count)
        for axis in PhaseReferenceLibrary.lowIndexZoneAxes {
            let base = PhaseReferenceLibrary.projectedVectors(
                reflections: reflections, crystal: crystal,
                zoneAxis: axis, settings: library.settings)
            guard !base.isEmpty else { continue }
            out.append(PhaseOrientationReference(
                phaseIndex: library.matrixPhaseIndex, zoneAxis: axis,
                inPlaneRotationRad: 0, vectors: base))
        }
        return out
    }

    /// Does any orientation of the matrix crystal explain `surviving` at least
    /// as well as the winning candidate does? Returns the orientation that
    /// does, or nil.
    ///
    /// THE RULE, and why it is what it is. A challenger must explain STRICTLY
    /// MORE of the pattern than the winning candidate, at a smaller mean
    /// distance, having cleared the same floors a candidate clears.
    ///
    /// "Strictly more" is not fussiness; it is what makes the guard unable to
    /// erase a precipitate. A candidate that already explains every surviving
    /// vector cannot be beaten on count by anything, so a genuine precipitate
    /// pattern is safe BY CONSTRUCTION. The first version required only "at
    /// least as many", and a review measured what that costs: with a sparse
    /// precipitate of 3 to 6 vectors and 0.004 Å⁻¹ of jitter — a third of a
    /// detector pixel on the demonstration scan, inside this file's own stated
    /// measurement budget — a 49-axis search found a coincidental matrix fit
    /// with the same count and a smaller distance, and took 26.5 % of
    /// three-vector precipitates. A second matrix grain is the opposite case:
    /// the candidate explains part of it and the matrix explains all of it.
    ///
    /// Rotations are enumerated, not searched on a grid: an orientation that
    /// explains the pattern must put some reference of matching length under
    /// each observed vector, so every rotation worth testing is
    /// `angle(u) − angle(v)` for an observed `u` and a reference `v` with
    /// `||u| − |v||` inside the pair radius, taken over every surviving
    /// vector. The enumerated angle is exact only for a zero-distance match —
    /// its error is about `asin(ρ/|u|)` — and that error biases toward keeping
    /// the candidate label, which is the safe direction.
    public static func challengeByMatrix(
        surviving: [SIMD2<Double>],
        bases: [PhaseOrientationReference],
        settings: PhaseVectorSettings,
        accessibleRadius: Double,
        beating winner: (matched: Int, score: Double),
        scratch: Scratch
    ) -> (entry: PhaseOrientationReference, score: Double, matched: Int)? {
        guard !surviving.isEmpty, !bases.isEmpty else { return nil }
        let radius = settings.pairRadiusInvAngstrom
        // EVERY surviving vector seeds a rotation, not a chosen few. The first
        // version took the three LONGEST, on the reasoning that an angle is
        // most precise far out — and a review measured that backwards: spurious
        // maxima are typically FARTHER out than any real reflection, so the
        // longest-three rule hands the seeds to exactly them. With three
        // spurious peaks beyond the outermost real one, the catch rate for a
        // second matrix grain went from 100 % to 0 % (2026-09-14, 200 trials
        // per scenario). Seeding on all of them removes the failure mode and
        // the arbitrary constant with it; the cost is linear in the vector
        // count, measured at ~1.6 µs per surviving vector per challenged
        // position.
        //
        // The BEST qualifying orientation is returned, not the first. The
        // numbers it carries reach the result and the user, so "whichever
        // zone axis came first in the list" is not good enough for them.
        var best: (entry: PhaseOrientationReference, score: Double, matched: Int)?
        var rotated = surviving
        for base in bases {
            let chance = base.chanceMatchFraction(
                pairRadius: radius, accessibleRadius: accessibleRadius)
                * Double(surviving.count)
            for u in surviving {
                let uLength = simd_length(u)
                guard uLength > 0 else { continue }
                let uAngle = atan2(u.y, u.x)
                for v in base.vectors where abs(v.length - uLength) <= radius {
                    // Rotating the references by θ puts v under u; scoring the
                    // observations rotated by −θ against the unrotated base is
                    // the same comparison and rotates 8 vectors instead of 48.
                    let theta = uAngle - atan2(v.q.y, v.q.x)
                    let c = cos(theta), s = sin(theta)
                    for i in surviving.indices {
                        rotated[i] = SIMD2(surviving[i].x * c + surviving[i].y * s,
                                           -surviving[i].x * s + surviving[i].y * c)
                    }
                    guard let candidate = score(vectors: rotated, against: base,
                                                pairRadius: radius, scratch: scratch)
                    else { continue }
                    guard candidate.matched > winner.matched,
                          candidate.score < winner.score,
                          candidate.matched >= settings.minimumMatchedVectors,
                          Double(candidate.matched) >= settings.chanceMatchMultiple * chance
                    else { continue }
                    if let current = best,
                       current.matched > candidate.matched
                        || (current.matched == candidate.matched
                            && current.score <= candidate.score) { continue }
                    // Rotated, because the entry claims `theta`: a
                    // `PhaseOrientationReference`'s vectors are defined to be
                    // the ones AFTER its in-plane rotation, and handing back
                    // the unrotated base with a non-zero angle would be a trap
                    // for the next caller.
                    best = (PhaseOrientationReference(
                        phaseIndex: base.phaseIndex, zoneAxis: base.zoneAxis,
                        inPlaneRotationRad: theta,
                        vectors: PhaseReferenceLibrary.rotate(base.vectors, by: theta)),
                        candidate.score, candidate.matched)
                }
            }
        }
        return best
    }

    // MARK: The orientation relationship

    /// The in-plane azimuth, radians, of a Cartesian vector projected into
    /// an entry's detector frame (`ZoneAxisFrame.detectorBasis` of its
    /// zone axis, then the entry's own in-plane rotation) — nil when the
    /// vector is not in the zone (|v·n̂| > 1e-6·|v|, `n̂` the unit zone axis)
    /// or is zero.
    public static func projectedAzimuth(of v: SIMD3<Double>, zoneAxis n: SIMD3<Double>,
                                         inPlaneRotationRad: Double) -> Double? {
        let vLen = simd_length(v)
        guard vLen.isFinite, vLen > 1e-12 else { return nil }
        let nLen = simd_length(n)
        guard nLen.isFinite, nLen > 1e-12 else { return nil }
        let nHat = n / nLen
        guard abs(simd_dot(v, nHat)) <= 1e-6 * vLen else { return nil }

        let basis = ZoneAxisFrame.detectorBasis(zoneAxis: nHat)
        let x0 = simd_dot(v, basis.columns.0)
        let y0 = simd_dot(v, basis.columns.1)
        let c = cos(inPlaneRotationRad), s = sin(inPlaneRotationRad)
        // The library's own rotate convention (`PhaseReferenceLibrary.rotate`):
        // x' = c·x − s·y, y' = s·x + c·y, i.e. azimuth + θ.
        let x = c * x0 - s * y0
        let y = s * x0 + c * y0
        return atan2(y, x)
    }

    /// Is a candidate entry consistent with any listed relationship, given
    /// the fitted matrix entry?
    ///
    /// For each pair: the candidate's azimuth (its vector in the candidate
    /// crystal, projected through the candidate entry's own zone and
    /// rotation) minus the matrix's azimuth (likewise, through the matrix
    /// entry's own zone and rotation) must be within `toleranceDeg` of 0
    /// modulo 180° (see `PhaseDefinition.orientationRelationships` for why
    /// 180 and not 360). A pair whose vector is not in the respective zone
    /// does not apply to this (candidate, matrix) pair of entries; if no
    /// listed pair applies, the entry is free (true). An empty list is free.
    public static func orientationConsistent(candidate: PhaseOrientationReference,
                                               candidateCrystal: CrystalStructure,
                                               matrix: PhaseOrientationReference,
                                               matrixCrystal: CrystalStructure,
                                               relationships: [OrientationRelationship],
                                               toleranceDeg: Double) -> Bool {
        guard !relationships.isEmpty else { return true }
        guard let candidateZone = PhaseReferenceLibrary.cartesianZoneAxis(
                candidate.zoneAxis, crystal: candidateCrystal),
              let matrixZone = PhaseReferenceLibrary.cartesianZoneAxis(
                matrix.zoneAxis, crystal: matrixCrystal)
        else { return true } // a degenerate zone axis never reaches an entry; defensive.

        var anyApplicable = false
        for relationship in relationships {
            let candidateVector = relationship.candidate.cartesian(in: candidateCrystal)
            let matrixVector = relationship.matrix.cartesian(in: matrixCrystal)
            guard let candidateAz = projectedAzimuth(
                    of: candidateVector, zoneAxis: candidateZone,
                    inPlaneRotationRad: candidate.inPlaneRotationRad),
                  let matrixAz = projectedAzimuth(
                    of: matrixVector, zoneAxis: matrixZone,
                    inPlaneRotationRad: matrix.inPlaneRotationRad)
            else { continue }
            anyApplicable = true
            var d = (candidateAz - matrixAz) * 180 / .pi
            d = d.truncatingRemainder(dividingBy: 180)
            if d < 0 { d += 180 }
            let distance = min(d, 180 - d)
            if distance < toleranceDeg { return true }
        }
        return !anyApplicable
    }

    // MARK: The whole scan

    /// The phase map of a scan: one `PhaseVectorResult` per position, row-major, or nil when the run was refused or
    /// cancelled (never a partial map).
    ///
    /// `scan` is row-major (`scanWidth` × `scanHeight`), one array of peak positions per position; a malformed scan
    /// (width or height not positive, or a count that is not width · height) is refused (nil), as is a
    /// `matrixEntryIndex` outside the library, and a library whose candidate set is empty after the orientation
    /// relationships. `maximumWorkers` bounds the threads (nil: the active processor count); each position's result
    /// depends only on its own vectors, so the worker count moves no number. `progress` receives 0...1 from several
    /// threads; `isCancelled` is polled between positions.
    public static func map(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                            library: PhaseReferenceLibrary,
                            settings: PhaseVectorSettings,
                            originX: Float, originY: Float,
                            invAngstromPerPixel: Double,
                            matrixEntryIndex: Int? = nil,
                            maximumWorkers: Int? = nil,
                            progress: (@Sendable (Double) -> Void)? = nil,
                            isCancelled: (@Sendable () -> Bool)? = nil) -> PhaseMap? {
        map(scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: .detectorPixels(originX: originX, originY: originY,
                                         invAngstromPerPixel: invAngstromPerPixel),
            matrixEntryIndex: matrixEntryIndex, maximumWorkers: maximumWorkers,
            progress: progress, isCancelled: isCancelled)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func map(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                            library: PhaseReferenceLibrary,
                            settings: PhaseVectorSettings,
                            coordinates: PeakCoordinates,
                            matrixEntryIndex: Int? = nil,
                            maximumWorkers: Int? = nil,
                            progress: (@Sendable (Double) -> Void)? = nil,
                            isCancelled: (@Sendable () -> Bool)? = nil) -> PhaseMap? {
        guard scanIsWellFormed(scan, scanWidth: scanWidth, scanHeight: scanHeight) else { return nil }
        guard library.isConsistent else { return nil }
        if let matrixEntryIndex, !library.entries.indices.contains(matrixEntryIndex) { return nil }

        let fitted = matrixEntryIndex ?? fitMatrixOrientation(
            scan: scan, scanWidth: scanWidth, scanHeight: scanHeight, library: library, settings: settings,
            coordinates: coordinates, isCancelled: isCancelled
        )?.entryIndex
        guard isCancelled?() != true else { return nil }

        let matrixEntry = fitted.map { library.entries[$0] }
        // The orientation relationship (2026-09-15): once the matrix is
        // fitted, a candidate entry whose phase lists an orientation
        // relationship is scored only when the relationship's derived
        // azimuth agrees with the fitted matrix entry's — this is what keeps
        // a θ′ edge-on entry, rotated 45° so its (002) at 0.345 Å⁻¹ sits on
        // face-on's (110) at 0.350, from being offered as face-on. A phase
        // with an empty list stays free, so the default (no matrix, or no
        // list anywhere) moves nothing.
        let candidates: [Int]
        if let matrixEntry {
            let matrixCrystal = library.phases[library.matrixPhaseIndex].crystal
            candidates = library.candidateEntryIndices.filter { index in
                let entry = library.entries[index]
                let relationships = library.phases[entry.phaseIndex].orientationRelationships
                // An empty list means free: it would otherwise empty the
                // candidate set and refuse the whole map.
                guard !relationships.isEmpty else { return true }
                return Self.orientationConsistent(
                    candidate: entry, candidateCrystal: library.phases[entry.phaseIndex].crystal,
                    matrix: matrixEntry, matrixCrystal: matrixCrystal,
                    relationships: relationships, toleranceDeg: settings.orientationRelationshipToleranceDeg)
            }
        } else {
            candidates = library.candidateEntryIndices
        }
        guard !candidates.isEmpty else { return nil }

        // Built once for the scan, scored only at positions that would
        // otherwise be labelled — see `classify` step 5 for what it is for.
        let challenge = matrixChallengeBases(library: library)

        var map = PhaseMap(width: scanWidth, height: scanHeight,
                           matrixEntryIndex: fitted ?? -1,
                           phaseNames: library.phases.map(\.displayName),
                           matrixPhaseIndex: library.matrixPhaseIndex)

        let capacity = library.entries.map(\.vectors.count).max() ?? 1
        let total = scan.count
        let counter = ProgressCounter()

        map.results.withUnsafeMutableBufferPointer { out in
            let ptr = SendableBox(out.baseAddress!)
            let workers = max(1, min(maximumWorkers ?? ProcessInfo.processInfo.activeProcessorCount, total))
            let progressBatch = max(1, total / (workers * 100))
            DispatchQueue.concurrentPerform(iterations: workers) { worker in
                guard isCancelled?() != true else { return }
                let scratch = Scratch(capacity: capacity)
                var pending = 0
                for position in Swift.stride(from: worker, to: total, by: workers) {
                    if isCancelled?() == true { break }
                    let vectors = experimentalVectors(
                        peaks: scan[position], coordinates: coordinates,
                        directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                        maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom
                    )
                    ptr.value[position] = classify(
                        vectors: vectors, library: library, settings: settings,
                        matrixEntry: matrixEntry, candidateEntryIndices: candidates,
                        scratch: scratch, matrixChallenge: challenge
                    )
                    pending += 1
                    if progress != nil && pending >= progressBatch {
                        let completed = counter.add(pending)
                        pending = 0
                        progress?(Double(completed) / Double(total))
                    }
                }
                if progress != nil && pending > 0 {
                    let completed = counter.add(pending)
                    progress?(Double(completed) / Double(total))
                }
            }
        }
        guard isCancelled?() != true else { return nil }
        return map
    }

    /// Minimal wrapper to pass a mutable pointer into the concurrent closure.
    private struct SendableBox: @unchecked Sendable {
        let value: UnsafeMutablePointer<PhaseVectorResult>
        init(_ v: UnsafeMutablePointer<PhaseVectorResult>) { value = v }
    }

    /// The positions done so far, shared by the workers: a lock around a captured counter, as in the apps' original.
    private final class ProgressCounter: @unchecked Sendable {
        private let lock = NSLock()
        private var done = 0
        func add(_ n: Int) -> Int {
            lock.lock(); done += n; let completed = done; lock.unlock()
            return completed
        }
    }
}

// MARK: - Which detected disk did which phase claim?

/// What the phase map made of one detected Bragg disk at one scan position.
public enum PhaseDiskClaim: Sendable, Equatable {
    /// Never offered to the matcher: inside the direct-beam radius, past the
    /// maximum vector, or not a finite position.
    case notMatched
    /// Removed as the fitted matrix orientation's (step 1 of `classify`).
    case matrix
    /// Paired, within the pair radius, with a reference vector of the winning
    /// entry of this phase (its index in the library's `phases`).
    case phase(Int)
    /// Offered to the matcher, and claimed by nothing: it survived matrix
    /// removal and the winning entry has no reference vector within the pair
    /// radius — or the position has no winning entry (matrix by exclusion, not
    /// indexed), so nothing beyond the matrix could claim it.
    case unexplained
}

extension PhaseVectorMatcher {

    /// The claim of every peak of one scan position, in `peaks` order.
    ///
    /// A REPLAY OF THE MATCHER'S OWN PAIRING, not a second opinion: usable
    /// vectors from `experimentalVector`, matrix removal by the same
    /// `nearest(_:in:radius:)` at `matrixToleranceInvAngstrom` that
    /// `classify` step 1 (and `classifyKnownVariants` step a) uses, and the
    /// phase claim by the same `nearest` at `pairRadiusInvAngstrom` against
    /// the entry the map recorded for this position
    /// (`result.entryIndex`) — the pairing `score(vectors:against:…)` counts
    /// as `matchedCount`. So, for an `.indexed` position, the number of
    /// `.phase` claims equals `result.matchedCount`, the `.matrix` claims
    /// `removedCount`, and matrix + phase + unexplained = the usable vectors.
    ///
    /// Nothing is re-decided here: the verdict, the winning phase and entry
    /// are read from `result`, so a claim can name only the phase the map
    /// drew. A position the map did not index gives no `.phase` claim.
    ///
    /// `peaks` must be the CALIBRATED peaks the run matched, and `matrixEntry`
    /// the run's fitted matrix entry, nil when none was fitted.
    ///
    /// LIMIT: where the verdict is `.matrix` because the matrix challenge
    /// (`classify` step 5) found the matrix on another zone axis, that axis
    /// is not recorded (`entryIndex` is cleared), so the disks it explains
    /// are reported as unexplained rather than matrix.
    public static func claims(peaks: [SIMD2<Double>],
                               originX: Float, originY: Float,
                               invAngstromPerPixel: Double,
                               settings: PhaseVectorSettings,
                               matrixEntry: PhaseOrientationReference?,
                               result: PhaseVectorResult,
                               library: PhaseReferenceLibrary) -> [PhaseDiskClaim] {
        claims(peaks: peaks,
               coordinates: .detectorPixels(originX: originX, originY: originY,
                                            invAngstromPerPixel: invAngstromPerPixel),
               settings: settings, matrixEntry: matrixEntry, result: result, library: library)
    }

    /// The same, for either form of input (`PeakCoordinates`).
    public static func claims(peaks: [SIMD2<Double>],
                               coordinates: PeakCoordinates,
                               settings: PhaseVectorSettings,
                               matrixEntry: PhaseOrientationReference?,
                               result: PhaseVectorResult,
                               library: PhaseReferenceLibrary) -> [PhaseDiskClaim] {
        var winning: PhaseOrientationReference?
        if result.verdict == .indexed, result.entryIndex >= 0,
           Int(result.entryIndex) < library.entries.count {
            winning = library.entries[Int(result.entryIndex)]
        }
        return peaks.map { peak in
            guard let u = experimentalVector(
                peak: peak, coordinates: coordinates,
                directBeamRadiusInvAngstrom: settings.directBeamRadiusInvAngstrom,
                maximumVectorInvAngstrom: settings.maximumVectorInvAngstrom)
            else { return .notMatched }
            if let matrixEntry,
               nearest(u, in: matrixEntry.vectors,
                       radius: settings.matrixToleranceInvAngstrom) != nil {
                return .matrix
            }
            if let winning,
               nearest(u, in: winning.vectors, radius: settings.pairRadiusInvAngstrom) != nil {
                return .phase(Int(result.phaseIndex))
            }
            return .unexplained
        }
    }
}
