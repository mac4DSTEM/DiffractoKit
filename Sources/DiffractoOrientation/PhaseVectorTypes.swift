//
//  PhaseVectorTypes.swift
//  Role: the value types of vector-matched phase mapping: settings, resolution, verdict, result, map and input form.
//  Split out of PhaseVectorMatching.swift, which carries the method, its sources and its limits.
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
