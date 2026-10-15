//
//  EllipseCalibration+Policy.swift
//  Role: the choices the ellipse fit makes that physics does not fix (DiffractoKit 0.2.0, docs/registration-0.2.0.md §1):
//        the azimuthal coverage rule, the residual caps, the one-ring ratio, the robust-refit cut, the ring-width floors
//        and the thresholds of the pre-fit checks. Each is a field with its documented default; `Policy.default` is the
//        set every earlier release used, and the fits run bit for bit as before with it.
//

import Foundation

extension EllipseCalibration {

    /// The policy of the ellipse fit: every number below was chosen by the suite (measured on its fixtures and on the
    /// owner's data, or derived and rounded), none is a physical constant. `Policy.default` is what the fits used before
    /// 0.2.0, and a fit given `.default` is bit identical to one given no policy.
    ///
    /// An invalid policy (see `isValid`) is refused by every fit that takes it with `FitError.invalidInput`, and gives
    /// empty ring evidence; it is never clamped silently.
    public struct Policy: Sendable, Equatable {
        /// Azimuthal bins the coverage checks divide the annulus into. Default 36 (10° each).
        public var angularBinCount: Int
        /// The fewest occupied bins an annulus needs to be fitted without the caller's assertion of one ring. Default
        /// `angularBinCount * 5 / 6` (30 of 36): a gap above about 60° refuses.
        public var degeneracyBoundBins: Int
        /// The fewest occupied bins even the `acceptSparseCoverage` path accepts. Default `angularBinCount / 3` (12 of 36).
        public var sparseFloorBins: Int
        /// A sample is "strong" (counts as occupying its bin) from `minimum + strongSignalFraction * range`. Default 0.2.
        public var strongSignalFraction: Double
        /// On the sparse path, the largest mean-radius ratio across bins before the annulus is called more than one ring.
        /// Default 1.10 (fcc {111}/{200} are 15.5 % apart).
        public var oneRingRatio: Double
        /// The fewest finite annulus pixels the conic fit needs. Default 20.
        public var minimumConicSamples: Int
        /// The fewest finite annulus pixels the profile fit needs. Default 32.
        public var minimumProfileSamples: Int
        /// The signal test: the annulus's range must exceed `|maximum| * minimumRelativeSignal`. Default 1e-6.
        public var minimumRelativeSignal: Double
        /// A pixel is an isolated impulse when it lies outside its neighbours' range by more than this many robust ranges
        /// of the annulus (P99.5 - P0.5). Default 2; 0 turns the impulse exclusion off.
        public var impulseRangeFactor: Double
        /// The conic fit's residual cap when the caller gives none. Default 0.2.
        public var maximumConicResidual: Double
        /// The profile fit's residual cap when the caller gives none. Default 0.22.
        public var maximumProfileResidual: Double
        /// The conic that seeds the profile fit is accepted up to this residual or the conic cap, whichever is larger: a
        /// broad ring may exceed the conic's publication residual and still seed a profile. Default 0.4.
        public var profileSeedResidual: Double
        /// The largest axis ratio a/b of an accepted ellipse; the conic's smallest b/a is its reciprocal. Default 5 (0.2).
        public var maximumAxisRatio: Double
        /// The robust refit keeps residuals within this multiple of the median absolute deviation. Default 4 * 1.4826
        /// (four sigma of a Gaussian).
        public var robustCutMultiplier: Double
        /// The robust refit runs only when at least this percentage of the pixels are inliers. Default 70.
        public var minimumInlierPercent: Int
        /// The smallest ring sigma the profile fit may take, pixels. Default 0.25.
        public var ringWidthFloor: Double
        /// A ring sigma within this of the floor counts as at the floor. Default 0.05 (so 0.30 px).
        public var ringWidthFloorTolerance: Double
        /// The floor of the one-sided refit (one edge sharper than `ringWidthFloor`). Default 0.01.
        public var loweredRingWidthFloor: Double
        /// The ring-evidence bounds the profile fit reads when it keeps a ring with one sharp edge (its evidence must
        /// pass), and with which `fitBestAvailable(requireRingEvidence: true)` measures the evidence it attaches.
        public var evidence: RingEvidence.Policy

        /// `degeneracyBoundBins` and `sparseFloorBins` follow the bin count (5/6 and 1/3 of it) unless given.
        public init(
            angularBinCount: Int = 36,
            degeneracyBoundBins: Int? = nil,
            sparseFloorBins: Int? = nil,
            strongSignalFraction: Double = 0.2,
            oneRingRatio: Double = 1.10,
            minimumConicSamples: Int = 20,
            minimumProfileSamples: Int = 32,
            minimumRelativeSignal: Double = 1e-6,
            impulseRangeFactor: Double = 2,
            maximumConicResidual: Double = 0.2,
            maximumProfileResidual: Double = 0.22,
            profileSeedResidual: Double = 0.4,
            maximumAxisRatio: Double = 5,
            robustCutMultiplier: Double = 4 * 1.4826,
            minimumInlierPercent: Int = 70,
            ringWidthFloor: Double = 0.25,
            ringWidthFloorTolerance: Double = 0.05,
            loweredRingWidthFloor: Double = 0.01,
            evidence: RingEvidence.Policy = .default
        ) {
            self.angularBinCount = angularBinCount
            self.degeneracyBoundBins = degeneracyBoundBins ?? angularBinCount * 5 / 6
            self.sparseFloorBins = sparseFloorBins ?? angularBinCount / 3
            self.strongSignalFraction = strongSignalFraction
            self.oneRingRatio = oneRingRatio
            self.minimumConicSamples = minimumConicSamples
            self.minimumProfileSamples = minimumProfileSamples
            self.minimumRelativeSignal = minimumRelativeSignal
            self.impulseRangeFactor = impulseRangeFactor
            self.maximumConicResidual = maximumConicResidual
            self.maximumProfileResidual = maximumProfileResidual
            self.profileSeedResidual = profileSeedResidual
            self.maximumAxisRatio = maximumAxisRatio
            self.robustCutMultiplier = robustCutMultiplier
            self.minimumInlierPercent = minimumInlierPercent
            self.ringWidthFloor = ringWidthFloor
            self.ringWidthFloorTolerance = ringWidthFloorTolerance
            self.loweredRingWidthFloor = loweredRingWidthFloor
            self.evidence = evidence
        }

        /// The suite's values, as every release before 0.2.0 used them.
        public static let `default` = Policy()

        /// Whether every field is in its domain: a bin count of 6 to 3600 with 1 <= sparse floor <= degeneracy bound <=
        /// the bin count; fractions in [0, 1); positive finite caps and floors (the lowered floor not above the floor);
        /// a ratio and an axis ratio of at least 1; an inlier percentage in 0...100; sample minimums that leave the
        /// fit determined (at least 6 for the five-parameter conic, 12 for the eleven-parameter profile).
        public var isValid: Bool {
            let bins: Bool = angularBinCount >= 6 && angularBinCount <= 3600
                && sparseFloorBins >= 1 && sparseFloorBins <= degeneracyBoundBins
                && degeneracyBoundBins <= angularBinCount
            let fractions: Bool = strongSignalFraction.isFinite && strongSignalFraction >= 0 && strongSignalFraction < 1
                && minimumRelativeSignal.isFinite && minimumRelativeSignal >= 0
                && impulseRangeFactor.isFinite && impulseRangeFactor >= 0
            let ratios: Bool = oneRingRatio.isFinite && oneRingRatio >= 1
                && maximumAxisRatio.isFinite && maximumAxisRatio >= 1
                && robustCutMultiplier.isFinite && robustCutMultiplier > 0
            let caps: Bool = maximumConicResidual.isFinite && maximumConicResidual > 0
                && maximumProfileResidual.isFinite && maximumProfileResidual > 0
                && profileSeedResidual.isFinite && profileSeedResidual > 0
            let counts: Bool = minimumConicSamples >= 6 && minimumProfileSamples >= 12
                && minimumInlierPercent >= 0 && minimumInlierPercent <= 100
            let floors: Bool = ringWidthFloor.isFinite && ringWidthFloor > 0
                && ringWidthFloorTolerance.isFinite && ringWidthFloorTolerance >= 0
                && loweredRingWidthFloor.isFinite && loweredRingWidthFloor > 0
                && loweredRingWidthFloor <= ringWidthFloor
            return bins && fractions && ratios && caps && counts && floors && evidence.isValid
        }
    }

    static let invalidPolicyDetail = "the ellipse-fit policy has a field outside its domain"
}
