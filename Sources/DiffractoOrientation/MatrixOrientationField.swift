//
//  MatrixOrientationField.swift
//  Role: the matrix orientation at every scan position, found with the ACOM templates, and the removal of the spots that
//        orientation explains. Step 1 of "find phase orientations from the structure alone" (DiffractoKit#81).
//
//  Why: `PhaseVectorMatcher.map` fits ONE matrix entry for the whole scan, at the zone axis the phase list states. A scan with
//  a second matrix grain keeps every spot of that grain, and the candidate phases are then asked to explain a matrix pattern;
//  the least bad of them wins (on the synthetic Al-Mg-Si demo an Al [011] grain is labelled β″ [001] at 0.007 Å⁻¹, better
//  than a true needle scores). No threshold separates the two: the matrix has to explain its spots first, at its own
//  orientation, position by position. ASTAR, pyxem and py4DSTEM find that orientation by template matching over the
//  fundamental zone; so does this, with the kit's `ACOMMatcher`.
//
//  What it is not: a phase decision. Only the matrix crystal is matched, and only its spots are removed. What is left is
//  the caller's to classify (`PhaseVectorMatcher.map` on the residual scan). No Python reference does this step as one call;
//  it is checked against planted truth (`MatrixOrientationFieldTests`).
//

import Foundation
import simd
import DiffractoCore

/// The matrix orientation of every scan position, row-major.
public struct MatrixOrientationField: Sendable, Equatable {
    public let width: Int
    public let height: Int
    /// The template match of every position (`OrientationResult.empty` where the position has fewer than
    /// ``ACOMMatcher/minimumPeakCount`` peaks): the score is the template correlation, a measure of how well the matrix
    /// templates describe the pattern, not a phase verdict.
    public let results: [OrientationResult]
    /// The orientation matrix of every position (columns: lab x, y, z in crystal coordinates, the convention of
    /// ``ACOMTemplateLibrary/orientationMatrix(templateIndex:inPlaneAngle:mirrored:)``), before symmetry reduction; nil where
    /// there is no match.
    public let orientations: [simd_double3x3?]

    public init(width: Int, height: Int, results: [OrientationResult], orientations: [simd_double3x3?]) {
        self.width = width
        self.height = height
        self.results = results
        self.orientations = orientations
    }
}

/// A scan with the matrix's spots taken out, and what was taken.
public struct MatrixResidualScan: Sendable, Equatable {
    /// The peaks of every position that the matrix orientation of that position does not explain, in the input's
    /// coordinates and order. The direct beam and every peak the matcher would not use stay: they are not the matrix's.
    public let scan: [[SIMD2<Double>]]
    /// How many peaks were removed at each position.
    public let removedCounts: [Int]
    /// The in-plane correction fitted at each position before the removal, radians (0 where none was fitted).
    public let inPlaneCorrectionsRad: [Double]

    public init(scan: [[SIMD2<Double>]], removedCounts: [Int], inPlaneCorrectionsRad: [Double]) {
        self.scan = scan
        self.removedCounts = removedCounts
        self.inPlaneCorrectionsRad = inPlaneCorrectionsRad
    }
}

extension MatrixResidualScan {
    /// The phase map of the residual scan, read as the map of the whole scan.
    ///
    /// `PhaseVectorMatcher.map` on ``scan`` sees only what the matrix did not explain: a position whose spots were all the
    /// matrix's reaches it with no usable vector and is called "no data", and no result counts the spots removed here. This
    /// puts both back: every position's `removedCount` gains the spots removed here, and a "no data" position that lost at
    /// least one spot to the matrix is the matrix (`phaseIndex` the map's matrix phase). Nothing else is changed.
    ///
    /// This is the map's matrix-by-exclusion rule at its defaults. It does not re-apply
    /// `PhaseVectorSettings.minimumMatrixEvidence` or `matrixFallbackExplainedFraction`: the residual map judged those
    /// against its own removed count, which no longer holds the spots removed here. With either above 0 the verdicts of
    /// this function are not the settings' and the caller has to decide them from the summed `removedCount`.
    ///
    /// - Returns: nil when the map is not of this scan (another number of positions).
    public func completing(_ map: PhaseMap) -> PhaseMap? {
        guard map.results.count == removedCounts.count else { return nil }
        var out = map
        for i in out.results.indices where removedCounts[i] > 0 {
            out.results[i].removedCount = Int32(clamping: Int(out.results[i].removedCount) + removedCounts[i])
            if out.results[i].verdict == .noData {
                out.results[i].verdict = .matrix
                out.results[i].phaseIndex = Int32(clamping: map.matrixPhaseIndex)
            }
        }
        return out
    }
}

extension PhaseVectorMatcher {

    /// The matrix orientation at every position of a scan.
    ///
    /// Every position's peaks are matched against the matrix crystal's templates (`matcher.library`, built by
    /// ``ACOMTemplateLibrary/generate(crystal:kMax:settings:isCancelled:progress:)`` over the crystal's fundamental zone).
    /// The peaks carry no intensity here: each counts 1, so the match is by position, as the vector matcher's is.
    /// Only the peaks the vector matcher would use are matched (``experimentalVector(peak:originX:originY:invAngstromPerPixel:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:)``):
    /// the direct beam, a peak past `maximumVectorInvAngstrom` and a non-finite position are left out, so a position needs
    /// ``ACOMMatcher/minimumPeakCount`` usable peaks to have an orientation. There is no score floor: where the matrix is
    /// absent the best template still answers, and ``MatrixOrientationField/results`` carries the score for the caller.
    ///
    /// - Parameters:
    ///   - scan: row-major, `scanWidth` × `scanHeight`, the peak positions of every position in detector pixels.
    ///   - matcher: a matcher on the MATRIX crystal's library.
    ///   - originX: the pattern centre in pixels, the same for every position (as for ``map(scan:scanWidth:scanHeight:library:settings:originX:originY:invAngstromPerPixel:matrixEntryIndex:maximumWorkers:progress:isCancelled:)``).
    ///   - originY: likewise.
    ///   - invAngstromPerPixel: the Q pixel size.
    ///   - directBeamRadiusInvAngstrom: peaks inside it are the direct beam (`PhaseVectorSettings`' default).
    ///   - maximumVectorInvAngstrom: peaks at or beyond it are not used; 0 (the default) uses all.
    ///   - maximumWorkers: the most threads used; nil uses the active processor count. It moves no number.
    ///   - backend: `.cpu` (the reference) or `.metal`.
    ///   - isCancelled: polled between positions; a cancelled run returns nil.
    ///   - progress: the fraction of positions done, from worker threads.
    /// - Returns: nil for a malformed scan (width or height not positive, or a count that is not width · height), when the
    ///   Metal backend is asked for and unavailable, or when the run was cancelled.
    public static func matrixOrientations(scan: [[SIMD2<Double>]], scanWidth: Int, scanHeight: Int,
                                          matcher: ACOMMatcher,
                                          originX: Float, originY: Float, invAngstromPerPixel: Double,
                                          directBeamRadiusInvAngstrom: Double = 0.15, maximumVectorInvAngstrom: Double = 0,
                                          maximumWorkers: Int? = nil, backend: ACOMBackend = .cpu,
                                          isCancelled: (@Sendable () -> Bool)? = nil,
                                          progress: (@Sendable (Double) -> Void)? = nil) -> MatrixOrientationField? {
        let (positions, overflow) = scanWidth.multipliedReportingOverflow(by: scanHeight)
        guard scanWidth > 0, scanHeight > 0, !overflow, scan.count == positions else { return nil }
        // The usable peaks of every position, once: the same rule the removal and the vector matcher apply.
        let usable: [[ACOMPeak]] = scan.map { peaks in
            peaks.compactMap { p in
                experimentalVector(peak: p, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                   directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                                   maximumVectorInvAngstrom: maximumVectorInvAngstrom) == nil
                    ? nil : ACOMPeak(x: Float(p.x), y: Float(p.y), intensity: 1)
            }
        }
        let matched = matcher.matchAll(
            width: scanWidth, height: scanHeight, originX: originX, originY: originY,
            invAngstromPerPixel: invAngstromPerPixel, selection: .full, maxWorkers: maximumWorkers, backend: backend,
            isCancelled: isCancelled, progress: progress,
            peaksAt: { usable[$0] })
        guard let matched else { return nil }
        let library = matcher.library
        let orientations: [simd_double3x3?] = matched.results.map { r in
            guard r.templateIndex >= 0, r.templateIndex < library.count else { return nil }
            return library.orientationMatrix(templateIndex: r.templateIndex, inPlaneAngle: r.inPlaneAngle, mirrored: r.mirrored)
        }
        return MatrixOrientationField(width: scanWidth, height: scanHeight, results: matched.results, orientations: orientations)
    }

    /// The reflections of a crystal that an orientation brings into the pattern, as in-plane vectors (Å⁻¹) in the frame of
    /// ``experimentalVectors(peaks:originX:originY:invAngstromPerPixel:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:)``:
    /// g·x and g·y of the orientation's lab axes, for every reflection whose component along the beam, |g·z|, is within the slab.
    ///
    /// `reflections` must hold −g for every g (a kinematic list does): a mirrored match is the same spot set only then. A
    /// slab wider than the spacing of the reciprocal layers along the beam admits spots of the next layer (for aluminium
    /// down [122] the first layer is 0.082 Å⁻¹ away), which the removal would then treat as the matrix's.
    ///
    /// - Parameters:
    ///   - orientation: columns lab x, y, z in crystal coordinates (``MatrixOrientationField/orientations``).
    ///   - excitationSlabInvAngstrom: half the thickness of the slab through the origin of reciprocal space, Å⁻¹.
    public static func matrixSpots(reflections: [ACOMReflection], orientation: simd_double3x3,
                                   excitationSlabInvAngstrom: Double) -> [SIMD2<Double>] {
        var spots: [SIMD2<Double>] = []
        for r in reflections where abs(simd_dot(r.g, orientation.columns.2)) <= excitationSlabInvAngstrom {
            let p = SIMD2(simd_dot(r.g, orientation.columns.0), simd_dot(r.g, orientation.columns.1))
            if simd_length_squared(p) > 0 { spots.append(p) }
        }
        return spots
    }

    /// The default slab of ``removingMatrixSpots(scan:field:reflections:originX:originY:invAngstromPerPixel:toleranceInvAngstrom:excitationSlabInvAngstrom:directBeamRadiusInvAngstrom:)``:
    /// the templates' own reach along the beam (`ACOMTemplateSettings.sgMax`'s default). The template nearest to the true
    /// beam direction is up to a few degrees off it, which moves a zero-order reflection of length g by g·sin(tilt) along
    /// the beam (0.04 Å⁻¹ at 0.75 Å⁻¹ and 3°); the slab has to hold that.
    public static let matrixRemovalSlabInvAngstrom = 0.1

    /// The most the in-plane angle is corrected before the removal, radians (5°). The template match gives the angle to half
    /// an azimuthal bin (1.4° at 128 bins) and the nearest template's beam direction adds to it; at 0.7 Å⁻¹ one bin is
    /// 0.034 Å⁻¹, more than the removal tolerance, so the angle is fitted to the spots first.
    public static let matrixInPlaneCorrectionLimitRad = 5 * Double.pi / 180

    /// The scan without the spots the matrix explains.
    ///
    /// At every position with an orientation: the matrix's spots at that orientation (``matrixSpots(reflections:orientation:excitationSlabInvAngstrom:)``)
    /// are turned in the plane by the rotation that best lays them on the measured vectors (least squares over the pairs
    /// within three tolerances, twice, at most ``matrixInPlaneCorrectionLimitRad``; none with fewer than two pairs), and
    /// every usable measured vector within `toleranceInvAngstrom` of one of them is removed. A position without an
    /// orientation loses nothing.
    ///
    /// - Parameters:
    ///   - scan: the scan the field was made from, in detector pixels.
    ///   - reflections: the MATRIX crystal's reflections (the list its template library was built from).
    ///   - toleranceInvAngstrom: how close a measured vector must be to a matrix spot to be the matrix's
    ///     (`PhaseVectorSettings.matrixToleranceInvAngstrom`'s default, 0.02).
    ///   - excitationSlabInvAngstrom: see ``matrixSpots(reflections:orientation:excitationSlabInvAngstrom:)``.
    ///   - directBeamRadiusInvAngstrom: vectors inside it are the direct beam and stay (`PhaseVectorSettings`' default).
    ///   - maximumVectorInvAngstrom: vectors at or beyond it are not the matcher's and stay; 0 (the default) uses all.
    /// - Returns: nil when the scan and the field do not have the same number of positions.
    public static func removingMatrixSpots(scan: [[SIMD2<Double>]], field: MatrixOrientationField,
                                           reflections: [ACOMReflection],
                                           originX: Float, originY: Float, invAngstromPerPixel: Double,
                                           toleranceInvAngstrom: Double = 0.02,
                                           excitationSlabInvAngstrom: Double = matrixRemovalSlabInvAngstrom,
                                           directBeamRadiusInvAngstrom: Double = 0.15,
                                           maximumVectorInvAngstrom: Double = 0) -> MatrixResidualScan? {
        guard scan.count == field.orientations.count, scan.count == field.width * field.height else { return nil }
        var out: [[SIMD2<Double>]] = []
        out.reserveCapacity(scan.count)
        var removedCounts = [Int](repeating: 0, count: scan.count)
        var corrections = [Double](repeating: 0, count: scan.count)
        for (position, peaks) in scan.enumerated() {
            guard let orientation = field.orientations[position] else { out.append(peaks); continue }
            var spots = matrixSpots(reflections: reflections, orientation: orientation,
                                    excitationSlabInvAngstrom: excitationSlabInvAngstrom)
            guard !spots.isEmpty else { out.append(peaks); continue }
            // The usable vector of every peak, nil for the direct beam, a peak past the maximum and a non-finite position (they stay).
            let vectors: [SIMD2<Double>?] = peaks.map {
                experimentalVector(peak: $0, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                   directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                                   maximumVectorInvAngstrom: maximumVectorInvAngstrom)
            }
            func nearest(_ u: SIMD2<Double>) -> (spot: SIMD2<Double>, distance: Double)? {
                var best: (SIMD2<Double>, Double)?
                for s in spots {
                    let d = simd_length(u - s)
                    if d < (best?.1 ?? .infinity) { best = (s, d) }
                }
                return best.map { (spot: $0.0, distance: $0.1) }
            }
            // The in-plane rotation that lays the matrix spots on the measured vectors (`fittingInPlane`).
            let fitted = fittingInPlane(spots: spots, vectors: vectors.compactMap { $0 }, toleranceInvAngstrom: toleranceInvAngstrom)
            spots = fitted.spots
            corrections[position] = fitted.rotationRad
            var kept: [SIMD2<Double>] = []
            kept.reserveCapacity(peaks.count)
            for (peak, vector) in zip(peaks, vectors) {
                if let u = vector, let n = nearest(u), n.distance <= toleranceInvAngstrom {
                    removedCounts[position] += 1
                } else {
                    kept.append(peak)
                }
            }
            out.append(kept)
        }
        return MatrixResidualScan(scan: out, removedCounts: removedCounts, inPlaneCorrectionsRad: corrections)
    }
}

// MARK: - One orientation per grain

/// A matrix orientation field after ``PhaseVectorMatcher/regularizingMatrixOrientations(_:scan:matcher:originX:originY:invAngstromPerPixel:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:neighbourhoodRadius:keepFraction:)``.
public struct MatrixOrientationRegularization: Sendable, Equatable {
    public let field: MatrixOrientationField
    /// The positions (row-major indices, ascending) whose orientation was replaced by their neighbourhood's.
    public let changedPositions: [Int]

    public init(field: MatrixOrientationField, changedPositions: [Int]) {
        self.field = field
        self.changedPositions = changedPositions
    }
}

extension PhaseVectorMatcher {

    /// The usable peaks of one position for the template match: the rule of
    /// ``experimentalVector(peak:originX:originY:invAngstromPerPixel:directBeamRadiusInvAngstrom:maximumVectorInvAngstrom:)``,
    /// each peak counting 1.
    static func usableACOMPeaks(_ peaks: [SIMD2<Double>], originX: Float, originY: Float, invAngstromPerPixel: Double,
                                directBeamRadiusInvAngstrom: Double, maximumVectorInvAngstrom: Double) -> [ACOMPeak] {
        peaks.compactMap { p in
            experimentalVector(peak: p, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                               directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                               maximumVectorInvAngstrom: maximumVectorInvAngstrom) == nil
                ? nil : ACOMPeak(x: Float(p.x), y: Float(p.y), intensity: 1)
        }
    }

    /// The field with the matrix orientation held to the grain where the pattern allows it.
    ///
    /// A grain has one orientation, but a pattern full of precipitate spots can pull the template match of a single position
    /// to another zone (on Thronsen et al. dataset A about 1 % of the positions of a one-grain scan land 40° from the
    /// matrix zone). For every position whose zone disagrees with its neighbourhood:
    ///
    /// 1. the neighbourhood's zone is the zone that most of the other positions in the (2r + 1)² window share (beam
    ///    directions within the library's distinct angle of one another); it must be shared by more than half of ALL the
    ///    other positions of the window, with or without an orientation;
    /// 2. the position's OWN pattern is matched against the templates around that zone
    ///    (``ACOMMatcher/match(peaks:originX:originY:invAngstromPerPixel:near:withinRadians:)``);
    /// 3. the position takes that match only if its score is at least `keepFraction` of the position's own best score.
    ///
    /// So nothing is decided by the vote alone: a position whose own best is clearly better keeps it. A grain boundary
    /// stays where it is, and a grain of a few positions inside another stays, because their patterns do not fit the
    /// neighbours' zone. The neighbourhood is read from the input field, so the result does not depend on the order of
    /// the positions. A replaced position carries the template's bin angle (the matcher's in-plane refinement is not applied
    /// to it); the removal fits the angle anyway.
    ///
    /// - Parameters:
    ///   - scan: the scan the field was made from, in detector pixels.
    ///   - matcher: the matcher the field was made with.
    ///   - neighbourhoodRadius: r, in positions; 2 is a 5 × 5 window.
    ///   - keepFraction: see step 3; in 0...1.
    /// - Returns: nil when the scan and the field do not have the same number of positions, or a result names a template
    ///   the matcher's library does not hold.
    public static func regularizingMatrixOrientations(_ field: MatrixOrientationField, scan: [[SIMD2<Double>]],
                                                      matcher: ACOMMatcher,
                                                      originX: Float, originY: Float, invAngstromPerPixel: Double,
                                                      directBeamRadiusInvAngstrom: Double = 0.15,
                                                      maximumVectorInvAngstrom: Double = 0,
                                                      neighbourhoodRadius: Int = 2,
                                                      keepFraction: Float = 0.9) -> MatrixOrientationRegularization? {
        let width = field.width, height = field.height
        guard scan.count == field.results.count, scan.count == field.orientations.count, neighbourhoodRadius >= 1 else { return nil }
        let (positions, overflow) = width.multipliedReportingOverflow(by: height)
        guard width > 0, height > 0, !overflow, positions == scan.count else { return nil }
        let library = matcher.library
        guard field.results.allSatisfy({ $0.templateIndex < library.count }) else { return nil }
        let cosDistinct = cos(library.distinctOrientationRadians)
        func zone(_ i: Int) -> SIMD3<Double>? {
            let t = field.results[i].templateIndex
            return field.orientations[i] != nil && t >= 0 ? library.zoneAxes[t] : nil
        }
        var results = field.results, orientations = field.orientations
        var changed: [Int] = []
        for row in 0..<height {
            for column in 0..<width {
                let i = row * width + column
                guard let own = zone(i) else { continue }
                // The other positions of the window that have an orientation.
                var neighbours: [SIMD3<Double>] = []
                var others = 0
                for r in max(0, row - min(neighbourhoodRadius, height))...min(height - 1, row + min(neighbourhoodRadius, height)) {
                    for c in max(0, column - min(neighbourhoodRadius, width))...min(width - 1, column + min(neighbourhoodRadius, width)) where !(r == row && c == column) {
                        others += 1
                        if let z = zone(r * width + c) { neighbours.append(z) }
                    }
                }
                guard !neighbours.isEmpty else { continue }
                var dominant = neighbours[0], dominantCount = 0
                for z in neighbours {
                    let count = neighbours.reduce(0) { $0 + (simd_dot($1, z) > cosDistinct ? 1 : 0) }
                    if count > dominantCount { dominant = z; dominantCount = count }
                }
                guard 2 * dominantCount > others, simd_dot(own, dominant) <= cosDistinct else { continue }
                let peaks = usableACOMPeaks(scan[i], originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                            directBeamRadiusInvAngstrom: directBeamRadiusInvAngstrom,
                                            maximumVectorInvAngstrom: maximumVectorInvAngstrom)
                let candidate = matcher.match(peaks: peaks, originX: originX, originY: originY, invAngstromPerPixel: invAngstromPerPixel,
                                              near: dominant, withinRadians: library.distinctOrientationRadians)
                guard candidate.templateIndex >= 0, candidate.score >= keepFraction * field.results[i].score else { continue }
                results[i] = candidate
                orientations[i] = library.orientationMatrix(templateIndex: candidate.templateIndex,
                                                            inPlaneAngle: candidate.inPlaneAngle, mirrored: candidate.mirrored)
                changed.append(i)
            }
        }
        return MatrixOrientationRegularization(
            field: MatrixOrientationField(width: width, height: height, results: results, orientations: orientations),
            changedPositions: changed)
    }
}
