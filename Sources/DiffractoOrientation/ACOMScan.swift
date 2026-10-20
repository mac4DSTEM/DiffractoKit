//
//  ACOMScan.swift
//  Role: ACOM matching over a whole scan on the CPU: which scan positions to match (`ACOMScanSelection`), the matched grid
//        (`ACOMScanResults`) and `ACOMMatcher.matchAll`, which spreads the positions over worker threads and reports progress
//        and honours cancellation. The per-pattern numbers are `ACOMMatcher.match`; this file only schedules them (the apps'
//        `OrientationMatching.matchAll` CPU path, `ACOMScanSelection`, `OrientationMap` container without its presentation).
//

import Foundation

/// Which scan positions an ACOM run matches. Positions that are not selected stay `OrientationResult.empty`, except in
/// `preview`, where each matched position is copied over its coarse block.
public enum ACOMScanSelection: Sendable, Equatable {
    /// Every position.
    case full
    /// One representative position (the block's top-left corner) per square block of the step that fits the longer scan side
    /// into `maxDimension` blocks; the result is expanded over the blocks by ``ACOMScanResults/expandedPreview(selection:)``.
    case preview(maxDimension: Int)
    /// The positions within `radius` of the centre in x and y (a square of side `2 * radius + 1`), clipped to the scan.
    case square(centerX: Int, centerY: Int, radius: Int)

    /// The selected positions, row-major (`y * width + x`), ascending. Empty for an empty scan or a square outside it.
    public func sourceIndices(width: Int, height: Int) -> [Int] {
        guard width > 0, height > 0 else { return [] }
        switch self {
        case .full:
            return Array(0..<(width * height))
        case .preview(let maxDimension):
            let step = previewStep(width: width, height: height, maxDimension: maxDimension)
            var indices: [Int] = []
            indices.reserveCapacity(((width + step - 1) / step) * ((height + step - 1) / step))
            for y in stride(from: 0, to: height, by: step) {
                for x in stride(from: 0, to: width, by: step) { indices.append(y * width + x) }
            }
            return indices
        case .square(let centerX, let centerY, let radius):
            // Clamped to +-2^40 so the sums cannot overflow; positions beyond that are far outside any scan anyway.
            let lim = 1 << 40
            let r = min(max(0, radius), lim), cx = min(max(centerX, -lim), lim), cy = min(max(centerY, -lim), lim)
            let x0 = max(0, cx - r), x1 = min(width - 1, cx + r)
            let y0 = max(0, cy - r), y1 = min(height - 1, cy + r)
            guard x0 <= x1, y0 <= y1 else { return [] }
            var indices: [Int] = []
            indices.reserveCapacity((x1 - x0 + 1) * (y1 - y0 + 1))
            for y in y0...y1 { for x in x0...x1 { indices.append(y * width + x) } }
            return indices
        }
    }

    /// How many patterns a run matches.
    public func positionCount(width: Int, height: Int) -> Int { sourceIndices(width: width, height: height).count }

    /// The block side of `preview` for this scan (1 for the other selections).
    public func previewStep(width: Int, height: Int) -> Int {
        guard case .preview(let maxDimension) = self else { return 1 }
        return previewStep(width: width, height: height, maxDimension: maxDimension)
    }

    private func previewStep(width: Int, height: Int, maxDimension: Int) -> Int {
        let target = max(1, maxDimension)
        let longer = max(width, height)
        return max(1, longer / target + (longer % target == 0 ? 0 : 1))
    }
}

/// The matched scan: one `OrientationResult` per position, row-major.
public struct ACOMScanResults: Sendable, Equatable {
    public let width: Int
    public let height: Int
    /// Templates in the library the run used.
    public let templateCount: Int
    public var results: [OrientationResult]

    /// All positions `OrientationResult.empty`.
    public init(width: Int, height: Int, templateCount: Int) {
        self.width = max(0, width)
        self.height = max(0, height)
        self.templateCount = templateCount
        self.results = Array(repeating: .empty, count: self.width * self.height)
    }

    public subscript(x: Int, y: Int) -> OrientationResult {
        get { results[y * width + x] }
        set { results[y * width + x] = newValue }
    }

    /// Copies each block's top-left result over the whole block, for a `preview` selection; any other selection returns the
    /// results unchanged.
    public func expandedPreview(selection: ACOMScanSelection) -> ACOMScanResults {
        guard case .preview = selection else { return self }
        let step = selection.previewStep(width: width, height: height)
        var expanded = self
        for y in 0..<height {
            let sourceY = (y / step) * step
            for x in 0..<width { expanded[x, y] = self[(x / step) * step, sourceY] }
        }
        return expanded
    }
}

extension ACOMMatcher {
    /// Matches a scan of patterns on worker threads.
    ///
    /// - Parameters:
    ///   - peaks: the detected peaks of every position, row-major, `width * height` entries.
    ///   - originX: the pattern centre, pixels, the same for every position.
    ///   - originY: likewise.
    ///   - invAngstromPerPixel: the Q pixel size, as in ``match(peaks:originX:originY:invAngstromPerPixel:)``.
    ///   - selection: the positions to match; the others stay empty (`preview` is expanded over its blocks).
    ///   - maxWorkers: the most threads used, at least 1; nil uses the active processor count. The positions are dealt to the
    ///     workers round-robin, so a spatially uneven peak count does not leave one worker with the slow rows.
    ///   - isCancelled: polled before every position; once it returns true no further position is started.
    ///   - progress: called from worker threads with the fraction of the selected positions done, in 1 % batches, ending at 1.
    /// - Returns: nil when `peaks.count != width * height`, the selection is empty, or the run was cancelled. A cancelled
    ///   run returns nil even if every position had finished by the time it was noticed.
    public func matchAll(peaks: [[ACOMPeak]], width: Int, height: Int,
                         originX: Float, originY: Float, invAngstromPerPixel: Double = 1,
                         selection: ACOMScanSelection = .full, maxWorkers: Int? = nil,
                         isCancelled: (@Sendable () -> Bool)? = nil,
                         progress: (@Sendable (Double) -> Void)? = nil) -> ACOMScanResults? {
        guard width > 0, height > 0, peaks.count == width * height else { return nil }
        if isCancelled?() == true { return nil }
        let selected = selection.sourceIndices(width: width, height: height)
        guard !selected.isEmpty else { return nil }

        let counter = Progress(total: selected.count)
        let workers = max(1, min(maxWorkers ?? ProcessInfo.processInfo.activeProcessorCount, selected.count))
        let batch = max(1, selected.count / (workers * 100))
        let results = SlotBuffer(count: width * height)
        DispatchQueue.concurrentPerform(iterations: workers) { worker in
            var pending = 0
            for k in stride(from: worker, to: selected.count, by: workers) {
                if isCancelled?() == true { break }
                let source = selected[k]
                results.write(match(peaks: peaks[source], originX: originX, originY: originY,
                                    invAngstromPerPixel: invAngstromPerPixel), at: source)
                pending += 1
                if progress != nil, pending >= batch {
                    counter.add(pending, report: progress)
                    pending = 0
                }
            }
            if progress != nil, pending > 0 { counter.add(pending, report: progress) }
        }
        if isCancelled?() == true { return nil }
        var map = ACOMScanResults(width: width, height: height, templateCount: library.count)
        map.results = results.take()
        return map.expandedPreview(selection: selection)
    }

    /// The finished-position count behind the progress callback.
    private final class Progress: @unchecked Sendable {
        private let lock = NSLock()
        private var done = 0
        private let total: Int
        init(total: Int) { self.total = total }
        /// Adds and reports under the lock, so the callback sees non-decreasing fractions and the last one is 1.
        func add(_ n: Int, report: (@Sendable (Double) -> Void)?) {
            lock.lock(); defer { lock.unlock() }
            done += n
            report?(Double(done) / Double(total))
        }
    }

    /// Disjoint-slot writes from several threads into one buffer (raw memory, so Swift's exclusivity checks on an array
    /// property do not fire for writers that never share an index).
    private final class SlotBuffer: @unchecked Sendable {
        private let base: UnsafeMutablePointer<OrientationResult>
        private let count: Int
        init(count: Int) {
            self.count = count
            base = UnsafeMutablePointer<OrientationResult>.allocate(capacity: max(count, 1))
            base.initialize(repeating: .empty, count: count)
        }
        deinit {
            base.deinitialize(count: count)
            base.deallocate()
        }
        /// Each index is written by exactly one worker.
        func write(_ value: OrientationResult, at index: Int) { base[index] = value }
        func take() -> [OrientationResult] { Array(UnsafeBufferPointer(start: base, count: count)) }
    }
}
