//
//  ACOMMetal.swift
//  Role: the Metal option of `ACOMMatcher.matchAll` (`ACOMBackend.metal`). The polar image and the forward ring FFTs of each
//        pattern stay on the CPU (`experimentalFFT`, the same code as the reference); the GPU computes, for every position and
//        template, the correlation spectrum Σ_ring conj(E)·T and its inverse transform over the azimuth, and returns the best
//        bin and its score (the forward and the mirrored pass); the CPU then applies the reference's selection of the winner and
//        the distinct runner-up and its symmetry reduction. Ported in structure from the mac4DSTEM Metal matcher
//        (`MetalACOMMatcher`, `ACOMMatching.metal`), with two changes: the inverse transform reads a precomputed Float cos/sin
//        table indexed by `(f * a) mod n` instead of calling sin/cos on an argument up to 2π·n (which loses 3 to 4 digits in
//        Float), and the shader is compiled at run time from the source below with safe (non-fast) math, so the package needs no
//        Metal build step and compiles unchanged where Metal is not used.
//

import Foundation
import DiffractoCore
#if canImport(Metal)
import Metal
#endif

extension ACOMMatcher {
    /// Whether a Metal device is present and this library's azimuthal length (a power of two up to 1024) can run on it.
    public var isMetalAvailable: Bool {
        #if canImport(Metal)
        return MTLCreateSystemDefaultDevice() != nil && library.geometry.nAzimuthal <= 1024
            && library.geometry.nAzimuthal & (library.geometry.nAzimuthal - 1) == 0
        #else
        return false
        #endif
    }

    /// Matches `selected` positions with the GPU correlation; results in the order of `selected`. nil when Metal is unavailable, a
    /// GPU step fails or the run is cancelled.
    func matchSelectedWithMetal(peaks: [[ACOMPeak]], selected: [Int], originX: Float, originY: Float,
                                invAngstromPerPixel: Double, maxWorkers: Int?,
                                isCancelled: (@Sendable () -> Bool)?,
                                progress: (@Sendable (Double) -> Void)?) -> [OrientationResult]? {
        #if canImport(Metal)
        guard isMetalAvailable, library.count > 0, let gpu = ACOMMetalContext(library: library) else { return nil }
        let na = library.geometry.nAzimuthal, nr = library.geometry.nRadial, grid = nr * na
        let workers = max(1, min(maxWorkers ?? ProcessInfo.processInfo.activeProcessorCount, 64))
        var results = [OrientationResult](repeating: .empty, count: selected.count)
        var done = 0
        for lower in stride(from: 0, to: selected.count, by: gpu.batchSize) {
            if isCancelled?() == true { return nil }
            let upper = min(selected.count, lower + gpu.batchSize), count = upper - lower
            // CPU preparation, the reference's own code, spread over the workers (disjoint slices).
            let slices = Slices(count: count, grid: grid)
            DispatchQueue.concurrentPerform(iterations: min(workers, count)) { w in
                for local in stride(from: w, to: count, by: min(workers, count)) {
                    guard let f = experimentalFFT(peaks: peaks[selected[lower + local]], originX: originX, originY: originY,
                                                  invAngstromPerPixel: invAngstromPerPixel) else { continue }
                    slices.write(local: local, re: f.re, im: f.im)
                }
            }
            if isCancelled?() == true { return nil }
            let valid = slices.isValid
            guard let out = gpu.run(re: slices.real, im: slices.imaginary, positions: count) else { return nil }
            for local in 0..<count where valid[local] {
                let o = local * library.count
                var scores = Array(out.forwardScores[o..<(o + library.count)])
                var bins = Array(out.forwardBins[o..<(o + library.count)])
                var mirrored = [Bool](repeating: false, count: library.count)
                for t in 0..<library.count where out.mirrorScores[o + t] > scores[t] {   // strictly larger, as the reference
                    scores[t] = out.mirrorScores[o + t]; bins[t] = out.mirrorBins[o + t]; mirrored[t] = true
                }
                let pick = Self.selectOrientation(zoneAxes: library.zoneAxes, scores: scores, bins: bins,
                                                  distinctOrientationRadians: library.distinctOrientationRadians)
                guard pick.template >= 0 else { continue }
                results[lower + local] = makeResult(bestTemplate: pick.template, bestBin: pick.bin, bestScore: pick.score,
                                                    secondScore: pick.secondScore, mirrored: mirrored[pick.template])
            }
            done += count
            progress?(Double(done) / Double(selected.count))
        }
        return isCancelled?() == true ? nil : results
        #else
        return nil
        #endif
    }
}

#if canImport(Metal)
/// A batch's experimental ring spectra, written from several threads into disjoint slices of memory it owns.
private final class Slices: @unchecked Sendable {
    private let re: UnsafeMutablePointer<Float>, im: UnsafeMutablePointer<Float>, valid: UnsafeMutablePointer<Bool>
    private let count: Int, grid: Int
    init(count: Int, grid: Int) {
        self.count = count; self.grid = grid
        re = .allocate(capacity: max(1, count * grid)); re.initialize(repeating: 0, count: count * grid)
        im = .allocate(capacity: max(1, count * grid)); im.initialize(repeating: 0, count: count * grid)
        valid = .allocate(capacity: max(1, count)); valid.initialize(repeating: false, count: count)
    }
    deinit { re.deallocate(); im.deallocate(); valid.deallocate() }
    /// Each `local` is written by exactly one worker.
    func write(local: Int, re source: [Float], im sourceIm: [Float]) {
        for i in 0..<grid { re[local * grid + i] = source[i]; im[local * grid + i] = sourceIm[i] }
        valid[local] = true
    }
    var real: [Float] { Array(UnsafeBufferPointer(start: re, count: count * grid)) }
    var imaginary: [Float] { Array(UnsafeBufferPointer(start: im, count: count * grid)) }
    var isValid: [Bool] { Array(UnsafeBufferPointer(start: valid, count: count)) }
}

/// One matchAll's Metal state: device, pipelines, the template spectra, the twiddle tables.
private final class ACOMMetalContext {
    let device: MTLDevice, queue: MTLCommandQueue
    let spectrum: MTLComputePipelineState, argmax: MTLComputePipelineState
    let templateRe: MTLBuffer, templateIm: MTLBuffer, cosTable: MTLBuffer, sinTable: MTLBuffer
    let templates: Int, radial: Int, azimuthal: Int
    let batchSize: Int

    struct Params { var positions: UInt32, templates: UInt32, radial: UInt32, azimuthal: UInt32, sign: Float }

    init?(library: ACOMTemplateLibrary) {
        guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else { return nil }
        let options = MTLCompileOptions()
        options.mathMode = .safe
        guard let lib = try? device.makeLibrary(source: Self.source, options: options),
              let f1 = lib.makeFunction(name: "acomCorrelationSpectrum"), let f2 = lib.makeFunction(name: "acomInverseArgmax"),
              let p1 = try? device.makeComputePipelineState(function: f1), let p2 = try? device.makeComputePipelineState(function: f2)
        else { return nil }
        let n = library.geometry.nAzimuthal
        guard n <= p2.maxTotalThreadsPerThreadgroup, n & (n - 1) == 0 else { return nil }
        let bytes = library.templateFFTRe.count * MemoryLayout<Float>.stride
        // The shader indexes with 32-bit unsigned integers.
        guard bytes > 0, bytes <= device.maxBufferLength, library.templateFFTRe.count < 1 << 32,
              let tr = device.makeBuffer(bytes: library.templateFFTRe, length: bytes, options: .storageModeShared),
              let ti = device.makeBuffer(bytes: library.templateFFTIm, length: bytes, options: .storageModeShared)
        else { return nil }
        // Float tables of the Double cos/sin of 2π k / n, k = 0..<n (the inverse DFT's twiddles, rounded once).
        let cosines = (0..<n).map { Float(cos(2 * Double.pi * Double($0) / Double(n))) }
        let sines = (0..<n).map { Float(sin(2 * Double.pi * Double($0) / Double(n))) }
        guard let ct = device.makeBuffer(bytes: cosines, length: n * 4, options: .storageModeShared),
              let st = device.makeBuffer(bytes: sines, length: n * 4, options: .storageModeShared) else { return nil }
        self.device = device; self.queue = queue; self.spectrum = p1; self.argmax = p2
        templateRe = tr; templateIm = ti; cosTable = ct; sinTable = st
        templates = library.count; radial = library.geometry.nRadial; azimuthal = n
        // Per position: both spectra, the scores and bins of two passes, and the experimental grid (the reference app's budget).
        let perPosition = templates * n * 2 * 4 + 2 * templates * 8 + radial * n * 2 * 4
        let budget = min(Int(device.recommendedMaxWorkingSetSize / 8), 128 << 20)
        batchSize = max(1, min(64, budget / max(1, perPosition)))
    }

    /// Both passes (forward = conj(E)·T, mirrored = E·T) for `positions` experimental ring spectra.
    func run(re: [Float], im: [Float], positions: Int)
        -> (forwardScores: [Float], forwardBins: [UInt32], mirrorScores: [Float], mirrorBins: [UInt32])? {
        autoreleasepool { runUnpooled(re: re, im: im, positions: positions) }
    }

    private func runUnpooled(re: [Float], im: [Float], positions: Int)
        -> (forwardScores: [Float], forwardBins: [UInt32], mirrorScores: [Float], mirrorBins: [UInt32])? {
        let spectrumCount = positions * templates * azimuthal, resultCount = positions * templates
        guard let inRe = device.makeBuffer(bytes: re, length: re.count * 4, options: .storageModeShared),
              let inIm = device.makeBuffer(bytes: im, length: im.count * 4, options: .storageModeShared),
              let specRe = device.makeBuffer(length: spectrumCount * 4, options: .storageModePrivate),
              let specIm = device.makeBuffer(length: spectrumCount * 4, options: .storageModePrivate),
              let scores0 = device.makeBuffer(length: resultCount * 4, options: .storageModeShared),
              let bins0 = device.makeBuffer(length: resultCount * 4, options: .storageModeShared),
              let scores1 = device.makeBuffer(length: resultCount * 4, options: .storageModeShared),
              let bins1 = device.makeBuffer(length: resultCount * 4, options: .storageModeShared),
              let command = queue.makeCommandBuffer() else { return nil }
        for (pass, scores, bins) in [(0, scores0, bins0), (1, scores1, bins1)] {
            var params = Params(positions: UInt32(positions), templates: UInt32(templates), radial: UInt32(radial),
                                azimuthal: UInt32(azimuthal), sign: pass == 0 ? 1 : -1)
            guard let e1 = command.makeComputeCommandEncoder() else { return nil }
            e1.setComputePipelineState(spectrum)
            e1.setBuffer(inRe, offset: 0, index: 0); e1.setBuffer(inIm, offset: 0, index: 1)
            e1.setBuffer(templateRe, offset: 0, index: 2); e1.setBuffer(templateIm, offset: 0, index: 3)
            e1.setBuffer(specRe, offset: 0, index: 4); e1.setBuffer(specIm, offset: 0, index: 5)
            e1.setBytes(&params, length: MemoryLayout<Params>.stride, index: 6)
            e1.dispatchThreads(MTLSize(width: spectrumCount, height: 1, depth: 1),
                               threadsPerThreadgroup: MTLSize(width: min(spectrum.maxTotalThreadsPerThreadgroup, 256), height: 1, depth: 1))
            e1.endEncoding()
            guard let e2 = command.makeComputeCommandEncoder() else { return nil }
            e2.setComputePipelineState(argmax)
            e2.setBuffer(specRe, offset: 0, index: 0); e2.setBuffer(specIm, offset: 0, index: 1)
            e2.setBuffer(scores, offset: 0, index: 2); e2.setBuffer(bins, offset: 0, index: 3)
            e2.setBuffer(cosTable, offset: 0, index: 4); e2.setBuffer(sinTable, offset: 0, index: 5)
            e2.setBytes(&params, length: MemoryLayout<Params>.stride, index: 6)
            e2.setThreadgroupMemoryLength(azimuthal * 4, index: 0)
            e2.setThreadgroupMemoryLength(azimuthal * 4, index: 1)
            e2.dispatchThreadgroups(MTLSize(width: 1, height: templates, depth: positions),
                                    threadsPerThreadgroup: MTLSize(width: azimuthal, height: 1, depth: 1))
            e2.endEncoding()
        }
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else { return nil }
        func floats(_ b: MTLBuffer) -> [Float] { Array(UnsafeBufferPointer(start: b.contents().bindMemory(to: Float.self, capacity: resultCount), count: resultCount)) }
        func uints(_ b: MTLBuffer) -> [UInt32] { Array(UnsafeBufferPointer(start: b.contents().bindMemory(to: UInt32.self, capacity: resultCount), count: resultCount)) }
        return (floats(scores0), uints(bins0), floats(scores1), uints(bins1))
    }

    /// The two kernels. Pass 0 (`sign` +1) is conj(E)·T, pass 1 (`sign` -1) is E·T, as ACOMMatcher.correlate. The inverse DFT
    /// is exact over the azimuth (`n` power of two); ties between bins go to the lowest bin, as the reference's first-maximum scan.
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    struct ACOMParams { uint positions; uint templates; uint radial; uint azimuthal; float sign; };
    kernel void acomCorrelationSpectrum(
        const device float *er [[buffer(0)]], const device float *ei [[buffer(1)]],
        const device float *tr [[buffer(2)]], const device float *ti [[buffer(3)]],
        device float *cr [[buffer(4)]], device float *ci [[buffer(5)]],
        constant ACOMParams &p [[buffer(6)]], uint index [[thread_position_in_grid]])
    {
        if (index >= p.positions * p.templates * p.azimuthal) return;
        const uint f = index % p.azimuthal;
        const uint t = (index / p.azimuthal) % p.templates;
        const uint pos = index / (p.azimuthal * p.templates);
        const uint eBase = pos * p.radial * p.azimuthal, tBase = t * p.radial * p.azimuthal;
        float re = 0.0f, im = 0.0f;
        for (uint r = 0; r < p.radial; ++r) {
            const uint o = r * p.azimuthal + f;
            const float a = er[eBase + o], b = ei[eBase + o], c = tr[tBase + o], d = ti[tBase + o];
            re += a * c + p.sign * b * d;
            im += a * d - p.sign * b * c;
        }
        cr[index] = re; ci[index] = im;
    }
    kernel void acomInverseArgmax(
        const device float *cr [[buffer(0)]], const device float *ci [[buffer(1)]],
        device float *bestScores [[buffer(2)]], device uint *bestBins [[buffer(3)]],
        const device float *cosT [[buffer(4)]], const device float *sinT [[buffer(5)]],
        constant ACOMParams &p [[buffer(6)]],
        threadgroup float *scores [[threadgroup(0)]], threadgroup uint *bins [[threadgroup(1)]],
        uint a [[thread_index_in_threadgroup]], uint3 group [[threadgroup_position_in_grid]])
    {
        const uint pos = group.z, t = group.y;
        if (pos >= p.positions || t >= p.templates || a >= p.azimuthal) return;
        const uint base = (pos * p.templates + t) * p.azimuthal;
        const uint mask = p.azimuthal - 1;
        float value = 0.0f;
        for (uint f = 0; f < p.azimuthal; ++f) {
            const uint k = (f * a) & mask;
            value += cr[base + f] * cosT[k] - ci[base + f] * sinT[k];
        }
        scores[a] = value / float(p.azimuthal);
        bins[a] = a;
        threadgroup_barrier(mem_flags::mem_threadgroup);
        for (uint stride = p.azimuthal / 2; stride > 0; stride >>= 1) {
            if (a < stride && scores[a + stride] > scores[a]) { scores[a] = scores[a + stride]; bins[a] = bins[a + stride]; }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
        if (a == 0) { bestScores[pos * p.templates + t] = scores[0]; bestBins[pos * p.templates + t] = bins[0]; }
    }
    """
}
#endif
