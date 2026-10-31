//
//  ACOMMetalArgmaxTieTests.swift
//  Role: the Metal inverse-transform argmax keeps the reference rule of `ACOMMatcher.correlate` (NumPy argmax, as py4DSTEM):
//        on an exact tie between in-plane bins the lowest bin wins, and NaN never wins (an all-NaN row gives bin 0 and the
//        CPU start value). The correlation is planted: one ring, an experimental spectrum of ones (so both passes equal the
//        template spectrum) and a template spectrum whose inverse transform over 4 bins is exactly 16 x [0, 1, 1, 0] in Float
//        arithmetic, a tie between bins 1 and 2. Skipped where Metal is unavailable.
//

import Foundation
import Testing
import DiffractoOrientation

#if canImport(Metal)
import Metal
private let hasMetal = MTLCreateSystemDefaultDevice() != nil
#else
private let hasMetal = false
#endif

/// The reference rule on the CPU path (`ACOMMatcher.correlate`): start below every finite value, a strictly larger value wins.
private func firstMaximum(_ values: [Float]) -> (score: Float, bin: Int) {
    var best: Float = -.greatestFiniteMagnitude, bin = 0
    for a in 0..<values.count where values[a] > best { best = values[a]; bin = a }
    return (best, bin)
}

struct ACOMMetalArgmaxTieTests {
    @Test(.enabled(if: hasMetal, "no Metal device"))
    func anExactTieBetweenBinsGoesToTheLowestBinOnMetalAsOnTheCPU() throws {
        // x = [0, 4, 4, 0] over 4 bins; its forward DFT X_k = sum x_a exp(-2 pi i a k / 4).
        let x: [Float] = [0, 4, 4, 0]
        let templateRe: [Float] = [8, -4, 0, -4], templateIm: [Float] = [0, -4, 0, 4]
        let cpu = firstMaximum(x)    // the inverse transform divided by n: 16 x [0, 1, 1, 0] / 4
        #expect(cpu.bin == 1 && cpu.score == 4)
        let gpu = try #require(ACOMMetalKernels.bestBins(experimentalRe: [1, 1, 1, 1], experimentalIm: [0, 0, 0, 0],
                                                         templateRe: templateRe, templateIm: templateIm,
                                                         templates: 1, radial: 1, azimuthal: 4, positions: 1))
        #expect(gpu.forwardScores == [cpu.score], "the planted inverse transform is exact")
        #expect(gpu.forwardBins == [UInt32(cpu.bin)], "forward pass: the lowest of the tied bins")
        #expect(gpu.mirrorScores == [cpu.score])
        #expect(gpu.mirrorBins == [UInt32(cpu.bin)], "mirrored pass: the lowest of the tied bins")
    }

    @Test(.enabled(if: hasMetal, "no Metal device"))
    func anAllNaNCorrelationGivesBinZeroAndTheCPUStartValueOnMetal() throws {
        let cpu = firstMaximum([Float](repeating: .nan, count: 4))
        let gpu = try #require(ACOMMetalKernels.bestBins(experimentalRe: [1, 1, 1, 1], experimentalIm: [0, 0, 0, 0],
                                                         templateRe: [.nan, 0, 0, 0], templateIm: [0, 0, 0, 0],
                                                         templates: 1, radial: 1, azimuthal: 4, positions: 1))
        #expect(gpu.forwardBins == [UInt32(cpu.bin)] && gpu.mirrorBins == [UInt32(cpu.bin)])
        // The CPU path scales its start value by 1/n after the scan (-greatestFiniteMagnitude / 4 is exact).
        #expect(gpu.forwardScores == [cpu.score / 4] && gpu.mirrorScores == [cpu.score / 4])
    }
}
