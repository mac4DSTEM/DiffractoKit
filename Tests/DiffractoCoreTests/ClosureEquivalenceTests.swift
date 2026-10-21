//
//  ClosureEquivalenceTests.swift
//  0.9.1: `SymmetryDetection.closedSubset` (the closure check that took 95 % of a 32-atom CIF import) uses a rotation product
//  table and a bucketed translation search. It must return exactly what composing every pair and searching the product's
//  rotation linearly returns. The naive version is written out here, and compared on a closed set (all 192 operations of
//  fcc in the cube), on sets with operations removed and a stray added (the iteration has to drop the operations that
//  no longer compose), and on translations that sit on the bucket edges.
//

import Foundation
import Testing
import simd
@testable import DiffractoCore

struct ClosureEquivalenceTests {

    /// The pre-0.9.1 implementation, verbatim.
    private func naive(_ operations: [SymmetryOperation], tolerance: Double) -> [SymmetryOperation] {
        var set = operations
        while true {
            var byRotation: [IntegerMatrix: [SIMD3<Double>]] = [:]
            for op in set { byRotation[op.rotation, default: []].append(op.translation) }
            func contains(_ op: SymmetryOperation) -> Bool {
                byRotation[op.rotation]?.contains { SymmetryDetection.translationsEqual($0, op.translation, tolerance) } ?? false
            }
            let kept = set.filter { g in g.rotation.isIdentity && SymmetryOperation.isWholeNumber(g.translation)
                || set.allSatisfy { h in contains(g.composed(with: h)) } }
            if kept.count == set.count { return set }
            set = kept
        }
    }

    private func fccOperations() -> [SymmetryOperation] {
        let cell = UnitCell.cubic(a: 4.05)
        let sites = [SIMD3(0.0, 0, 0), SIMD3(0.5, 0.5, 0), SIMD3(0.5, 0, 0.5), SIMD3(0, 0.5, 0.5)].map { AtomSite(z: 13, fractional: $0) }
        return SymmetryDetection.detect(sites: sites, cell: cell).operations
    }

    @Test func closedSetIsReturnedUnchanged() {
        let ops = fccOperations()
        #expect(ops.count == 192)
        #expect(SymmetryDetection.closedSubset(ops, tolerance: 5e-3) == naive(ops, tolerance: 5e-3))
        #expect(SymmetryDetection.closedSubset(ops, tolerance: 5e-3).count == 192)
    }

    @Test func brokenSetsLoseTheSameOperations() {
        let ops = fccOperations()
        var rng = SplitMix(seed: 20261021)
        for trial in 0..<12 {
            var set = ops
            for _ in 0..<(1 + trial % 4) { set.remove(at: Int.random(in: 1..<set.count, using: &rng)) }
            let source = ops[Int.random(in: 0..<ops.count, using: &rng)]
            // A stray: a rotation of the set with a translation no operation has, some on a bucket edge (k / 128, 0.5, 1 - 1e-12).
            let edge: [Double] = [0.13, 0.5 - 1e-12, 0.5 + 1e-12, 1 - 1e-12, 63.0 / 128, 0.2]
            set.append(SymmetryOperation(rotation: source.rotation, translation: SIMD3(edge[trial % 6], 0.25, edge[(trial + 1) % 6])))
            for tolerance in [5e-3, 1e-9, 0.2] {
                #expect(SymmetryDetection.closedSubset(set, tolerance: tolerance) == naive(set, tolerance: tolerance),
                        "trial \(trial), tolerance \(tolerance)")
            }
        }
    }

    /// Two translations of one rotation that are equal within tolerance but sit either side of a bucket boundary (63 / 128), and a pair
    /// straddling 0 / 1, each reached also as a product: the neighbour-bucket and wrap paths.
    @Test func translationsEitherSideOfABucketBoundaryAndOfTheWrap() {
        let ops = fccOperations()
        let r = ops[1].rotation
        for centre in [63.0 / 128, 0.0, 1.0] {
            var set = ops
            set.append(SymmetryOperation(rotation: r, translation: SIMD3(centre - 1e-4, 0.5, 0.25)))
            set.append(SymmetryOperation(rotation: r, translation: SIMD3(centre + 1e-4, 0.5, 0.25)))
            for tolerance in [5e-3, 1e-5, 1e-3] {
                #expect(SymmetryDetection.closedSubset(set, tolerance: tolerance) == naive(set, tolerance: tolerance),
                        "centre \(centre), tolerance \(tolerance)")
            }
        }
    }
}

private struct SplitMix: RandomNumberGenerator {
    var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9E37_79B9_7F4A_7C15
        var z = state
        z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
        z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
        return z ^ (z >> 31)
    }
}
