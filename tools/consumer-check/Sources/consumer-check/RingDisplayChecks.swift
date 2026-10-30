//
//  RingDisplayChecks.swift
//  Role: the public API of the macTEM pieces (CrystalLoading, RingOverlay, RingScale, IndexingDisplay) used from outside
//        the package, against values a reader verifies by hand. `runRingDisplayChecks` is called from main.swift.
//

import Foundation
import simd
import DiffractoCore
import DiffractoGeometry

func runRingDisplayChecks(_ report: inout Report) {
    report.section("macTEM pieces: CrystalLoading") { r in
        // SHA-256 of "abc" (FIPS 180-2 test vector)
        r.equal("fingerprint of abc", CrystalLoading.fingerprint(of: Data("abc".utf8)),
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
        var refused = false
        do { _ = try CrystalLoading.load(bytes: Data(count: CrystalLoading.maximumBytes + 1), fileName: "x.cif") }
        catch is CrystalLoading.Refusal { refused = true } catch {}
        r.equal("one byte over the cap is refused", refused, true)
    }
    report.section("macTEM pieces: RingOverlay and RingScale") { r in
        let circle = RingEllipse(x: 1, y: 2, a: 3, b: 3, theta: 0)
        r.check("a point on the circle is at distance 0", abs(RingOverlay.distance(fromX: 1, y: 5, to: circle)) < 1e-9)
        r.equal("the polyline has the default point count", RingOverlay.polyline(circle).count, RingOverlay.polylinePoints)
        r.equal("linear: R 100 px at dq 0.01", RingScale.linear(radius: 100, deltaQ: 0.01), 1.0)
        r.check("fit versus file is the relative difference", abs(RingScale.fitVersusFile(fitDeltaQ: 0.0105, fileDeltaQ: 0.01) - 0.05) < 1e-12)
    }
}
