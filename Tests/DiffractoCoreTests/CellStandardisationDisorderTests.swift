//
//  CellStandardisationDisorderTests.swift
//  Role: `CellStandardisation` with a disordered site (several species at one position, as `MaterialsProjectConversion` expands
//        them): the occupancies and labels survive the conversion and the per-element count rule holds, so the app's adapter
//        needs no expansion of its own. Exact (a fixed matrix on a fixed cell).
//

import DiffractoKitTestSupport
import Foundation
import Testing
import simd
import DiffractoCore

struct CellStandardisationDisorderTests {
    @Test func aDisorderedFaceCentredPrimitiveCellKeepsOccupanciesAndLabels() throws {
        let a = 4.0
        // fcc primitive vectors of a conventional cell of edge a.
        let rows = [SIMD3(0.0, a / 2, a / 2), SIMD3(a / 2, 0, a / 2), SIMD3(a / 2, a / 2, 0)]
        // One site, 60 % Fe (Z 26) and 40 % Co (Z 27) at the origin.
        let sites = [AtomSite(z: 26, fractional: .zero, occupancy: 0.6, label: "Fe0"),
                     AtomSite(z: 27, fractional: .zero, occupancy: 0.4, label: "Co0")]
        let result = try CellStandardisation.standardise(latticeRows: rows, sites: sites, symbol: "Fm-3m", crystalSystem: "cubic")
        #expect(result.outcome == .conventional(centring: "F"))
        #expect(result.sites.count == 8)
        #expect(result.sites.filter { $0.z == 26 }.allSatisfy { $0.occupancy == 0.6 && $0.label == "Fe0" })
        #expect(result.sites.filter { $0.z == 27 }.allSatisfy { $0.occupancy == 0.4 && $0.label == "Co0" })
        #expect(result.sites.filter { $0.z == 26 }.count == 4)
        // Total occupancy per position stays 1: the species of one position share it.
        let positions = Dictionary(grouping: result.sites, by: { "\($0.fractional.x),\($0.fractional.y),\($0.fractional.z)" })
        #expect(positions.count == 4)
        #expect(positions.values.allSatisfy { abs($0.reduce(0) { $0 + $1.occupancy } - 1) < 1e-12 })
    }
}
