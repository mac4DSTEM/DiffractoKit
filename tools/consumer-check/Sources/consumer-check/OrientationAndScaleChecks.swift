//
//  OrientationAndScaleChecks.swift
//  Role: the public surface 0.10.0 adds, used from outside the package (only public API is named): `Misorientation`, `Orientation`,
//        `OrientationVariantSet` with `OrientationRelationship.rotation` / `.variants`, `LinearScale.assign` / `.fit`, the pattern measurements and the
//        public `ReflectionFamily` initialiser. Planted values only; the pins are the package's own tests.
//

import Foundation
import simd
import DiffractoCore
import DiffractoGeometry
import DiffractoOrientation

func runOrientationAndScaleChecks(_ r: inout Report) throws {
    // beta'' in Al: 12 variants
    let rotation = try OrientationRelationship.rotation(
        OrientationRelationship(candidate: .direction(SIMD3(1, 0, 0)), matrix: .direction(SIMD3(3, 2, 0))),
        OrientationRelationship(candidate: .direction(SIMD3(0, 1, 0)), matrix: .direction(SIMD3(0, 0, 1))),
        precipitate: CrystalPresets.betaDoublePrime, matrix: CrystalPresets.aluminum)
    let cubic = CrystalOrientationSymmetry(laueClass: .mThreeBarM), monoclinic = CrystalOrientationSymmetry(laueClass: .twoOverM)
    let set = try OrientationVariantSet(matrix: rotation, matrixSymmetry: cubic, precipitateSymmetry: monoclinic)
    r.equal("beta'' in Al has 12 variants", set.variants.count, 12)
    let hit = set.nearestVariant(of: set.variants[5])
    r.check("a variant is its own nearest variant", hit?.index == 5 && (hit?.disorientationRad ?? 1) < 1e-12)
    let reduced = try Misorientation(rotation, source: monoclinic, target: cubic).reduce()
    r.near("the relation's reduced disorientation is arctan(2/3)", reduced.angleRad, atan2(2.0, 3.0), 1e-12)
    r.equal("the equivalents are |target| x |source| matrices", try Misorientation(rotation, source: monoclinic, target: cubic).equivalent().count, 48)
    let o = try Orientation(matrix_identity_double3x3, symmetry: cubic)
    r.check("the Euler region of the identity is the origin", o.inEulerFundamentalRegion() == SIMD3(0, 0, 0))
    r.check("two groups cannot be compared", (try? o.angle(with: Orientation(matrix_identity_double3x3, symmetry: monoclinic))) == nil)

    // ring-sequence assignment on Al [001]
    let al = [3.0, 4.0, 8.0, 11.0].map { $0.squareRoot() / 4.0495 }
    let radii = [al[1], al[2]].map { $0 / 0.01 }
    let positions = [[Double]](repeating: radii + radii, count: 9)
    do {
        let a = try LinearScale.assign(radiiPerPosition: positions, shellLengths: al, labels: ["{111}", "{200}", "{220}", "{311}"])
        r.check("the first ring of Al [001] is {200} or the call refuses", a?.shellIndex == 1 || a == nil)
    } catch let error as LinearScaleAssignmentError {
        r.check("a refusal is a LocalizedError with a sentence", !error.localizedDescription.isEmpty)
    }
    let q = LinearScale.fit(ringRadiiPixels: al.map { $0 / 0.0125 }, indexedInvAngstrom: al)
    r.near("delta q from four rings", q?.invAngstromPerPixel ?? .nan, 0.0125, 1e-12)

    // pattern measurements
    let n = 128
    var px = [Float](repeating: 10, count: n * n)
    for y in 0 ..< n { for x in 0 ..< n {
        let dy = Double(y) - 63.4, dx = Double(x) - 64.6
        px[y * n + x] += Float(5000 * exp(-(dy * dy + dx * dx) / 18))
    } }
    let pattern = DiffractionPattern(width: n, height: n, pixels: px)
    let centre = PatternCentre.find(pattern: pattern)
    r.near("the beam centre (row)", centre?.row ?? .nan, 63.4, 0.05)
    let spot = PatternSpot.refine(pattern: pattern, row: 64, column: 65, half: 10)
    r.near("the spot centre (column)", spot?.column ?? .nan, 64.6, 0.05)
    r.check("the radial profile has bins", (PatternRadialProfile.about(pattern: pattern, centerRow: 63.4, centerColumn: 64.6)?.mean.count ?? 0) > 50)
    let family = ReflectionFamily(g: 0.5, members: [MillerIndex(1, 1, 1)], orbits: [MillerIndex(1, 1, 1)], name: "{1 1 1}")
    r.check("a ReflectionFamily can be built", family.d == 2)
}
