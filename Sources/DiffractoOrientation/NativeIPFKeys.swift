// Native IPF keys and direction samplers, ported from GPL-3.0 source c809679c.
// Oracle: tools/pins/orientation-symmetry/pins.py compiles the original git-show source.
import Foundation
import simd
/// Cubic crystal-symmetry operations and the m-3m inverse-pole-figure color
/// key. Its arithmetic is pinned to the original Swift implementation.
internal enum NativeCubicSymmetry {
    /// Fold a direction into the cubic m-3m IPF sector z ≥ x ≥ y ≥ 0.
    static func reduceDirection(_ direction: SIMD3<Double>) -> SIMD3<Double> {
        let length = simd_length(direction)
        guard length.isFinite, length > 1e-15 else { return [0, 0, 1] }
        let values = [abs(direction.x), abs(direction.y), abs(direction.z)].sorted()
        return simd_normalize(SIMD3(values[1], values[0], values[2]))
    }

    /// Deterministic, approximately uniform sampling of the cubic m-3m
    /// direction fundamental zone [001]-[101]-[111]. A dense Fibonacci sphere
    /// is folded by symmetry, then spherical farthest-point selection avoids
    /// spending templates on symmetry-equivalent beam directions.
    static func sampleFundamentalZone(count: Int) -> [SIMD3<Double>] {
        guard count > 0 else { return [] }
        let vertices = [
            SIMD3<Double>(0, 0, 1),
            simd_normalize(SIMD3<Double>(1, 0, 1)),
            simd_normalize(SIMD3<Double>(1, 1, 1)),
        ]
        if count <= vertices.count { return Array(vertices.prefix(count)) }

        let candidateCount = max(128, count * 12)
        let golden = Double.pi * (3 - 5.0.squareRoot())
        var candidates: [SIMD3<Double>] = []
        candidates.reserveCapacity(candidateCount)
        for index in 0..<candidateCount {
            let z = 1 - 2 * (Double(index) + 0.5) / Double(candidateCount)
            let radius = max(0, 1 - z * z).squareRoot()
            let azimuth = golden * Double(index)
            candidates.append(reduceDirection(SIMD3(
                radius * cos(azimuth), radius * sin(azimuth), z
            )))
        }

        var selected = vertices
        var minimumChord = candidates.map { candidate in
            vertices.map { 1 - simd_dot(candidate, $0) }.min() ?? 0
        }
        while selected.count < count {
            guard let best = minimumChord.indices.max(by: {
                minimumChord[$0] < minimumChord[$1]
            }) else { break }
            let next = candidates[best]
            selected.append(next)
            for index in candidates.indices {
                minimumChord[index] = min(
                    minimumChord[index], 1 - simd_dot(candidates[index], next)
                )
            }
            minimumChord[best] = -.infinity
        }
        return selected
    }

    /// EDAX/TSL-style cubic IPF color, numerically aligned with orix's m-3m
    /// direction color key. The original numerical residuals at the vertices are preserved.
    static func ipfColor(direction: SIMD3<Double>) -> SIMD3<Float> {
        let v = reduceDirection(direction)
        let azimuth = correctedAzimuth(rawAzimuth(v))

        let vCenterNormal = unit(simd_cross(v, sectorCenter))
        var radial = Double.infinity
        for normal in sectorNormals {
            let boundary = unit(simd_cross(vCenterNormal, normal))
            let numerator = vectorAngle(-v, boundary)
            let denominator = vectorAngle(-sectorCenter, boundary)
            let ratio = numerator / denominator
            radial = min(radial, ratio.isFinite ? ratio : 1)
        }
        if !radial.isFinite { radial = 1 }
        return hslToRGB(hue: azimuth / (2 * .pi), lightness: 0.5 + 0.5 * radial)
    }

    private static let sectorCenter = simd_normalize(SIMD3(
        0.4281523501333333, 0.19245008973333333, 0.7614856834666667
    ))
    private static let sectorVertices = [
        simd_normalize(SIMD3(1.0, 1.0, 1.0)),
        simd_normalize(SIMD3(1.0, 0.0, 1.0)),
        SIMD3(0.0, 0.0, 1.0),
    ]
    private static let sectorNormals = [
        SIMD3(1.0, -1.0, 0.0),
        SIMD3(-1.0, 0.0, 1.0),
        SIMD3(0.0, 1.0, 0.0),
    ]

    /// A one-time equalized angular lookup around the asymmetric spherical
    /// triangle. This is the polar-coordinate construction used by MTEX/TSL.
    private static let azimuthLookup: [Double] = {
        let count = 1000
        let rawAngles = (0..<count).map { 2 * Double.pi * Double($0) / Double(count - 1) }
        let rx = SIMD3(0.0, 0.0, 1.0) - sectorCenter
        let tangent = unit(simd_cross(unit(rx - simd_dot(rx, sectorCenter) * sectorCenter),
                                      sectorCenter))
        var distances = [Double]()
        distances.reserveCapacity(count - 1)
        for angle in rawAngles.dropLast() {
            let rotated = rodrigues(tangent, axis: sectorCenter, angle: angle)
            let boundaryDistances = sectorNormals.map {
                vectorAngle(simd_cross($0, rotated), sectorCenter)
            }
            distances.append(boundaryDistances.min() ?? 0)
        }

        var cuts = sectorVertices.map { rawAzimuth($0) }.sorted().map {
            Int((Double(count) * $0 / (2 * .pi)).rounded())
        }
        if cuts.first.map({ $0 < 10 }) == true { cuts.removeFirst() }
        if cuts.count == 2 {
            for range in [0..<cuts[0], cuts[0]..<cuts[1], cuts[1]..<distances.count] {
                let sum = range.reduce(0.0) { $0 + distances[$1] }
                if sum > 0 {
                    let scale = 3 / sum
                    for index in range { distances[index] *= scale }
                }
            }
        }
        let total = distances.reduce(0, +)
        var cumulative = [0.0]
        cumulative.reserveCapacity(count)
        var sum = 0.0
        for distance in distances {
            sum += distance
            cumulative.append(2 * .pi * sum / total)
        }
        return cumulative
    }()

    private static func rawAzimuth(_ v: SIMD3<Double>) -> Double {
        let initial = SIMD3(0.0, 0.0, 1.0) - sectorCenter
        let x = unit(initial - simd_dot(initial, sectorCenter) * sectorCenter)
        let y = unit(simd_cross(sectorCenter, x))
        let distance = unit(v - sectorCenter)
        guard simd_length_squared(distance) > 0 else { return 0 }
        let angle = atan2(simd_dot(y, distance), simd_dot(x, distance))
        return angle < 0 ? angle + 2 * .pi : angle
    }

    private static func correctedAzimuth(_ angle: Double) -> Double {
        let position = max(0, min(Double(azimuthLookup.count - 1),
                                  angle / (2 * .pi) * Double(azimuthLookup.count - 1)))
        let lower = Int(position.rounded(.down))
        let upper = min(lower + 1, azimuthLookup.count - 1)
        let fraction = position - Double(lower)
        return azimuthLookup[lower] * (1 - fraction) + azimuthLookup[upper] * fraction
    }

    private static func unit(_ value: SIMD3<Double>) -> SIMD3<Double> {
        let length = simd_length(value)
        return length > 1e-15 && length.isFinite ? value / length : .zero
    }

    private static func vectorAngle(_ lhs: SIMD3<Double>, _ rhs: SIMD3<Double>) -> Double {
        let a = unit(lhs), b = unit(rhs)
        guard simd_length_squared(a) > 0, simd_length_squared(b) > 0 else { return .nan }
        return acos(min(1, max(-1, simd_dot(a, b))))
    }

    private static func rodrigues(_ vector: SIMD3<Double>, axis: SIMD3<Double>,
                                  angle: Double) -> SIMD3<Double> {
        let c = cos(angle), s = sin(angle)
        return vector * c + simd_cross(axis, vector) * s
            + axis * simd_dot(axis, vector) * (1 - c)
    }

    private static func hslToRGB(hue: Double, lightness: Double) -> SIMD3<Float> {
        let h = hue - floor(hue)
        let l2 = 2 * lightness
        let s2 = l2 <= 1 ? l2 : 2 - l2
        let saturation = l2 + s2 > 0 ? 2 * s2 / (l2 + s2) : 0
        let value = (l2 + s2) / 2
        let sector = h * 6
        let index = Int(floor(sector)) % 6
        let f = sector - floor(sector)
        let p = value * (1 - saturation)
        let q = value * (1 - saturation * f)
        let t = value * (1 - saturation * (1 - f))
        let rgb: SIMD3<Double>
        switch index {
        case 0: rgb = [value, t, p]
        case 1: rgb = [q, value, p]
        case 2: rgb = [p, value, t]
        case 3: rgb = [p, q, value]
        case 4: rgb = [t, p, value]
        default: rgb = [value, p, q]
        }
        return SIMD3(Float(rgb.x), Float(rgb.y), Float(rgb.z))
    }
}

/// Proper rotational subgroup of the hexagonal 6/mmm Laue class (D6): six
/// rotations about c plus six twofold rotations about basal-plane axes. The
/// diffraction/Friedel direction key folds the resulting sector to
/// `[0001]-[2-1-10]-[10-10]`.
internal enum NativeHexagonalSymmetry {
    static func reduceDirection(_ direction: SIMD3<Double>) -> SIMD3<Double> {
        let length = simd_length(direction)
        guard length.isFinite, length > 1e-15 else { return [0, 0, 1] }
        var value = direction / length
        if value.z < 0 { value = -value }
        let basal = hypot(value.x, value.y)
        var azimuth = atan2(value.y, value.x).truncatingRemainder(dividingBy: .pi / 3)
        if azimuth < 0 { azimuth += .pi / 3 }
        if azimuth > .pi / 6 { azimuth = .pi / 3 - azimuth }
        return simd_normalize(SIMD3(
            basal * cos(azimuth), basal * sin(azimuth), value.z
        ))
    }

    static func sampleFundamentalZone(count: Int) -> [SIMD3<Double>] {
        guard count > 0 else { return [] }
        let vertices = [
            SIMD3<Double>(0, 0, 1),
            SIMD3<Double>(1, 0, 0),
            SIMD3<Double>(cos(.pi / 6), sin(.pi / 6), 0),
        ]
        if count <= vertices.count { return Array(vertices.prefix(count)) }
        let candidateCount = max(192, count * 16)
        let golden = Double.pi * (3 - 5.0.squareRoot())
        var candidates: [SIMD3<Double>] = []
        for index in 0..<candidateCount {
            let z = 1 - 2 * (Double(index) + 0.5) / Double(candidateCount)
            let radius = max(0, 1 - z * z).squareRoot()
            candidates.append(reduceDirection(SIMD3(
                radius * cos(golden * Double(index)),
                radius * sin(golden * Double(index)), z
            )))
        }
        var selected = vertices
        var minimumChord = candidates.map { candidate in
            vertices.map { 1 - simd_dot(candidate, $0) }.min() ?? 0
        }
        while selected.count < count {
            guard let best = minimumChord.indices.max(by: {
                minimumChord[$0] < minimumChord[$1]
            }) else { break }
            let next = candidates[best]
            selected.append(next)
            for index in candidates.indices {
                minimumChord[index] = min(
                    minimumChord[index], 1 - simd_dot(candidates[index], next)
                )
            }
            minimumChord[best] = -.infinity
        }
        return selected
    }

    /// One corner of the 6/mmm key triangle: its crystal-Cartesian direction
    /// (py4DSTEM frame: a along +x, c along +z), the four-index label the
    /// legend prints, and the colour `ipfColor` gives it.
    /// Native 6/mmm key: c red, basal azimuth 0 green, basal azimuth 30° blue.
    /// Square-root intensity and maximum normalization are part of this convention.
    static func ipfColor(direction: SIMD3<Double>) -> SIMD3<Float> {
        let value = reduceDirection(direction)
        let tilt = min(1, max(0, acos(min(1, max(0, value.z))) / (.pi / 2)))
        let azimuth = min(.pi / 6, max(0, atan2(value.y, value.x)))
        let fraction = azimuth / (.pi / 6)
        var rgb = SIMD3<Double>(
            sqrt(max(0, 1 - tilt)),
            sqrt(max(0, tilt * (1 - fraction))),
            sqrt(max(0, tilt * fraction))
        )
        let maximum = max(rgb.x, rgb.y, rgb.z)
        if maximum > 0 { rgb /= maximum }
        return SIMD3(Float(rgb.x), Float(rgb.y), Float(rgb.z))
    }
}

/// Native b-unique monoclinic direction key.
internal enum NativeMonoclinicSymmetry {
    static func reduceDirection(_ direction: SIMD3<Double>) -> SIMD3<Double> {
        let length = simd_length(direction)
        guard length.isFinite, length > 1e-15 else { return [0, 1, 0] }
        var value = direction / length
        if value.y < 0 { value = -value }
        let basal = hypot(value.x, value.z)
        var azimuth = atan2(value.z, value.x).truncatingRemainder(dividingBy: .pi)
        if azimuth < 0 { azimuth += .pi }
        // A rounding hair below pi is the same direction as azimuth 0.
        if azimuth > .pi - 1e-15 { azimuth = 0 }
        return simd_normalize(SIMD3(basal * cos(azimuth), value.y, basal * sin(azimuth)))
    }

    /// Farthest-point sampling of the quarter sphere, seeded with its three corners
    /// (b, a-direction, z-direction). Distance to a point is the smaller of the distances to it
    /// and to its half-turn image, because the azimuth wraps at 180 degrees.
    static func sampleFundamentalZone(count: Int) -> [SIMD3<Double>] {
        guard count > 0 else { return [] }
        let vertices = [
            SIMD3<Double>(0, 1, 0),
            SIMD3<Double>(1, 0, 0),
            SIMD3<Double>(0, 0, 1),
        ]
        if count <= vertices.count { return Array(vertices.prefix(count)) }
        func chord(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double {
            1 - max(simd_dot(a, b), simd_dot(a, SIMD3(-b.x, b.y, -b.z)))
        }
        let candidateCount = max(192, count * 16)
        let golden = Double.pi * (3 - 5.0.squareRoot())
        var candidates: [SIMD3<Double>] = []
        for index in 0..<candidateCount {
            let y = 1 - 2 * (Double(index) + 0.5) / Double(candidateCount)
            let radius = max(0, 1 - y * y).squareRoot()
            candidates.append(reduceDirection(SIMD3(
                radius * cos(golden * Double(index)), y,
                radius * sin(golden * Double(index))
            )))
        }
        var selected = vertices
        var minimumChord = candidates.map { candidate in
            vertices.map { chord(candidate, $0) }.min() ?? 0
        }
        while selected.count < count {
            guard let best = minimumChord.indices.max(by: {
                minimumChord[$0] < minimumChord[$1]
            }) else { break }
            let next = candidates[best]
            selected.append(next)
            for index in candidates.indices {
                minimumChord[index] = min(minimumChord[index], chord(candidates[index], next))
            }
            minimumChord[best] = -.infinity
        }
        return selected
    }

    /// 2/m IPF key for the beam direction. The zone is the quarter sphere: polar angle from b in
    /// [0, 90 degrees], azimuth phi about b in [0, 180 degrees), with phi and phi + 180 degrees
    /// the same direction. A key must be injective on that zone and seamless across phi = 0 = 180,
    /// so it is the HSV disk: hue = 120 degrees + 2*phi (one full turn of the wheel over the
    /// period, so phi = 0 and 180 meet; the a-direction is green as before, phi = 60 degrees is
    /// blue), saturation = polar angle / 90 degrees (the pole b is the single achromatic colour,
    /// white), value fixed at 1. The earlier (1 - cos 2*phi)/2 weight was symmetric under
    /// phi -> 180 - phi, so (x, y, z) and (-x, y, z), which are inequivalent under 2/m, shared a
    /// colour (refuter RB). Unlike the cubic and hexagonal keys (three primaries, sqrt
    /// intensity) this one is a wheel, because the zone is a disk with no corner structure to
    /// carry three primaries; the cost is the usual one of hue keys: red-green and blue-yellow
    /// pairs a colour-blind reader cannot separate are told apart by saturation and by the
    /// legend's position on the disk, not by hue alone. No fraction or threshold is tuned.
    static func ipfColor(direction: SIMD3<Double>) -> SIMD3<Float> {
        let value = reduceDirection(direction)
        let tilt = min(1, max(0, acos(min(1, max(0, value.y))) / (.pi / 2)))
        var azimuth = atan2(value.z, value.x)
        if azimuth < 0 { azimuth += .pi }
        let hue = (120 + 2 * azimuth * 180 / .pi).truncatingRemainder(dividingBy: 360) / 60
        let sector = Int(hue) % 6
        let f = hue - Double(Int(hue))
        // HSV -> RGB with value 1 and saturation tilt.
        let p = 1 - tilt, q = 1 - tilt * f, t = 1 - tilt * (1 - f)
        let rgb: SIMD3<Double>
        switch sector {
        case 0: rgb = SIMD3(1, t, p)
        case 1: rgb = SIMD3(q, 1, p)
        case 2: rgb = SIMD3(p, 1, t)
        case 3: rgb = SIMD3(p, q, 1)
        case 4: rgb = SIMD3(t, p, 1)
        default: rgb = SIMD3(1, p, q)
        }
        return SIMD3(Float(rgb.x), Float(rgb.y), Float(rgb.z))
    }
}
