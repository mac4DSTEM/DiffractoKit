//
//  StructureFactors.swift
//  Role: kinematic structure factors F(hkl) = (1/V) Σ_sites f(Z, |g|²) · occ · exp(2πi (h x + k y + l z)) in Å⁻² (f in Å
//        over V in Å³; py4DSTEM's and the reference implementations' "A" label names f's unit, not F's), and the
//        reflection list after the |F| cut with its rings.
//        The loop is reference commit 7d6fdbd1's `Crystal.reflections(kMax:tolerance:)`, in the same order of operations, so
//        the numbers can be bit identical to the reference implementations'; the enumeration and the ring rule are DiffractoCore's.
//

import Foundation
import simd
import DiffractoCore

/// A form factor: f(Z, g²) for g² = |g|² in Å⁻², or nil when the element has none.
public typealias FormFactor = (_ z: Int, _ gSquared: Double) -> Double?

/// F(hkl) of one reflection.
public struct StructureFactor: Sendable, Equatable {
    public let hkl: MillerIndex
    /// g = h b₁ + k b₂ + l b₃ (Å⁻¹, Cartesian).
    public let g: SIMD3<Double>
    /// |g| (Å⁻¹).
    public let gLength: Double
    /// Re F and Im F (Å⁻²: f_e in Å over V in Å³).
    public let real: Double
    public let imaginary: Double
    /// |F|² = Re² + Im² (Å⁻⁴).
    public let intensity: Double

    /// |F| (Å⁻²).
    public var magnitude: Double { intensity.squareRoot() }
    /// d = 1/|g| (Å).
    public var d: Double { 1 / gLength }
}

public enum StructureFactors {

    /// The reference implementations' |F| cut-off (Å⁻²): a reflection is kept when |F| > cut. A fixed number against a quantity that scales as
    /// 1/V: a large cell keeps fewer reflections at the same cut.
    public static let defaultCut = 1e-4

    /// F(hkl) for `structure`, with the electron form factor unless another is given. nil when a site's element has
    /// no form factor.
    public static func compute(_ hkl: MillerIndex, in structure: CrystalStructure,
                               formFactor: FormFactor = ElectronScatteringFactors.f) -> StructureFactor? {
        compute(hkl, g: structure.cell.gVector(hkl), sites: structure.sites, volume: structure.cell.volume,
                formFactor: formFactor)
    }

    /// The loop, in the reference implementations' order: |g| = length(g), g² = |g|·|g|; per site f, phase = 2π(hx + ky + lz),
    /// Re += f·occ·cos, Im += f·occ·sin; then ÷ V.
    static func compute(_ hkl: MillerIndex, g: SIMD3<Double>, sites: [AtomSite], volume: Double,
                        formFactor: FormFactor) -> StructureFactor? {
        let gLen = length(g)
        let gSq = gLen * gLen
        var re = 0.0, im = 0.0
        for site in sites {
            guard let fe = formFactor(site.z, gSq) else { return nil }
            let phase = 2 * Double.pi * (Double(hkl.h) * site.fractional.x
                                       + Double(hkl.k) * site.fractional.y
                                       + Double(hkl.l) * site.fractional.z)
            re += fe * site.occupancy * cos(phase)
            im += fe * site.occupancy * sin(phase)
        }
        re /= volume; im /= volume
        return StructureFactor(hkl: hkl, g: g, gLength: gLen, real: re, imaginary: im, intensity: re * re + im * im)
    }
}

public enum KinematicError: Error, Equatable, Sendable, CustomStringConvertible {
    /// Atomic numbers the form factor has no entry for; checked before any F is computed.
    case unsupportedElements([Int])
    /// DiffractoCore's refusal of the reflection list (k_max, or a cell too large or invalid).
    case reflectionList(ReflectionListError)
    /// The cut is not a finite number ≥ 0.
    case invalidCut(Double)
    /// The structure fails `CrystalStructure.validationIssues` (non-finite or out-of-range occupancy or coordinate, a
    /// degenerate cell, a close contact): refused before any F, never silently scored as nothing.
    case invalidStructure([String])
    /// The form factor returned a non-finite value, or F came out non-finite, at this reflection.
    case nonFiniteStructureFactor(MillerIndex)

    public var description: String {
        switch self {
        case .unsupportedElements(let z):
            return "The scattering table has no entry for Z = \(z.map(String.init).joined(separator: ", ")), so this "
                + "crystal's structure factors would be wrong."
        case .reflectionList(let e): return e.description
        case .invalidCut(let c): return "The |F| cut must be a finite number of Å⁻² at least 0; \(c) is not."
        case .invalidStructure(let issues):
            return "The structure is not valid, so its structure factors would be meaningless: \(issues.joined(separator: "; "))."
        case .nonFiniteStructureFactor(let hkl):
            return "The structure factor of \(hkl) is not a finite number (the form factor returned NaN or infinity)."
        }
    }
}

/// `localizedDescription` (what an app shows) is the same text as `description`.
extension KinematicError: LocalizedError {
    public var errorDescription: String? { description }
}

/// A reflection of DiffractoCore's list with its structure factor.
public struct KinematicReflection: Sendable, Equatable {
    public let reflection: Reflection
    public let structureFactor: StructureFactor
}

/// A ring after the cut: DiffractoCore's family of the kept reflections, with Σ|F|² over its members.
public struct KinematicFamily: Sendable, Equatable {
    public let family: ReflectionFamily
    /// Σ|F|² over the members (Å⁻⁴).
    public let intensity: Double
}

/// Every reflection with 0 < |g| ≤ k_max and its F; the ones with |F| > cut; their rings.
public struct KinematicReflectionList: Sendable, Equatable {
    public let kMax: Double
    public let cut: Double
    /// Every reflection of DiffractoCore's list (allowed or absent), in its order, with F. EMPTY when the list was built with
    /// `retainingAll: false` (`CrystalStructure.kinematicReflections`): `kept` and `families` are then the same values.
    public let all: [KinematicReflection]
    /// |F| > cut, in the same order.
    public let kept: [KinematicReflection]
    /// DiffractoCore's ring rule on `kept`.
    public let families: [KinematicFamily]

    /// Kept reflections that the detected space-group operations forbid: their |F| is non-zero only because the sites'
    /// written coordinates break the symmetry slightly (⅓ written as 0.333 gives |F| up to 3e-4 in magnesium). The reference implementations
    /// keep them by |F| alone and so does the list; this names them so a caller can tell.
    public var keptButOperationAbsent: [KinematicReflection] { kept.filter { !$0.reflection.isAllowed } }
}

extension CrystalStructure {

    /// The kinematic list: DiffractoCore's reflections, F by `formFactor`, the cut, the rings.
    ///
    /// **Memory (`retainingAll`).** By default the list keeps every reflection before the cut in `all` (and the Core list
    /// it came from is held while it is built): a 30 Å cell at k_max 2.5 Å⁻¹ peaks at some hundreds of megabytes, and the
    /// 1e7-reflection limit at several gigabytes. With `retainingAll: false` the structure factor of each reflection is
    /// computed as it is enumerated and only the reflections above the cut are kept: `all` is EMPTY, the pre-cut list is
    /// never held, and `kept` and `families` are bit for bit the values the default gives. Everything else is the same,
    /// the refusals included (for a form factor that fails at several reflections, the one named is the first in the
    /// list's order, as in the default).
    public func kinematicReflections(kMax: Double, cut: Double = StructureFactors.defaultCut,
                                     formFactor: FormFactor = ElectronScatteringFactors.f,
                                     policy: ReflectionList.Policy = .default,
                                     validation: CrystalStructure.ValidationPolicy = .default,
                                     retainingAll: Bool = true)
        throws(KinematicError) -> KinematicReflectionList {
        guard cut.isFinite, cut >= 0 else { throw .invalidCut(cut) }
        let missing = Array(Set(sites.map(\.z).filter { formFactor($0, 0) == nil })).sorted()
        guard missing.isEmpty else { throw .unsupportedElements(missing) }
        let issues = validationIssues(using: validation)
        guard issues.isEmpty else { throw .invalidStructure(issues.map(\.message)) }
        if !retainingAll {
            return try keptKinematicReflections(kMax: kMax, cut: cut, formFactor: formFactor, policy: policy)
        }
        let list: ReflectionList
        do { list = try reflections(kMax: kMax, policy: policy) } catch { throw .reflectionList(error) }
        var all: [KinematicReflection] = []
        all.reserveCapacity(list.reflections.count)
        for r in list.reflections {
            // Unreachable nil: every element was checked above.
            guard let f = StructureFactors.compute(r.hkl, g: r.g, sites: sites, volume: cell.volume,
                                                   formFactor: formFactor)
            else {
                // A form factor that has no value at this g² for some element: name the elements.
                let gSq = r.gLength * r.gLength
                throw .unsupportedElements(Array(Set(sites.map(\.z).filter { formFactor($0, gSq) == nil })).sorted())
            }
            guard f.real.isFinite, f.imaginary.isFinite else { throw .nonFiniteStructureFactor(r.hkl) }
            all.append(KinematicReflection(reflection: r, structureFactor: f))
        }
        let kept = all.filter { $0.structureFactor.magnitude > cut }
        let rings = ReflectionList.families(of: kept.map(\.reflection), laueGroup: list.laueGroup,
                                            millerBravais: list.usesMillerBravais,
                                            familyTolerance: policy.familyTolerance)
        var byHKL: [MillerIndex: Double] = [:]
        for k in kept { byHKL[k.reflection.hkl] = k.structureFactor.intensity }
        let families = rings.map { family in
            KinematicFamily(family: family, intensity: family.members.reduce(0) { $0 + (byHKL[$1] ?? 0) })
        }
        return KinematicReflectionList(kMax: kMax, cut: cut, all: all, kept: kept, families: families)
    }

    /// `kinematicReflections(retainingAll: false)`: the same reflections enumerated by `ReflectionList.enumerate`, F computed
    /// as each one is met, only those above the cut kept, then put in the list's order (`ReflectionList.ordered`).
    private func keptKinematicReflections(kMax: Double, cut: Double, formFactor: FormFactor,
                                          policy: ReflectionList.Policy) throws(KinematicError) -> KinematicReflectionList {
        do { try ReflectionList.checkBounds(cell: cell, kMax: kMax, policy: policy) } catch { throw .reflectionList(error) }
        let absenceTolerance = policy.effectiveAbsenceTolerance
        let volume = cell.volume
        let sites = self.sites
        var kept: [KinematicReflection] = []
        // The default path throws at the first failing reflection in the list's order; this one meets them in the
        // enumeration's order, so it keeps the earliest in the list's order and throws it at the end.
        var failure: (gLength: Double, hkl: MillerIndex, error: KinematicError)?
        ReflectionList.enumerate(cell: cell, kMax: kMax) { hkl, g, gLen in
            func record(_ error: KinematicError) {
                if let current = failure {
                    // the list's order: by |g|, and within a tie (the list's `orderTolerance`) by descending (h, k, l)
                    let tied = abs(gLen - current.gLength) <= ReflectionList.orderTolerance * max(gLen, current.gLength)
                    let earlier = tied ? hkl > current.hkl : gLen < current.gLength
                    if !earlier { return }
                }
                failure = (gLen, hkl, error)
            }
            guard let f = StructureFactors.compute(hkl, g: g, sites: sites, volume: volume, formFactor: formFactor) else {
                let gSq = gLen * gLen
                record(.unsupportedElements(Array(Set(sites.map(\.z).filter { formFactor($0, gSq) == nil })).sorted()))
                return
            }
            guard f.real.isFinite, f.imaginary.isFinite else { record(.nonFiniteStructureFactor(hkl)); return }
            guard f.magnitude > cut else { return }
            let reflection = Reflection(hkl: hkl, g: g, gLength: gLen,
                                        isAllowed: ReflectionList.isAllowed(hkl, by: operations, tolerance: absenceTolerance))
            kept.append(KinematicReflection(reflection: reflection, structureFactor: f))
        }
        if let failure { throw failure.error }
        let ordered = ReflectionList.ordered(kept, by: { $0.reflection })
        let rings = ReflectionList.families(
            of: ordered.map(\.reflection), laueGroup: laueGroup,
            millerBravais: ReflectionList.usesMillerBravais(laueClass: laueClass, cell: cell),
            familyTolerance: policy.familyTolerance)
        var byHKL: [MillerIndex: Double] = [:]
        for k in ordered { byHKL[k.reflection.hkl] = k.structureFactor.intensity }
        let families = rings.map { family in
            KinematicFamily(family: family, intensity: family.members.reduce(0) { $0 + (byHKL[$1] ?? 0) })
        }
        return KinematicReflectionList(kMax: kMax, cut: cut, all: [], kept: ordered, families: families)
    }
}
