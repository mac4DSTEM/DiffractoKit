//
//  PolicyTests.swift
//  DiffractoKit 0.2.0 (docs/registration-0.2.0.md §1, §2): the policies of DiffractoScattering. The |F| cut is a parameter
//  since slice 2; the list's refusal bounds and ring rule and the structure validation reach the kinematic list through
//  `policy:` and `validation:`. Tier 2: a non-default value moves the outcome in the stated direction, an invalid one is
//  refused or taken as documented. Every test names under "RED UNDER" the production mutation that turns a default-path test
//  red (tools/kit-0.2.0/mutants-policy.py).
//

import DiffractoKitTestSupport
import Testing
import Foundation
import simd
import DiffractoCore
import DiffractoScattering

struct ScatteringPolicyTests {

    static func aluminium(_ cell: UnitCell = .cubic(a: 4.0495)) -> CrystalStructure {
        let fcc: [SIMD3<Double>] = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        let sites: [AtomSite] = fcc.map { AtomSite(z: 13, fractional: $0) }
        return CrystalStructure.detectingSymmetry(cell: cell, sites: sites)
    }

    /// The |F| cut is a policy with the default 1e-4 Å⁻²: a lower cut keeps at least as many reflections, a higher one at most
    /// as many, and a cut that is negative or NaN is `invalidCut`.
    /// RED UNDER: the default cut changed; the cut ignored.
    @Test func theCutMovesWhatIsKept() throws {
        #expect(StructureFactors.defaultCut == 1e-4)
        let structure = Self.aluminium()
        let ordinary = try structure.kinematicReflections(kMax: 1.5)
        let lower = try structure.kinematicReflections(kMax: 1.5, cut: 0)
        let magnitudes: [Double] = ordinary.kept.map { $0.structureFactor.magnitude }.sorted()
        let median: Double = magnitudes[magnitudes.count / 2]
        let higher = try structure.kinematicReflections(kMax: 1.5, cut: median)
        #expect(lower.kept.count >= ordinary.kept.count)
        #expect(higher.kept.count < ordinary.kept.count)
        #expect(lower.kept.count > higher.kept.count)
        for bad in [-1.0, Double.nan] {
            var refused = false
            do { _ = try structure.kinematicReflections(kMax: 1.5, cut: bad) } catch { guard case .invalidCut = error else { throw error }; refused = true }
            #expect(refused)
        }
    }

    /// The list's policy reaches the kinematic list: the k_max cap, and the ring rule (a tetragonal cell with c/a = 1 + 1e-7
    /// has more rings at the default tolerance than at 1e-6).
    /// RED UNDER: `policy:` not passed to the reflection list; the family tolerance not passed to the ring rule.
    @Test func theListPolicyReachesTheKinematicList() throws {
        let structure = Self.aluminium(.cubic(a: 1.2))
        var refused = false
        do { _ = try structure.kinematicReflections(kMax: 4.5) }
        catch { guard case .reflectionList(let inner) = error else { throw error }; refused = inner == .invalidKMax(4.5) }
        #expect(refused)
        let raised = try structure.kinematicReflections(kMax: 4.5, policy: ReflectionList.Policy(maximumKMax: 5))
        #expect(raised.kept.count > 0)

        let split = Self.aluminium(.tetragonal(a: 3, c: 3 * (1 + 1e-7)))
        let separate = try split.kinematicReflections(kMax: 1.5)
        let joined = try split.kinematicReflections(kMax: 1.5, policy: ReflectionList.Policy(familyTolerance: 1e-6))
        #expect(joined.families.count < separate.families.count)
    }

    /// The validation policy reaches the kinematic list: two sites 0.3 Å apart are `invalidStructure` at the default and a
    /// list at a close-contact limit of 0.2.
    /// RED UNDER: `validation:` not passed to the structure check.
    @Test func theValidationPolicyReachesTheKinematicList() throws {
        let close = CrystalStructure(cell: .cubic(a: 4),
                                     sites: [AtomSite(z: 13, fractional: [0, 0, 0]), AtomSite(z: 13, fractional: [0.075, 0, 0])],
                                     operations: [.identity])
        var refused = false
        do { _ = try close.kinematicReflections(kMax: 1) } catch { guard case .invalidStructure = error else { throw error }; refused = true }
        #expect(refused)
        let relaxed = try close.kinematicReflections(kMax: 1, validation: CrystalStructure.ValidationPolicy(closeContactLimit: 0.2))
        #expect(relaxed.all.count > 0)
    }

    /// A list given the default policies equals a list given none.
    @Test func aListGivenTheDefaultsEqualsAListGivenNone() throws {
        let structure = Self.aluminium()
        let plain = try structure.kinematicReflections(kMax: 1.5)
        let explicit = try structure.kinematicReflections(kMax: 1.5, policy: .default, validation: .default)
        #expect(plain == explicit)
    }
}
