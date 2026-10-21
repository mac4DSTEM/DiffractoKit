//
//  SymmetryOperationTests.swift
//  The x,y,z grammar (registration 3.3): what parses, what it means, what is refused, and the round trip.
//

import DiffractoKitTestSupport
import Testing
import simd
import DiffractoCore

struct SymmetryOperationTests {

    @Test func testIdentityAndTranslatedIdentity() throws {
        let identity = try SymmetryOperation(xyz: "x,y,z")
        #expect(identity.rotation == .identity)
        #expect(identity.translation == SIMD3(0, 0, 0))
        #expect(identity.isIdentity)
        #expect(try SymmetryOperation(xyz: "x+1,y,z").isIdentity)
        #expect(!(try SymmetryOperation(xyz: "x+1/2,y,z").isIdentity))
        #expect(!(try SymmetryOperation(xyz: "-x,-y,-z").isIdentity))
    }

    /// `-x,y+1/2,-z`: the 2₁ screw along b. Mutant: "1/2+x" read as "x".
    @Test func testScrewAxisParsesRotationAndTranslation() throws {
        let op = try SymmetryOperation(xyz: "-x,y+1/2,-z")
        #expect(op.rotation == IntegerMatrix([-1, 0, 0], [0, 1, 0], [0, 0, -1]))
        #expect(op.translation == SIMD3(0, 0.5, 0))
        let centring = try SymmetryOperation(xyz: "1/2+x,1/2-y,z")
        #expect(centring.rotation == IntegerMatrix([1, 0, 0], [0, -1, 0], [0, 0, 1]))
        #expect(centring.translation == SIMD3(0.5, 0.5, 0))
        let hex = try SymmetryOperation(xyz: "x-y,x,-z+1/3")
        #expect(hex.rotation == IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, -1]))
        #expect(abs(hex.translation.z - (1.0 / 3)) <= 1e-15)
        let spaced = try SymmetryOperation(xyz: " -X , Y + 1/2 , -Z ")
        #expect(spaced.rotation == op.rotation)
        #expect(spaced.translation == op.translation)
        let decimal = try SymmetryOperation(xyz: "x+0.5,y,z")
        #expect(decimal.translation.x == 0.5)
    }

    @Test func testApplyAndCompose() throws {
        let op = try SymmetryOperation(xyz: "-x,y+1/2,-z")
        let p = SIMD3(0.1, 0.2, 0.3)
        #expect(op.apply(p) == SIMD3(-0.1, 0.7, -0.3))
        // op ∘ op = (x, y + 1, z): a lattice translation, the identity mod 1.
        let twice = op.composed(with: op)
        #expect(twice.isIdentity)
        #expect(twice.translationMod1 == SIMD3(0, 0, 0))
        // Composition order: self(other(p)).
        let a = try SymmetryOperation(xyz: "-y,x,z"), b = try SymmetryOperation(xyz: "x,y,z+1/2")
        #expect(a.composed(with: b).apply(p) == a.apply(b.apply(p)))
    }

    @Test func testRefusals() {
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "x,y") }
            #expect(error == .unparsable("x,y"))
        }
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "a,b,c") }
            #expect(error == .unparsable("a,b,c"))
        }
        #expect(throws: (any Error).self) { try SymmetryOperation(xyz: "x,,z") }
        #expect(throws: (any Error).self) { try SymmetryOperation(xyz: "x+,y,z") }
        // "x/2" is "x" then the fraction "/2" with an empty numerator: unparsable (the reference implementations' parser read it as x + 2).
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "x/2,y,z") }
            #expect(error == .unparsable("x/2,y,z"))
        }
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "0.5x,y,z") }
            #expect(error == .nonIntegerRotation("0.5x,y,z"))
        }
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "1/3x,y,z") }
            #expect(error == .nonIntegerRotation("1/3x,y,z"))
        }
        // A garbled coefficient or translation is refused, not trapped (Int(1e21) would trap).
        do {
            let error = #expect(throws: SymmetryOperationError.self) { try SymmetryOperation(xyz: "1000000000000000000000x,y,z") }
            #expect(error == .unparsable("1000000000000000000000x,y,z"))
        }
        #expect(throws: (any Error).self) { try SymmetryOperation(xyz: "x+1000000000000000000000,y,z") }
        #expect(throws: (any Error).self) { try SymmetryOperation(xyz: "x+inf,y,z") }
        #expect(throws: (any Error).self) { try SymmetryOperation(xyz: "x/0,y,z") }
        // reference commit 7d6fdbd1's rule: empty fraction pieces, and a term that does not end at a sign or the end.
        for raw in ["x+2/,y,z", "x+1//2,y,z", "xy,y,z", "x2,y,z", "x.5,y,z", "-y1/2,x,z", "1/2x,y,z"] {
            #expect(throws: (any Error).self, "\(raw)") { try SymmetryOperation(xyz: raw) }
        }
        #expect(throws: Never.self) { try SymmetryOperation(xyz: "-y+1/2,x-y,z+1/3") }
        #expect(throws: Never.self) { try SymmetryOperation(xyz: "1/2+x,1/2-y,-z") }
    }

    /// `xyzString` parses back to an equivalent operation, for every operation of C2/m, P6₃/mmc's screw and a 1/3.
    @Test func testRoundTrip() throws {
        for text in ["x,y,z", "-x,y,-z", "1/2+x,1/2+y,z", "1/2-x,1/2-y,-z", "x-y,x,z+1/2", "-y,x-y,z+1/3",
                     "y,x,-z+3/4", "-x+y,y,-z"] {
            let op = try SymmetryOperation(xyz: text)
            let back = try SymmetryOperation(xyz: op.xyzString)
            #expect(back.isEquivalent(to: op), "\(text) -> \(op.xyzString)")
        }
        #expect(try SymmetryOperation(xyz: "1/2+x,1/2+y,z").xyzString == "x+1/2,y+1/2,z")
        #expect(try SymmetryOperation(xyz: "-x,y,-z").xyzString == "-x,y,-z")
    }

    /// Distinct operations: rotation exact, translation mod 1 to 1e-3; a garbled coefficient counts, not traps.
    @Test func testDistinctCount() throws {
        let ops = try ["x,y,z", "x+1,y,z", "-x,-y,-z", "-x,-y,-z", "x+1/2,y,z", "x+0.5001,y,z"].map { try SymmetryOperation(xyz: $0) }
        #expect(SymmetryOperation.distinctCount(ops) == 3)
        // A whole-number translation, however large, is the identity mod 1 (1e21 is an integer-valued double) and
        // must not trap in the quantisation; a non-finite one counts as its own operation.
        let huge = SymmetryOperation(rotation: .identity, translation: SIMD3(1e21, 0, 0))
        #expect(SymmetryOperation.distinctCount([huge, .identity]) == 1)
        let infinite = SymmetryOperation(rotation: .identity, translation: SIMD3(.infinity, 0, 0))
        #expect(SymmetryOperation.distinctCount([infinite, .identity]) == 2)
    }

    @Test func testOrdersAndDeterminants() {
        let c6 = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        #expect(c6.order == 6)
        #expect(c6.determinant == 1)
        #expect(c6.multiplied(by: c6).order == 3)
        #expect(IntegerMatrix.inversion.order == 2)
        #expect(IntegerMatrix.inversion.determinant == -1)
        #expect(IntegerMatrix.identity.order == 1)
        #expect(IntegerMatrix([1, 1, 0], [0, 1, 0], [0, 0, 1]).order == nil, "a shear has no finite order")
        #expect(c6.transposed.multiplied(by: c6).trace == 3 + 1, "the transpose is not the inverse in the hexagonal basis")
    }
}
