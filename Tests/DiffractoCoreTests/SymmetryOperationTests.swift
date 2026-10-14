//
//  SymmetryOperationTests.swift
//  The x,y,z grammar (registration 3.3): what parses, what it means, what is refused, and the round trip.
//

import XCTest
import simd
import DiffractoCore

final class SymmetryOperationTests: XCTestCase {

    func testIdentityAndTranslatedIdentity() throws {
        let identity = try SymmetryOperation(xyz: "x,y,z")
        XCTAssertEqual(identity.rotation, .identity)
        XCTAssertEqual(identity.translation, SIMD3(0, 0, 0))
        XCTAssertTrue(identity.isIdentity)
        XCTAssertTrue(try SymmetryOperation(xyz: "x+1,y,z").isIdentity)
        XCTAssertFalse(try SymmetryOperation(xyz: "x+1/2,y,z").isIdentity)
        XCTAssertFalse(try SymmetryOperation(xyz: "-x,-y,-z").isIdentity)
    }

    /// `-x,y+1/2,-z`: the 2₁ screw along b. Mutant: "1/2+x" read as "x".
    func testScrewAxisParsesRotationAndTranslation() throws {
        let op = try SymmetryOperation(xyz: "-x,y+1/2,-z")
        XCTAssertEqual(op.rotation, IntegerMatrix([-1, 0, 0], [0, 1, 0], [0, 0, -1]))
        XCTAssertEqual(op.translation, SIMD3(0, 0.5, 0))
        let centring = try SymmetryOperation(xyz: "1/2+x,1/2-y,z")
        XCTAssertEqual(centring.rotation, IntegerMatrix([1, 0, 0], [0, -1, 0], [0, 0, 1]))
        XCTAssertEqual(centring.translation, SIMD3(0.5, 0.5, 0))
        let hex = try SymmetryOperation(xyz: "x-y,x,-z+1/3")
        XCTAssertEqual(hex.rotation, IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, -1]))
        XCTAssertEqual(hex.translation.z, 1.0 / 3, accuracy: 1e-15)
        let spaced = try SymmetryOperation(xyz: " -X , Y + 1/2 , -Z ")
        XCTAssertEqual(spaced.rotation, op.rotation)
        XCTAssertEqual(spaced.translation, op.translation)
        let decimal = try SymmetryOperation(xyz: "x+0.5,y,z")
        XCTAssertEqual(decimal.translation.x, 0.5)
    }

    func testApplyAndCompose() throws {
        let op = try SymmetryOperation(xyz: "-x,y+1/2,-z")
        let p = SIMD3(0.1, 0.2, 0.3)
        XCTAssertEqual(op.apply(p), SIMD3(-0.1, 0.7, -0.3))
        // op ∘ op = (x, y + 1, z): a lattice translation, the identity mod 1.
        let twice = op.composed(with: op)
        XCTAssertTrue(twice.isIdentity)
        XCTAssertEqual(twice.translationMod1, SIMD3(0, 0, 0))
        // Composition order: self(other(p)).
        let a = try SymmetryOperation(xyz: "-y,x,z"), b = try SymmetryOperation(xyz: "x,y,z+1/2")
        XCTAssertEqual(a.composed(with: b).apply(p), a.apply(b.apply(p)))
    }

    func testRefusals() {
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x,y")) { XCTAssertEqual($0 as? SymmetryOperationError, .unparsable("x,y")) }
        XCTAssertThrowsError(try SymmetryOperation(xyz: "a,b,c")) { XCTAssertEqual($0 as? SymmetryOperationError, .unparsable("a,b,c")) }
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x,,z"))
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x+,y,z"))
        // "x/2" is "x" then the fraction "/2" with an empty numerator: unparsable (the apps' parser read it as x + 2).
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x/2,y,z")) {
            XCTAssertEqual($0 as? SymmetryOperationError, .unparsable("x/2,y,z"))
        }
        XCTAssertThrowsError(try SymmetryOperation(xyz: "0.5x,y,z")) {
            XCTAssertEqual($0 as? SymmetryOperationError, .nonIntegerRotation("0.5x,y,z"))
        }
        XCTAssertThrowsError(try SymmetryOperation(xyz: "1/3x,y,z")) {
            XCTAssertEqual($0 as? SymmetryOperationError, .nonIntegerRotation("1/3x,y,z"))
        }
        // A garbled coefficient or translation is refused, not trapped (Int(1e21) would trap).
        XCTAssertThrowsError(try SymmetryOperation(xyz: "1000000000000000000000x,y,z")) {
            XCTAssertEqual($0 as? SymmetryOperationError, .unparsable("1000000000000000000000x,y,z"))
        }
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x+1000000000000000000000,y,z"))
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x+inf,y,z"))
        XCTAssertThrowsError(try SymmetryOperation(xyz: "x/0,y,z"))
    }

    /// `xyzString` parses back to an equivalent operation, for every operation of C2/m, P6₃/mmc's screw and a 1/3.
    func testRoundTrip() throws {
        for text in ["x,y,z", "-x,y,-z", "1/2+x,1/2+y,z", "1/2-x,1/2-y,-z", "x-y,x,z+1/2", "-y,x-y,z+1/3",
                     "y,x,-z+3/4", "-x+y,y,-z"] {
            let op = try SymmetryOperation(xyz: text)
            let back = try SymmetryOperation(xyz: op.xyzString)
            XCTAssertTrue(back.isEquivalent(to: op), "\(text) -> \(op.xyzString)")
        }
        XCTAssertEqual(try SymmetryOperation(xyz: "1/2+x,1/2+y,z").xyzString, "x+1/2,y+1/2,z")
        XCTAssertEqual(try SymmetryOperation(xyz: "-x,y,-z").xyzString, "-x,y,-z")
    }

    /// Distinct operations: rotation exact, translation mod 1 to 1e-3; a garbled coefficient counts, not traps.
    func testDistinctCount() throws {
        let ops = try ["x,y,z", "x+1,y,z", "-x,-y,-z", "-x,-y,-z", "x+1/2,y,z", "x+0.5001,y,z"].map { try SymmetryOperation(xyz: $0) }
        XCTAssertEqual(SymmetryOperation.distinctCount(ops), 3)
        // A whole-number translation, however large, is the identity mod 1 (1e21 is an integer-valued double) and
        // must not trap in the quantisation; a non-finite one counts as its own operation.
        let huge = SymmetryOperation(rotation: .identity, translation: SIMD3(1e21, 0, 0))
        XCTAssertEqual(SymmetryOperation.distinctCount([huge, .identity]), 1)
        let infinite = SymmetryOperation(rotation: .identity, translation: SIMD3(.infinity, 0, 0))
        XCTAssertEqual(SymmetryOperation.distinctCount([infinite, .identity]), 2)
    }

    func testOrdersAndDeterminants() {
        let c6 = IntegerMatrix([1, -1, 0], [1, 0, 0], [0, 0, 1])
        XCTAssertEqual(c6.order, 6)
        XCTAssertEqual(c6.determinant, 1)
        XCTAssertEqual(c6.multiplied(by: c6).order, 3)
        XCTAssertEqual(IntegerMatrix.inversion.order, 2)
        XCTAssertEqual(IntegerMatrix.inversion.determinant, -1)
        XCTAssertEqual(IntegerMatrix.identity.order, 1)
        XCTAssertNil(IntegerMatrix([1, 1, 0], [0, 1, 0], [0, 0, 1]).order, "a shear has no finite order")
        XCTAssertEqual(c6.transposed.multiplied(by: c6).trace, 3 + 1, "the transpose is not the inverse in the hexagonal basis")
    }
}
