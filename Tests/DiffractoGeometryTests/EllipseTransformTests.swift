//
//  EllipseTransformTests.swift
//  DiffractoKit slice 3b: `EllipseTransform` (reference commit 7d6fdbd1 `Calibration.ellipseTransform` / `ellipseCorrectedOffset` /
//  `ellipseUncorrectedOffset` / `acceptedEllipse`) against hand arithmetic, and the four `CalibrationEllipseRow29Tests` of the
//  seed re-expressed on the kit's type (same names, same numbers; the seed's `Calibration` becomes the failable `EllipseTransform`
//  and its static `correctedOffset`, which leaves the offset alone where the seed's `hasEllipse` was false).
//
//  Hand arithmetic (e = b/a, s = sin(theta - pi/2), c = cos(theta - pi/2), t00 = e s^2 + c^2, t01 = s c (1 - e),
//  t11 = s^2 + e c^2, qx' = t00 dy + t01 dx, qy' = t01 dy + t11 dx, returned (x: qy', y: qx')):
//    theta = pi/2: s = 0, c = 1: t00 = 1, t01 = 0, t11 = e, so (dx, dy) -> (e dx, dy).
//    theta = 0:    s = -1, c ~ 0: t00 = e, t01 ~ 0, t11 = 1, so (dx, dy) -> (dx, e dy).
//  The ellipse (a, b, theta) lands on a circle of radius b: a point at distance a along (cos theta, sin theta) in (row, column)
//  keeps its length times b/a, a point at distance b across keeps b.
//
//  RED UNDER: e = a/b for e = b/a; the axis swap of the returned pair dropped; the inverse's signs; `accepted` letting b <= 0 pass.
//

import DiffractoKitTestSupport
import Testing
import Foundation
import DiffractoGeometry

struct EllipseTransformTests {

    @Test func testTheMatrixIsHandArithmeticAtThetaZeroAndAtHalfPi() throws {
        // a = 20, b = 10: e = 0.5. Offsets are Float, the matrix is Double: exact on these values.
        let atHalfPi = try #require(EllipseTransform(a: 20, b: 10, theta: .pi / 2))
        let first = atHalfPi.correctedOffset(dx: 4, dy: 3)
        #expect(first.x == 2 && first.y == 3, "theta pi/2 halves the column offset: \(first)")
        let atZero = try #require(EllipseTransform(a: 20, b: 10, theta: 0))
        let second = atZero.correctedOffset(dx: 4, dy: 3)
        #expect(abs(second.x - 4) < 1e-6 && abs(second.y - 1.5) < 1e-6, "theta 0 halves the row offset: \(second)")
    }

    @Test func testAnEllipseIsCorrectedOntoACircleOfTheMinorRadius() throws {
        // A point of the ellipse a = 21, b = 19, theta = 0.6 (a along (cos theta, sin theta) in (row, column)).
        let a = 21.0, b = 19.0, theta = 0.6
        let transform = try #require(EllipseTransform(a: a, b: b, theta: theta))
        for phi in stride(from: 0.0, to: 6.3, by: 0.7) {
            let row = a * cos(phi) * cos(theta) - b * sin(phi) * sin(theta)
            let column = a * cos(phi) * sin(theta) + b * sin(phi) * cos(theta)
            let corrected = transform.correctedOffset(dx: Float(column), dy: Float(row))
            let radius = Double(corrected.x * corrected.x + corrected.y * corrected.y).squareRoot()
            #expect(abs(radius - b) < 1e-4, "phi \(phi): radius \(radius) against \(b)")
        }
    }

    @Test func testTheInverseUndoesTheCorrectionAndKeepsAnSemiMinorFirstEllipse() throws {
        for (a, b, theta) in [(19.5, 17.97, 0.4), (17.97, 19.5, 0.4), (30.0, 12.0, 2.2)] {
            let transform = try #require(EllipseTransform(a: a, b: b, theta: theta))
            for (dx, dy) in [(3.0, 4.0), (-7.5, 2.25), (0.0, -9.0)] as [(Float, Float)] {
                let there = transform.correctedOffset(dx: dx, dy: dy)
                let back = transform.uncorrectedOffset(dx: there.x, dy: there.y)
                #expect(abs(back.x - dx) < 1e-4 && abs(back.y - dy) < 1e-4, "a \(a) b \(b): \(back) against (\(dx), \(dy))")
            }
        }
    }

    // MARK: The seed's CalibrationEllipseRow29Tests

    @Test func testRow29NonPositiveAxesAreNotAnEllipseAndCorrectNothing() {
        for (a, b) in [(19.5, 0.0), (19.5, -17.9), (-19.5, 17.9), (-19.5, -17.9), (0.0, 5.0)] {
            #expect(EllipseTransform(a: a, b: b, theta: 0.4) == nil, "a \(a), b \(b)")
            let out = EllipseTransform.correctedOffset(dx: 3, dy: 4, a: a, b: b, theta: 0.4)
            #expect(out.x == 3, "a \(a), b \(b) must leave the vector alone")
            #expect(out.y == 4, "a \(a), b \(b) must leave the vector alone")
        }
    }

    @Test func testRow29ValidEllipseStillCorrects() {
        #expect(EllipseTransform(a: 19.5, b: 17.97, theta: 0.4) != nil)
        let out = EllipseTransform.correctedOffset(dx: 3, dy: 4, a: 19.5, b: 17.97, theta: 0.4)
        #expect(out.x != 3)
    }

    /// An a < b ellipse is a valid transform (py4DSTEM's `_transform` uses e = b/a with no ordering), so an imported one is kept.
    @Test func testRow29SemiMinorFirstFileEllipseIsKept() {
        #expect(EllipseTransform(a: 17.97, b: 19.5, theta: 0.4) != nil)
    }

    @Test func testRow29AcceptedEllipseDropsTheWholeTripleOnABadAxis() {
        #expect(EllipseTransform.accepted(a: 19.5, b: -1, theta: 0.3) == nil)
        #expect(EllipseTransform.accepted(a: nil, b: 1, theta: 0.3) == nil)
        #expect(EllipseTransform.accepted(a: 19.5, b: 1, theta: .nan) == nil)
        #expect(EllipseTransform.accepted(a: .infinity, b: 1, theta: 0.3) == nil)
        let ok = EllipseTransform.accepted(a: 19.5, b: 17.9, theta: 0.3)
        #expect(ok?.a == 19.5)
        #expect(ok?.b == 17.9)
        #expect(ok?.theta == 0.3)
    }
}
