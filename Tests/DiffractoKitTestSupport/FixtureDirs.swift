//
//  FixtureDirs.swift
//  Role: where the committed fixtures are. Tests read them by #filePath, not as resources; the fixture folders stay
//        beside the test targets they were written for, so a parity test names the target it reads from.
//

import Foundation

package enum FixtureDirs {
    /// The package's Tests directory.
    package static let tests: URL = URL(fileURLWithPath: #filePath)
        .deletingLastPathComponent().deletingLastPathComponent()

    /// `Tests/<target>`.
    package static func target(_ name: String) -> URL { tests.appendingPathComponent(name) }

    /// `Tests/<target>/Fixtures`.
    package static func fixtures(of name: String) -> URL { target(name).appendingPathComponent("Fixtures") }
}
