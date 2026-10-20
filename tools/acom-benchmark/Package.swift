// swift-tools-version: 6.2
// DiffractoKit ACOM benchmark: a SEPARATE package that depends on DiffractoKit by path and times `ACOMMatcher.matchAll`
// (the speed number in docs/parity.md P25). Run it with `run.sh`; it prints, never asserts.
import Foundation
import PackageDescription

let kit = URL(fileURLWithPath: Context.packageDirectory).deletingLastPathComponent().deletingLastPathComponent()
    .lastPathComponent

let package = Package(
    name: "ACOMBenchmark",
    platforms: [.macOS("15.0")],
    dependencies: [.package(path: "../..")],
    targets: [
        .executableTarget(name: "acom-benchmark", dependencies: [
            .product(name: "DiffractoOrientation", package: kit),
            .product(name: "DiffractoCore", package: kit),
        ]),
    ],
    swiftLanguageModes: [.v6]
)
