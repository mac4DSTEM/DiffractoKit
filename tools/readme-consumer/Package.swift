// swift-tools-version: 6.2
// DiffractoKit readme-consumer: a SEPARATE package that depends on DiffractoKit by path and compiles, from outside the
// package, exactly the program in the README's "First ten minutes" (Sources/readme-consumer/main.swift is that code and
// nothing else; run.sh checks it). See README.md.
import Foundation
import PackageDescription

// A path dependency's identity is its directory's name (lower-cased): "diffractokit" in the repository, but not in a
// scratch copy; `.product` needs that name, so it is taken from the directory two levels above this manifest.
let kit = URL(fileURLWithPath: Context.packageDirectory).deletingLastPathComponent().deletingLastPathComponent()
    .lastPathComponent

let package = Package(
    name: "ReadmeConsumer",
    platforms: [.macOS("27.0")],
    dependencies: [
        .package(path: "../.."),
    ],
    targets: [
        .executableTarget(name: "readme-consumer", dependencies: [
            .product(name: "DiffractoCore", package: kit),
            .product(name: "DiffractoScattering", package: kit),
            .product(name: "DiffractoGeometry", package: kit),
        ]),
    ],
    swiftLanguageModes: [.v6]
)
