//
//  CrystalLoadingTests.swift
//  `CrystalLoading` (ported from macTEM's CrystalLoading.swift; the app's source is the oracle, there is no Python
//  reference). The SHA-256 values below are NOT from the Swift code: they were made with `shasum -a 256` on the same bytes
//  (the NIST vectors "abc" and the empty string, the 5 MiB of zeros, and the gold CIF below written with `printf`, one
//  line per row, a final newline). In-memory data only.
//

import DiffractoKitTestSupport
import Foundation
import Testing
import DiffractoCore

struct CrystalLoadingTests {

    static let goldLines = [
        "data_gold", "_cell_length_a 4.0782", "_cell_length_b 4.0782", "_cell_length_c 4.0782",
        "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_",
        "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "Au1 Au 0 0 0", "Au2 Au 0.5 0.5 0", "Au3 Au 0.5 0 0.5", "Au4 Au 0 0.5 0.5",
    ]
    static var goldBytes: Data { Data((goldLines.joined(separator: "\n") + "\n").utf8) }
    static let goldSHA = "96ca667484cbf141583df74ecbf0af043c3907f53250f014c96d2a1a9326e67a"

    @Test func fingerprintMatchesShasum() {
        #expect(CrystalLoading.fingerprint(of: Data("abc".utf8)) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
        #expect(CrystalLoading.fingerprint(of: Data()) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
        #expect(CrystalLoading.fingerprint(of: Data(count: CrystalLoading.maximumBytes))
                == "c036cbb7553a909f8b8877d4461924307f27ecb66cff928eeeafd569c3887e29")
    }

    @Test func loadsAGoldCIFWithItsFingerprint() throws {
        let loaded = try CrystalLoading.load(bytes: Self.goldBytes, fileName: "ignored/other.cif")
        #expect(loaded.sha256 == Self.goldSHA)
        #expect(loaded.byteCount == Self.goldBytes.count)
        #expect(loaded.displayName == "gold")                  // the data_ block's name wins
        #expect(loaded.structure.sites.count == 4)
    }

    @Test func fileNameIsTheFallbackDisplayName() throws {
        let anonymous = Data(Self.goldBytes.dropFirst("data_gold\n".utf8.count))
        let loaded = try CrystalLoading.load(bytes: Data("data_\n".utf8) + anonymous, fileName: "/x/y/my crystal.cif")
        #expect(loaded.displayName == "my crystal")
    }

    @Test func theCapIsExactlyFiveMiB() throws {
        #expect(CrystalLoading.maximumBytes == 5_242_880)
        // A comment line pads the CIF to exactly the cap: read; one byte more: refused before the importer.
        var padded = Self.goldBytes + Data("#".utf8)
        padded.append(Data(repeating: UInt8(ascii: " "), count: CrystalLoading.maximumBytes - padded.count))
        #expect(padded.count == CrystalLoading.maximumBytes)
        let loaded = try CrystalLoading.load(bytes: padded, fileName: "a.cif")
        #expect(loaded.byteCount == CrystalLoading.maximumBytes)
        padded.append(UInt8(ascii: " "))
        #expect(throws: CrystalLoading.Refusal.tooLarge(bytes: CrystalLoading.maximumBytes + 1)) {
            try CrystalLoading.load(bytes: padded, fileName: "a.cif")
        }
    }

    @Test func sizeIsCheckedBeforeEncoding() {
        let big = Data(repeating: 0xFF, count: CrystalLoading.maximumBytes + 1)
        #expect(throws: CrystalLoading.Refusal.tooLarge(bytes: big.count)) { try CrystalLoading.load(bytes: big, fileName: "a.cif") }
    }

    @Test func invalidUTF8IsRefused() {
        let bad = Self.goldBytes + Data([0xFF, 0xFE])
        #expect(throws: CrystalLoading.Refusal.notUTF8) { try CrystalLoading.load(bytes: bad, fileName: "a.cif") }
        // Overlong encoding of "/" and a lone continuation byte are not UTF-8 either.
        #expect(throws: CrystalLoading.Refusal.notUTF8) { try CrystalLoading.load(bytes: Data([0xC0, 0xAF]), fileName: "a.cif") }
        #expect(throws: CrystalLoading.Refusal.notUTF8) { try CrystalLoading.load(bytes: Data([0x80]), fileName: "a.cif") }
    }

    @Test func importerErrorsPropagateUnchanged() {
        #expect(throws: CIFImportError.noAtomSites) { try CrystalLoading.load(bytes: Data("data_x\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n".utf8), fileName: "x.cif") }
    }

    @Test func refusalSentencesAndNames() {
        let tooLarge = CrystalLoading.Refusal.tooLarge(bytes: 6_000_000)
        #expect(tooLarge.name == "tooLarge")
        #expect(tooLarge.errorDescription == "This file is 6000000 bytes; a structure file larger than 5242880 bytes (5 MiB) is not read.")
        #expect(CrystalLoading.Refusal.notUTF8.name == "notUTF8")
        #expect(CrystalLoading.Refusal.notUTF8.errorDescription == "This file is not valid UTF-8 text, so it cannot be a CIF.")
    }
}
