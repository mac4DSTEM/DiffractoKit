//
//  CrystalLoading.swift
//  Role: a crystal from the bytes of a CIF a person chose: the size cap, strict UTF-8, `CIFImport`, and a SHA-256 of the
//        bytes as the file's fingerprint. Pure: it reads no file; the caller hands over the bytes and the file's name.
//
//  Ported from macTEM `TEMCore/Diffraction/CrystalLoading.swift` (the app's source is the oracle; no Python reference
//  exists). `CIFImportError` propagates unchanged; this adds only the two refusals that happen before the importer sees text.
//

import Foundation
import CryptoKit

/// A crystal read from the bytes of a CIF, with the fingerprint of those bytes.
public struct LoadedCrystal: Sendable, Equatable {
    /// The importer's document.
    public let document: CIFDocument
    /// The name to show: the `data_` block's, else the file's name without its extension.
    public var displayName: String { document.displayName }
    /// The structure of the document.
    public var structure: CrystalStructure { document.structure }
    /// Lower-case hex SHA-256 of the file's bytes exactly as given (not of any decoded or re-encoded text).
    public let sha256: String
    /// The number of bytes read.
    public let byteCount: Int
}

/// Loads a CIF from bytes: size cap, strict UTF-8, import, fingerprint.
public enum CrystalLoading {

    /// The most bytes a CIF may have: 5 MiB. One byte more is refused.
    public static let maximumBytes = 5 * 1024 * 1024

    /// Why a file is refused before the importer reads it.
    public enum Refusal: Error, Equatable, Sendable, LocalizedError {
        /// The file has more than `CrystalLoading.maximumBytes` bytes.
        case tooLarge(bytes: Int)
        /// The bytes are not valid UTF-8.
        case notUTF8

        /// The case's name without its numbers.
        public var name: String {
            switch self {
            case .tooLarge: return "tooLarge"
            case .notUTF8: return "notUTF8"
            }
        }

        /// One sentence for the person reading it.
        public var errorDescription: String? {
            switch self {
            case .tooLarge(let bytes):
                return "This file is \(bytes) bytes; a structure file larger than \(CrystalLoading.maximumBytes) bytes (5 MiB) is not read."
            case .notUTF8:
                return "This file is not valid UTF-8 text, so it cannot be a CIF."
            }
        }
    }

    /// Reads a CIF. Throws `Refusal` (size first, then encoding), then the importer's `CIFImportError` unchanged.
    /// `fileName` supplies the fallback display name (its last component without extension).
    public static func load(bytes: Data, fileName: String) throws -> LoadedCrystal {
        guard bytes.count <= maximumBytes else { throw Refusal.tooLarge(bytes: bytes.count) }
        guard let text = String(validating: bytes, as: UTF8.self) else { throw Refusal.notUTF8 }
        let base = URL(fileURLWithPath: fileName).deletingPathExtension().lastPathComponent
        let document = try CIFImport.read(text, fileBaseName: base)
        return LoadedCrystal(document: document, sha256: fingerprint(of: bytes), byteCount: bytes.count)
    }

    /// Lower-case hex SHA-256 of `bytes`.
    public static func fingerprint(of bytes: Data) -> String {
        SHA256.hash(data: bytes).map { String(format: "%02x", $0) }.joined()
    }
}
