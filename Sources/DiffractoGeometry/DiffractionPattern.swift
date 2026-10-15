//
//  DiffractionPattern.swift
//  Role: the three-field pattern the ring and ellipse fit reads: width, height and row-major Float32 pixels. Slice 3b of
//        DiffractoKit (docs/registration-slice-3.md section 2).
//
//  SEED: reference commit 7d6fdbd1 `Core/Data/DiffractionPattern.swift`, only the 3-field struct (its display helpers
//  `FloatImage`, `RGBAImage`, the log contrast and the -1 masked-pixel sentinel are presentation and stay in the reference implementations). The
//  seed names its fields `qy` (the ROW count) and `qx` (the COLUMN count); here they are `height` and `width`. The fit's
//  coordinates keep py4DSTEM's order: `qx` is the ROW index and `qy` the COLUMN index, so the pixel at row r and column c is
//  `pixels[r * width + c]` and a seed `DiffractionPattern(qy: rows, qx: columns, pixels:)` is `DiffractionPattern(width:
//  columns, height: rows, pixels:)`. Pixels are Float32 although the fits compute in Double.
//

/// A diffraction pattern: `height` rows by `width` columns of Float32, rows first (row-major, columns fastest).
public struct DiffractionPattern: Sendable, Equatable {
    /// The number of columns.
    public let width: Int
    /// The number of rows.
    public let height: Int
    public let pixels: [Float]

    // Explicit so the memberwise initializer is `public` (synthesized ones are internal).
    public init(width: Int, height: Int, pixels: [Float]) {
        self.width = width
        self.height = height
        self.pixels = pixels
    }
}
