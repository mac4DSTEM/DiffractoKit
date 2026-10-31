//
//  FriedelCentre.swift
//  Role: the centre of one pattern whose unscattered beam is blocked (a beam stop, a saturated or masked centre), from
//        Friedel symmetry. A port of py4DSTEM f050d207 `get_origin_friedel` (process/calibration/origin.py, line 332)
//        for one pattern with upsample_factor 1, CPU Double. A kinematic pattern is centrosymmetric about the
//        unscattered beam, so the correlation of the zero-padded pattern with its point reflection,
//        Re(IFFT2(FFT2(im)^2)), peaks at twice the centre whether or not the beam itself is visible. With a mask (true in
//        the pattern, false under the beam stop) it is the reference's three-term masked correlation, transcribed as
//        written (the padding of the mask is ones, the image itself is not masked in the first two terms).
//        DEVIATIONS: the third term, Re(IFFT2(FFT2(im mask_pad))), is taken as im mask_pad exactly (the reference carries
//        about 1e-16 relative rounding there, so where the denominator vanishes the port gets 0/0 and refuses); and the
//        refusals listed on `find` (an all-false mask, a non-finite pixel or result) return nil where the reference
//        returns a number, inf or NaN.
//
//  Checked against the reference on synthetic blocked patterns at assert_allclose float64 defaults (rtol 1e-7; the
//  reference has no test of the function): `FriedelCentreTests`, fixture written by tools/pins/beam-centre.
//  Not ported: upsample_factor above 1 (the DFT-upsampled peak) and the GPU path. `validation: "none"`.
//

import Foundation
import ImagoCore

/// The centre of a pattern found from its Friedel symmetry.
public struct FriedelCentreFit: Sendable, Equatable {
    public let row: Double
    public let column: Double
    /// "none": no standard has validated the finder.
    public var validation: String { "none" }
}

public enum FriedelCentre {
    /// The centre of `pattern` from Friedel symmetry: py4DSTEM `get_origin_friedel` on one pattern, upsample_factor 1 (the
    /// correlation peak refined by a three-point parabola on each axis). Use it where `PatternCentre.find` refuses: a
    /// beam stop, a saturated or masked centre. `mask` (row-major, `width * height`) is true in the pattern and false under the
    /// beam stop; without one the blocked region simply counts as dark, which the reference says may fail with a large stop.
    ///
    /// Returns nil where the reference would give no usable number: an empty pattern or one under 2 x 2, a pixel count or
    /// mask length that does not match, a mask that keeps no pixel, a pixel that is not finite (the reference would return
    /// the first NaN's position), a correlation peak on the last row or column of the padded grid (the reference raises
    /// IndexError there), or a parabola with a zero denominator (the reference returns inf or NaN). Coordinates are
    /// (row, column) in pixels, wrapped by NumPy's `% shape` as the reference does (so in [0, size], the upper end only for a
    /// tiny negative value).
    public static func find(pattern: DiffractionPattern, mask: [Bool]? = nil) -> FriedelCentreFit? {
        let rows = pattern.height, columns = pattern.width
        guard rows >= 2, columns >= 2, pattern.pixels.count == rows * columns, pattern.pixels.allSatisfy(\.isFinite) else { return nil }
        if let mask {
            guard mask.count == rows * columns, mask.contains(true) else { return nil }
        }
        let ny = 2 * rows, nx = 2 * columns, total = ny * nx
        guard let fft = FFT2DDouble(nx: nx, ny: ny) else { return nil }
        let zeros = [Double](repeating: 0, count: total)

        // im: the pattern at the top-left of a grid twice its size, zeros elsewhere
        var image = zeros
        for r in 0 ..< rows { for c in 0 ..< columns { image[r * nx + c] = Double(pattern.pixels[r * columns + c]) } }
        // IFFT2(FFT2(im)^2), complex; NumPy's ifft2 divides by the size, as scaleInverse does
        var gRe = image, gIm = zeros
        fft.transform(re: &gRe, im: &gIm, forward: true)
        square(&gRe, &gIm)
        fft.transform(re: &gRe, im: &gIm, forward: false, scaleInverse: true)

        var cc: [Double]
        if let mask {
            // mask_pad: the mask at the top-left, ones in the padding (np.pad constant_values 1)
            var maskPad = [Double](repeating: 1, count: total)
            for r in 0 ..< rows { for c in 0 ..< columns { maskPad[r * nx + c] = mask[r * columns + c] ? 1 : 0 } }
            var mRe = maskPad, mIm = zeros
            fft.transform(re: &mRe, im: &mIm, forward: true)                     // M = FFT2(mask_pad)
            var m2Re = mRe, m2Im = mIm
            square(&m2Re, &m2Im)
            fft.transform(re: &m2Re, im: &m2Im, forward: false, scaleInverse: true)   // IFFT2(M^2)
            var sRe = image.map { $0 * $0 }, sIm = zeros
            fft.transform(re: &sRe, im: &sIm, forward: true)                     // FFT2(im^2)
            for k in 0 ..< total {
                let a = sRe[k], b = sIm[k]
                sRe[k] = a * mRe[k] - b * mIm[k]; sIm[k] = a * mIm[k] + b * mRe[k]
            }
            fft.transform(re: &sRe, im: &sIm, forward: false, scaleInverse: true)    // IFFT2(FFT2(im^2) M)
            cc = zeros
            for k in 0 ..< total {
                let term1 = gRe[k] * m2Re[k] - gIm[k] * m2Im[k]                   // Re(IFFT2(G^2) IFFT2(M^2))
                let term2 = sRe[k]
                let term3 = image[k] * maskPad[k]                                // Re(IFFT2(FFT2(im mask_pad))) = im mask_pad
                cc[k] = (term1 - term3) / (term2 - term3)
            }
        } else {
            cc = gRe
        }

        let peak = argmax(cc)
        let x = peak / nx, y = peak % nx
        guard x + 1 < ny, y + 1 < nx else { return nil }                        // cc[x + 1, y] raises IndexError in the reference
        // index -1 wraps to the last row or column, as NumPy's negative index does
        let at = { (r: Int, c: Int) in cc[((r + ny) % ny) * nx + (c + nx) % nx] }
        let dx = parabola(cc[peak], at(x + 1, y), at(x - 1, y))
        let dy = parabola(cc[peak], at(x, y + 1), at(x, y - 1))
        let row = remainder((Double(x) + dx) / 2, Double(rows))
        let column = remainder((Double(y) + dy) / 2, Double(columns))
        guard row.isFinite, column.isFinite else { return nil }
        return FriedelCentreFit(row: row, column: column)
    }

    /// (re, im) := (re, im)^2, elementwise.
    static func square(_ re: inout [Double], _ im: inout [Double]) {
        for k in re.indices { let a = re[k], b = im[k]; re[k] = a * a - b * b; im[k] = 2 * a * b }
    }

    /// NumPy's argmax: the first NaN if there is one, else the first maximum in row-major order.
    static func argmax(_ v: [Double]) -> Int {
        if let nan = v.firstIndex(where: \.isNaN) { return nan }
        var best = 0
        for i in 1 ..< v.count where v[i] > v[best] { best = i }
        return best
    }

    /// The reference's `(c+ - c-) / (4 c0 - 2 c+ - 2 c-)`, in that order of operations.
    static func parabola(_ c0: Double, _ plus: Double, _ minus: Double) -> Double {
        (plus - minus) / (4.0 * c0 - 2.0 * plus - 2.0 * minus)
    }

    /// NumPy's `%` on float64 (np.remainder): the result takes the divisor's sign.
    static func remainder(_ a: Double, _ b: Double) -> Double {
        var mod = fmod(a, b)
        if mod != 0 { if (b < 0) != (mod < 0) { mod += b } } else { mod = copysign(0, b) }
        return mod
    }
}
