// Compiled with the KIT's EllipseCalibration.swift, DiffractionPattern.swift, RingCalibration.swift and ElectronWavelength.swift
// (run.sh, kit-check step) so that main.swift, written against the apps' `DiffractionPattern(qy: rows, qx: columns, pixels:)`,
// runs unchanged on the kit's pattern type. The kit's `height` is the apps' `qy`, its `width` the apps' `qx`.
extension DiffractionPattern {
    init(qy: Int, qx: Int, pixels: [Float]) { self.init(width: qx, height: qy, pixels: pixels) }
}
