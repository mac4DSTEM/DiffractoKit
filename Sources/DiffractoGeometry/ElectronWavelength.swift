//
//  ElectronWavelength.swift
//  Role: the one electron wavelength of the suite, the relativistic de Broglie wavelength at an accelerating voltage.
//        Slice 3a (docs/registration-slice-3.md section 1). Ported from macTEM 9482831 `RingCalibration.wavelengthAngstroms`
//        with the same arithmetic order, so the results are bit identical.
//
//  λ = hc / √(eV (2 m₀c² + eV)),  V in volts, λ in Å.
//
//  CONSTANTS, CODATA 2018 (Tiesinga, Mohr, Newell, Taylor, Rev. Mod. Phys. 93, 025010 (2021)):
//    h  = 6.626 070 15e-34 J s       exact (the 2019 SI definition)
//    c  = 299 792 458 m/s            exact
//    e  = 1.602 176 634e-19 C        exact
//    m₀ = 9.109 383 7015(28)e-31 kg  recommended value (relative standard uncertainty 3.0e-10)
//  hc/e (12398.41984332003 eV·Å) and m₀c²/e (510998.9499961642 eV) are derived from these in code, never typed; the
//  second agrees with CODATA's own electron rest energy, 0.510 998 950 00(15) MeV, to 7.5e-12 relative.
//
//  Outside its domain (V not finite, or not positive) the function returns NaN, never a plausible number.
//

import Foundation

/// The relativistic electron wavelength (CODATA 2018 constants). The one wavelength of the suite.
public enum ElectronWavelength {
    /// Planck constant, J s (exact).
    public static let planckConstant = 6.626_070_15e-34
    /// Speed of light in vacuum, m/s (exact).
    public static let speedOfLight = 299_792_458.0
    /// Elementary charge, C (exact).
    public static let elementaryCharge = 1.602_176_634e-19
    /// Electron rest mass, kg (CODATA 2018 recommended value).
    public static let electronMass = 9.109_383_7015e-31
    /// hc/e in eV·Å: 12398.41984332003.
    public static let planckTimesLightEVAngstroms = planckConstant * speedOfLight / elementaryCharge * 1e10
    /// m₀c²/e in eV: 510998.9499961642.
    public static let electronRestEnergyEV = electronMass * speedOfLight * speedOfLight / elementaryCharge

    /// λ in Å at an accelerating voltage in VOLTS (200 kV is 200_000). NaN unless volts is finite and positive.
    public static func angstroms(volts: Double) -> Double {
        guard volts.isFinite, volts > 0 else { return .nan }
        return planckTimesLightEVAngstroms / (volts * (2 * electronRestEnergyEV + volts)).squareRoot()
    }
}
