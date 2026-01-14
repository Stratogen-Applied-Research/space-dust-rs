//! Physical and mathematical constants used throughout the library.
//!
//! This module provides commonly used constants for astrodynamics calculations,
//! including angular conversions, time constants, and physical parameters.

use std::f64::consts::PI;

/// Collection of physical and mathematical constants.
pub struct Constants;

impl Constants {
    // ============================================================================
    // Mathematical Constants
    // ============================================================================

    /// Two times pi (2π)
    pub const TWO_PI: f64 = 2.0 * PI;

    /// Pi
    pub const PI: f64 = PI;

    // ============================================================================
    // Angular Conversion Constants
    // ============================================================================

    /// Degrees to radians conversion factor
    pub const DEG_TO_RAD: f64 = PI / 180.0;

    /// Radians to degrees conversion factor
    pub const RAD_TO_DEG: f64 = 180.0 / PI;

    /// Arcseconds to radians conversion factor
    pub const ARCSEC_TO_RAD: f64 = PI / 648_000.0;

    /// Radians to arcseconds conversion factor
    pub const RAD_TO_ARCSEC: f64 = 648_000.0 / PI;

    /// Thousandths of arcseconds to radians (for TT nutation calculations)
    pub const TT_ARCSEC_TO_RAD: f64 = Self::ARCSEC_TO_RAD / 1000.0;

    // ============================================================================
    // Time Constants
    // ============================================================================

    /// Seconds per day
    pub const SECONDS_PER_DAY: f64 = 86400.0;

    /// Seconds per hour
    pub const SECONDS_PER_HOUR: f64 = 3600.0;

    /// Seconds per minute
    pub const SECONDS_PER_MINUTE: f64 = 60.0;

    /// Days per Julian century
    pub const JULIAN_CENTURY: f64 = 36525.0;

    /// Julian Date of J2000.0 epoch (2000-01-01 12:00:00 TT)
    pub const J2000_JD: f64 = 2_451_545.0;

    /// Julian Date of Unix epoch (1970-01-01 00:00:00 UTC)
    pub const UNIX_EPOCH_JD: f64 = 2_440_587.5;

    /// Modified Julian Date offset
    pub const MJD_OFFSET: f64 = 2_400_000.5;

    /// GPS epoch in Unix seconds (1980-01-06 00:00:00 UTC)
    pub const GPS_EPOCH_UNIX: i64 = 315_964_800;

    /// Seconds per GPS week
    pub const SECONDS_PER_GPS_WEEK: f64 = 604_800.0;

    /// TT - TAI offset in seconds
    pub const TT_TAI_OFFSET: f64 = 32.184;

    /// GPS - TAI offset in seconds
    pub const GPS_TAI_OFFSET: f64 = -19.0;

    // ============================================================================
    // Earth Parameters (WGS84)
    // ============================================================================

    /// Earth gravitational parameter (μ) in km³/s²
    pub const EARTH_MU_KM: f64 = 398600.4418;

    /// Earth gravitational parameter (μ) in m³/s²
    pub const EARTH_MU_M: f64 = 3.986004418e14;

    /// Earth equatorial radius in km (WGS84)
    pub const EARTH_RADIUS_EQ_KM: f64 = 6378.137;

    /// Earth equatorial radius in m (WGS84)
    pub const EARTH_RADIUS_EQ_M: f64 = 6_378_137.0;

    /// Earth polar radius in km (WGS84)
    pub const EARTH_RADIUS_POLAR_KM: f64 = 6356.752314245;

    /// Earth polar radius in m (WGS84)
    pub const EARTH_RADIUS_POLAR_M: f64 = 6_356_752.314245;

    /// Earth flattening (WGS84)
    pub const EARTH_FLATTENING: f64 = 1.0 / 298.257223563;

    /// Earth first eccentricity squared (WGS84)
    pub const EARTH_E2: f64 =
        2.0 * Self::EARTH_FLATTENING - Self::EARTH_FLATTENING * Self::EARTH_FLATTENING;

    /// Earth sidereal rotation rate in rad/s
    pub const EARTH_ROTATION_RATE: f64 = 7.292115e-5;

    /// Earth sidereal rotation period in seconds
    pub const EARTH_SIDEREAL_DAY: f64 = Self::TWO_PI / Self::EARTH_ROTATION_RATE;

    /// Earth J2 zonal harmonic coefficient
    pub const EARTH_J2: f64 = 1.08262668355315e-3;

    /// Earth J3 zonal harmonic coefficient
    pub const EARTH_J3: f64 = -2.53265648533224e-6;

    /// Earth J4 zonal harmonic coefficient
    pub const EARTH_J4: f64 = -1.619621591367e-6;

    /// Earth J5 zonal harmonic coefficient
    pub const EARTH_J5: f64 = -2.27296082868698e-7;

    /// Earth J6 zonal harmonic coefficient
    pub const EARTH_J6: f64 = 5.40681239107085e-7;

    // ============================================================================
    // Sun Parameters
    // ============================================================================

    /// Sun gravitational parameter (μ) in m³/s²
    pub const SUN_MU: f64 = 1.32712440018e20;

    /// Astronomical Unit in meters
    pub const AU_M: f64 = 149_597_870_700.0;

    /// Astronomical Unit in kilometers
    pub const AU_KM: f64 = 149_597_870.7;

    /// Solar radius in meters
    pub const SUN_RADIUS_M: f64 = 6.96e8;

    /// Solar radius in kilometers
    pub const SUN_RADIUS_KM: f64 = 696_000.0;

    // ============================================================================
    // Moon Parameters
    // ============================================================================

    /// Moon gravitational parameter (μ) in m³/s²
    pub const MOON_MU: f64 = 4.902800066e12;

    /// Mean Earth-Moon distance in meters
    pub const MOON_MEAN_DISTANCE_M: f64 = 384_400_000.0;

    /// Mean Earth-Moon distance in kilometers
    pub const MOON_MEAN_DISTANCE_KM: f64 = 384_400.0;

    /// Lunar radius in meters
    pub const MOON_RADIUS_M: f64 = 1_737_400.0;

    /// Lunar radius in kilometers
    pub const MOON_RADIUS_KM: f64 = 1737.4;

    // ============================================================================
    // Physical Constants
    // ============================================================================

    /// Speed of light in m/s
    pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

    /// Gravitational constant in m³/(kg·s²)
    pub const G: f64 = 6.67430e-11;
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Convert degrees to radians
#[inline]
pub fn deg_to_rad(degrees: f64) -> f64 {
    degrees * Constants::DEG_TO_RAD
}

/// Convert radians to degrees
#[inline]
pub fn rad_to_deg(radians: f64) -> f64 {
    radians * Constants::RAD_TO_DEG
}

/// Convert arcseconds to radians
#[inline]
pub fn arcsec_to_rad(arcseconds: f64) -> f64 {
    arcseconds * Constants::ARCSEC_TO_RAD
}

/// Convert radians to arcseconds
#[inline]
pub fn rad_to_arcsec(radians: f64) -> f64 {
    radians * Constants::RAD_TO_ARCSEC
}

/// Normalize an angle to the range [0, 2π)
#[inline]
pub fn normalize_angle(angle: f64) -> f64 {
    let result = angle % Constants::TWO_PI;
    if result < 0.0 {
        result + Constants::TWO_PI
    } else {
        result
    }
}

/// Normalize an angle to the range [-π, π)
#[inline]
pub fn normalize_angle_signed(angle: f64) -> f64 {
    let result = normalize_angle(angle);
    if result >= PI {
        result - Constants::TWO_PI
    } else {
        result
    }
}

/// Normalize degrees to the range [0, 360)
#[inline]
pub fn normalize_degrees(degrees: f64) -> f64 {
    let result = degrees % 360.0;
    if result < 0.0 {
        result + 360.0
    } else {
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deg_to_rad() {
        assert!((deg_to_rad(180.0) - PI).abs() < 1e-15);
        assert!((deg_to_rad(90.0) - PI / 2.0).abs() < 1e-15);
        assert!((deg_to_rad(360.0) - Constants::TWO_PI).abs() < 1e-15);
    }

    #[test]
    fn test_rad_to_deg() {
        assert!((rad_to_deg(PI) - 180.0).abs() < 1e-12);
        assert!((rad_to_deg(PI / 2.0) - 90.0).abs() < 1e-12);
        assert!((rad_to_deg(Constants::TWO_PI) - 360.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_angle() {
        assert!((normalize_angle(0.0) - 0.0).abs() < 1e-15);
        assert!((normalize_angle(Constants::TWO_PI) - 0.0).abs() < 1e-15);
        assert!((normalize_angle(-PI) - PI).abs() < 1e-15);
        assert!((normalize_angle(3.0 * PI) - PI).abs() < 1e-15);
    }

    #[test]
    fn test_normalize_degrees() {
        assert!((normalize_degrees(0.0) - 0.0).abs() < 1e-12);
        assert!((normalize_degrees(360.0) - 0.0).abs() < 1e-12);
        assert!((normalize_degrees(-90.0) - 270.0).abs() < 1e-12);
        assert!((normalize_degrees(450.0) - 90.0).abs() < 1e-12);
    }

    #[test]
    fn test_earth_parameters() {
        // Verify WGS84 eccentricity squared calculation
        let expected_e2 = 0.00669437999014;
        assert!((Constants::EARTH_E2 - expected_e2).abs() < 1e-10);
    }
}
