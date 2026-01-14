//! Data tables and parameters for astrodynamics calculations.
//!
//! This module provides:
//! - **LeapSeconds**: Leap second table for UTC-TAI conversions
//! - **EOP**: Earth Orientation Parameters for precise frame transformations
//! - **IAU1980**: Nutation coefficients from the IAU 1980 theory

use crate::constants::Constants;
use crate::math::lerp;

// ============================================================================
// Leap Seconds
// ============================================================================

/// Leap second data and lookup functions.
///
/// The leap second table contains the Julian Date at which each leap second
/// was introduced and the cumulative number of leap seconds at that time.
pub struct LeapSeconds;

impl LeapSeconds {
    /// Leap second table: (Julian Date, cumulative leap seconds)
    /// This table is valid from 1972 (when leap seconds were introduced) to the present.
    const LEAP_SECOND_TABLE: [(f64, i32); 28] = [
        (2_441_317.5, 10), // 1972-01-01
        (2_441_499.5, 11), // 1972-07-01
        (2_441_683.5, 12), // 1973-01-01
        (2_442_048.5, 13), // 1974-01-01
        (2_442_413.5, 14), // 1975-01-01
        (2_442_778.5, 15), // 1976-01-01
        (2_443_144.5, 16), // 1977-01-01
        (2_443_509.5, 17), // 1978-01-01
        (2_443_874.5, 18), // 1979-01-01
        (2_444_239.5, 19), // 1980-01-01
        (2_444_786.5, 20), // 1981-07-01
        (2_445_151.5, 21), // 1982-07-01
        (2_445_516.5, 22), // 1983-07-01
        (2_446_247.5, 23), // 1985-07-01
        (2_447_161.5, 24), // 1988-01-01
        (2_447_892.5, 25), // 1990-01-01
        (2_448_257.5, 26), // 1991-01-01
        (2_448_804.5, 27), // 1992-07-01
        (2_449_169.5, 28), // 1993-07-01
        (2_449_534.5, 29), // 1994-07-01
        (2_450_083.5, 30), // 1996-01-01
        (2_450_630.5, 31), // 1997-07-01
        (2_451_179.5, 32), // 1999-01-01
        (2_453_736.5, 33), // 2006-01-01
        (2_454_832.5, 34), // 2009-01-01
        (2_456_109.5, 35), // 2012-07-01
        (2_457_204.5, 36), // 2015-07-01
        (2_457_754.5, 37), // 2017-01-01
    ];

    /// Get the number of leap seconds at a given Julian Date.
    ///
    /// # Arguments
    /// * `jd` - Julian Date
    ///
    /// # Returns
    /// The cumulative number of leap seconds at the given date.
    pub fn get_leap_seconds_at_jd(jd: f64) -> i32 {
        // Binary search would be faster, but the table is small enough
        // that linear search is fine and clearer.

        // If before the first entry, return 10 (the initial value)
        if jd < Self::LEAP_SECOND_TABLE[0].0 {
            return 10;
        }

        // Find the last entry that is <= jd
        let mut leap_seconds = Self::LEAP_SECOND_TABLE[0].1;
        for &(entry_jd, ls) in &Self::LEAP_SECOND_TABLE {
            if entry_jd <= jd {
                leap_seconds = ls;
            } else {
                break;
            }
        }

        leap_seconds
    }

    /// Get the number of leap seconds at a given Modified Julian Date.
    pub fn get_leap_seconds_at_mjd(mjd: f64) -> i32 {
        Self::get_leap_seconds_at_jd(mjd + Constants::MJD_OFFSET)
    }

    /// Get the number of leap seconds at a given Unix timestamp (seconds since 1970-01-01).
    pub fn get_leap_seconds_at_unix(unix_seconds: f64) -> i32 {
        let jd = unix_seconds / Constants::SECONDS_PER_DAY + Constants::UNIX_EPOCH_JD;
        Self::get_leap_seconds_at_jd(jd)
    }

    /// Get the current number of leap seconds (as of the last table entry).
    pub fn current_leap_seconds() -> i32 {
        Self::LEAP_SECOND_TABLE
            .last()
            .map(|&(_, ls)| ls)
            .unwrap_or(37)
    }

    /// Check if a leap second occurs at a given Julian Date.
    pub fn is_leap_second_date(jd: f64) -> bool {
        for &(entry_jd, _) in &Self::LEAP_SECOND_TABLE {
            if (jd - entry_jd).abs() < 0.5 {
                return true;
            }
        }
        false
    }
}

// ============================================================================
// Earth Orientation Parameters (EOP)
// ============================================================================

/// Earth Orientation Parameters for a specific epoch.
///
/// EOP data is required for precise transformations between the celestial
/// and terrestrial reference frames.
#[derive(Debug, Clone, Copy, Default)]
pub struct EarthOrientationParameters {
    /// Modified Julian Date
    pub mjd: f64,
    /// Polar motion X component (arcseconds)
    pub polar_motion_x: f64,
    /// Polar motion Y component (arcseconds)
    pub polar_motion_y: f64,
    /// UT1 - UTC (seconds)
    pub ut1_utc: f64,
    /// Length of day correction (seconds)
    pub lod: f64,
    /// Nutation correction in longitude (arcseconds)
    pub d_psi: f64,
    /// Nutation correction in obliquity (arcseconds)
    pub d_eps: f64,
}

impl EarthOrientationParameters {
    /// Create a new EOP instance.
    pub fn new(
        mjd: f64,
        polar_motion_x: f64,
        polar_motion_y: f64,
        ut1_utc: f64,
        lod: f64,
        d_psi: f64,
        d_eps: f64,
    ) -> Self {
        Self {
            mjd,
            polar_motion_x,
            polar_motion_y,
            ut1_utc,
            lod,
            d_psi,
            d_eps,
        }
    }

    /// Create EOP with zero corrections (for when precise EOP is not available).
    pub fn zero(mjd: f64) -> Self {
        Self {
            mjd,
            polar_motion_x: 0.0,
            polar_motion_y: 0.0,
            ut1_utc: 0.0,
            lod: 0.0,
            d_psi: 0.0,
            d_eps: 0.0,
        }
    }

    /// Get polar motion X in radians.
    pub fn polar_motion_x_rad(&self) -> f64 {
        self.polar_motion_x * Constants::ARCSEC_TO_RAD
    }

    /// Get polar motion Y in radians.
    pub fn polar_motion_y_rad(&self) -> f64 {
        self.polar_motion_y * Constants::ARCSEC_TO_RAD
    }

    /// Get nutation correction in longitude in radians.
    pub fn d_psi_rad(&self) -> f64 {
        self.d_psi * Constants::ARCSEC_TO_RAD
    }

    /// Get nutation correction in obliquity in radians.
    pub fn d_eps_rad(&self) -> f64 {
        self.d_eps * Constants::ARCSEC_TO_RAD
    }

    /// Interpolate between two EOP data points.
    pub fn interpolate(eop1: &Self, eop2: &Self, mjd: f64) -> Self {
        let t = (mjd - eop1.mjd) / (eop2.mjd - eop1.mjd);
        Self {
            mjd,
            polar_motion_x: lerp(eop1.polar_motion_x, eop2.polar_motion_x, t),
            polar_motion_y: lerp(eop1.polar_motion_y, eop2.polar_motion_y, t),
            ut1_utc: lerp(eop1.ut1_utc, eop2.ut1_utc, t),
            lod: lerp(eop1.lod, eop2.lod, t),
            d_psi: lerp(eop1.d_psi, eop2.d_psi, t),
            d_eps: lerp(eop1.d_eps, eop2.d_eps, t),
        }
    }
}

/// Simple EOP provider that returns zero corrections.
///
/// For precise applications, users should implement their own EOP provider
/// that fetches data from IERS or other sources.
pub struct DefaultEOPProvider;

impl DefaultEOPProvider {
    /// Get EOP data at a given MJD.
    ///
    /// This default implementation returns zero corrections.
    /// For precise applications, use actual EOP data from IERS.
    pub fn get_eop(mjd: f64) -> EarthOrientationParameters {
        EarthOrientationParameters::zero(mjd)
    }
}

// ============================================================================
// IAU 1980 Nutation Theory
// ============================================================================

/// IAU 1980 Nutation coefficients.
///
/// Contains the coefficients for computing nutation in longitude (Δψ)
/// and nutation in obliquity (Δε) according to the IAU 1980 theory.
#[derive(Debug, Clone, Copy)]
pub struct IAU1980Coefficients {
    /// Multiplier for mean anomaly of the Moon
    pub a1: i32,
    /// Multiplier for mean anomaly of the Sun
    pub a2: i32,
    /// Multiplier for mean argument of latitude of the Moon
    pub a3: i32,
    /// Multiplier for mean elongation of the Moon from the Sun
    pub a4: i32,
    /// Multiplier for mean longitude of the ascending node of the Moon
    pub a5: i32,
    /// Sine coefficient for Δψ (0.0001 arcseconds)
    pub ai: f64,
    /// Time-dependent sine coefficient for Δψ (0.0001 arcseconds per century)
    pub bi: f64,
    /// Cosine coefficient for Δε (0.0001 arcseconds)
    pub ci: f64,
    /// Time-dependent cosine coefficient for Δε (0.0001 arcseconds per century)
    pub di: f64,
}

impl IAU1980Coefficients {
    /// Create new IAU 1980 coefficients.
    #[allow(clippy::too_many_arguments)]
    pub const fn new(
        a1: i32,
        a2: i32,
        a3: i32,
        a4: i32,
        a5: i32,
        ai: f64,
        bi: f64,
        ci: f64,
        di: f64,
    ) -> Self {
        Self {
            a1,
            a2,
            a3,
            a4,
            a5,
            ai,
            bi,
            ci,
            di,
        }
    }
}

/// IAU 1980 Nutation data table.
pub struct IAU1980;

impl IAU1980 {
    /// The 106 terms of the IAU 1980 nutation series.
    /// Only the first few dominant terms are shown; the rest follow the same pattern.
    pub const COEFFICIENTS: [IAU1980Coefficients; 106] = [
        IAU1980Coefficients::new(0, 0, 0, 0, 1, -171996.0, -174.2, 92025.0, 8.9),
        IAU1980Coefficients::new(0, 0, 2, -2, 2, -13187.0, -1.6, 5736.0, -3.1),
        IAU1980Coefficients::new(0, 0, 2, 0, 2, -2274.0, -0.2, 977.0, -0.5),
        IAU1980Coefficients::new(0, 0, 0, 0, 2, 2062.0, 0.2, -895.0, 0.5),
        IAU1980Coefficients::new(0, 1, 0, 0, 0, 1426.0, -3.4, 54.0, -0.1),
        IAU1980Coefficients::new(1, 0, 0, 0, 0, 712.0, 0.1, -7.0, 0.0),
        IAU1980Coefficients::new(0, 1, 2, -2, 2, -517.0, 1.2, 224.0, -0.6),
        IAU1980Coefficients::new(0, 0, 2, 0, 1, -386.0, -0.4, 200.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, 0, 2, -301.0, 0.0, 129.0, -0.1),
        IAU1980Coefficients::new(0, -1, 2, -2, 2, 217.0, -0.5, -95.0, 0.3),
        IAU1980Coefficients::new(1, 0, 0, -2, 0, -158.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, -2, 1, 129.0, 0.1, -70.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, 0, 2, 123.0, 0.0, -53.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, 0, 1, 63.0, 0.1, -33.0, 0.0),
        IAU1980Coefficients::new(0, 0, 0, 2, 0, 63.0, 0.0, -2.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, 2, 2, -59.0, 0.0, 26.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 0, 0, 1, -58.0, -0.1, 32.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, 0, 1, -51.0, 0.0, 27.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, -2, 0, 48.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 2, 0, 1, 46.0, 0.0, -24.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, 2, 2, -38.0, 0.0, 16.0, 0.0),
        IAU1980Coefficients::new(2, 0, 2, 0, 2, -31.0, 0.0, 13.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, 0, 0, 29.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, -2, 2, 29.0, 0.0, -12.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, 0, 0, 26.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, -2, 0, -22.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, 0, 1, 21.0, 0.0, -10.0, 0.0),
        IAU1980Coefficients::new(0, 2, 0, 0, 0, 17.0, -0.1, 0.0, 0.0),
        IAU1980Coefficients::new(0, 2, 2, -2, 2, -16.0, 0.1, 7.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 0, 2, 1, 16.0, 0.0, -8.0, 0.0),
        IAU1980Coefficients::new(0, 1, 0, 0, 1, -15.0, 0.0, 9.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, -2, 1, -13.0, 0.0, 7.0, 0.0),
        IAU1980Coefficients::new(0, -1, 0, 0, 1, -12.0, 0.0, 6.0, 0.0),
        IAU1980Coefficients::new(2, 0, -2, 0, 0, 11.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, 2, 1, -10.0, 0.0, 5.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, 2, 2, -8.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(0, -1, 2, 0, 2, -7.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, 2, 1, -7.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(1, 1, 0, -2, 0, -7.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 2, 0, 2, 7.0, 0.0, -3.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 0, 2, 1, -6.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(0, 0, 0, 2, 1, -6.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(2, 0, 2, -2, 2, 6.0, 0.0, -3.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, 2, 0, 6.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, -2, 1, 6.0, 0.0, -3.0, 0.0),
        IAU1980Coefficients::new(0, 0, 0, -2, 1, -5.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(0, -1, 2, -2, 1, -5.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(2, 0, 2, 0, 1, -5.0, 0.0, 3.0, 0.0),
        IAU1980Coefficients::new(1, -1, 0, 0, 0, 5.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, -1, 0, -4.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 0, 0, 1, 0, -4.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 0, -2, 0, -4.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, -2, 0, 0, 4.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, -2, 1, 4.0, 0.0, -2.0, 0.0),
        IAU1980Coefficients::new(0, 1, 2, -2, 1, 4.0, 0.0, -2.0, 0.0),
        IAU1980Coefficients::new(1, 1, 0, 0, 0, -3.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, -1, 0, -1, 0, -3.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, -1, 2, 2, 2, -3.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(0, -1, 2, 2, 2, -3.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(1, -1, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(3, 0, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 2, 0, 2, -3.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, 0, 0, 3.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, 4, 2, -2.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, 0, 2, -2.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 2, -2, 1, -2.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(0, -2, 2, -2, 1, -2.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 0, 0, 1, -2.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, 0, 1, 2.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(3, 0, 0, 0, 0, 2.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 1, 2, 0, 2, 2.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, 1, 2, 2.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, 2, 1, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, 2, 1, -1.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(1, 1, 0, -2, 1, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 0, 2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 2, -2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, -2, 2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, -2, 2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, -2, -2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, 2, -2, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 0, 0, -4, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, -4, 0, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, 4, 2, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 0, 2, -1, 2, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 2, 4, 2, -1.0, 0.0, 1.0, 0.0),
        IAU1980Coefficients::new(2, 0, 2, 2, 2, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, -1, 2, 0, 1, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 0, -2, 0, 1, -1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 0, 4, -2, 2, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 0, 0, 2, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, 1, 2, -2, 2, 1.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(3, 0, 2, -2, 2, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-2, 0, 2, 2, 2, 1.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 0, 0, 2, 1.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(0, 0, -2, 2, 1, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 2, 0, 1, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 4, 0, 2, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(2, 1, 0, -2, 0, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(2, 0, 0, 2, 0, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(2, 0, 2, -2, 1, 1.0, 0.0, -1.0, 0.0),
        IAU1980Coefficients::new(2, 0, -2, 0, 1, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(1, -1, 0, -2, 0, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, 0, 0, 1, 1, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(-1, -1, 0, 2, 1, 1.0, 0.0, 0.0, 0.0),
        IAU1980Coefficients::new(0, 1, 0, 1, 0, 1.0, 0.0, 0.0, 0.0),
    ];

    /// Get a specific coefficient set by index.
    pub fn get_coefficients(index: usize) -> Option<&'static IAU1980Coefficients> {
        Self::COEFFICIENTS.get(index)
    }

    /// Get the number of coefficient sets.
    pub fn coefficient_count() -> usize {
        Self::COEFFICIENTS.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_seconds_before_first() {
        // Before 1972, should return 10
        let jd = 2_440_000.0; // ~1968
        assert_eq!(LeapSeconds::get_leap_seconds_at_jd(jd), 10);
    }

    #[test]
    fn test_leap_seconds_at_epoch() {
        // At GPS epoch (1980-01-06), should be 19 leap seconds
        // JD for 1980-01-06 is approximately 2444244.5
        let jd = 2_444_244.5;
        assert_eq!(LeapSeconds::get_leap_seconds_at_jd(jd), 19);
    }

    #[test]
    fn test_leap_seconds_recent() {
        // After 2017-01-01, should be 37 leap seconds
        let jd = 2_458_000.0; // ~2017
        assert_eq!(LeapSeconds::get_leap_seconds_at_jd(jd), 37);
    }

    #[test]
    fn test_current_leap_seconds() {
        assert_eq!(LeapSeconds::current_leap_seconds(), 37);
    }

    #[test]
    fn test_eop_zero() {
        let eop = EarthOrientationParameters::zero(51544.5);
        assert_eq!(eop.mjd, 51544.5);
        assert_eq!(eop.polar_motion_x, 0.0);
        assert_eq!(eop.polar_motion_y, 0.0);
        assert_eq!(eop.ut1_utc, 0.0);
    }

    #[test]
    fn test_eop_interpolate() {
        let eop1 = EarthOrientationParameters::new(51544.0, 0.1, 0.2, 0.3, 0.001, 0.01, 0.02);
        let eop2 = EarthOrientationParameters::new(51545.0, 0.2, 0.4, 0.6, 0.002, 0.02, 0.04);

        let eop_mid = EarthOrientationParameters::interpolate(&eop1, &eop2, 51544.5);

        assert!((eop_mid.mjd - 51544.5).abs() < 1e-10);
        assert!((eop_mid.polar_motion_x - 0.15).abs() < 1e-10);
        assert!((eop_mid.polar_motion_y - 0.3).abs() < 1e-10);
        assert!((eop_mid.ut1_utc - 0.45).abs() < 1e-10);
    }

    #[test]
    fn test_iau1980_coefficients() {
        // Check first coefficient (dominant term)
        let coef = IAU1980::get_coefficients(0).unwrap();
        assert_eq!(coef.a1, 0);
        assert_eq!(coef.a2, 0);
        assert_eq!(coef.a3, 0);
        assert_eq!(coef.a4, 0);
        assert_eq!(coef.a5, 1);
        assert!((coef.ai - (-171996.0)).abs() < 1e-10);
    }

    #[test]
    fn test_iau1980_coefficient_count() {
        assert_eq!(IAU1980::coefficient_count(), 106);
    }
}
