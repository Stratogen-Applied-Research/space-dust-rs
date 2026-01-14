//! Time systems and transformations for astrodynamics.
//!
//! This module provides representations and conversions between different
//! astronomical time scales:
//!
//! - **UTC** - Coordinated Universal Time (civil time standard)
//! - **TAI** - International Atomic Time (continuous, no leap seconds)
//! - **TT** - Terrestrial Time (modern astronomical standard)
//! - **GPS** - GPS Time
//! - **JulianDate** - Julian Date and Modified Julian Date
//! - **GMST** - Greenwich Mean Sidereal Time
//!
//! ## Time Scale Relationships
//!
//! ```text
//! UTC <---> TAI <---> TT
//!   |         |
//!   |         +--> GPS
//!   |
//!   +--> JD/MJD
//!   |
//!   +--> GMST
//! ```
//!
//! - TAI = UTC + leap_seconds
//! - TT = TAI + 32.184s
//! - GPS = TAI - 19s

use crate::constants::Constants;
use crate::data::LeapSeconds;
use crate::math::poly_eval;
use chrono::{DateTime, NaiveDate, NaiveDateTime, NaiveTime, TimeZone, Utc};
use std::f64::consts::PI;

// ============================================================================
// UTC - Coordinated Universal Time
// ============================================================================

/// Coordinated Universal Time (UTC).
///
/// UTC is the primary civil time standard. It is kept within 0.9 seconds of UT1
/// (astronomical time based on Earth's rotation) by the insertion of leap seconds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UTC {
    /// Seconds since Unix epoch (1970-01-01 00:00:00 UTC)
    unix_seconds: f64,
}

impl UTC {
    /// Create a new UTC time from Unix seconds.
    pub fn new(unix_seconds: f64) -> Self {
        Self { unix_seconds }
    }

    /// Create a UTC time from a chrono DateTime.
    pub fn from_datetime(dt: &DateTime<Utc>) -> Self {
        let secs = dt.timestamp() as f64;
        let nanos = dt.timestamp_subsec_nanos() as f64 / 1_000_000_000.0;
        Self::new(secs + nanos)
    }

    /// Convert to a chrono DateTime.
    pub fn to_datetime(&self) -> DateTime<Utc> {
        let secs = self.unix_seconds.floor() as i64;
        let nanos = ((self.unix_seconds - secs as f64) * 1_000_000_000.0) as u32;
        Utc.timestamp_opt(secs, nanos).unwrap()
    }

    /// Create from year, month, day, hour, minute, second components.
    pub fn from_components(
        year: i32,
        month: u32,
        day: u32,
        hour: u32,
        minute: u32,
        second: f64,
    ) -> Self {
        let whole_second = second.floor() as u32;
        let nanos = ((second - whole_second as f64) * 1_000_000_000.0) as u32;

        let date = NaiveDate::from_ymd_opt(year, month, day).unwrap();
        let time = NaiveTime::from_hms_nano_opt(hour, minute, whole_second, nanos).unwrap();
        let naive_dt = NaiveDateTime::new(date, time);
        let dt = DateTime::<Utc>::from_naive_utc_and_offset(naive_dt, Utc);

        Self::from_datetime(&dt)
    }

    /// Get the Unix seconds value.
    pub fn unix_seconds(&self) -> f64 {
        self.unix_seconds
    }

    /// Convert to Julian Date.
    pub fn to_jd(&self) -> f64 {
        self.unix_seconds / Constants::SECONDS_PER_DAY + Constants::UNIX_EPOCH_JD
    }

    /// Create from Julian Date.
    pub fn from_jd(jd: f64) -> Self {
        let unix_seconds = (jd - Constants::UNIX_EPOCH_JD) * Constants::SECONDS_PER_DAY;
        Self::new(unix_seconds)
    }

    /// Convert to Modified Julian Date.
    pub fn to_mjd(&self) -> f64 {
        self.to_jd() - Constants::MJD_OFFSET
    }

    /// Create from Modified Julian Date.
    pub fn from_mjd(mjd: f64) -> Self {
        Self::from_jd(mjd + Constants::MJD_OFFSET)
    }

    /// Add seconds to the UTC time.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self::new(self.unix_seconds + seconds)
    }

    /// Difference between two UTC times in seconds.
    pub fn diff(&self, other: &UTC) -> f64 {
        self.unix_seconds - other.unix_seconds
    }
}

impl From<DateTime<Utc>> for UTC {
    fn from(dt: DateTime<Utc>) -> Self {
        Self::from_datetime(&dt)
    }
}

impl From<UTC> for DateTime<Utc> {
    fn from(utc: UTC) -> Self {
        utc.to_datetime()
    }
}

impl Default for UTC {
    fn default() -> Self {
        Self::new(0.0)
    }
}

// ============================================================================
// TAI - International Atomic Time
// ============================================================================

/// International Atomic Time (TAI).
///
/// TAI is a continuous time scale based on atomic clocks. Unlike UTC, TAI does not
/// include leap seconds, making it ideal for precise scientific calculations.
///
/// TAI = UTC + leap_seconds
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TAI {
    /// TAI seconds since the epoch corresponding to Unix epoch
    tai_seconds: f64,
}

impl TAI {
    /// Create a new TAI time from TAI seconds.
    pub fn new(tai_seconds: f64) -> Self {
        Self { tai_seconds }
    }

    /// Get the TAI seconds value.
    pub fn seconds(&self) -> f64 {
        self.tai_seconds
    }

    /// Add seconds to the TAI time.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self::new(self.tai_seconds + seconds)
    }

    /// Difference between two TAI times in seconds.
    pub fn diff(&self, other: &TAI) -> f64 {
        self.tai_seconds - other.tai_seconds
    }

    /// Convert to Julian Date (in TAI scale).
    pub fn to_jd(&self) -> f64 {
        self.tai_seconds / Constants::SECONDS_PER_DAY + Constants::UNIX_EPOCH_JD
    }

    /// Create from Julian Date (in TAI scale).
    pub fn from_jd(jd: f64) -> Self {
        let tai_seconds = (jd - Constants::UNIX_EPOCH_JD) * Constants::SECONDS_PER_DAY;
        Self::new(tai_seconds)
    }
}

impl Default for TAI {
    fn default() -> Self {
        Self::new(0.0)
    }
}

// ============================================================================
// TT - Terrestrial Time
// ============================================================================

/// Terrestrial Time (TT).
///
/// TT is the modern astronomical time standard for geocentric ephemerides.
/// It provides a uniform time scale for planetary and lunar ephemeris calculations.
///
/// TT = TAI + 32.184 seconds
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TT {
    /// TT seconds since the epoch corresponding to Unix epoch
    tt_seconds: f64,
}

impl TT {
    /// Create a new TT time from TT seconds.
    pub fn new(tt_seconds: f64) -> Self {
        Self { tt_seconds }
    }

    /// Get the TT seconds value.
    pub fn seconds(&self) -> f64 {
        self.tt_seconds
    }

    /// Add seconds to the TT time.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self::new(self.tt_seconds + seconds)
    }

    /// Difference between two TT times in seconds.
    pub fn diff(&self, other: &TT) -> f64 {
        self.tt_seconds - other.tt_seconds
    }

    /// Convert to Julian Date (in TT scale).
    pub fn to_jd(&self) -> f64 {
        self.tt_seconds / Constants::SECONDS_PER_DAY + Constants::UNIX_EPOCH_JD
    }

    /// Create from Julian Date (in TT scale).
    pub fn from_jd(jd: f64) -> Self {
        let tt_seconds = (jd - Constants::UNIX_EPOCH_JD) * Constants::SECONDS_PER_DAY;
        Self::new(tt_seconds)
    }

    /// Calculate Julian centuries since J2000.0.
    /// This is commonly used for precession/nutation calculations.
    pub fn julian_centuries_j2000(&self) -> f64 {
        (self.to_jd() - Constants::J2000_JD) / Constants::JULIAN_CENTURY
    }
}

impl Default for TT {
    fn default() -> Self {
        Self::new(0.0)
    }
}

// ============================================================================
// GPS - GPS Time
// ============================================================================

/// GPS Time.
///
/// GPS Time is a continuous time scale used by the Global Positioning System.
/// It started at midnight UTC on January 6, 1980 (the GPS epoch) and does not
/// include leap seconds.
///
/// GPS = TAI - 19 seconds
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GPS {
    /// GPS seconds since GPS epoch (1980-01-06 00:00:00 UTC)
    gps_seconds: f64,
}

impl GPS {
    /// Create a new GPS time from GPS seconds (since GPS epoch).
    pub fn new(gps_seconds: f64) -> Self {
        Self { gps_seconds }
    }

    /// Get the GPS seconds value (since GPS epoch).
    pub fn seconds(&self) -> f64 {
        self.gps_seconds
    }

    /// Create from GPS week number and seconds of week.
    pub fn from_week_and_seconds(week: u32, seconds_of_week: f64) -> Self {
        let gps_seconds = week as f64 * Constants::SECONDS_PER_GPS_WEEK + seconds_of_week;
        Self::new(gps_seconds)
    }

    /// Convert to GPS week number and seconds of week.
    pub fn to_week_and_seconds(&self) -> (u32, f64) {
        let week = (self.gps_seconds / Constants::SECONDS_PER_GPS_WEEK).floor() as u32;
        let seconds_of_week = self.gps_seconds - week as f64 * Constants::SECONDS_PER_GPS_WEEK;
        (week, seconds_of_week)
    }

    /// Get GPS week number.
    pub fn week_number(&self) -> u32 {
        (self.gps_seconds / Constants::SECONDS_PER_GPS_WEEK).floor() as u32
    }

    /// Get seconds of week (0 to 604800).
    pub fn seconds_of_week(&self) -> f64 {
        let week = self.week_number();
        self.gps_seconds - week as f64 * Constants::SECONDS_PER_GPS_WEEK
    }

    /// Get day of week (0 = Sunday).
    pub fn day_of_week(&self) -> u32 {
        (self.seconds_of_week() / Constants::SECONDS_PER_DAY).floor() as u32
    }

    /// Add seconds to the GPS time.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self::new(self.gps_seconds + seconds)
    }

    /// Difference between two GPS times in seconds.
    pub fn diff(&self, other: &GPS) -> f64 {
        self.gps_seconds - other.gps_seconds
    }

    /// GPS epoch in Unix seconds.
    pub const GPS_EPOCH_UNIX: i64 = Constants::GPS_EPOCH_UNIX;
}

impl Default for GPS {
    fn default() -> Self {
        Self::new(0.0)
    }
}

// ============================================================================
// JulianDate - Julian Date
// ============================================================================

/// Julian Date (JD).
///
/// The Julian Date is a continuous count of days since the beginning of the
/// Julian Period on January 1, 4713 BC (proleptic Julian calendar) at noon
/// Universal Time.
///
/// Common reference epochs:
/// - J2000.0: JD 2451545.0 (2000-01-01 12:00:00 TT)
/// - Unix epoch: JD 2440587.5 (1970-01-01 00:00:00 UTC)
///
/// Modified Julian Date (MJD) = JD - 2400000.5
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDate {
    /// Julian Date value
    jd: f64,
}

impl JulianDate {
    /// Create a new Julian Date.
    pub fn new(jd: f64) -> Self {
        Self { jd }
    }

    /// Get the Julian Date value.
    pub fn value(&self) -> f64 {
        self.jd
    }

    /// Create from Modified Julian Date.
    pub fn from_mjd(mjd: f64) -> Self {
        Self::new(mjd + Constants::MJD_OFFSET)
    }

    /// Convert to Modified Julian Date.
    pub fn to_mjd(&self) -> f64 {
        self.jd - Constants::MJD_OFFSET
    }

    /// Get J2000.0 epoch as Julian Date.
    pub fn j2000() -> Self {
        Self::new(Constants::J2000_JD)
    }

    /// Calculate Julian centuries since J2000.0.
    pub fn julian_centuries_j2000(&self) -> f64 {
        (self.jd - Constants::J2000_JD) / Constants::JULIAN_CENTURY
    }

    /// Calculate Julian millennia since J2000.0.
    pub fn julian_millennia_j2000(&self) -> f64 {
        (self.jd - Constants::J2000_JD) / (Constants::JULIAN_CENTURY * 10.0)
    }

    /// Add days to the Julian Date.
    pub fn add_days(&self, days: f64) -> Self {
        Self::new(self.jd + days)
    }

    /// Add seconds to the Julian Date.
    pub fn add_seconds(&self, seconds: f64) -> Self {
        Self::new(self.jd + seconds / Constants::SECONDS_PER_DAY)
    }

    /// Difference between two Julian Dates in days.
    pub fn diff_days(&self, other: &JulianDate) -> f64 {
        self.jd - other.jd
    }

    /// Difference between two Julian Dates in seconds.
    pub fn diff_seconds(&self, other: &JulianDate) -> f64 {
        (self.jd - other.jd) * Constants::SECONDS_PER_DAY
    }

    /// Convert to year and fractional day of year.
    /// Returns (year, day_of_year) where day_of_year is 1-based and fractional.
    pub fn to_year_and_day(&self) -> (i32, f64) {
        let jd = self.jd;
        let z = (jd + 0.5).floor() as i64;
        let f = jd + 0.5 - z as f64;

        let a = if z < 2_299_161 {
            z
        } else {
            let alpha = ((z as f64 - 1_867_216.25) / 36524.25).floor() as i64;
            z + 1 + alpha - alpha / 4
        };

        let b = a + 1524;
        let c = ((b as f64 - 122.1) / 365.25).floor() as i64;
        let d = (365.25 * c as f64).floor() as i64;
        let e = ((b - d) as f64 / 30.6001).floor() as i64;

        let _day = (b - d) as f64 - (30.6001 * e as f64).floor() + f;

        let month = if e < 14 { e - 1 } else { e - 13 };
        let year = if month > 2 { c - 4716 } else { c - 4715 };

        // Calculate day of year
        let jan1_jd = Self::from_calendar(year as i32, 1, 1.0);
        let day_of_year = jd - jan1_jd.jd + 1.0;

        (year as i32, day_of_year)
    }

    /// Create Julian Date from calendar date.
    pub fn from_calendar(year: i32, month: u32, day: f64) -> Self {
        let (y, m) = if month <= 2 {
            (year - 1, month + 12)
        } else {
            (year, month)
        };

        let a = y / 100;
        let b = 2 - a + a / 4;

        let jd = (365.25 * (y + 4716) as f64).floor()
            + (30.6001 * (m + 1) as f64).floor()
            + day
            + b as f64
            - 1524.5;

        Self::new(jd)
    }
}

impl Default for JulianDate {
    fn default() -> Self {
        Self::j2000()
    }
}

// ============================================================================
// GMST - Greenwich Mean Sidereal Time
// ============================================================================

/// Greenwich Mean Sidereal Time (GMST).
///
/// GMST is the hour angle of the mean vernal equinox measured from the Greenwich
/// meridian. It represents the Earth's rotation angle and is essential for
/// converting between Earth-fixed and inertial reference frames.
///
/// GMST differs from Greenwich Apparent Sidereal Time (GAST) by the equation
/// of the equinoxes (nutation in right ascension).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GMST {
    /// GMST angle in radians
    radians: f64,
}

impl GMST {
    /// Create a new GMST from radians.
    pub fn new(radians: f64) -> Self {
        Self { radians }
    }

    /// Create from degrees.
    pub fn from_degrees(degrees: f64) -> Self {
        Self::new(degrees * Constants::DEG_TO_RAD)
    }

    /// Create from hours.
    pub fn from_hours(hours: f64) -> Self {
        Self::new(hours * PI / 12.0)
    }

    /// Get the GMST angle in radians.
    pub fn to_radians(&self) -> f64 {
        self.radians
    }

    /// Get the GMST angle in degrees.
    pub fn to_degrees(&self) -> f64 {
        self.radians * Constants::RAD_TO_DEG
    }

    /// Get the GMST angle in hours.
    pub fn to_hours(&self) -> f64 {
        self.radians * 12.0 / PI
    }

    /// Normalize to range [0, 2π).
    pub fn normalize(&self) -> Self {
        let mut normalized = self.radians % Constants::TWO_PI;
        if normalized < 0.0 {
            normalized += Constants::TWO_PI;
        }
        Self::new(normalized)
    }
}

impl Default for GMST {
    fn default() -> Self {
        Self::new(0.0)
    }
}

// ============================================================================
// TimeTransforms - Time System Conversions
// ============================================================================

/// Time system transformations.
///
/// Provides conversions between different astronomical time scales.
pub struct TimeTransforms;

impl TimeTransforms {
    // GMST polynomial coefficients (IAU 1982)
    const GMST_POLY: [f64; 4] = [
        -6.2e-6,        // T³ coefficient (seconds per century³)
        0.093104,       // T² coefficient (seconds per century²)
        8640184.812866, // T coefficient (seconds per century)
        24110.54841,    // constant (seconds)
    ];

    // ========================================================================
    // UTC <-> TAI conversions
    // ========================================================================

    /// Convert UTC to TAI.
    ///
    /// TAI = UTC + leap_seconds
    pub fn utc_to_tai(utc: &UTC) -> TAI {
        let leap_seconds = LeapSeconds::get_leap_seconds_at_jd(utc.to_jd());
        TAI::new(utc.unix_seconds() + leap_seconds as f64)
    }

    /// Convert TAI to UTC.
    ///
    /// UTC = TAI - leap_seconds
    ///
    /// Note: This requires knowing the leap seconds at the TAI time, which
    /// technically requires iterating since leap seconds are defined in UTC.
    /// We use an approximation that works for normal use cases.
    pub fn tai_to_utc(tai: &TAI) -> UTC {
        // First approximation: use current TAI time to estimate UTC
        let approx_utc_jd = tai.tai_seconds / Constants::SECONDS_PER_DAY + Constants::UNIX_EPOCH_JD;
        let leap_seconds = LeapSeconds::get_leap_seconds_at_jd(approx_utc_jd);
        UTC::new(tai.tai_seconds - leap_seconds as f64)
    }

    // ========================================================================
    // TAI <-> TT conversions
    // ========================================================================

    /// Convert TAI to TT.
    ///
    /// TT = TAI + 32.184 seconds
    pub fn tai_to_tt(tai: &TAI) -> TT {
        TT::new(tai.tai_seconds + Constants::TT_TAI_OFFSET)
    }

    /// Convert TT to TAI.
    ///
    /// TAI = TT - 32.184 seconds
    pub fn tt_to_tai(tt: &TT) -> TAI {
        TAI::new(tt.tt_seconds - Constants::TT_TAI_OFFSET)
    }

    // ========================================================================
    // TAI <-> GPS conversions
    // ========================================================================

    /// Convert TAI to GPS time.
    ///
    /// GPS = TAI - 19 seconds (offset since GPS epoch)
    pub fn tai_to_gps(tai: &TAI) -> GPS {
        // GPS epoch in TAI = Unix epoch + leap seconds at GPS epoch (19)
        let gps_epoch_tai = Constants::GPS_EPOCH_UNIX as f64 + 19.0;
        GPS::new(tai.tai_seconds - gps_epoch_tai + Constants::GPS_TAI_OFFSET)
    }

    /// Convert GPS time to TAI.
    pub fn gps_to_tai(gps: &GPS) -> TAI {
        let gps_epoch_tai = Constants::GPS_EPOCH_UNIX as f64 + 19.0;
        TAI::new(gps.gps_seconds + gps_epoch_tai - Constants::GPS_TAI_OFFSET)
    }

    // ========================================================================
    // UTC <-> TT convenience conversions
    // ========================================================================

    /// Convert UTC directly to TT.
    pub fn utc_to_tt(utc: &UTC) -> TT {
        let tai = Self::utc_to_tai(utc);
        Self::tai_to_tt(&tai)
    }

    /// Convert TT directly to UTC.
    pub fn tt_to_utc(tt: &TT) -> UTC {
        let tai = Self::tt_to_tai(tt);
        Self::tai_to_utc(&tai)
    }

    // ========================================================================
    // UTC <-> GPS convenience conversions
    // ========================================================================

    /// Convert UTC directly to GPS time.
    pub fn utc_to_gps(utc: &UTC) -> GPS {
        let tai = Self::utc_to_tai(utc);
        Self::tai_to_gps(&tai)
    }

    /// Convert GPS time directly to UTC.
    pub fn gps_to_utc(gps: &GPS) -> UTC {
        let tai = Self::gps_to_tai(gps);
        Self::tai_to_utc(&tai)
    }

    // ========================================================================
    // Julian Date conversions
    // ========================================================================

    /// Convert UTC to Julian Date.
    pub fn utc_to_jd(utc: &UTC) -> JulianDate {
        JulianDate::new(utc.to_jd())
    }

    /// Convert Julian Date to UTC.
    pub fn jd_to_utc(jd: &JulianDate) -> UTC {
        UTC::from_jd(jd.value())
    }

    // ========================================================================
    // GMST calculations
    // ========================================================================

    /// Calculate GMST from UTC time using the IAU 1982 expression.
    pub fn utc_to_gmst(utc: &UTC) -> GMST {
        let jd = utc.to_jd();
        Self::calculate_gmst(jd)
    }

    /// Calculate GMST from Julian Date.
    pub fn jd_to_gmst(jd: &JulianDate) -> GMST {
        Self::calculate_gmst(jd.value())
    }

    /// Internal GMST calculation using IAU 1982 formula.
    fn calculate_gmst(jd: f64) -> GMST {
        // Julian centuries from J2000.0
        let t = (jd - Constants::J2000_JD) / Constants::JULIAN_CENTURY;

        // GMST at 0h UT in seconds (polynomial evaluation)
        let gmst_0h = poly_eval(&Self::GMST_POLY, t);

        // Add the rotation for the fractional day
        // Earth rotates 360.98564736629 degrees per day (sidereal)
        let frac_day = (jd + 0.5).fract();
        let rotation_seconds = frac_day * Constants::SECONDS_PER_DAY * 1.00273790935;

        let gmst_seconds = gmst_0h + rotation_seconds;

        // Convert to radians (24 hours = 2π radians)
        let gmst_radians = gmst_seconds / Constants::SECONDS_PER_DAY * Constants::TWO_PI;

        GMST::new(gmst_radians).normalize()
    }

    // ========================================================================
    // DateTime conversions
    // ========================================================================

    /// Convert a chrono `DateTime<Utc>` to UTC.
    pub fn datetime_to_utc(dt: &DateTime<Utc>) -> UTC {
        UTC::from_datetime(dt)
    }

    /// Convert UTC to a chrono `DateTime<Utc>`.
    pub fn utc_to_datetime(utc: &UTC) -> DateTime<Utc> {
        utc.to_datetime()
    }

    /// Convert a chrono `DateTime<Utc>` directly to TT.
    pub fn datetime_to_tt(dt: &DateTime<Utc>) -> TT {
        let utc = UTC::from_datetime(dt);
        Self::utc_to_tt(&utc)
    }

    /// Convert a chrono `DateTime<Utc>` directly to Julian Date.
    pub fn datetime_to_jd(dt: &DateTime<Utc>) -> JulianDate {
        let utc = UTC::from_datetime(dt);
        Self::utc_to_jd(&utc)
    }

    /// Get Julian centuries since J2000.0 for a DateTime.
    pub fn julian_centuries_j2000(dt: &DateTime<Utc>) -> f64 {
        let tt = Self::datetime_to_tt(dt);
        tt.julian_centuries_j2000()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_utc_from_datetime() {
        let dt = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
        let utc = UTC::from_datetime(&dt);
        let back = utc.to_datetime();
        assert_eq!(dt, back);
    }

    #[test]
    fn test_utc_to_jd() {
        // J2000.0 epoch: 2000-01-01 12:00:00 TT
        // In UTC this is approximately 2000-01-01 11:58:55.816 UTC
        // For simplicity, test a known date
        let utc = UTC::from_components(2000, 1, 1, 12, 0, 0.0);
        let jd = utc.to_jd();
        // The JD should be close to J2000 (2451545.0)
        assert!((jd - 2451545.0).abs() < 1.0);
    }

    #[test]
    fn test_utc_to_mjd() {
        let utc = UTC::from_jd(2451545.0);
        let mjd = utc.to_mjd();
        assert!((mjd - 51544.5).abs() < EPSILON);
    }

    #[test]
    fn test_tai_conversion() {
        let utc = UTC::from_components(2020, 1, 1, 0, 0, 0.0);
        let tai = TimeTransforms::utc_to_tai(&utc);
        let utc_back = TimeTransforms::tai_to_utc(&tai);

        // Should be approximately equal (within leap second precision)
        assert!((utc.unix_seconds() - utc_back.unix_seconds()).abs() < 1.0);
    }

    #[test]
    fn test_tt_conversion() {
        let tai = TAI::new(1000000.0);
        let tt = TimeTransforms::tai_to_tt(&tai);
        let tai_back = TimeTransforms::tt_to_tai(&tt);

        assert!((tai.seconds() - tai_back.seconds()).abs() < EPSILON);
    }

    #[test]
    fn test_gps_week_conversion() {
        let gps = GPS::from_week_and_seconds(2000, 345600.0);
        let (week, sow) = gps.to_week_and_seconds();
        assert_eq!(week, 2000);
        assert!((sow - 345600.0).abs() < EPSILON);
    }

    #[test]
    fn test_julian_date_mjd() {
        let jd = JulianDate::new(2451545.0);
        let mjd = jd.to_mjd();
        let jd_back = JulianDate::from_mjd(mjd);
        assert!((jd.value() - jd_back.value()).abs() < EPSILON);
    }

    #[test]
    fn test_julian_centuries() {
        let jd = JulianDate::j2000();
        let t = jd.julian_centuries_j2000();
        assert!(t.abs() < EPSILON);

        // One century later
        let jd2 = JulianDate::new(Constants::J2000_JD + Constants::JULIAN_CENTURY);
        let t2 = jd2.julian_centuries_j2000();
        assert!((t2 - 1.0).abs() < EPSILON);
    }

    #[test]
    fn test_gmst_normalize() {
        let gmst = GMST::new(3.0 * PI);
        let normalized = gmst.normalize();
        assert!((normalized.to_radians() - PI).abs() < EPSILON);

        let gmst_neg = GMST::new(-PI / 2.0);
        let normalized_neg = gmst_neg.normalize();
        assert!((normalized_neg.to_radians() - 3.0 * PI / 2.0).abs() < EPSILON);
    }

    #[test]
    fn test_gmst_conversions() {
        let gmst = GMST::from_hours(12.0);
        assert!((gmst.to_radians() - PI).abs() < EPSILON);
        assert!((gmst.to_degrees() - 180.0).abs() < EPSILON);
        assert!((gmst.to_hours() - 12.0).abs() < EPSILON);
    }

    #[test]
    fn test_julian_date_from_calendar() {
        // Test J2000.0: 2000-01-01.5 (noon)
        let jd = JulianDate::from_calendar(2000, 1, 1.5);
        assert!((jd.value() - 2451545.0).abs() < EPSILON);
    }

    #[test]
    fn test_utc_add_seconds() {
        let utc1 = UTC::new(1000.0);
        let utc2 = utc1.add_seconds(500.0);
        assert!((utc2.unix_seconds() - 1500.0).abs() < EPSILON);
    }

    #[test]
    fn test_utc_diff() {
        let utc1 = UTC::new(1500.0);
        let utc2 = UTC::new(1000.0);
        assert!((utc1.diff(&utc2) - 500.0).abs() < EPSILON);
    }
}
