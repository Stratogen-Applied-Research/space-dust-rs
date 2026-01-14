//! Two-Line Element Set (TLE) parsing and SGP4 propagation.
//!
//! This module provides:
//! - TLE parsing from standard two-line format
//! - SGP4 satellite propagation using the `sgp4` crate
//! - Conversion to state vectors in TEME frame
//!
//! ## Example
//!
//! ```rust,no_run
//! use space_dust::tle::Tle;
//! use chrono::Utc;
//!
//! let line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9002";
//! let line2 = "2 25544  51.6400 208.1200 0001234  85.0000 275.0000 15.48919100123456";
//!
//! let tle = Tle::parse(line1, line2).unwrap();
//! println!("Satellite: {}", tle.catalog_number());
//! println!("Epoch: {:?}", tle.epoch());
//!
//! // Propagate to current time
//! let state = tle.propagate(&Utc::now()).unwrap();
//! println!("Position: {:?}", state.position);
//! ```

use crate::math::Vector3;
use crate::state::TEMEState;
use chrono::{DateTime, TimeZone, Utc};
use sgp4::{Constants as Sgp4Constants, Elements, Prediction};
use std::error::Error;
use std::fmt;

// ============================================================================
// Error Types
// ============================================================================

/// Errors that can occur during TLE parsing or propagation.
#[derive(Debug)]
pub enum TleError {
    /// TLE line has incorrect length
    InvalidLineLength {
        line: u8,
        expected: usize,
        actual: usize,
    },
    /// Failed to parse a field in the TLE
    ParseError { field: String, message: String },
    /// SGP4 propagation error
    PropagationError { message: String },
    /// Elements are invalid for SGP4
    InvalidElements { message: String },
}

impl fmt::Display for TleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TleError::InvalidLineLength {
                line,
                expected,
                actual,
            } => {
                write!(
                    f,
                    "TLE line {} has invalid length: expected {}, got {}",
                    line, expected, actual
                )
            }
            TleError::ParseError { field, message } => {
                write!(f, "Failed to parse TLE field '{}': {}", field, message)
            }
            TleError::PropagationError { message } => {
                write!(f, "SGP4 propagation error: {}", message)
            }
            TleError::InvalidElements { message } => {
                write!(f, "Invalid TLE elements: {}", message)
            }
        }
    }
}

impl Error for TleError {}

/// Result type for TLE operations.
pub type TleResult<T> = Result<T, TleError>;

// ============================================================================
// TLE Structure
// ============================================================================

/// Two-Line Element Set (TLE) data structure.
///
/// A TLE is a standardized format for describing the orbit of an Earth-orbiting
/// object. TLEs are used with the SGP4/SDP4 propagators to predict satellite
/// positions at future times.
#[derive(Debug, Clone)]
pub struct Tle {
    /// Raw TLE line 1
    line1: String,
    /// Raw TLE line 2
    line2: String,
    /// NORAD catalog number
    catalog_number: String,
    /// Security classification (U = Unclassified)
    classification: char,
    /// International designator
    international_designator: String,
    /// TLE epoch as DateTime<Utc>
    epoch: DateTime<Utc>,
    /// First derivative of mean motion (rev/day²)
    mean_motion_dot: f64,
    /// Second derivative of mean motion (rev/day³)
    mean_motion_double_dot: f64,
    /// BSTAR drag coefficient (1/Earth radii)
    bstar: f64,
    /// Ephemeris type (usually 0)
    ephemeris_type: u8,
    /// Element set number
    element_set_number: u32,
    /// Orbital inclination (degrees)
    inclination_deg: f64,
    /// Right Ascension of Ascending Node (degrees)
    raan_deg: f64,
    /// Orbital eccentricity
    eccentricity: f64,
    /// Argument of perigee (degrees)
    arg_perigee_deg: f64,
    /// Mean anomaly at epoch (degrees)
    mean_anomaly_deg: f64,
    /// Mean motion (revolutions per day)
    mean_motion: f64,
    /// Revolution number at epoch
    rev_number: u32,
    /// Parsed SGP4 elements (kept for potential future use/debugging)
    #[allow(dead_code)]
    sgp4_elements: Elements,
    /// SGP4 propagator constants
    sgp4_constants: Sgp4Constants,
}

impl Tle {
    /// Parse a TLE from two lines.
    ///
    /// # Arguments
    /// * `line1` - First line of the TLE (69 characters)
    /// * `line2` - Second line of the TLE (69 characters)
    ///
    /// # Returns
    /// * `Ok(Tle)` - Successfully parsed TLE
    /// * `Err(TleError)` - Parsing failed
    pub fn parse(line1: &str, line2: &str) -> TleResult<Self> {
        // Validate line lengths
        if line1.len() != 69 {
            return Err(TleError::InvalidLineLength {
                line: 1,
                expected: 69,
                actual: line1.len(),
            });
        }
        if line2.len() != 69 {
            return Err(TleError::InvalidLineLength {
                line: 2,
                expected: 69,
                actual: line2.len(),
            });
        }

        // Parse epoch
        let epoch_year = parse_field::<i32>(&line1[18..20], "epoch_year")?;
        let epoch_day = parse_field::<f64>(&line1[20..32], "epoch_day")?;

        let full_year = if epoch_year > 56 {
            1900 + epoch_year
        } else {
            2000 + epoch_year
        };

        let epoch = epoch_from_year_and_day(full_year, epoch_day)?;

        // Parse mean motion derivatives
        let mean_motion_dot = parse_mean_motion_dot(&line1[33..43])?;
        let mean_motion_double_dot = parse_exponential(&line1[44..52])?;
        let bstar = parse_exponential(&line1[53..61])?;

        // Parse other line 1 fields
        let catalog_number = line1[2..7].trim().to_string();
        let classification = line1.chars().nth(7).unwrap_or('U');
        let international_designator = line1[9..17].trim().to_string();
        let ephemeris_type = parse_field::<u8>(&line1[62..63], "ephemeris_type").unwrap_or(0);
        let element_set_number = parse_field::<u32>(line1[64..68].trim(), "element_set_number")?;

        // Parse line 2 fields
        let inclination_deg = parse_field::<f64>(line2[8..16].trim(), "inclination")?;
        let raan_deg = parse_field::<f64>(line2[17..25].trim(), "raan")?;
        let eccentricity =
            parse_field::<f64>(&format!("0.{}", &line2[26..33].trim()), "eccentricity")?;
        let arg_perigee_deg = parse_field::<f64>(line2[34..42].trim(), "arg_perigee")?;
        let mean_anomaly_deg = parse_field::<f64>(line2[43..51].trim(), "mean_anomaly")?;
        let mean_motion = parse_field::<f64>(line2[52..63].trim(), "mean_motion")?;
        let rev_number = parse_field::<u32>(line2[63..68].trim(), "rev_number")?;

        // Parse using sgp4 crate for propagation
        let sgp4_elements = sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes())
            .map_err(|e| TleError::ParseError {
                field: "sgp4_elements".to_string(),
                message: format!("{:?}", e),
            })?;

        // Create SGP4 constants
        let sgp4_constants = Sgp4Constants::from_elements(&sgp4_elements).map_err(|e| {
            TleError::InvalidElements {
                message: format!("{:?}", e),
            }
        })?;

        Ok(Self {
            line1: line1.to_string(),
            line2: line2.to_string(),
            catalog_number,
            classification,
            international_designator,
            epoch,
            mean_motion_dot,
            mean_motion_double_dot,
            bstar,
            ephemeris_type,
            element_set_number,
            inclination_deg,
            raan_deg,
            eccentricity,
            arg_perigee_deg,
            mean_anomaly_deg,
            mean_motion,
            rev_number,
            sgp4_elements,
            sgp4_constants,
        })
    }

    /// Parse a TLE from three lines (name + two TLE lines).
    ///
    /// # Arguments
    /// * `name` - Satellite name (ignored)
    /// * `line1` - First line of the TLE
    /// * `line2` - Second line of the TLE
    pub fn parse_3le(_name: &str, line1: &str, line2: &str) -> TleResult<Self> {
        Self::parse(line1, line2)
    }

    // ========================================================================
    // Accessors
    // ========================================================================

    /// Get the raw TLE line 1.
    pub fn line1(&self) -> &str {
        &self.line1
    }

    /// Get the raw TLE line 2.
    pub fn line2(&self) -> &str {
        &self.line2
    }

    /// Get the NORAD catalog number.
    pub fn catalog_number(&self) -> &str {
        &self.catalog_number
    }

    /// Get the security classification.
    pub fn classification(&self) -> char {
        self.classification
    }

    /// Get the international designator.
    pub fn international_designator(&self) -> &str {
        &self.international_designator
    }

    /// Get the TLE epoch.
    pub fn epoch(&self) -> DateTime<Utc> {
        self.epoch
    }

    /// Get the first derivative of mean motion (rev/day²).
    pub fn mean_motion_dot(&self) -> f64 {
        self.mean_motion_dot
    }

    /// Get the second derivative of mean motion (rev/day³).
    pub fn mean_motion_double_dot(&self) -> f64 {
        self.mean_motion_double_dot
    }

    /// Get the BSTAR drag coefficient.
    pub fn bstar(&self) -> f64 {
        self.bstar
    }

    /// Get the ephemeris type.
    pub fn ephemeris_type(&self) -> u8 {
        self.ephemeris_type
    }

    /// Get the element set number.
    pub fn element_set_number(&self) -> u32 {
        self.element_set_number
    }

    /// Get the orbital inclination in degrees.
    pub fn inclination_deg(&self) -> f64 {
        self.inclination_deg
    }

    /// Get the Right Ascension of Ascending Node in degrees.
    pub fn raan_deg(&self) -> f64 {
        self.raan_deg
    }

    /// Get the orbital eccentricity.
    pub fn eccentricity(&self) -> f64 {
        self.eccentricity
    }

    /// Get the argument of perigee in degrees.
    pub fn arg_perigee_deg(&self) -> f64 {
        self.arg_perigee_deg
    }

    /// Get the mean anomaly at epoch in degrees.
    pub fn mean_anomaly_deg(&self) -> f64 {
        self.mean_anomaly_deg
    }

    /// Get the mean motion in revolutions per day.
    pub fn mean_motion(&self) -> f64 {
        self.mean_motion
    }

    /// Get the revolution number at epoch.
    pub fn rev_number(&self) -> u32 {
        self.rev_number
    }

    /// Get the orbital period in minutes.
    pub fn period_minutes(&self) -> f64 {
        1440.0 / self.mean_motion
    }

    /// Get the orbital period in seconds.
    pub fn period_seconds(&self) -> f64 {
        86400.0 / self.mean_motion
    }

    // ========================================================================
    // Propagation
    // ========================================================================

    /// Propagate the TLE to a specific epoch using SGP4.
    ///
    /// Returns the satellite state in the TEME (True Equator Mean Equinox) frame.
    ///
    /// # Arguments
    /// * `epoch` - Target epoch for propagation
    ///
    /// # Returns
    /// * `Ok(TEMEState)` - Position and velocity in TEME frame
    /// * `Err(TleError)` - Propagation failed
    pub fn propagate(&self, epoch: &DateTime<Utc>) -> TleResult<TEMEState> {
        // Calculate minutes since TLE epoch
        let minutes_since_epoch = self.minutes_since_epoch(epoch);

        // Propagate using SGP4
        let prediction = self
            .sgp4_constants
            .propagate(sgp4::MinutesSinceEpoch(minutes_since_epoch))
            .map_err(|e| TleError::PropagationError {
                message: format!("{:?}", e),
            })?;

        // Convert to TEMEState
        Ok(prediction_to_teme_state(&prediction, *epoch))
    }

    /// Propagate the TLE by a number of minutes from the TLE epoch.
    ///
    /// # Arguments
    /// * `minutes` - Minutes since TLE epoch (can be negative)
    pub fn propagate_minutes(&self, minutes: f64) -> TleResult<TEMEState> {
        let prediction = self
            .sgp4_constants
            .propagate(sgp4::MinutesSinceEpoch(minutes))
            .map_err(|e| TleError::PropagationError {
                message: format!("{:?}", e),
            })?;

        // Calculate the target epoch
        let target_epoch =
            self.epoch + chrono::Duration::milliseconds((minutes * 60.0 * 1000.0) as i64);

        Ok(prediction_to_teme_state(&prediction, target_epoch))
    }

    /// Get position and velocity at a specific epoch.
    ///
    /// Returns ((x, y, z), (vx, vy, vz)) in km and km/s.
    pub fn get_rv_at_time(
        &self,
        epoch: &DateTime<Utc>,
    ) -> TleResult<((f64, f64, f64), (f64, f64, f64))> {
        let state = self.propagate(epoch)?;
        Ok((state.position.to_tuple(), state.velocity.to_tuple()))
    }

    /// Calculate minutes since TLE epoch.
    pub fn minutes_since_epoch(&self, target: &DateTime<Utc>) -> f64 {
        let duration = target.signed_duration_since(self.epoch);
        duration.num_milliseconds() as f64 / 60_000.0
    }

    /// Check if propagation to a given epoch is likely to be accurate.
    ///
    /// TLEs are typically valid for a few days before/after their epoch.
    /// This function returns a warning level based on the time difference.
    ///
    /// # Returns
    /// * `0` - Within ±1 day (good accuracy)
    /// * `1` - Within ±3 days (acceptable)
    /// * `2` - Within ±7 days (degraded)
    /// * `3` - Beyond ±7 days (poor accuracy)
    pub fn accuracy_warning(&self, epoch: &DateTime<Utc>) -> u8 {
        let days = self.minutes_since_epoch(epoch).abs() / 1440.0;
        if days <= 1.0 {
            0
        } else if days <= 3.0 {
            1
        } else if days <= 7.0 {
            2
        } else {
            3
        }
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Parse a field from a TLE string.
fn parse_field<T: std::str::FromStr>(s: &str, field_name: &str) -> TleResult<T> {
    s.trim().parse().map_err(|_| TleError::ParseError {
        field: field_name.to_string(),
        message: format!("Could not parse '{}'", s.trim()),
    })
}

/// Parse the mean motion derivative field (special format).
fn parse_mean_motion_dot(s: &str) -> TleResult<f64> {
    let trimmed = s.trim();
    if trimmed.is_empty() || trimmed == "." {
        return Ok(0.0);
    }

    // Handle sign
    let (sign, value_str) = if trimmed.starts_with('-') {
        (-1.0, &trimmed[1..])
    } else if trimmed.starts_with('+') {
        (1.0, &trimmed[1..])
    } else {
        (1.0, trimmed)
    };

    // Parse as decimal
    let value: f64 = if value_str.starts_with('.') {
        format!("0{}", value_str)
    } else {
        value_str.to_string()
    }
    .trim()
    .parse()
    .map_err(|_| TleError::ParseError {
        field: "mean_motion_dot".to_string(),
        message: format!("Could not parse '{}'", s),
    })?;

    Ok(sign * value)
}

/// Parse an exponential notation field (like BSTAR).
fn parse_exponential(s: &str) -> TleResult<f64> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Ok(0.0);
    }

    // Handle sign
    let (sign, rest) = if trimmed.starts_with('-') {
        (-1.0, &trimmed[1..])
    } else if trimmed.starts_with('+') {
        (1.0, &trimmed[1..])
    } else if trimmed.starts_with(' ') {
        (1.0, &trimmed[1..])
    } else {
        (1.0, trimmed)
    };

    // TLE format: mantissa followed by exponent (e.g., "12345-4" means 0.12345e-4)
    if rest.len() < 2 {
        return Ok(0.0);
    }

    // Find the exponent (last 2 characters including sign)
    let exp_start = rest.len() - 2;
    let mantissa_str = &rest[..exp_start];
    let exp_str = &rest[exp_start..];

    // Parse mantissa as 0.xxxxx
    let mantissa: f64 = format!("0.{}", mantissa_str.trim()).parse().unwrap_or(0.0);

    // Parse exponent
    let exponent: i32 = exp_str.trim().parse().unwrap_or(0);

    Ok(sign * mantissa * 10.0_f64.powi(exponent))
}

/// Convert year and day of year to DateTime<Utc>.
fn epoch_from_year_and_day(year: i32, day_of_year: f64) -> TleResult<DateTime<Utc>> {
    // Day 1 is January 1st, so we subtract 1 to get days to add
    let days_to_add = (day_of_year - 1.0).floor() as i64;
    let fraction = day_of_year - day_of_year.floor();
    let microseconds = (fraction * 86_400_000_000.0) as i64;

    let start_of_year = Utc
        .with_ymd_and_hms(year, 1, 1, 0, 0, 0)
        .single()
        .ok_or_else(|| TleError::ParseError {
            field: "epoch".to_string(),
            message: format!("Invalid year: {}", year),
        })?;

    let epoch = start_of_year
        + chrono::Duration::days(days_to_add)
        + chrono::Duration::microseconds(microseconds);

    Ok(epoch)
}

/// Convert SGP4 prediction to TEMEState.
fn prediction_to_teme_state(prediction: &Prediction, epoch: DateTime<Utc>) -> TEMEState {
    let position = Vector3::new(
        prediction.position[0],
        prediction.position[1],
        prediction.position[2],
    );
    let velocity = Vector3::new(
        prediction.velocity[0],
        prediction.velocity[1],
        prediction.velocity[2],
    );
    TEMEState::new(epoch, position, velocity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Datelike;

    // Real ISS TLE from Celestrak with correct checksums
    const ISS_LINE1: &str = "1 25544U 98067A   26014.62805721  .00006818  00000+0  13044-3 0  9991";
    const ISS_LINE2: &str = "2 25544  51.6333 339.6562 0007763  17.9854 342.1408 15.49289811547943";

    #[test]
    fn test_tle_parse() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2);
        assert!(tle.is_ok());

        let tle = tle.unwrap();
        assert_eq!(tle.catalog_number(), "25544");
        assert_eq!(tle.classification(), 'U');
        assert!(tle.inclination_deg() > 51.0 && tle.inclination_deg() < 52.0);
        assert!(tle.eccentricity() < 0.01); // Nearly circular
    }

    #[test]
    fn test_tle_invalid_length() {
        let short_line = "1 25544U";
        let result = Tle::parse(short_line, ISS_LINE2);
        assert!(result.is_err());

        if let Err(TleError::InvalidLineLength { line, .. }) = result {
            assert_eq!(line, 1);
        } else {
            panic!("Expected InvalidLineLength error");
        }
    }

    #[test]
    fn test_tle_epoch() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        let epoch = tle.epoch();

        // Should be in 2026
        assert_eq!(epoch.year(), 2026);
        assert_eq!(epoch.month(), 1);
    }

    #[test]
    fn test_tle_period() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();

        // ISS period is about 92 minutes
        let period = tle.period_minutes();
        assert!(period > 90.0 && period < 95.0);
    }

    #[test]
    fn test_tle_propagate() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();

        // Propagate to TLE epoch (0 minutes)
        let state = tle.propagate_minutes(0.0);
        assert!(state.is_ok());

        let state = state.unwrap();

        // Position should be roughly in LEO (300-450 km altitude, so ~6700-6800 km from center)
        let radius = state.position.magnitude();
        assert!(radius > 6500.0 && radius < 7000.0);

        // Velocity should be roughly 7.5 km/s for LEO
        let speed = state.velocity.magnitude();
        assert!(speed > 7.0 && speed < 8.0);
    }

    #[test]
    fn test_tle_propagate_forward() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();

        // Propagate forward 1 orbit (about 92 minutes)
        let state0 = tle.propagate_minutes(0.0).unwrap();
        let state1 = tle.propagate_minutes(92.0).unwrap();

        // After one orbit, should be approximately back to same position
        // (not exact due to J2 perturbations and precession)
        let diff = state0.position - state1.position;
        let distance = diff.magnitude();

        // Should be within ~200 km after one orbit
        assert!(distance < 500.0);
    }

    #[test]
    fn test_tle_accuracy_warning() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();

        // At TLE epoch
        assert_eq!(tle.accuracy_warning(&tle.epoch()), 0);

        // 2 days later
        let later = tle.epoch() + chrono::Duration::days(2);
        assert_eq!(tle.accuracy_warning(&later), 1);

        // 5 days later
        let much_later = tle.epoch() + chrono::Duration::days(5);
        assert_eq!(tle.accuracy_warning(&much_later), 2);

        // 10 days later
        let way_later = tle.epoch() + chrono::Duration::days(10);
        assert_eq!(tle.accuracy_warning(&way_later), 3);
    }

    #[test]
    fn test_parse_exponential() {
        // Test BSTAR-style exponential parsing
        assert!((parse_exponential(" 10270-3").unwrap() - 0.0001027).abs() < 1e-10);
        assert!((parse_exponential("-10270-3").unwrap() - (-0.0001027)).abs() < 1e-10);
        assert!((parse_exponential(" 00000+0").unwrap()).abs() < 1e-15);
    }

    #[test]
    fn test_get_rv_at_time() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        let result = tle.get_rv_at_time(&tle.epoch());

        assert!(result.is_ok());
        let ((x, y, z), (vx, vy, vz)) = result.unwrap();

        // Position magnitude check
        let r = (x * x + y * y + z * z).sqrt();
        assert!(r > 6500.0 && r < 7000.0);

        // Velocity magnitude check
        let v = (vx * vx + vy * vy + vz * vz).sqrt();
        assert!(v > 7.0 && v < 8.0);
    }
}
