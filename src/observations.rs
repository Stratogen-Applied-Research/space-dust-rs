//! Angular observation calculations for ground-based tracking.
//!
//! This module provides:
//! - **AzEl** - Azimuth/Elevation observations
//! - **RaDec** - Right Ascension/Declination observations
//! - **Observations** - Functions for computing observations from observer and target states

use crate::constants::Constants;
use crate::math::Vector3;
use crate::state::{ECIState, GeodeticState, StateTransforms};
use chrono::{DateTime, Utc};
use std::f64::consts::PI;

// ============================================================================
// AzEl - Azimuth/Elevation Observation
// ============================================================================

/// Azimuth and Elevation observation.
///
/// Represents an angular observation in the local horizontal coordinate system:
/// - Azimuth (Az): Angle measured clockwise from true North, in radians [0, 2π)
/// - Elevation (El): Angle measured above the local horizon, in radians [-π/2, +π/2]
#[derive(Debug, Clone, Copy)]
pub struct AzEl {
    /// UTC epoch of the observation
    pub epoch: DateTime<Utc>,
    /// Azimuth in radians [0, 2π), measured clockwise from North
    pub azimuth: f64,
    /// Elevation in radians [-π/2, +π/2], measured above horizon
    pub elevation: f64,
    /// Slant range in kilometers (optional)
    pub range: Option<f64>,
    /// Range rate in km/s (optional)
    pub range_rate: Option<f64>,
    /// Azimuth rate in rad/s (optional)
    pub azimuth_rate: Option<f64>,
    /// Elevation rate in rad/s (optional)
    pub elevation_rate: Option<f64>,
}

impl AzEl {
    /// Create a new Az/El observation.
    ///
    /// # Arguments
    /// * `epoch` - UTC epoch of the observation
    /// * `azimuth` - Azimuth in radians [0, 2π), measured clockwise from North
    /// * `elevation` - Elevation in radians [-π/2, +π/2], measured above horizon
    pub fn new(epoch: DateTime<Utc>, azimuth: f64, elevation: f64) -> Self {
        Self {
            epoch,
            azimuth: normalize_azimuth(azimuth),
            elevation,
            range: None,
            range_rate: None,
            azimuth_rate: None,
            elevation_rate: None,
        }
    }

    /// Create Az/El from degrees.
    pub fn from_degrees(epoch: DateTime<Utc>, az_deg: f64, el_deg: f64) -> Self {
        Self::new(
            epoch,
            az_deg * Constants::DEG_TO_RAD,
            el_deg * Constants::DEG_TO_RAD,
        )
    }

    /// Set the range.
    pub fn with_range(mut self, range: f64) -> Self {
        self.range = Some(range);
        self
    }

    /// Set the range rate.
    pub fn with_range_rate(mut self, range_rate: f64) -> Self {
        self.range_rate = Some(range_rate);
        self
    }

    /// Set the azimuth rate.
    pub fn with_azimuth_rate(mut self, azimuth_rate: f64) -> Self {
        self.azimuth_rate = Some(azimuth_rate);
        self
    }

    /// Set the elevation rate.
    pub fn with_elevation_rate(mut self, elevation_rate: f64) -> Self {
        self.elevation_rate = Some(elevation_rate);
        self
    }

    /// Get azimuth in degrees.
    pub fn azimuth_deg(&self) -> f64 {
        self.azimuth * Constants::RAD_TO_DEG
    }

    /// Get elevation in degrees.
    pub fn elevation_deg(&self) -> f64 {
        self.elevation * Constants::RAD_TO_DEG
    }

    /// Convert to degrees, returning (azimuth_deg, elevation_deg).
    pub fn to_degrees(&self) -> (f64, f64) {
        (self.azimuth_deg(), self.elevation_deg())
    }

    /// Check if the target is above the horizon.
    pub fn above_horizon(&self) -> bool {
        self.elevation > 0.0
    }

    /// Check if the target is above a minimum elevation.
    ///
    /// # Arguments
    /// * `min_elevation_deg` - Minimum elevation in degrees
    pub fn above_elevation(&self, min_elevation_deg: f64) -> bool {
        self.elevation >= min_elevation_deg * Constants::DEG_TO_RAD
    }

    /// Get compass direction from azimuth.
    pub fn compass_direction(&self) -> &'static str {
        let az_deg = self.azimuth_deg();
        if az_deg < 22.5 || az_deg >= 337.5 {
            "N"
        } else if az_deg < 67.5 {
            "NE"
        } else if az_deg < 112.5 {
            "E"
        } else if az_deg < 157.5 {
            "SE"
        } else if az_deg < 202.5 {
            "S"
        } else if az_deg < 247.5 {
            "SW"
        } else if az_deg < 292.5 {
            "W"
        } else {
            "NW"
        }
    }

    /// Format as a human-readable string.
    pub fn to_string(&self) -> String {
        format!(
            "Az: {:.2}° ({}), El: {:.2}°",
            self.azimuth_deg(),
            self.compass_direction(),
            self.elevation_deg()
        )
    }
}

// ============================================================================
// RaDec - Right Ascension/Declination Observation
// ============================================================================

/// Right Ascension and Declination observation.
///
/// Represents an angular observation in the equatorial coordinate system:
/// - Right Ascension (RA): Angle measured eastward along the celestial equator
///   from the vernal equinox, in radians [0, 2π)
/// - Declination (Dec): Angle measured north (+) or south (-) from the
///   celestial equator, in radians [-π/2, +π/2]
#[derive(Debug, Clone, Copy)]
pub struct RaDec {
    /// UTC epoch of the observation
    pub epoch: DateTime<Utc>,
    /// Right ascension in radians [0, 2π)
    pub right_ascension: f64,
    /// Declination in radians [-π/2, +π/2]
    pub declination: f64,
    /// Distance in kilometers (optional)
    pub range: Option<f64>,
    /// Range rate in km/s (optional)
    pub range_rate: Option<f64>,
    /// RA rate in rad/s (optional)
    pub right_ascension_rate: Option<f64>,
    /// Dec rate in rad/s (optional)
    pub declination_rate: Option<f64>,
}

impl RaDec {
    /// Create a new RA/Dec observation.
    ///
    /// # Arguments
    /// * `epoch` - UTC epoch of the observation
    /// * `right_ascension` - Right ascension in radians [0, 2π)
    /// * `declination` - Declination in radians [-π/2, +π/2]
    pub fn new(epoch: DateTime<Utc>, right_ascension: f64, declination: f64) -> Self {
        Self {
            epoch,
            right_ascension: normalize_ra(right_ascension),
            declination,
            range: None,
            range_rate: None,
            right_ascension_rate: None,
            declination_rate: None,
        }
    }

    /// Create RA/Dec from degrees.
    pub fn from_degrees(epoch: DateTime<Utc>, ra_deg: f64, dec_deg: f64) -> Self {
        Self::new(
            epoch,
            ra_deg * Constants::DEG_TO_RAD,
            dec_deg * Constants::DEG_TO_RAD,
        )
    }

    /// Create RA/Dec from hours/minutes/seconds (RA) and degrees/arcmin/arcsec (Dec).
    pub fn from_hms_dms(
        epoch: DateTime<Utc>,
        ra_h: i32,
        ra_m: i32,
        ra_s: f64,
        dec_d: i32,
        dec_am: i32,
        dec_as: f64,
    ) -> Self {
        // Convert RA from HMS to radians (1 hour = 15 degrees)
        let ra_hours = ra_h as f64 + ra_m as f64 / 60.0 + ra_s / 3600.0;
        let ra_rad = ra_hours * 15.0 * Constants::DEG_TO_RAD;

        // Convert Dec from DMS to radians
        let sign = if dec_d < 0 { -1.0 } else { 1.0 };
        let dec_deg = dec_d.abs() as f64 + dec_am as f64 / 60.0 + dec_as / 3600.0;
        let dec_rad = sign * dec_deg * Constants::DEG_TO_RAD;

        Self::new(epoch, ra_rad, dec_rad)
    }

    /// Set the range.
    pub fn with_range(mut self, range: f64) -> Self {
        self.range = Some(range);
        self
    }

    /// Set the range rate.
    pub fn with_range_rate(mut self, range_rate: f64) -> Self {
        self.range_rate = Some(range_rate);
        self
    }

    /// Set the RA rate.
    pub fn with_ra_rate(mut self, ra_rate: f64) -> Self {
        self.right_ascension_rate = Some(ra_rate);
        self
    }

    /// Set the Dec rate.
    pub fn with_dec_rate(mut self, dec_rate: f64) -> Self {
        self.declination_rate = Some(dec_rate);
        self
    }

    /// Get right ascension in degrees.
    pub fn ra_deg(&self) -> f64 {
        self.right_ascension * Constants::RAD_TO_DEG
    }

    /// Get declination in degrees.
    pub fn dec_deg(&self) -> f64 {
        self.declination * Constants::RAD_TO_DEG
    }

    /// Convert to degrees, returning (ra_deg, dec_deg).
    pub fn to_degrees(&self) -> (f64, f64) {
        (self.ra_deg(), self.dec_deg())
    }

    /// Get right ascension in hours.
    pub fn ra_hours(&self) -> f64 {
        self.right_ascension * 12.0 / PI
    }

    /// Convert RA to hours/minutes/seconds string format.
    pub fn ra_to_hms(&self) -> String {
        let hours_total = self.ra_hours();
        let h = hours_total.floor() as i32;
        let m_total = (hours_total - h as f64) * 60.0;
        let m = m_total.floor() as i32;
        let s = (m_total - m as f64) * 60.0;
        format!("{}h {}m {:.2}s", h, m, s)
    }

    /// Convert Dec to degrees/arcmin/arcsec string format.
    pub fn dec_to_dms(&self) -> String {
        let dec_deg = self.dec_deg();
        let sign = if dec_deg < 0.0 { "-" } else { "+" };
        let dec_deg = dec_deg.abs();
        let d = dec_deg.floor() as i32;
        let m_total = (dec_deg - d as f64) * 60.0;
        let m = m_total.floor() as i32;
        let s = (m_total - m as f64) * 60.0;
        format!("{}{}° {}' {:.2}\"", sign, d, m, s)
    }

    /// Format as a human-readable string.
    pub fn to_string(&self) -> String {
        format!("RA: {}, Dec: {}", self.ra_to_hms(), self.dec_to_dms())
    }
}

// ============================================================================
// Observations - Observation Calculations
// ============================================================================

/// Functions for computing angular observations from observer and target states.
pub struct Observations;

impl Observations {
    // ========================================================================
    // RA/Dec Observations
    // ========================================================================

    /// Compute RA/Dec observation of a target from an observer.
    ///
    /// # Arguments
    /// * `observer` - GeodeticState of the ground observer
    /// * `target` - ECIState of the target
    ///
    /// # Returns
    /// RaDec observation struct
    pub fn compute_ra_dec(observer: &GeodeticState, target: &ECIState) -> RaDec {
        // Get observer position in ECI
        let observer_eci = observer.to_eci(target.epoch);

        // Compute relative position vector (target - observer)
        let rel_pos = target.position - observer_eci.position;
        let rel_vel = target.velocity - observer_eci.velocity;

        let (dx, dy, dz) = (rel_pos.x, rel_pos.y, rel_pos.z);
        let (_dvx, _dvy, _dvz) = (rel_vel.x, rel_vel.y, rel_vel.z);

        // Compute range
        let range = rel_pos.magnitude();

        // Compute RA and Dec
        let mut ra = dy.atan2(dx);
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }
        let dec = (dz / range).asin();

        RaDec::new(target.epoch, ra, dec).with_range(range)
    }

    /// Compute RA/Dec observation with angular rates.
    pub fn compute_ra_dec_with_rates(observer: &GeodeticState, target: &ECIState) -> RaDec {
        let observer_eci = observer.to_eci(target.epoch);
        let rel_pos = target.position - observer_eci.position;
        let rel_vel = target.velocity - observer_eci.velocity;

        let (dx, dy, dz) = (rel_pos.x, rel_pos.y, rel_pos.z);
        let (dvx, dvy, dvz) = (rel_vel.x, rel_vel.y, rel_vel.z);

        let range = rel_pos.magnitude();

        let mut ra = dy.atan2(dx);
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }
        let dec = (dz / range).asin();

        // Range rate: (r · v) / |r|
        let range_rate = (dx * dvx + dy * dvy + dz * dvz) / range;

        // RA rate: d/dt[atan2(y, x)] = (x*vy - y*vx) / (x² + y²)
        let xy_sq = dx * dx + dy * dy;
        let ra_rate = if xy_sq > 1e-10 {
            (dx * dvy - dy * dvx) / xy_sq
        } else {
            0.0
        };

        // Dec rate: d/dt[asin(z/r)] = (vz*r - z*dr/dt) / (r² * sqrt(1 - (z/r)²))
        let cos_dec = dec.cos();
        let dec_rate = if cos_dec > 1e-10 {
            (dvz * range - dz * range_rate) / (range * range * cos_dec)
        } else {
            0.0
        };

        RaDec::new(target.epoch, ra, dec)
            .with_range(range)
            .with_range_rate(range_rate)
            .with_ra_rate(ra_rate)
            .with_dec_rate(dec_rate)
    }

    /// Compute geocentric RA/Dec directly from ECI position.
    ///
    /// This gives the geocentric RA/Dec (as seen from Earth's center),
    /// not topocentric (as seen from a ground observer).
    pub fn geocentric_ra_dec(target: &ECIState) -> RaDec {
        let (x, y, z) = (target.position.x, target.position.y, target.position.z);
        let range = target.position.magnitude();

        let mut ra = y.atan2(x);
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }
        let dec = (z / range).asin();

        RaDec::new(target.epoch, ra, dec).with_range(range)
    }

    // ========================================================================
    // Az/El Observations
    // ========================================================================

    /// Compute Az/El observation of a target from an observer.
    ///
    /// # Arguments
    /// * `observer` - GeodeticState of the ground observer
    /// * `target` - ECIState of the target
    ///
    /// # Returns
    /// AzEl observation struct
    pub fn compute_az_el(observer: &GeodeticState, target: &ECIState) -> AzEl {
        // Get observer ECEF position
        let observer_ecef = observer.to_ecef(target.epoch);

        // Convert target to ECEF
        let target_ecef = StateTransforms::eci_to_ecef(target);

        // Compute relative position in ECEF
        let rel_ecef = target_ecef.position - observer_ecef.position;

        // Transform to SEZ (South-East-Zenith) local coordinates
        let (south, east, zenith) = ecef_to_sez(&rel_ecef, observer.latitude, observer.longitude);

        // Compute range
        let range = (south * south + east * east + zenith * zenith).sqrt();

        // Compute Az and El
        // Azimuth is measured clockwise from North
        // In SEZ: North = -South, so Az = atan2(East, -South)
        let mut az = east.atan2(-south);
        if az < 0.0 {
            az += Constants::TWO_PI;
        }

        // Elevation is angle above horizon
        let el = (zenith / range).asin();

        AzEl::new(target.epoch, az, el).with_range(range)
    }

    /// Compute Az/El observation with angular rates.
    pub fn compute_az_el_with_rates(observer: &GeodeticState, target: &ECIState) -> AzEl {
        let observer_ecef = observer.to_ecef(target.epoch);
        let target_ecef = StateTransforms::eci_to_ecef(target);

        let rel_ecef = target_ecef.position - observer_ecef.position;
        let (south, east, zenith) = ecef_to_sez(&rel_ecef, observer.latitude, observer.longitude);

        let range = (south * south + east * east + zenith * zenith).sqrt();

        let mut az = east.atan2(-south);
        if az < 0.0 {
            az += Constants::TWO_PI;
        }
        let el = (zenith / range).asin();

        // Relative velocity in ECEF and transform to SEZ
        let rel_vel_ecef = target_ecef.velocity - observer_ecef.velocity;
        let (ds, de, dz_vel) = ecef_to_sez(&rel_vel_ecef, observer.latitude, observer.longitude);

        // Range rate
        let range_rate = (south * ds + east * de + zenith * dz_vel) / range;

        // Azimuth rate: d/dt[atan2(E, -S)] = (-S*dE - E*(-dS)) / (S² + E²)
        let horiz_sq = south * south + east * east;
        let az_rate = if horiz_sq > 1e-10 {
            (south * de - east * ds) / horiz_sq
        } else {
            0.0
        };

        // Elevation rate: d/dt[asin(Z/r)]
        let cos_el = el.cos();
        let el_rate = if cos_el > 1e-10 {
            (dz_vel * range - zenith * range_rate) / (range * range * cos_el)
        } else {
            0.0
        };

        AzEl::new(target.epoch, az, el)
            .with_range(range)
            .with_range_rate(range_rate)
            .with_azimuth_rate(az_rate)
            .with_elevation_rate(el_rate)
    }

    // ========================================================================
    // Coordinate Conversions
    // ========================================================================

    /// Convert RA/Dec to Az/El for a given observer and time.
    pub fn ra_dec_to_az_el(ra_dec: &RaDec, observer: &GeodeticState) -> AzEl {
        let lat_rad = observer.latitude_rad();
        let lst = observer.local_sidereal_time(&ra_dec.epoch);

        // Hour angle = LST - RA
        let ha = lst - ra_dec.right_ascension;

        let sin_dec = ra_dec.declination.sin();
        let cos_dec = ra_dec.declination.cos();
        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_ha = ha.sin();
        let cos_ha = ha.cos();

        // Elevation
        let sin_el = sin_dec * sin_lat + cos_dec * cos_lat * cos_ha;
        let el = sin_el.asin();

        // Azimuth
        let cos_el = el.cos();
        let sin_az = -cos_dec * sin_ha / cos_el;
        let cos_az = (sin_dec - sin_lat * sin_el) / (cos_lat * cos_el);

        let mut az = sin_az.atan2(cos_az);
        if az < 0.0 {
            az += Constants::TWO_PI;
        }

        let mut az_el = AzEl::new(ra_dec.epoch, az, el);
        if let Some(range) = ra_dec.range {
            az_el = az_el.with_range(range);
        }
        az_el
    }

    /// Convert Az/El to RA/Dec for a given observer and time.
    pub fn az_el_to_ra_dec(az_el: &AzEl, observer: &GeodeticState) -> RaDec {
        let lat_rad = observer.latitude_rad();
        let lst = observer.local_sidereal_time(&az_el.epoch);

        let sin_el = az_el.elevation.sin();
        let cos_el = az_el.elevation.cos();
        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_az = az_el.azimuth.sin();
        let cos_az = az_el.azimuth.cos();

        // Declination
        let sin_dec = sin_el * sin_lat + cos_el * cos_lat * cos_az;
        let dec = sin_dec.asin();

        // Hour angle
        let cos_dec = dec.cos();
        let sin_ha = -cos_el * sin_az / cos_dec;
        let cos_ha = (sin_el - sin_lat * sin_dec) / (cos_lat * cos_dec);
        let ha = sin_ha.atan2(cos_ha);

        // RA = LST - HA
        let mut ra = lst - ha;
        ra = ra % Constants::TWO_PI;
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }

        let mut ra_dec = RaDec::new(az_el.epoch, ra, dec);
        if let Some(range) = az_el.range {
            ra_dec = ra_dec.with_range(range);
        }
        ra_dec
    }

    // ========================================================================
    // Line of Sight Vectors
    // ========================================================================

    /// Compute unit vector in ECI frame from RA/Dec.
    pub fn ra_dec_to_eci_direction(ra_dec: &RaDec) -> Vector3 {
        let cos_dec = ra_dec.declination.cos();
        Vector3::new(
            cos_dec * ra_dec.right_ascension.cos(),
            cos_dec * ra_dec.right_ascension.sin(),
            ra_dec.declination.sin(),
        )
    }

    /// Compute unit vector in local SEZ frame from Az/El.
    pub fn az_el_to_sez_direction(az_el: &AzEl) -> Vector3 {
        let cos_el = az_el.elevation.cos();
        // SEZ: South-East-Zenith
        // Az measured from North, clockwise
        Vector3::new(
            -cos_el * az_el.azimuth.cos(), // South = -North
            cos_el * az_el.azimuth.sin(),  // East
            az_el.elevation.sin(),         // Zenith
        )
    }

    // ========================================================================
    // Visibility Checks
    // ========================================================================

    /// Check if a target is visible from an observer (above horizon).
    ///
    /// # Arguments
    /// * `observer` - Ground observer location
    /// * `target` - Target state in ECI
    /// * `min_elevation_deg` - Minimum elevation in degrees (default 0.0)
    pub fn is_visible(observer: &GeodeticState, target: &ECIState, min_elevation_deg: f64) -> bool {
        let az_el = Self::compute_az_el(observer, target);
        az_el.above_elevation(min_elevation_deg)
    }

    /// Compute the look angles and range from observer to target.
    ///
    /// Returns (az_deg, el_deg, range_km).
    pub fn look_angles(observer: &GeodeticState, target: &ECIState) -> (f64, f64, f64) {
        let az_el = Self::compute_az_el(observer, target);
        (
            az_el.azimuth_deg(),
            az_el.elevation_deg(),
            az_el.range.unwrap_or(0.0),
        )
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Normalize azimuth to [0, 2π) range.
fn normalize_azimuth(az: f64) -> f64 {
    let mut normalized = az % Constants::TWO_PI;
    if normalized < 0.0 {
        normalized += Constants::TWO_PI;
    }
    normalized
}

/// Normalize right ascension to [0, 2π) range.
fn normalize_ra(ra: f64) -> f64 {
    let mut normalized = ra % Constants::TWO_PI;
    if normalized < 0.0 {
        normalized += Constants::TWO_PI;
    }
    normalized
}

/// Transform ECEF vector to SEZ (South-East-Zenith) local frame.
fn ecef_to_sez(ecef: &Vector3, lat_deg: f64, lon_deg: f64) -> (f64, f64, f64) {
    let lat_rad = lat_deg * Constants::DEG_TO_RAD;
    let lon_rad = lon_deg * Constants::DEG_TO_RAD;

    let sin_lat = lat_rad.sin();
    let cos_lat = lat_rad.cos();
    let sin_lon = lon_rad.sin();
    let cos_lon = lon_rad.cos();

    let (x, y, z) = (ecef.x, ecef.y, ecef.z);

    // Rotation matrix ECEF -> SEZ
    let south = sin_lat * cos_lon * x + sin_lat * sin_lon * y - cos_lat * z;
    let east = -sin_lon * x + cos_lon * y;
    let zenith = cos_lat * cos_lon * x + cos_lat * sin_lon * y + sin_lat * z;

    (south, east, zenith)
}

/// Transform SEZ vector to ECEF.
#[allow(dead_code)]
fn sez_to_ecef(south: f64, east: f64, zenith: f64, lat_deg: f64, lon_deg: f64) -> Vector3 {
    let lat_rad = lat_deg * Constants::DEG_TO_RAD;
    let lon_rad = lon_deg * Constants::DEG_TO_RAD;

    let sin_lat = lat_rad.sin();
    let cos_lat = lat_rad.cos();
    let sin_lon = lon_rad.sin();
    let cos_lon = lon_rad.cos();

    // Inverse rotation matrix SEZ -> ECEF
    let x = sin_lat * cos_lon * south - sin_lon * east + cos_lat * cos_lon * zenith;
    let y = sin_lat * sin_lon * south + cos_lon * east + cos_lat * sin_lon * zenith;
    let z = -cos_lat * south + sin_lat * zenith;

    Vector3::new(x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    const EPSILON: f64 = 1e-6;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn test_az_el_creation() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let az_el = AzEl::from_degrees(epoch, 90.0, 45.0);

        assert!(approx_eq(az_el.azimuth_deg(), 90.0, EPSILON));
        assert!(approx_eq(az_el.elevation_deg(), 45.0, EPSILON));
    }

    #[test]
    fn test_az_el_above_horizon() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        let above = AzEl::from_degrees(epoch, 180.0, 10.0);
        let below = AzEl::from_degrees(epoch, 180.0, -10.0);

        assert!(above.above_horizon());
        assert!(!below.above_horizon());
    }

    #[test]
    fn test_az_el_compass_direction() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        assert_eq!(
            AzEl::from_degrees(epoch, 0.0, 45.0).compass_direction(),
            "N"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 45.0, 45.0).compass_direction(),
            "NE"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 90.0, 45.0).compass_direction(),
            "E"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 135.0, 45.0).compass_direction(),
            "SE"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 180.0, 45.0).compass_direction(),
            "S"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 225.0, 45.0).compass_direction(),
            "SW"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 270.0, 45.0).compass_direction(),
            "W"
        );
        assert_eq!(
            AzEl::from_degrees(epoch, 315.0, 45.0).compass_direction(),
            "NW"
        );
    }

    #[test]
    fn test_ra_dec_creation() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let ra_dec = RaDec::from_degrees(epoch, 180.0, 45.0);

        assert!(approx_eq(ra_dec.ra_deg(), 180.0, EPSILON));
        assert!(approx_eq(ra_dec.dec_deg(), 45.0, EPSILON));
    }

    #[test]
    fn test_ra_dec_hms() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        // 6 hours RA = 90 degrees
        let ra_dec = RaDec::from_hms_dms(epoch, 6, 0, 0.0, 45, 0, 0.0);

        assert!(approx_eq(ra_dec.ra_deg(), 90.0, 0.01));
        assert!(approx_eq(ra_dec.dec_deg(), 45.0, 0.01));
        assert!(approx_eq(ra_dec.ra_hours(), 6.0, 0.01));
    }

    #[test]
    fn test_compute_az_el() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        // Observer at mid-latitudes
        let observer = GeodeticState::new(40.0, -105.0, 1.6);

        // Target in a high orbit position that should be visible
        // Put satellite at a position that's definitely above the horizon
        let target_pos = Vector3::new(20000.0, 20000.0, 15000.0);
        let target = ECIState::new(epoch, target_pos, Vector3::zero());

        let az_el = Observations::compute_az_el(&observer, &target);

        // Range should be positive
        assert!(az_el.range.unwrap() > 0.0);

        // Azimuth should be in valid range
        assert!(az_el.azimuth >= 0.0 && az_el.azimuth < Constants::TWO_PI);
    }

    #[test]
    fn test_is_visible() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        // Observer in Denver
        let observer = GeodeticState::new(39.7392, -104.9903, 1.6);

        // Target in orbit
        let target_pos = Vector3::new(7000.0, 1000.0, 500.0);
        let target_vel = Vector3::new(-1.0, 7.0, 0.5);
        let target = ECIState::new(epoch, target_pos, target_vel);

        // Just check that the function runs without error
        let _visible = Observations::is_visible(&observer, &target, 0.0);
    }

    #[test]
    fn test_look_angles() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        let observer = GeodeticState::new(40.0, -105.0, 1.6);
        let target_pos = Vector3::new(7000.0, 0.0, 3000.0);
        let target = ECIState::new(epoch, target_pos, Vector3::zero());

        let (az, _el, range) = Observations::look_angles(&observer, &target);

        // Azimuth should be in valid range
        assert!(az >= 0.0 && az < 360.0);
        // Range should be positive
        assert!(range > 0.0);
    }

    #[test]
    fn test_normalize_azimuth() {
        assert!(approx_eq(normalize_azimuth(0.0), 0.0, EPSILON));
        assert!(approx_eq(
            normalize_azimuth(Constants::TWO_PI),
            0.0,
            EPSILON
        ));
        assert!(approx_eq(
            normalize_azimuth(-PI / 2.0),
            3.0 * PI / 2.0,
            EPSILON
        ));
    }

    #[test]
    fn test_ecef_to_sez() {
        // At equator, prime meridian, a point directly above should be in zenith direction
        let ecef = Vector3::new(1000.0, 0.0, 0.0);
        let (south, east, zenith) = ecef_to_sez(&ecef, 0.0, 0.0);

        // For a point on the X-axis at equator/prime meridian, it should be mostly zenith
        assert!(zenith > south.abs());
        assert!(zenith > east.abs());
    }
}
