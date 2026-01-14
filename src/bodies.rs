//! Celestial body parameters and position calculations.
//!
//! This module provides:
//! - **Earth**: Parameters, precession, and nutation calculations
//! - **Sun**: Solar position calculations
//! - **Moon**: Lunar position calculations

use crate::constants::Constants;
use crate::data::{DefaultEOPProvider, IAU1980};
use crate::math::{poly_eval, Vector3};
use crate::time::{TimeTransforms, UTC};
use chrono::{DateTime, Utc};

// ============================================================================
// Earth
// ============================================================================

/// Precession angles for the Earth (in radians).
#[derive(Debug, Clone, Copy)]
pub struct PrecessionAngles {
    /// Zeta angle
    pub zeta: f64,
    /// Theta angle
    pub theta: f64,
    /// Z angle
    pub z: f64,
}

/// Nutation angles for the Earth (in radians).
#[derive(Debug, Clone, Copy)]
pub struct NutationAngles {
    /// Nutation in longitude (Δψ)
    pub d_psi: f64,
    /// Nutation in obliquity (Δε)
    pub d_eps: f64,
    /// Mean obliquity of the ecliptic
    pub m_eps: f64,
    /// True obliquity of the ecliptic (m_eps + d_eps)
    pub eps: f64,
    /// Equation of the equinoxes
    pub eq_eq: f64,
    /// Greenwich Apparent Sidereal Time
    pub gast: f64,
}

/// Earth parameters and orientation calculations.
pub struct Earth;

impl Earth {
    // Precession polynomial coefficients (arcseconds)
    const ZETA_POLY: [f64; 4] = [
        0.017988 * Constants::ARCSEC_TO_RAD,
        0.30188 * Constants::ARCSEC_TO_RAD,
        2306.2181 * Constants::ARCSEC_TO_RAD,
        0.0,
    ];

    const THETA_POLY: [f64; 4] = [
        -0.041833 * Constants::ARCSEC_TO_RAD,
        -0.42665 * Constants::ARCSEC_TO_RAD,
        2004.3109 * Constants::ARCSEC_TO_RAD,
        0.0,
    ];

    const Z_POLY: [f64; 4] = [
        0.018203 * Constants::ARCSEC_TO_RAD,
        1.09468 * Constants::ARCSEC_TO_RAD,
        2306.2181 * Constants::ARCSEC_TO_RAD,
        0.0,
    ];

    // Nutation polynomial coefficients (degrees -> radians)
    const LUNAR_ANOMALY_POLY: [f64; 4] = [
        1.78e-5 * Constants::DEG_TO_RAD,
        0.0086972 * Constants::DEG_TO_RAD,
        (1325.0 * 360.0 + 198.8673981) * Constants::DEG_TO_RAD,
        134.96298139 * Constants::DEG_TO_RAD,
    ];

    const SOLAR_ANOMALY_POLY: [f64; 4] = [
        -3.3e-6 * Constants::DEG_TO_RAD,
        -0.0001603 * Constants::DEG_TO_RAD,
        (99.0 * 360.0 + 359.05034) * Constants::DEG_TO_RAD,
        357.52772333 * Constants::DEG_TO_RAD,
    ];

    const LUNAR_LATITUDE_POLY: [f64; 4] = [
        3.1e-6 * Constants::DEG_TO_RAD,
        -0.0036825 * Constants::DEG_TO_RAD,
        (1342.0 * 360.0 + 82.0175381) * Constants::DEG_TO_RAD,
        93.27191028 * Constants::DEG_TO_RAD,
    ];

    const SUN_ELONGATION_POLY: [f64; 4] = [
        5.3e-6 * Constants::DEG_TO_RAD,
        -0.0019142 * Constants::DEG_TO_RAD,
        (1236.0 * 360.0 + 307.1114800) * Constants::DEG_TO_RAD,
        297.85036306 * Constants::DEG_TO_RAD,
    ];

    const LUNAR_RAAN_POLY: [f64; 4] = [
        2.2e-6 * Constants::DEG_TO_RAD,
        0.0020708 * Constants::DEG_TO_RAD,
        -(5.0 * 360.0 + 134.1362608) * Constants::DEG_TO_RAD,
        125.04452222 * Constants::DEG_TO_RAD,
    ];

    const MEAN_EPSILON_POLY: [f64; 4] = [
        5.04e-7 * Constants::DEG_TO_RAD,
        -1.64e-7 * Constants::DEG_TO_RAD,
        -0.0130042 * Constants::DEG_TO_RAD,
        23.439291 * Constants::DEG_TO_RAD,
    ];

    /// Earth gravitational parameter in m³/s²
    pub const MU: f64 = Constants::EARTH_MU_M;

    /// Earth gravitational parameter in km³/s²
    pub const MU_KM: f64 = Constants::EARTH_MU_KM;

    /// Earth equatorial radius in meters
    pub const EQUATORIAL_RADIUS: f64 = Constants::EARTH_RADIUS_EQ_M;

    /// Earth equatorial radius in km
    pub const EQUATORIAL_RADIUS_KM: f64 = Constants::EARTH_RADIUS_EQ_KM;

    /// Earth flattening
    pub const FLATTENING: f64 = Constants::EARTH_FLATTENING;

    /// Earth polar radius in meters
    pub fn polar_radius() -> f64 {
        Self::EQUATORIAL_RADIUS * (1.0 - Self::FLATTENING)
    }

    /// Earth mean radius in meters
    pub fn mean_radius() -> f64 {
        (2.0 * Self::EQUATORIAL_RADIUS + Self::polar_radius()) / 3.0
    }

    /// Earth sidereal rotation rate in rad/s
    pub const ROTATION_RATE: f64 = Constants::EARTH_ROTATION_RATE;

    /// Earth J2 coefficient
    pub const J2: f64 = Constants::EARTH_J2;

    /// Calculate precession angles for a given epoch.
    ///
    /// Returns the precession angles (zeta, theta, z) in radians.
    pub fn precession_angles(epoch: &DateTime<Utc>) -> PrecessionAngles {
        let utc = UTC::from_datetime(epoch);
        let jd = utc.to_jd();
        let t = (jd - Constants::J2000_JD) / Constants::JULIAN_CENTURY;

        PrecessionAngles {
            zeta: poly_eval(&Self::ZETA_POLY, t),
            theta: poly_eval(&Self::THETA_POLY, t),
            z: poly_eval(&Self::Z_POLY, t),
        }
    }

    /// Calculate nutation angles for a given epoch.
    ///
    /// Uses the IAU 1980 nutation theory with optional EOP corrections.
    /// The `num_coeffs` parameter controls how many terms to use (default 4).
    pub fn nutation_angles(epoch: &DateTime<Utc>) -> NutationAngles {
        Self::nutation_angles_with_coeffs(epoch, 4, true)
    }

    /// Calculate nutation angles with configurable precision.
    ///
    /// # Arguments
    /// * `epoch` - UTC epoch
    /// * `num_coeffs` - Number of IAU 1980 coefficients to use (more = more precise)
    /// * `use_eop` - Whether to apply EOP corrections
    pub fn nutation_angles_with_coeffs(
        epoch: &DateTime<Utc>,
        num_coeffs: usize,
        use_eop: bool,
    ) -> NutationAngles {
        let utc = UTC::from_datetime(epoch);
        let tt = TimeTransforms::utc_to_tt(&utc);
        let t = tt.julian_centuries_j2000();

        // Compute fundamental arguments
        let lunar_anomaly = poly_eval(&Self::LUNAR_ANOMALY_POLY, t);
        let solar_anomaly = poly_eval(&Self::SOLAR_ANOMALY_POLY, t);
        let lunar_latitude = poly_eval(&Self::LUNAR_LATITUDE_POLY, t);
        let sun_elongation = poly_eval(&Self::SUN_ELONGATION_POLY, t);
        let lunar_raan = poly_eval(&Self::LUNAR_RAAN_POLY, t);

        // Sum nutation terms
        let mut delta_psi = 0.0;
        let mut delta_eps = 0.0;

        let count = num_coeffs.min(IAU1980::coefficient_count());
        for i in 0..count {
            if let Some(coef) = IAU1980::get_coefficients(i) {
                let arg = coef.a1 as f64 * lunar_anomaly
                    + coef.a2 as f64 * solar_anomaly
                    + coef.a3 as f64 * lunar_latitude
                    + coef.a4 as f64 * sun_elongation
                    + coef.a5 as f64 * lunar_raan;

                let sin_c = coef.ai + coef.bi * t;
                let cos_c = coef.ci + coef.di * t;

                delta_psi += sin_c * arg.sin();
                delta_eps += cos_c * arg.cos();
            }
        }

        // Convert from 0.0001 arcseconds to radians
        let mut d_psi = delta_psi * Constants::TT_ARCSEC_TO_RAD;
        let mut d_eps = delta_eps * Constants::TT_ARCSEC_TO_RAD;

        // Apply EOP corrections if requested
        if use_eop {
            let eop = DefaultEOPProvider::get_eop(utc.to_mjd());
            d_psi += eop.d_psi_rad();
            d_eps += eop.d_eps_rad();
        }

        // Mean and true obliquity
        let m_eps = poly_eval(&Self::MEAN_EPSILON_POLY, t);
        let eps = m_eps + d_eps;

        // Greenwich Mean Sidereal Time
        let gmst = TimeTransforms::utc_to_gmst(&utc);

        // Equation of the equinoxes
        let eq_eq = d_psi * eps.cos()
            + 0.00264 * Constants::ARCSEC_TO_RAD * lunar_raan.sin()
            + 0.000063 * Constants::ARCSEC_TO_RAD * (2.0 * lunar_raan).sin();

        // Greenwich Apparent Sidereal Time
        let gast = gmst.to_radians() + eq_eq;

        NutationAngles {
            d_psi,
            d_eps,
            m_eps,
            eps,
            eq_eq,
            gast,
        }
    }
}

// ============================================================================
// Sun
// ============================================================================

/// Solar position calculations.
///
/// Provides functions to compute the Sun's position in various coordinate frames.
/// Uses the low-precision solar coordinates from the Astronomical Almanac
/// (approximate accuracy ~0.01 degrees).
pub struct Sun;

impl Sun {
    // Mean longitude polynomial coefficients (degrees)
    const MEAN_LONGITUDE_POLY: [f64; 2] = [36000.77005361, 280.4606184];

    // Mean anomaly polynomial coefficients (degrees)
    const MEAN_ANOMALY_POLY: [f64; 2] = [35999.05034, 357.5277233];

    // Obliquity polynomial coefficients (degrees)
    const OBLIQUITY_POLY: [f64; 2] = [-0.0130042, 23.439291];

    /// Sun's gravitational parameter in m³/s²
    pub const MU: f64 = Constants::SUN_MU;

    /// Astronomical Unit in meters
    pub const AU: f64 = Constants::AU_M;

    /// Solar radius in meters
    pub const RADIUS: f64 = Constants::SUN_RADIUS_M;

    /// Calculate the Sun's position in ECI frame (meters).
    ///
    /// Uses the low-precision solar coordinates from the Astronomical Almanac.
    pub fn eci_position(epoch: &DateTime<Utc>) -> Vector3 {
        let utc = UTC::from_datetime(epoch);
        let tt = TimeTransforms::utc_to_tt(&utc);
        let t = tt.julian_centuries_j2000();

        // Mean longitude (degrees)
        let l0 = poly_eval(&Self::MEAN_LONGITUDE_POLY, t);

        // Mean anomaly (degrees -> radians)
        let m = poly_eval(&Self::MEAN_ANOMALY_POLY, t);
        let m_rad = normalize_degrees(m) * Constants::DEG_TO_RAD;

        // Ecliptic longitude (degrees)
        let lambda = l0 + 1.9146 * m_rad.sin() + 0.0200 * (2.0 * m_rad).sin();
        let lambda_rad = lambda * Constants::DEG_TO_RAD;

        // Obliquity (degrees -> radians)
        let epsilon = poly_eval(&Self::OBLIQUITY_POLY, t);
        let epsilon_rad = epsilon * Constants::DEG_TO_RAD;

        // Distance from Earth to Sun in AU
        let r_au = 1.00014 - 0.01671 * m_rad.cos() - 0.00014 * (2.0 * m_rad).cos();
        let r_m = r_au * Self::AU;

        // Convert to ECI coordinates
        let cos_lambda = lambda_rad.cos();
        let sin_lambda = lambda_rad.sin();
        let cos_epsilon = epsilon_rad.cos();
        let sin_epsilon = epsilon_rad.sin();

        Vector3::new(
            r_m * cos_lambda,
            r_m * sin_lambda * cos_epsilon,
            r_m * sin_lambda * sin_epsilon,
        )
    }

    /// Calculate the Sun's position in ECI frame (kilometers).
    pub fn eci_position_km(epoch: &DateTime<Utc>) -> Vector3 {
        let pos_m = Self::eci_position(epoch);
        Vector3::new(pos_m.x / 1000.0, pos_m.y / 1000.0, pos_m.z / 1000.0)
    }

    /// Calculate the Sun's apparent right ascension and declination.
    ///
    /// Returns (right_ascension, declination) in radians.
    pub fn apparent_position(epoch: &DateTime<Utc>) -> (f64, f64) {
        let pos = Self::eci_position(epoch);
        let r = pos.magnitude();

        // Declination
        let dec = (pos.z / r).asin();

        // Right ascension
        let mut ra = pos.y.atan2(pos.x);
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }

        (ra, dec)
    }

    /// Calculate the distance from Earth to the Sun in meters.
    pub fn distance(epoch: &DateTime<Utc>) -> f64 {
        Self::eci_position(epoch).magnitude()
    }

    /// Calculate the Sun's unit vector from Earth in ECI frame.
    pub fn direction(epoch: &DateTime<Utc>) -> Vector3 {
        Self::eci_position(epoch).normalize()
    }
}

// ============================================================================
// Moon
// ============================================================================

/// Lunar position calculations.
///
/// Provides functions to compute the Moon's position in ECI frame.
/// Uses the low-precision lunar coordinates from the Astronomical Almanac
/// (approximate accuracy ~0.3 degrees in longitude, ~0.2 degrees in latitude).
pub struct Moon;

impl Moon {
    // Mean longitude polynomial coefficients (degrees)
    const MEAN_LONGITUDE_POLY: [f64; 5] = [
        -1.52e-6,
        1.0 / 538841.0,
        -0.0015786,
        481267.88123421,
        218.3164477,
    ];

    // Mean elongation polynomial coefficients (degrees)
    const MEAN_ELONGATION_POLY: [f64; 5] = [
        -8.78e-6,
        -1.0 / 113065.0,
        0.0018819,
        445267.1114034,
        297.8501921,
    ];

    // Sun's mean anomaly polynomial coefficients (degrees)
    const SUN_MEAN_ANOMALY_POLY: [f64; 5] = [
        1.59e-5,
        -1.0 / 24490000.0,
        -0.0001536,
        35999.0502909,
        357.5291092,
    ];

    // Moon's mean anomaly polynomial coefficients (degrees)
    const MOON_MEAN_ANOMALY_POLY: [f64; 5] = [
        3.239e-5,
        -1.0 / 69699.0,
        0.0087414,
        477198.8675055,
        134.9633964,
    ];

    // Moon's argument of latitude polynomial coefficients (degrees)
    const ARG_LATITUDE_POLY: [f64; 5] = [
        3.17e-6,
        -1.0 / 3526000.0,
        -0.0036539,
        483202.0175233,
        93.2720950,
    ];

    // Obliquity polynomial coefficients (degrees)
    const OBLIQUITY_POLY: [f64; 2] = [-0.0130042, 23.439291];

    /// Moon's gravitational parameter in m³/s²
    pub const MU: f64 = Constants::MOON_MU;

    /// Mean Earth-Moon distance in meters
    pub const MEAN_DISTANCE: f64 = Constants::MOON_MEAN_DISTANCE_M;

    /// Lunar radius in meters
    pub const RADIUS: f64 = Constants::MOON_RADIUS_M;

    /// Calculate the Moon's position in ECI frame (meters).
    ///
    /// Uses the low-precision lunar coordinates from the Astronomical Almanac.
    pub fn eci_position(epoch: &DateTime<Utc>) -> Vector3 {
        let utc = UTC::from_datetime(epoch);
        let tt = TimeTransforms::utc_to_tt(&utc);
        let t = tt.julian_centuries_j2000();

        // Compute fundamental arguments (degrees -> radians)
        let l_prime = poly_eval(&Self::MEAN_LONGITUDE_POLY, t);
        let d =
            normalize_degrees(poly_eval(&Self::MEAN_ELONGATION_POLY, t)) * Constants::DEG_TO_RAD;
        let m =
            normalize_degrees(poly_eval(&Self::SUN_MEAN_ANOMALY_POLY, t)) * Constants::DEG_TO_RAD;
        let m_prime =
            normalize_degrees(poly_eval(&Self::MOON_MEAN_ANOMALY_POLY, t)) * Constants::DEG_TO_RAD;
        let f = normalize_degrees(poly_eval(&Self::ARG_LATITUDE_POLY, t)) * Constants::DEG_TO_RAD;

        // Longitude corrections (degrees)
        let delta_lambda = 6.288774 * m_prime.sin()
            + 1.274027 * (2.0 * d - m_prime).sin()
            + 0.658314 * (2.0 * d).sin()
            + 0.213618 * (2.0 * m_prime).sin()
            - 0.185116 * m.sin()
            - 0.114332 * (2.0 * f).sin()
            + 0.058793 * (2.0 * (d - m_prime)).sin()
            + 0.057066 * (2.0 * d - m - m_prime).sin()
            + 0.053322 * (2.0 * d + m_prime).sin()
            + 0.045758 * (2.0 * d - m).sin()
            - 0.040923 * (m - m_prime).sin()
            - 0.034720 * d.sin()
            - 0.030383 * (m + m_prime).sin()
            + 0.015327 * (2.0 * (d - f)).sin()
            - 0.012528 * (m_prime + 2.0 * f).sin()
            + 0.010980 * (m_prime - 2.0 * f).sin();

        // Ecliptic longitude
        let lambda = l_prime + delta_lambda;
        let lambda_rad = lambda * Constants::DEG_TO_RAD;

        // Latitude corrections (degrees)
        let delta_beta = 5.128122 * f.sin()
            + 0.280602 * (m_prime + f).sin()
            + 0.277693 * (m_prime - f).sin()
            + 0.173237 * (2.0 * d - f).sin()
            + 0.055413 * (2.0 * d - m_prime + f).sin()
            + 0.046271 * (2.0 * d - m_prime - f).sin()
            + 0.032573 * (2.0 * d + f).sin()
            + 0.017198 * (2.0 * m_prime + f).sin()
            + 0.009266 * (2.0 * d + m_prime - f).sin();

        let beta_rad = delta_beta * Constants::DEG_TO_RAD;

        // Distance (in Earth radii)
        let delta_r = -0.40720 * m_prime.cos()
            - 0.18603 * (2.0 * d - m_prime).cos()
            - 0.01462 * (2.0 * d).cos()
            - 0.00122 * (2.0 * m_prime).cos()
            + 0.00079 * m.cos();

        // Mean distance in Earth radii is ~60.2666
        let r_earth_radii = 60.2666 + delta_r;
        let r_m = r_earth_radii * Constants::EARTH_RADIUS_EQ_M;

        // Obliquity
        let epsilon = poly_eval(&Self::OBLIQUITY_POLY, t);
        let epsilon_rad = epsilon * Constants::DEG_TO_RAD;

        // Convert ecliptic to equatorial (ECI)
        let cos_lambda = lambda_rad.cos();
        let sin_lambda = lambda_rad.sin();
        let cos_beta = beta_rad.cos();
        let sin_beta = beta_rad.sin();
        let cos_epsilon = epsilon_rad.cos();
        let sin_epsilon = epsilon_rad.sin();

        // Position in ecliptic coordinates
        let x_ecl = r_m * cos_beta * cos_lambda;
        let y_ecl = r_m * cos_beta * sin_lambda;
        let z_ecl = r_m * sin_beta;

        // Rotate to equatorial (ECI) frame
        Vector3::new(
            x_ecl,
            y_ecl * cos_epsilon - z_ecl * sin_epsilon,
            y_ecl * sin_epsilon + z_ecl * cos_epsilon,
        )
    }

    /// Calculate the Moon's position in ECI frame (kilometers).
    pub fn eci_position_km(epoch: &DateTime<Utc>) -> Vector3 {
        let pos_m = Self::eci_position(epoch);
        Vector3::new(pos_m.x / 1000.0, pos_m.y / 1000.0, pos_m.z / 1000.0)
    }

    /// Calculate the Moon's apparent right ascension and declination.
    ///
    /// Returns (right_ascension, declination) in radians.
    pub fn apparent_position(epoch: &DateTime<Utc>) -> (f64, f64) {
        let pos = Self::eci_position(epoch);
        let r = pos.magnitude();

        // Declination
        let dec = (pos.z / r).asin();

        // Right ascension
        let mut ra = pos.y.atan2(pos.x);
        if ra < 0.0 {
            ra += Constants::TWO_PI;
        }

        (ra, dec)
    }

    /// Calculate the distance from Earth to the Moon in meters.
    pub fn distance(epoch: &DateTime<Utc>) -> f64 {
        Self::eci_position(epoch).magnitude()
    }

    /// Calculate the Moon's unit vector from Earth in ECI frame.
    pub fn direction(epoch: &DateTime<Utc>) -> Vector3 {
        Self::eci_position(epoch).normalize()
    }

    /// Calculate the Moon's phase angle (angle between Sun and Moon as seen from Earth).
    ///
    /// Returns phase angle in radians (0 = new moon, π = full moon).
    pub fn phase_angle(epoch: &DateTime<Utc>) -> f64 {
        let moon_dir = Self::direction(epoch);
        let sun_dir = Sun::direction(epoch);

        // Phase angle is the angle between sun and moon vectors
        let cos_angle = moon_dir.dot(&sun_dir);
        cos_angle.clamp(-1.0, 1.0).acos()
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Normalize degrees to [0, 360) range.
fn normalize_degrees(deg: f64) -> f64 {
    let mut d = deg % 360.0;
    if d < 0.0 {
        d += 360.0;
    }
    d
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
    fn test_earth_parameters() {
        assert!(approx_eq(Earth::MU_KM, 398600.4418, 0.001));
        assert!(approx_eq(Earth::EQUATORIAL_RADIUS_KM, 6378.137, 0.001));
    }

    #[test]
    fn test_precession_angles_at_j2000() {
        // At J2000.0, precession angles should be approximately zero
        let epoch = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
        let precession = Earth::precession_angles(&epoch);

        // Should be close to zero at J2000
        assert!(precession.zeta.abs() < 0.01);
        assert!(precession.theta.abs() < 0.01);
        assert!(precession.z.abs() < 0.01);
    }

    #[test]
    fn test_nutation_angles() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let nutation = Earth::nutation_angles(&epoch);

        // Nutation in longitude should be small (typically < 20 arcsec = ~0.0001 rad)
        // With only 4 coefficients, the result may be larger, so use a more relaxed tolerance
        assert!(nutation.d_psi.abs() < 0.001); // ~200 arcsec tolerance

        // Mean obliquity should be around 23.4 degrees
        let m_eps_deg = nutation.m_eps * Constants::RAD_TO_DEG;
        assert!(m_eps_deg > 23.0 && m_eps_deg < 24.0);
    }

    #[test]
    fn test_sun_position() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let pos = Sun::eci_position(&epoch);

        // Distance should be approximately 1 AU
        let distance = pos.magnitude();
        assert!(distance > 0.98 * Sun::AU && distance < 1.02 * Sun::AU);
    }

    #[test]
    fn test_sun_direction() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let dir = Sun::direction(&epoch);

        // Direction should be a unit vector
        assert!(approx_eq(dir.magnitude(), 1.0, EPSILON));
    }

    #[test]
    fn test_moon_position() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let pos = Moon::eci_position(&epoch);

        // Distance should be approximately 384,400 km
        let distance = pos.magnitude();
        let distance_km = distance / 1000.0;
        assert!(distance_km > 350_000.0 && distance_km < 420_000.0);
    }

    #[test]
    fn test_moon_phase_angle() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let phase = Moon::phase_angle(&epoch);

        // Phase angle should be between 0 and π
        assert!(phase >= 0.0 && phase <= std::f64::consts::PI);
    }

    #[test]
    fn test_normalize_degrees() {
        assert!(approx_eq(normalize_degrees(0.0), 0.0, EPSILON));
        assert!(approx_eq(normalize_degrees(360.0), 0.0, EPSILON));
        assert!(approx_eq(normalize_degrees(-90.0), 270.0, EPSILON));
        assert!(approx_eq(normalize_degrees(450.0), 90.0, EPSILON));
    }
}
