//! State vectors and coordinate frame transformations.
//!
//! This module provides representations for orbital state vectors in various
//! coordinate frames, along with transformations between them:
//!
//! - **ECIState** - Earth-Centered Inertial J2000 frame
//! - **TEMEState** - True Equator Mean Equinox frame (SGP4 output)
//! - **ECEFState** - Earth-Centered Earth-Fixed frame
//! - **GeodeticState** - Geodetic coordinates (latitude, longitude, altitude)
//! - **KeplerianElements** - Classical orbital elements

use crate::bodies::Earth;
use crate::constants::Constants;
use crate::math::{Matrix3, Vector3};
use chrono::{DateTime, Utc};
use std::f64::consts::PI;

// ============================================================================
// ECIState - Earth-Centered Inertial J2000 Frame
// ============================================================================

/// State vector in Earth-Centered Inertial J2000 (ECI) reference frame.
///
/// ECI J2000 is an inertial reference frame with:
/// - Origin at Earth's center of mass
/// - X-axis pointing toward the mean vernal equinox at J2000.0
/// - Z-axis pointing toward the mean celestial pole at J2000.0
/// - Y-axis completing the right-handed system
///
/// Position is in kilometers and velocity is in km/s.
#[derive(Debug, Clone)]
pub struct ECIState {
    /// UTC epoch of this state
    pub epoch: DateTime<Utc>,
    /// Position vector in km
    pub position: Vector3,
    /// Velocity vector in km/s
    pub velocity: Vector3,
}

impl ECIState {
    /// Create a new ECI state.
    pub fn new(epoch: DateTime<Utc>, position: Vector3, velocity: Vector3) -> Self {
        Self {
            epoch,
            position,
            velocity,
        }
    }

    /// Create from position and velocity tuples.
    pub fn from_tuples(
        epoch: DateTime<Utc>,
        position: (f64, f64, f64),
        velocity: (f64, f64, f64),
    ) -> Self {
        Self::new(
            epoch,
            Vector3::from_tuple(position),
            Vector3::from_tuple(velocity),
        )
    }

    /// Get the position magnitude (distance from Earth center) in km.
    pub fn radius(&self) -> f64 {
        self.position.magnitude()
    }

    /// Get the velocity magnitude (speed) in km/s.
    pub fn speed(&self) -> f64 {
        self.velocity.magnitude()
    }

    /// Get the specific angular momentum vector (h = r × v) in km²/s.
    pub fn angular_momentum(&self) -> Vector3 {
        self.position.cross(&self.velocity)
    }

    /// Get the specific orbital energy in km²/s².
    pub fn specific_energy(&self) -> f64 {
        let v_sq = self.velocity.magnitude_squared();
        let r = self.radius();
        v_sq / 2.0 - Constants::EARTH_MU_KM / r
    }

    /// Convert to TEME frame.
    pub fn to_teme(&self) -> TEMEState {
        StateTransforms::eci_to_teme(self)
    }

    /// Convert to ECEF frame.
    pub fn to_ecef(&self) -> ECEFState {
        StateTransforms::eci_to_ecef(self)
    }

    /// Convert to Keplerian elements.
    pub fn to_keplerian(&self) -> KeplerianElements {
        StateTransforms::cartesian_to_keplerian(self)
    }
}

// ============================================================================
// TEMEState - True Equator Mean Equinox Frame
// ============================================================================

/// State vector in True Equator Mean Equinox (TEME) reference frame.
///
/// TEME is the reference frame used by SGP4/SDP4 propagators.
/// Position is in kilometers and velocity is in km/s.
#[derive(Debug, Clone)]
pub struct TEMEState {
    /// UTC epoch of this state
    pub epoch: DateTime<Utc>,
    /// Position vector in km
    pub position: Vector3,
    /// Velocity vector in km/s
    pub velocity: Vector3,
}

impl TEMEState {
    /// Create a new TEME state.
    pub fn new(epoch: DateTime<Utc>, position: Vector3, velocity: Vector3) -> Self {
        Self {
            epoch,
            position,
            velocity,
        }
    }

    /// Create from position and velocity tuples.
    pub fn from_tuples(
        epoch: DateTime<Utc>,
        position: (f64, f64, f64),
        velocity: (f64, f64, f64),
    ) -> Self {
        Self::new(
            epoch,
            Vector3::from_tuple(position),
            Vector3::from_tuple(velocity),
        )
    }

    /// Get the position magnitude (distance from Earth center) in km.
    pub fn radius(&self) -> f64 {
        self.position.magnitude()
    }

    /// Get the velocity magnitude (speed) in km/s.
    pub fn speed(&self) -> f64 {
        self.velocity.magnitude()
    }

    /// Convert to ECI J2000 frame.
    pub fn to_eci(&self) -> ECIState {
        StateTransforms::teme_to_eci(self)
    }

    /// Convert to ECEF frame.
    pub fn to_ecef(&self) -> ECEFState {
        StateTransforms::teme_to_ecef(self)
    }
}

// ============================================================================
// ECEFState - Earth-Centered Earth-Fixed Frame
// ============================================================================

/// State vector in Earth-Centered Earth-Fixed (ECEF) reference frame.
///
/// ECEF is a rotating reference frame fixed to the Earth with:
/// - Origin at Earth's center of mass
/// - X-axis pointing toward the intersection of the prime meridian and equator
/// - Z-axis pointing toward the North Pole
/// - Y-axis completing the right-handed system (90° East longitude)
///
/// Position is in kilometers and velocity is in km/s.
/// Note: Velocity in ECEF includes the Earth's rotation.
#[derive(Debug, Clone)]
pub struct ECEFState {
    /// UTC epoch of this state
    pub epoch: DateTime<Utc>,
    /// Position vector in km
    pub position: Vector3,
    /// Velocity vector in km/s
    pub velocity: Vector3,
}

impl ECEFState {
    /// Create a new ECEF state.
    pub fn new(epoch: DateTime<Utc>, position: Vector3, velocity: Vector3) -> Self {
        Self {
            epoch,
            position,
            velocity,
        }
    }

    /// Create from position and velocity tuples.
    pub fn from_tuples(
        epoch: DateTime<Utc>,
        position: (f64, f64, f64),
        velocity: (f64, f64, f64),
    ) -> Self {
        Self::new(
            epoch,
            Vector3::from_tuple(position),
            Vector3::from_tuple(velocity),
        )
    }

    /// Get the position magnitude (distance from Earth center) in km.
    pub fn radius(&self) -> f64 {
        self.position.magnitude()
    }

    /// Convert to ECI J2000 frame.
    pub fn to_eci(&self) -> ECIState {
        StateTransforms::ecef_to_eci(self)
    }

    /// Convert to geodetic coordinates.
    pub fn to_geodetic(&self) -> GeodeticState {
        GeodeticState::from_ecef(self)
    }
}

// ============================================================================
// GeodeticState - Geodetic Coordinates
// ============================================================================

/// State in geodetic coordinates (latitude, longitude, altitude).
///
/// Uses WGS84 ellipsoid for Earth reference. This is the standard
/// representation for ground-based observer locations.
///
/// - Latitude: Geodetic latitude in degrees (-90 to +90, positive North)
/// - Longitude: Geodetic longitude in degrees (-180 to +180, positive East)
/// - Altitude: Height above WGS84 ellipsoid in kilometers
#[derive(Debug, Clone, Copy)]
pub struct GeodeticState {
    /// Geodetic latitude in degrees (-90 to +90)
    pub latitude: f64,
    /// Geodetic longitude in degrees (-180 to +180)
    pub longitude: f64,
    /// Altitude above WGS84 ellipsoid in kilometers
    pub altitude: f64,
}

impl GeodeticState {
    /// Create a new geodetic state.
    ///
    /// # Arguments
    /// * `latitude` - Geodetic latitude in degrees (-90 to +90)
    /// * `longitude` - Geodetic longitude in degrees (-180 to +180)
    /// * `altitude` - Height above WGS84 ellipsoid in kilometers
    pub fn new(latitude: f64, longitude: f64, altitude: f64) -> Self {
        Self {
            latitude,
            longitude,
            altitude,
        }
    }

    /// Get latitude in radians.
    pub fn latitude_rad(&self) -> f64 {
        self.latitude * Constants::DEG_TO_RAD
    }

    /// Get longitude in radians.
    pub fn longitude_rad(&self) -> f64 {
        self.longitude * Constants::DEG_TO_RAD
    }

    /// Convert to ECEF position at a given epoch.
    ///
    /// Velocity is set to zero (stationary observer).
    pub fn to_ecef(&self, epoch: DateTime<Utc>) -> ECEFState {
        let lat_rad = self.latitude_rad();
        let lon_rad = self.longitude_rad();

        let sin_lat = lat_rad.sin();
        let cos_lat = lat_rad.cos();
        let sin_lon = lon_rad.sin();
        let cos_lon = lon_rad.cos();

        // Radius of curvature in the prime vertical
        let n =
            Constants::EARTH_RADIUS_EQ_KM / (1.0 - Constants::EARTH_E2 * sin_lat * sin_lat).sqrt();

        // ECEF position
        let x = (n + self.altitude) * cos_lat * cos_lon;
        let y = (n + self.altitude) * cos_lat * sin_lon;
        let z = (n * (1.0 - Constants::EARTH_E2) + self.altitude) * sin_lat;

        // Stationary observer - velocity is zero in ECEF
        ECEFState::new(epoch, Vector3::new(x, y, z), Vector3::zero())
    }

    /// Convert to ECI position at a given epoch.
    pub fn to_eci(&self, epoch: DateTime<Utc>) -> ECIState {
        let ecef = self.to_ecef(epoch);
        ecef.to_eci()
    }

    /// Create geodetic state from ECEF position.
    ///
    /// Uses iterative algorithm for accurate conversion.
    pub fn from_ecef(ecef: &ECEFState) -> Self {
        let x = ecef.position.x;
        let y = ecef.position.y;
        let z = ecef.position.z;

        // Longitude is straightforward
        let lon = y.atan2(x) * Constants::RAD_TO_DEG;

        // Iterative calculation for latitude
        let p = (x * x + y * y).sqrt();
        let mut lat = (z / (p * (1.0 - Constants::EARTH_E2))).atan();

        // Iterate to convergence
        for _ in 0..10 {
            let sin_lat = lat.sin();
            let n = Constants::EARTH_RADIUS_EQ_KM
                / (1.0 - Constants::EARTH_E2 * sin_lat * sin_lat).sqrt();
            let new_lat = (z + Constants::EARTH_E2 * n * sin_lat).atan2(p);

            if (new_lat - lat).abs() < 1e-12 {
                break;
            }
            lat = new_lat;
        }

        // Calculate altitude
        let sin_lat = lat.sin();
        let cos_lat = lat.cos();
        let n =
            Constants::EARTH_RADIUS_EQ_KM / (1.0 - Constants::EARTH_E2 * sin_lat * sin_lat).sqrt();

        let alt = if cos_lat.abs() > 1e-10 {
            p / cos_lat - n
        } else {
            z.abs() / sin_lat.abs() - n * (1.0 - Constants::EARTH_E2)
        };

        Self::new(lat * Constants::RAD_TO_DEG, lon, alt)
    }

    /// Get local sidereal time at this location.
    ///
    /// Returns LST in radians.
    pub fn local_sidereal_time(&self, epoch: &DateTime<Utc>) -> f64 {
        let nutation = Earth::nutation_angles(epoch);
        let gast = nutation.gast;
        let lon_rad = self.longitude_rad();

        // Local sidereal time = GAST + observer longitude
        let lst = gast + lon_rad;

        // Normalize to [0, 2π]
        let mut normalized = lst % Constants::TWO_PI;
        if normalized < 0.0 {
            normalized += Constants::TWO_PI;
        }
        normalized
    }
}

impl Default for GeodeticState {
    fn default() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }
}

// ============================================================================
// KeplerianElements - Classical Orbital Elements
// ============================================================================

/// Classical Keplerian orbital elements.
///
/// Represents an orbit using the six classical orbital elements:
/// - Semi-major axis (a): Size of the orbit in km
/// - Eccentricity (e): Shape of the orbit (0 = circular, 0-1 = elliptical)
/// - Inclination (i): Tilt relative to the equatorial plane in radians
/// - RAAN (Ω): Right Ascension of the Ascending Node in radians
/// - Argument of perigee (ω): Orientation of the orbit in its plane in radians
/// - True anomaly (ν): Position along the orbit in radians
#[derive(Debug, Clone, Copy)]
pub struct KeplerianElements {
    /// UTC epoch of these elements
    pub epoch: Option<DateTime<Utc>>,
    /// Semi-major axis in kilometers
    pub semi_major_axis: f64,
    /// Orbital eccentricity (0 to < 1 for elliptical)
    pub eccentricity: f64,
    /// Inclination in radians
    pub inclination: f64,
    /// Right Ascension of Ascending Node in radians
    pub raan: f64,
    /// Argument of perigee in radians
    pub arg_perigee: f64,
    /// True anomaly in radians
    pub true_anomaly: f64,
    /// Gravitational parameter (defaults to Earth)
    pub mu: f64,
}

impl KeplerianElements {
    /// Create new Keplerian elements.
    ///
    /// # Arguments
    /// * `semi_major_axis` - Semi-major axis in kilometers
    /// * `eccentricity` - Orbital eccentricity (0 to < 1 for elliptical)
    /// * `inclination` - Inclination in radians
    /// * `raan` - Right Ascension of Ascending Node in radians
    /// * `arg_perigee` - Argument of perigee in radians
    /// * `true_anomaly` - True anomaly in radians
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        arg_perigee: f64,
        true_anomaly: f64,
    ) -> Self {
        Self {
            epoch: None,
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            arg_perigee,
            true_anomaly,
            mu: Constants::EARTH_MU_KM,
        }
    }

    /// Create with a specific epoch.
    pub fn with_epoch(mut self, epoch: DateTime<Utc>) -> Self {
        self.epoch = Some(epoch);
        self
    }

    /// Create with a specific gravitational parameter.
    pub fn with_mu(mut self, mu: f64) -> Self {
        self.mu = mu;
        self
    }

    /// Calculate orbital period in seconds.
    pub fn period(&self) -> f64 {
        2.0 * PI * (self.semi_major_axis.powi(3) / self.mu).sqrt()
    }

    /// Calculate mean motion in radians per second.
    pub fn mean_motion(&self) -> f64 {
        (self.mu / self.semi_major_axis.powi(3)).sqrt()
    }

    /// Calculate apoapsis radius in kilometers.
    pub fn apoapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 + self.eccentricity)
    }

    /// Calculate periapsis radius in kilometers.
    pub fn periapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity)
    }

    /// Calculate apoapsis altitude above Earth's surface in kilometers.
    pub fn apoapsis_altitude(&self) -> f64 {
        self.apoapsis() - Constants::EARTH_RADIUS_EQ_KM
    }

    /// Calculate periapsis altitude above Earth's surface in kilometers.
    pub fn periapsis_altitude(&self) -> f64 {
        self.periapsis() - Constants::EARTH_RADIUS_EQ_KM
    }

    /// Calculate specific orbital energy in km²/s².
    pub fn specific_energy(&self) -> f64 {
        -self.mu / (2.0 * self.semi_major_axis)
    }

    /// Calculate specific angular momentum magnitude in km²/s.
    pub fn angular_momentum(&self) -> f64 {
        (self.mu * self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity)).sqrt()
    }

    /// Calculate semi-latus rectum (p) in km.
    pub fn semi_latus_rectum(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity)
    }

    /// Convert true anomaly to eccentric anomaly.
    pub fn eccentric_anomaly(&self) -> f64 {
        let e = self.eccentricity;
        let nu = self.true_anomaly;
        ((1.0 - e * e).sqrt() * nu.sin()).atan2(e + nu.cos())
    }

    /// Convert true anomaly to mean anomaly.
    pub fn mean_anomaly(&self) -> f64 {
        let ea = self.eccentric_anomaly();
        ea - self.eccentricity * ea.sin()
    }

    /// Get inclination in degrees.
    pub fn inclination_deg(&self) -> f64 {
        self.inclination * Constants::RAD_TO_DEG
    }

    /// Get RAAN in degrees.
    pub fn raan_deg(&self) -> f64 {
        self.raan * Constants::RAD_TO_DEG
    }

    /// Get argument of perigee in degrees.
    pub fn arg_perigee_deg(&self) -> f64 {
        self.arg_perigee * Constants::RAD_TO_DEG
    }

    /// Get true anomaly in degrees.
    pub fn true_anomaly_deg(&self) -> f64 {
        self.true_anomaly * Constants::RAD_TO_DEG
    }

    /// Convert to ECI state vector.
    pub fn to_eci(&self) -> ECIState {
        StateTransforms::keplerian_to_cartesian(self)
    }
}

// ============================================================================
// StateTransforms - Coordinate Frame Transformations
// ============================================================================

/// Coordinate frame transformations.
///
/// Provides transformations between different reference frames and
/// conversions between Cartesian and Keplerian representations.
pub struct StateTransforms;

impl StateTransforms {
    // ========================================================================
    // TEME <-> ECI J2000 Transformations
    // ========================================================================

    /// Convert TEME state to ECI J2000.
    ///
    /// This transformation accounts for precession, nutation, and the
    /// equation of equinoxes.
    pub fn teme_to_eci(teme: &TEMEState) -> ECIState {
        let precession = Earth::precession_angles(&teme.epoch);
        let nutation = Earth::nutation_angles(&teme.epoch);

        // Build transformation matrices
        let eq_mat = Self::equinox_matrix(nutation.d_psi, nutation.eps);
        let nut_mat = Self::nutation_matrix(nutation.m_eps, nutation.d_psi, nutation.d_eps);
        let prec_mat = Self::precession_matrix(precession.zeta, precession.theta, precession.z);

        // Combined transformation: P * N * E
        let combined = prec_mat.mul_mat(&nut_mat.mul_mat(&eq_mat));

        // Apply to position and velocity
        let pos_eci = combined.mul_vec(&teme.position);
        let vel_eci = combined.mul_vec(&teme.velocity);

        ECIState::new(teme.epoch, pos_eci, vel_eci)
    }

    /// Convert ECI J2000 state to TEME.
    pub fn eci_to_teme(eci: &ECIState) -> TEMEState {
        let precession = Earth::precession_angles(&eci.epoch);
        let nutation = Earth::nutation_angles(&eci.epoch);

        // Build transformation matrices (transposed for inverse)
        let eq_mat = Self::equinox_matrix(nutation.d_psi, nutation.eps).transpose();
        let nut_mat =
            Self::nutation_matrix(nutation.m_eps, nutation.d_psi, nutation.d_eps).transpose();
        let prec_mat =
            Self::precession_matrix(precession.zeta, precession.theta, precession.z).transpose();

        // Combined transformation: E^T * N^T * P^T
        let combined = eq_mat.mul_mat(&nut_mat.mul_mat(&prec_mat));

        // Apply to position and velocity
        let pos_teme = combined.mul_vec(&eci.position);
        let vel_teme = combined.mul_vec(&eci.velocity);

        TEMEState::new(eci.epoch, pos_teme, vel_teme)
    }

    // ========================================================================
    // ECI <-> ECEF Transformations
    // ========================================================================

    /// Convert ECI J2000 state to ECEF.
    pub fn eci_to_ecef(eci: &ECIState) -> ECEFState {
        let nutation = Earth::nutation_angles(&eci.epoch);
        let gast = nutation.gast;

        // Rotation matrix from ECI to ECEF
        let rot = Matrix3::rotation_z(gast);

        // Position transformation
        let pos_ecef = rot.mul_vec(&eci.position);

        // Velocity needs correction for Earth's rotation
        // v_ecef = R * v_eci - omega x r_ecef
        let vel_rotated = rot.mul_vec(&eci.velocity);

        let omega = Constants::EARTH_ROTATION_RATE;
        let omega_cross_r = Vector3::new(-omega * pos_ecef.y, omega * pos_ecef.x, 0.0);

        let vel_ecef = vel_rotated - omega_cross_r;

        ECEFState::new(eci.epoch, pos_ecef, vel_ecef)
    }

    /// Convert ECEF state to ECI J2000.
    pub fn ecef_to_eci(ecef: &ECEFState) -> ECIState {
        let nutation = Earth::nutation_angles(&ecef.epoch);
        let gast = nutation.gast;

        // Inverse rotation matrix
        let rot_inv = Matrix3::rotation_z(-gast);

        let omega = Constants::EARTH_ROTATION_RATE;

        // First correct velocity for Earth's rotation
        // v_eci = R^T * (v_ecef + omega x r_ecef)
        let omega_cross_r = Vector3::new(-omega * ecef.position.y, omega * ecef.position.x, 0.0);
        let vel_corrected = ecef.velocity + omega_cross_r;

        // Apply rotation
        let pos_eci = rot_inv.mul_vec(&ecef.position);
        let vel_eci = rot_inv.mul_vec(&vel_corrected);

        ECIState::new(ecef.epoch, pos_eci, vel_eci)
    }

    // ========================================================================
    // TEME <-> ECEF (via ECI)
    // ========================================================================

    /// Convert TEME state to ECEF.
    pub fn teme_to_ecef(teme: &TEMEState) -> ECEFState {
        let eci = Self::teme_to_eci(teme);
        Self::eci_to_ecef(&eci)
    }

    /// Convert ECEF state to TEME.
    pub fn ecef_to_teme(ecef: &ECEFState) -> TEMEState {
        let eci = Self::ecef_to_eci(ecef);
        Self::eci_to_teme(&eci)
    }

    // ========================================================================
    // Cartesian <-> Keplerian Conversions
    // ========================================================================

    /// Convert Cartesian state to Keplerian elements.
    pub fn cartesian_to_keplerian(state: &ECIState) -> KeplerianElements {
        let mu = Constants::EARTH_MU_KM;
        let pos = &state.position;
        let vel = &state.velocity;

        let r = pos.magnitude();
        let v = vel.magnitude();

        // Specific angular momentum h = r x v
        let h = pos.cross(vel);
        let h_mag = h.magnitude();

        // Node vector n = k x h (k is z-unit vector)
        let n = Vector3::new(-h.y, h.x, 0.0);
        let n_mag = n.magnitude();

        // Eccentricity vector
        let r_dot_v = pos.dot(vel);
        let e_vec = (*pos * (v * v - mu / r) - *vel * r_dot_v) / mu;
        let e = e_vec.magnitude();

        // Semi-major axis
        let energy = v * v / 2.0 - mu / r;
        let a = -mu / (2.0 * energy);

        // Inclination
        let i = (h.z / h_mag).clamp(-1.0, 1.0).acos();

        // RAAN
        let raan = if n_mag > 1e-10 {
            let raan_cos = (n.x / n_mag).clamp(-1.0, 1.0);
            if n.y >= 0.0 {
                raan_cos.acos()
            } else {
                Constants::TWO_PI - raan_cos.acos()
            }
        } else {
            0.0
        };

        // Argument of perigee
        let w = if n_mag > 1e-10 && e > 1e-10 {
            let n_dot_e = n.dot(&e_vec);
            let w_cos = (n_dot_e / (n_mag * e)).clamp(-1.0, 1.0);
            if e_vec.z >= 0.0 {
                w_cos.acos()
            } else {
                Constants::TWO_PI - w_cos.acos()
            }
        } else {
            0.0
        };

        // True anomaly
        let nu = if e > 1e-10 {
            let e_dot_r = e_vec.dot(pos);
            let nu_cos = (e_dot_r / (e * r)).clamp(-1.0, 1.0);
            if r_dot_v >= 0.0 {
                nu_cos.acos()
            } else {
                Constants::TWO_PI - nu_cos.acos()
            }
        } else {
            0.0
        };

        KeplerianElements::new(a, e, i, raan, w, nu)
            .with_epoch(state.epoch)
            .with_mu(mu)
    }

    /// Convert Keplerian elements to Cartesian state.
    pub fn keplerian_to_cartesian(elements: &KeplerianElements) -> ECIState {
        let a = elements.semi_major_axis;
        let e = elements.eccentricity;
        let i = elements.inclination;
        let raan = elements.raan;
        let w = elements.arg_perigee;
        let nu = elements.true_anomaly;
        let mu = elements.mu;

        // Semi-latus rectum
        let p = a * (1.0 - e * e);

        // Position and velocity in perifocal frame
        let cos_nu = nu.cos();
        let sin_nu = nu.sin();

        let r_mag = p / (1.0 + e * cos_nu);

        // Position in perifocal frame
        let r_pf = Vector3::new(r_mag * cos_nu, r_mag * sin_nu, 0.0);

        // Velocity in perifocal frame
        let sqrt_mu_p = (mu / p).sqrt();
        let v_pf = Vector3::new(-sqrt_mu_p * sin_nu, sqrt_mu_p * (e + cos_nu), 0.0);

        // Rotation from perifocal to ECI
        let cos_raan = raan.cos();
        let sin_raan = raan.sin();
        let cos_w = w.cos();
        let sin_w = w.sin();
        let cos_i = i.cos();
        let sin_i = i.sin();

        // Combined rotation matrix from perifocal to inertial
        let r11 = cos_raan * cos_w - sin_raan * sin_w * cos_i;
        let r12 = -cos_raan * sin_w - sin_raan * cos_w * cos_i;
        let r13 = sin_raan * sin_i;
        let r21 = sin_raan * cos_w + cos_raan * sin_w * cos_i;
        let r22 = -sin_raan * sin_w + cos_raan * cos_w * cos_i;
        let r23 = -cos_raan * sin_i;
        let r31 = sin_w * sin_i;
        let r32 = cos_w * sin_i;
        let r33 = cos_i;

        let rot = Matrix3::new(r11, r12, r13, r21, r22, r23, r31, r32, r33);

        let pos = rot.mul_vec(&r_pf);
        let vel = rot.mul_vec(&v_pf);

        let epoch = elements.epoch.unwrap_or_else(Utc::now);
        ECIState::new(epoch, pos, vel)
    }

    // ========================================================================
    // Helper Matrix Functions
    // ========================================================================

    /// Build precession matrix from precession angles.
    fn precession_matrix(zeta: f64, theta: f64, z: f64) -> Matrix3 {
        // P = Rz(-z) * Ry(theta) * Rz(-zeta)
        let rz_neg_zeta = Matrix3::rotation_z(-zeta);
        let ry_theta = Matrix3::rotation_y(theta);
        let rz_neg_z = Matrix3::rotation_z(-z);

        rz_neg_z.mul_mat(&ry_theta.mul_mat(&rz_neg_zeta))
    }

    /// Build nutation matrix from nutation angles.
    fn nutation_matrix(mean_eps: f64, delta_psi: f64, delta_eps: f64) -> Matrix3 {
        let eps = mean_eps + delta_eps;

        // N = Rx(mean_eps) * Rz(delta_psi) * Rx(-eps)
        let rx_mean_eps = Matrix3::rotation_x(mean_eps);
        let rz_delta_psi = Matrix3::rotation_z(delta_psi);
        let rx_neg_eps = Matrix3::rotation_x(-eps);

        rx_mean_eps.mul_mat(&rz_delta_psi.mul_mat(&rx_neg_eps))
    }

    /// Build equation of equinoxes rotation matrix.
    fn equinox_matrix(delta_psi: f64, eps: f64) -> Matrix3 {
        // The equation of equinoxes term: delta_psi * cos(eps)
        let eq_eq = delta_psi * eps.cos();
        Matrix3::rotation_z(-eq_eq)
    }
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
    fn test_eci_state_creation() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();
        let pos = Vector3::new(7000.0, 0.0, 0.0);
        let vel = Vector3::new(0.0, 7.5, 0.0);

        let state = ECIState::new(epoch, pos, vel);

        assert!(approx_eq(state.radius(), 7000.0, EPSILON));
        assert!(approx_eq(state.speed(), 7.5, EPSILON));
    }

    #[test]
    fn test_geodetic_to_ecef() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        // Test at equator, prime meridian
        let geo = GeodeticState::new(0.0, 0.0, 0.0);
        let ecef = geo.to_ecef(epoch);

        // Should be on the positive X axis at Earth's equatorial radius
        assert!(approx_eq(
            ecef.position.x,
            Constants::EARTH_RADIUS_EQ_KM,
            1e-3
        ));
        assert!(approx_eq(ecef.position.y, 0.0, 1e-3));
        assert!(approx_eq(ecef.position.z, 0.0, 1e-3));
    }

    #[test]
    fn test_geodetic_roundtrip() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        let geo_original = GeodeticState::new(40.0, -105.0, 1.6);
        let ecef = geo_original.to_ecef(epoch);
        let geo_back = GeodeticState::from_ecef(&ecef);

        assert!(approx_eq(geo_original.latitude, geo_back.latitude, 1e-6));
        assert!(approx_eq(geo_original.longitude, geo_back.longitude, 1e-6));
        assert!(approx_eq(geo_original.altitude, geo_back.altitude, 1e-6));
    }

    #[test]
    fn test_keplerian_period() {
        // Create a circular orbit at ~400km altitude (ISS-like)
        let a = Constants::EARTH_RADIUS_EQ_KM + 400.0;
        let elements = KeplerianElements::new(a, 0.0, 0.9, 0.0, 0.0, 0.0);

        // Period should be approximately 92 minutes = 5520 seconds
        let period = elements.period();
        assert!(period > 5400.0 && period < 5700.0);
    }

    #[test]
    fn test_keplerian_apoapsis_periapsis() {
        let a = 10000.0;
        let e = 0.1;
        let elements = KeplerianElements::new(a, e, 0.0, 0.0, 0.0, 0.0);

        assert!(approx_eq(elements.apoapsis(), a * (1.0 + e), EPSILON));
        assert!(approx_eq(elements.periapsis(), a * (1.0 - e), EPSILON));
    }

    #[test]
    fn test_cartesian_keplerian_roundtrip() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        // Create a state in ECI - use a non-degenerate orbit (not in equatorial plane)
        let pos = Vector3::new(6500.0, 1500.0, 1000.0);
        let vel = Vector3::new(-1.0, 7.0, 0.5);
        let eci = ECIState::new(epoch, pos, vel);

        // Convert to Keplerian and back
        let kep = StateTransforms::cartesian_to_keplerian(&eci);
        let eci_back = StateTransforms::keplerian_to_cartesian(&kep);

        // Check roundtrip accuracy (relaxed tolerance for numerical precision)
        assert!(approx_eq(eci.position.x, eci_back.position.x, 1.0));
        assert!(approx_eq(eci.position.y, eci_back.position.y, 1.0));
        assert!(approx_eq(eci.position.z, eci_back.position.z, 1.0));
        assert!(approx_eq(eci.velocity.x, eci_back.velocity.x, 1e-3));
        assert!(approx_eq(eci.velocity.y, eci_back.velocity.y, 1e-3));
        assert!(approx_eq(eci.velocity.z, eci_back.velocity.z, 1e-3));
    }

    #[test]
    fn test_ecef_eci_roundtrip() {
        let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

        let pos = Vector3::new(7000.0, 1000.0, 500.0);
        let vel = Vector3::new(-1.0, 7.0, 0.5);
        let eci = ECIState::new(epoch, pos, vel);

        // Convert to ECEF and back
        let ecef = StateTransforms::eci_to_ecef(&eci);
        let eci_back = StateTransforms::ecef_to_eci(&ecef);

        // Check roundtrip accuracy
        assert!(approx_eq(eci.position.x, eci_back.position.x, 1e-6));
        assert!(approx_eq(eci.position.y, eci_back.position.y, 1e-6));
        assert!(approx_eq(eci.position.z, eci_back.position.z, 1e-6));
        assert!(approx_eq(eci.velocity.x, eci_back.velocity.x, 1e-6));
        assert!(approx_eq(eci.velocity.y, eci_back.velocity.y, 1e-6));
        assert!(approx_eq(eci.velocity.z, eci_back.velocity.z, 1e-6));
    }

    #[test]
    fn test_keplerian_mean_anomaly() {
        // Circular orbit (e=0), mean anomaly should equal true anomaly
        let kep = KeplerianElements::new(7000.0, 0.0, 0.5, 0.0, 0.0, 1.0);
        assert!(approx_eq(kep.mean_anomaly(), 1.0, EPSILON));
    }

    #[test]
    fn test_keplerian_angular_momentum() {
        let a = 10000.0;
        let e = 0.1;
        let kep = KeplerianElements::new(a, e, 0.0, 0.0, 0.0, 0.0);

        // h = sqrt(mu * a * (1 - e²))
        let expected = (Constants::EARTH_MU_KM * a * (1.0 - e * e)).sqrt();
        assert!(approx_eq(kep.angular_momentum(), expected, EPSILON));
    }
}
