//! # SpaceDust
//!
//! A comprehensive astrodynamics library for Rust, providing tools for satellite tracking,
//! orbital mechanics, coordinate transformations, and ground-based observation calculations.
//!
//! ## Features
//!
//! - **Time Systems**: Conversions between UTC, TAI, TT, Julian Date, GPS, and GMST
//! - **Coordinate Frames**: Transformations between ECI J2000, TEME, ECEF, and Geodetic
//! - **Orbital Elements**: Conversions between Cartesian state vectors and Keplerian elements
//! - **TLE Parsing**: Parse and propagate Two-Line Element Sets using SGP4
//! - **Celestial Bodies**: Sun and Moon position calculations
//! - **Ground Observations**: Azimuth/Elevation and Right Ascension/Declination calculations
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use space_dust::tle::Tle;
//! use space_dust::state::{TEMEState, GeodeticState, StateTransforms};
//! use space_dust::observations::Observations;
//! use chrono::Utc;
//!
//! // Parse a TLE
//! let line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9002";
//! let line2 = "2 25544  51.6400 208.1200 0001234  85.0000 275.0000 15.48919100123456";
//! let tle = Tle::parse(line1, line2).unwrap();
//!
//! // Propagate to current time
//! let epoch = Utc::now();
//! let teme_state = tle.propagate(&epoch).unwrap();
//!
//! // Convert to ECI J2000
//! let eci_state = teme_state.to_eci();
//!
//! // Define a ground observer (Denver, CO)
//! let observer = GeodeticState::new(39.7392, -104.9903, 1.6);
//!
//! // Compute observation angles
//! let az_el = Observations::compute_az_el(&observer, &eci_state);
//! println!("Azimuth: {:.2}°, Elevation: {:.2}°", az_el.azimuth_deg(), az_el.elevation_deg());
//! ```

pub mod bodies;
pub mod constants;
pub mod data;
pub mod math;
pub mod observations;
pub mod state;
pub mod time;
pub mod tle;

// Re-export commonly used types
pub use bodies::{Earth, Moon, Sun};
pub use constants::Constants;
pub use math::{Matrix3, Vector3};
pub use observations::{AzEl, Observations, RaDec};
pub use state::{
    ECEFState, ECIState, GeodeticState, KeplerianElements, StateTransforms, TEMEState,
};
pub use time::{JulianDate, TimeTransforms, GMST, GPS, TAI, TT, UTC};
pub use tle::Tle;

/// Library version
pub const VERSION: &str = "0.1.0";

/// Returns the library version
pub fn version() -> &'static str {
    VERSION
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        assert_eq!(version(), "0.1.0");
    }
}
