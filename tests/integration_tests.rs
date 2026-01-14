//! Integration tests for the space-dust astrodynamics library.
//!
//! These tests verify the complete functionality of the library by
//! testing realistic use cases and end-to-end workflows.

use chrono::{TimeZone, Utc};
use space_dust::bodies::{Earth, Moon, Sun};
use space_dust::constants::Constants;
use space_dust::math::Vector3;
use space_dust::observations::{AzEl, Observations};
use space_dust::state::{ECIState, GeodeticState, KeplerianElements};
use space_dust::time::{JulianDate, TimeTransforms, UTC};
use space_dust::tle::Tle;

fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

// ============================================================================
// TLE and SGP4 Propagation Tests
// ============================================================================

#[test]
fn test_tle_full_workflow() {
    // Parse a real ISS TLE
    let line1 = "1 25544U 98067A   26014.62805721  .00006818  00000+0  13044-3 0  9991";
    let line2 = "2 25544  51.6333 339.6562 0007763  17.9854 342.1408 15.49289811547943";

    let tle = Tle::parse(line1, line2).expect("Failed to parse TLE");

    // Verify parsed values
    assert_eq!(tle.catalog_number(), "25544");
    assert!(tle.inclination_deg() > 51.0 && tle.inclination_deg() < 52.0);
    assert!(tle.eccentricity() < 0.01); // Nearly circular

    // Propagate to TLE epoch
    let state_at_epoch = tle.propagate_minutes(0.0).expect("Propagation failed");

    // ISS should be in LEO (~400km altitude, ~6700-6800km from Earth center)
    let radius = state_at_epoch.position.magnitude();
    assert!(
        radius > 6500.0 && radius < 7000.0,
        "Radius {} not in expected range",
        radius
    );

    // Velocity should be ~7.5 km/s for LEO
    let speed = state_at_epoch.velocity.magnitude();
    assert!(
        speed > 7.0 && speed < 8.0,
        "Speed {} not in expected range",
        speed
    );

    // Propagate forward one orbit
    let period_min = tle.period_minutes();
    let state_after_orbit = tle
        .propagate_minutes(period_min)
        .expect("Propagation failed");

    // After one orbit, radius should still be similar (LEO)
    let radius_after = state_after_orbit.position.magnitude();
    assert!(radius_after > 6500.0 && radius_after < 7000.0);
}

#[test]
fn test_tle_to_eci_conversion() {
    let line1 = "1 25544U 98067A   26014.62805721  .00006818  00000+0  13044-3 0  9991";
    let line2 = "2 25544  51.6333 339.6562 0007763  17.9854 342.1408 15.49289811547943";

    let tle = Tle::parse(line1, line2).unwrap();
    let teme_state = tle.propagate_minutes(0.0).unwrap();

    // Convert TEME to ECI J2000
    let eci_state = teme_state.to_eci();

    // Both should have similar magnitudes (rotation doesn't change magnitude)
    let teme_radius = teme_state.position.magnitude();
    let eci_radius = eci_state.position.magnitude();
    assert!(
        approx_eq(teme_radius, eci_radius, 1.0),
        "Radius changed from {} to {} during transformation",
        teme_radius,
        eci_radius
    );
}

// ============================================================================
// Ground Station Observation Tests
// ============================================================================

#[test]
fn test_satellite_visibility_from_ground() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 15, 12, 0, 0).unwrap();

    // Define a ground station (Boulder, CO)
    let observer = GeodeticState::new(40.0150, -105.2705, 1.655);

    // Create a satellite in a high orbit that should be visible
    // Geostationary satellite over the Americas
    let geo_distance = 42164.0; // km from Earth center
    let sat_position = Vector3::new(geo_distance, 0.0, 0.0);
    let sat_velocity = Vector3::new(0.0, 3.075, 0.0); // ~3 km/s for GEO

    let satellite = ECIState::new(epoch, sat_position, sat_velocity);

    // Compute observation angles
    let az_el = Observations::compute_az_el(&observer, &satellite);

    // Should be able to compute valid angles
    assert!(az_el.range.unwrap() > 0.0);
    assert!(az_el.azimuth >= 0.0 && az_el.azimuth < Constants::TWO_PI);
}

#[test]
fn test_az_el_ra_dec_conversion() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 15, 12, 0, 0).unwrap();

    // Observer at a specific location
    let observer = GeodeticState::new(45.0, -90.0, 0.0);

    // Create an Az/El observation
    let az_el = AzEl::from_degrees(epoch, 180.0, 45.0);

    // Convert to RA/Dec
    let ra_dec = Observations::az_el_to_ra_dec(&az_el, &observer);

    // Convert back to Az/El
    let az_el_back = Observations::ra_dec_to_az_el(&ra_dec, &observer);

    // Should be approximately equal (within small tolerance for numerical precision)
    assert!(
        approx_eq(az_el.azimuth, az_el_back.azimuth, 0.01),
        "Azimuth: {} vs {}",
        az_el.azimuth,
        az_el_back.azimuth
    );
    assert!(
        approx_eq(az_el.elevation, az_el_back.elevation, 0.01),
        "Elevation: {} vs {}",
        az_el.elevation,
        az_el_back.elevation
    );
}

// ============================================================================
// Coordinate Frame Transformation Tests
// ============================================================================

#[test]
fn test_eci_ecef_teme_roundtrip() {
    let epoch = Utc.with_ymd_and_hms(2024, 3, 20, 0, 0, 0).unwrap();

    // Create an ECI state
    let original_pos = Vector3::new(6500.0, 1500.0, 2000.0);
    let original_vel = Vector3::new(-1.5, 6.8, 0.8);
    let eci = ECIState::new(epoch, original_pos, original_vel);

    // ECI -> ECEF -> ECI roundtrip
    let ecef = eci.to_ecef();
    let eci_back = ecef.to_eci();

    assert!(
        approx_eq(eci.position.x, eci_back.position.x, 1e-6),
        "X position: {} vs {}",
        eci.position.x,
        eci_back.position.x
    );
    assert!(
        approx_eq(eci.position.y, eci_back.position.y, 1e-6),
        "Y position: {} vs {}",
        eci.position.y,
        eci_back.position.y
    );
    assert!(
        approx_eq(eci.position.z, eci_back.position.z, 1e-6),
        "Z position: {} vs {}",
        eci.position.z,
        eci_back.position.z
    );

    // ECI -> TEME -> ECI roundtrip
    let teme = eci.to_teme();
    let eci_from_teme = teme.to_eci();

    assert!(
        approx_eq(eci.position.x, eci_from_teme.position.x, 1e-3),
        "TEME roundtrip X: {} vs {}",
        eci.position.x,
        eci_from_teme.position.x
    );
}

#[test]
fn test_geodetic_ecef_roundtrip() {
    let epoch = Utc.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();

    // Test various locations
    let test_locations = vec![
        (0.0, 0.0, 0.0),           // Null Island (equator, prime meridian)
        (90.0, 0.0, 0.0),          // North Pole
        (-90.0, 0.0, 0.0),         // South Pole
        (40.7128, -74.0060, 0.01), // New York City
        (-33.8688, 151.2093, 0.0), // Sydney
        (35.6762, 139.6503, 0.0),  // Tokyo
    ];

    for (lat, lon, alt) in test_locations {
        let geo = GeodeticState::new(lat, lon, alt);
        let ecef = geo.to_ecef(epoch);
        let geo_back = GeodeticState::from_ecef(&ecef);

        assert!(
            approx_eq(lat, geo_back.latitude, 1e-6),
            "Lat: {} vs {} for ({}, {}, {})",
            lat,
            geo_back.latitude,
            lat,
            lon,
            alt
        );
        assert!(
            approx_eq(lon, geo_back.longitude, 1e-6),
            "Lon: {} vs {} for ({}, {}, {})",
            lon,
            geo_back.longitude,
            lat,
            lon,
            alt
        );
        assert!(
            approx_eq(alt, geo_back.altitude, 1e-4),
            "Alt: {} vs {} for ({}, {}, {})",
            alt,
            geo_back.altitude,
            lat,
            lon,
            alt
        );
    }
}

// ============================================================================
// Keplerian Elements Tests
// ============================================================================

#[test]
fn test_keplerian_circular_orbit() {
    // Create a circular orbit
    let a = 7000.0; // 7000 km semi-major axis (~620 km altitude)
    let elements = KeplerianElements::new(
        a, 0.0, // circular
        0.5, // 28.6 degree inclination
        1.0, // RAAN
        0.0, // argument of perigee (undefined for circular)
        0.0, // true anomaly
    );

    // Period should be calculable
    let period = elements.period();
    assert!(period > 0.0);

    // For circular orbit, apoapsis = periapsis = semi-major axis
    assert!(approx_eq(elements.apoapsis(), a, 1e-6));
    assert!(approx_eq(elements.periapsis(), a, 1e-6));
}

#[test]
fn test_keplerian_elliptical_orbit() {
    // Create an elliptical orbit (like a Molniya orbit)
    let a = 26600.0; // km
    let e = 0.74; // highly elliptical
    let i = 1.1; // ~63 degrees inclination
    let raan = 0.5; // radians
    let w = 4.7; // argument of perigee
    let nu = 0.0; // true anomaly

    let elements = KeplerianElements::new(a, e, i, raan, w, nu);

    // Check apoapsis and periapsis
    let expected_apoapsis = a * (1.0 + e);
    let expected_periapsis = a * (1.0 - e);

    assert!(approx_eq(elements.apoapsis(), expected_apoapsis, 1e-6));
    assert!(approx_eq(elements.periapsis(), expected_periapsis, 1e-6));

    // Verify orbital energy is negative (bound orbit)
    assert!(elements.specific_energy() < 0.0);
}

// ============================================================================
// Time System Tests
// ============================================================================

#[test]
fn test_time_conversions_chain() {
    // Start with a specific UTC time
    let utc = UTC::from_components(2024, 6, 21, 12, 0, 0.0);

    // Convert through all time systems
    let tai = TimeTransforms::utc_to_tai(&utc);
    let tt = TimeTransforms::tai_to_tt(&tai);
    let _gps = TimeTransforms::tai_to_gps(&tai);
    let jd = TimeTransforms::utc_to_jd(&utc);
    let gmst = TimeTransforms::utc_to_gmst(&utc);

    // TAI should be UTC + leap seconds (37 as of 2017)
    let leap_diff = tai.seconds() - utc.unix_seconds();
    assert!(leap_diff >= 37.0 && leap_diff < 40.0);

    // TT should be TAI + 32.184s
    let tt_tai_diff = tt.seconds() - tai.seconds();
    assert!(approx_eq(tt_tai_diff, 32.184, 1e-6));

    // Julian Date should be reasonable (around 2460000 for 2024)
    assert!(jd.value() > 2459000.0 && jd.value() < 2461000.0);

    // GMST should be in valid range
    assert!(gmst.to_radians() >= 0.0 && gmst.to_radians() < Constants::TWO_PI);
}

#[test]
fn test_julian_date_calculations() {
    // J2000.0 epoch: 2000-01-01 12:00:00 TT
    let j2000 = JulianDate::j2000();
    assert!(approx_eq(j2000.value(), 2451545.0, 1e-6));

    // Check Julian centuries calculation
    let t = j2000.julian_centuries_j2000();
    assert!(approx_eq(t, 0.0, 1e-10));

    // One century later
    let century_later = JulianDate::new(2451545.0 + 36525.0);
    let t_century = century_later.julian_centuries_j2000();
    assert!(approx_eq(t_century, 1.0, 1e-10));
}

// ============================================================================
// Celestial Body Tests
// ============================================================================

#[test]
fn test_sun_position() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    let sun_pos = Sun::eci_position(&epoch);
    let distance = sun_pos.magnitude();

    // Sun should be approximately 1 AU away
    let au = Constants::AU_M;
    assert!(
        distance > 0.98 * au && distance < 1.02 * au,
        "Sun distance {} not within expected range of 1 AU",
        distance
    );

    // Sun direction should be a unit vector
    let sun_dir = Sun::direction(&epoch);
    assert!(approx_eq(sun_dir.magnitude(), 1.0, 1e-10));
}

#[test]
fn test_moon_position() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    let moon_pos = Moon::eci_position(&epoch);
    let distance = moon_pos.magnitude();

    // Moon should be approximately 384,400 km away
    let expected_distance = 384_400_000.0; // meters
    assert!(
        distance > 0.9 * expected_distance && distance < 1.1 * expected_distance,
        "Moon distance {} not within expected range",
        distance
    );
}

#[test]
fn test_moon_phase() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    let phase = Moon::phase_angle(&epoch);

    // Phase angle should be between 0 and π
    assert!(
        phase >= 0.0 && phase <= std::f64::consts::PI,
        "Phase angle {} out of range",
        phase
    );
}

#[test]
fn test_earth_orientation() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    // Get precession angles
    let precession = Earth::precession_angles(&epoch);

    // Precession angles should be small (radians)
    assert!(precession.zeta.abs() < 0.1);
    assert!(precession.theta.abs() < 0.1);
    assert!(precession.z.abs() < 0.1);

    // Get nutation angles
    let nutation = Earth::nutation_angles(&epoch);

    // Mean obliquity should be around 23.4 degrees
    let obliquity_deg = nutation.m_eps * Constants::RAD_TO_DEG;
    assert!(
        obliquity_deg > 23.0 && obliquity_deg < 24.0,
        "Obliquity {} not in expected range",
        obliquity_deg
    );
}

// ============================================================================
// Complete Workflow Integration Tests
// ============================================================================

#[test]
fn test_satellite_tracking_workflow() {
    // This test simulates a complete satellite tracking workflow

    // 1. Parse TLE
    let line1 = "1 25544U 98067A   26014.62805721  .00006818  00000+0  13044-3 0  9991";
    let line2 = "2 25544  51.6333 339.6562 0007763  17.9854 342.1408 15.49289811547943";
    let tle = Tle::parse(line1, line2).expect("TLE parsing failed");

    // 2. Propagate to a specific time
    let observation_time = tle.epoch();
    let teme_state = tle
        .propagate(&observation_time)
        .expect("Propagation failed");

    // 3. Convert to ECI
    let eci_state = teme_state.to_eci();

    // 4. Define ground observer
    let observer = GeodeticState::new(40.0, -105.0, 1.6); // Denver, CO

    // 5. Compute observation angles
    let az_el = Observations::compute_az_el(&observer, &eci_state);
    let ra_dec = Observations::compute_ra_dec(&observer, &eci_state);

    // 6. Verify results are reasonable
    assert!(az_el.range.unwrap() > 0.0);
    assert!(ra_dec.range.unwrap() > 0.0);

    // 7. Convert ECI to Keplerian elements
    let keplerian = eci_state.to_keplerian();

    // 8. Verify Keplerian elements are reasonable for LEO
    assert!(keplerian.semi_major_axis > 6500.0 && keplerian.semi_major_axis < 7000.0);
    assert!(keplerian.eccentricity < 0.1);
}

#[test]
fn test_observation_geometry() {
    let epoch = Utc.with_ymd_and_hms(2024, 6, 15, 0, 0, 0).unwrap();

    // Multiple ground stations
    let stations = vec![
        ("Boulder", GeodeticState::new(40.0150, -105.2705, 1.655)),
        ("Tokyo", GeodeticState::new(35.6762, 139.6503, 0.04)),
        ("Sydney", GeodeticState::new(-33.8688, 151.2093, 0.058)),
    ];

    // High-altitude satellite (visible from multiple stations)
    let sat_pos = Vector3::new(30000.0, 10000.0, 5000.0);
    let sat_vel = Vector3::new(-0.5, 2.5, 0.1);
    let satellite = ECIState::new(epoch, sat_pos, sat_vel);

    for (name, station) in &stations {
        let az_el = Observations::compute_az_el(station, &satellite);
        let ra_dec = Observations::compute_ra_dec(station, &satellite);

        // All observations should have positive range
        assert!(az_el.range.unwrap() > 0.0, "{}: Invalid range", name);

        // RA should be in valid range [0, 2π)
        assert!(
            ra_dec.right_ascension >= 0.0 && ra_dec.right_ascension < Constants::TWO_PI,
            "{}: Invalid RA: {}",
            name,
            ra_dec.right_ascension
        );

        // Dec should be in valid range [-π/2, π/2]
        assert!(
            ra_dec.declination >= -std::f64::consts::FRAC_PI_2
                && ra_dec.declination <= std::f64::consts::FRAC_PI_2,
            "{}: Invalid Dec: {}",
            name,
            ra_dec.declination
        );
    }
}
