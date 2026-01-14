# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2024-01-15

### Added

- Initial release of SpaceDust for Rust
- **Time Systems Module** (`time`)
  - UTC, TAI, TT time system representations
  - Julian Date conversions
  - GPS time support
  - Greenwich Mean Sidereal Time (GMST) calculations
  - `TimeTransforms` trait for unified conversions
- **Coordinate Frames Module** (`state`)
  - ECI J2000 inertial frame
  - TEME (True Equator Mean Equinox) frame
  - ECEF (Earth-Centered Earth-Fixed) frame
  - Geodetic coordinates (latitude, longitude, altitude)
  - Keplerian orbital elements
  - `StateTransforms` trait for frame conversions
- **TLE Module** (`tle`)
  - Two-Line Element Set parsing
  - SGP4/SDP4 orbit propagation via `sgp4` crate
  - TLE validation and error handling
- **Observations Module** (`observations`)
  - Azimuth/Elevation calculations for ground observers
  - Right Ascension/Declination calculations
  - Topocentric coordinate transformations
- **Celestial Bodies Module** (`bodies`)
  - Sun position calculations (low-precision algorithm)
  - Moon position calculations (low-precision algorithm)
  - Earth model with WGS84 parameters
- **Constants Module** (`constants`)
  - Physical constants (speed of light, gravitational constant)
  - Earth parameters (radius, flattening, J2)
  - Astronomical constants (AU, solar parameters)
- **Math Module** (`math`)
  - 3D vector operations
  - 3x3 rotation matrices
  - Common trigonometric utilities
- **Data Module** (`data`)
  - Data handling utilities
- Optional `serde` feature for serialization support
- Optional `network` feature for TLE fetching
- Comprehensive documentation with examples
- Unit tests for all modules
- Integration tests for end-to-end workflows

### Notes

- Minimum supported Rust version: 1.70
- No `unsafe` code used
- Full documentation available on docs.rs

[Unreleased]: https://github.com/Stratogen-Applied-Research/space_dust/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Stratogen-Applied-Research/space_dust/releases/tag/v0.1.0