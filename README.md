# SpaceDust

[![Crates.io](https://img.shields.io/crates/v/space-dust.svg)](https://crates.io/crates/space-dust)
[![Documentation](https://docs.rs/space-dust/badge.svg)](https://docs.rs/space-dust)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive astrodynamics library for Rust, providing tools for satellite tracking, orbital mechanics, coordinate transformations, and ground-based observation calculations.

## Features

- **Time Systems**: Conversions between UTC, TAI, TT, Julian Date, GPS, and GMST
- **Coordinate Frames**: Transformations between ECI J2000, TEME, ECEF, and Geodetic coordinates
- **Orbital Elements**: Conversions between Cartesian state vectors and Keplerian elements
- **TLE Parsing**: Parse and propagate Two-Line Element Sets using SGP4
- **Celestial Bodies**: Sun and Moon position calculations
- **Ground Observations**: Azimuth/Elevation and Right Ascension/Declination calculations
- **High Performance**: Efficient numerical operations suitable for real-time applications

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
space-dust = "0.1"
```

### Optional Features

- `serde` - Enable serialization/deserialization support
- `network` - Enable network features for fetching TLE data
- `full` - Enable all optional features

```toml
[dependencies]
space-dust = { version = "0.1", features = ["full"] }
```

## Quick Start

### Satellite Tracking

```rust
use space_dust::tle::Tle;
use space_dust::state::{TEMEState, GeodeticState, StateTransforms};
use space_dust::observations::Observations;
use chrono::Utc;

// Parse a TLE (ISS example)
let line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9002";
let line2 = "2 25544  51.6400 208.1200 0001234  85.0000 275.0000 15.48919100123456";
let tle = Tle::parse(line1, line2).unwrap();

// Propagate to current time
let epoch = Utc::now();
let teme_state = tle.propagate(&epoch).unwrap();

// Convert to ECI J2000
let eci_state = teme_state.to_eci();

// Define a ground observer (Denver, CO)
let observer = GeodeticState::new(39.7392, -104.9903, 1.6);

// Compute observation angles
let az_el = Observations::compute_az_el(&observer, &eci_state);
println!("Azimuth: {:.2}°, Elevation: {:.2}°", az_el.azimuth_deg(), az_el.elevation_deg());
```

### Time System Conversions

```rust
use space_dust::time::{TimeTransforms, JulianDate, UTC, TAI, TT};
use chrono::Utc;

let now = Utc::now();

// Convert to Julian Date
let jd = JulianDate::from_utc(&now);
println!("Julian Date: {}", jd.value());

// Convert between time systems
let tai = TAI::from_utc(&now);
let tt = TT::from_tai(&tai);
```

### Coordinate Transformations

```rust
use space_dust::state::{ECIState, ECEFState, GeodeticState, StateTransforms};
use chrono::Utc;

// Create an ECI state (position in km, velocity in km/s)
let eci = ECIState::new(
    -6045.0, -3490.0, 2500.0,  // position
    -3.457, 6.618, 2.533,      // velocity
);

// Convert to ECEF
let epoch = Utc::now();
let ecef = eci.to_ecef(&epoch);

// Convert to Geodetic (lat/lon/alt)
let geodetic = ecef.to_geodetic();
println!("Lat: {:.4}°, Lon: {:.4}°, Alt: {:.2} km", 
    geodetic.latitude_deg(), 
    geodetic.longitude_deg(), 
    geodetic.altitude_km()
);
```

### Keplerian Elements

```rust
use space_dust::state::{ECIState, KeplerianElements};

// Create from Cartesian state
let eci = ECIState::new(
    -6045.0, -3490.0, 2500.0,
    -3.457, 6.618, 2.533,
);

let elements = KeplerianElements::from_cartesian(&eci);
println!("Semi-major axis: {:.2} km", elements.semi_major_axis());
println!("Eccentricity: {:.6}", elements.eccentricity());
println!("Inclination: {:.2}°", elements.inclination_deg());
```

### Celestial Body Positions

```rust
use space_dust::bodies::{Sun, Moon};
use chrono::Utc;

let epoch = Utc::now();

// Get Sun position in ECI coordinates
let sun_pos = Sun::position_eci(&epoch);
println!("Sun position: {:?} km", sun_pos);

// Get Moon position in ECI coordinates
let moon_pos = Moon::position_eci(&epoch);
println!("Moon position: {:?} km", moon_pos);
```

## Modules

| Module | Description |
|--------|-------------|
| `time` | Time system conversions (UTC, TAI, TT, JD, GPS, GMST) |
| `state` | Coordinate frames and state vectors (ECI, TEME, ECEF, Geodetic, Keplerian) |
| `tle` | Two-Line Element parsing and SGP4 propagation |
| `observations` | Ground-based observation calculations (Az/El, RA/Dec) |
| `bodies` | Celestial body positions (Sun, Moon, Earth) |
| `constants` | Physical and astronomical constants |
| `math` | Vector and matrix operations |
| `data` | Data handling utilities |

## Accuracy

- Time conversions: Sub-microsecond accuracy
- Coordinate transformations: Suitable for operational satellite tracking
- SGP4 propagation: Standard SGP4/SDP4 accuracy (typically < 1 km for well-maintained TLEs)
- Celestial positions: Low-precision algorithms suitable for most applications

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Based on algorithms from Vallado's "Fundamentals of Astrodynamics and Applications"
- SGP4 implementation provided by the [sgp4](https://crates.io/crates/sgp4) crate
- Translated to Rust from the original [SpaceDust](https://github.com/Stratogen-Applied-Research/space_dust) Elixir library
