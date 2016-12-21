extern crate euclid;
#[macro_use] extern crate approx;

use euclid::{Degrees, Radians};
use approx::ApproxEq;

pub mod ellipsoid;
use ellipsoid::*;

use std::f64;
use std::marker::PhantomData;
use std::fmt::{self, Debug, Formatter};

trait DegreesExt<T> {
    fn to_radians(&self) -> Radians<T>;
}

impl DegreesExt<f64> for Degrees<f64> {
    #[inline]
    fn to_radians(&self) -> Radians<f64> {
        Radians::new(self.get() * f64::consts::PI / 180.0)
    }
}

trait RadiansExt<T> {
    fn to_degrees(&self) -> Degrees<T>;
}

impl RadiansExt<f64> for Radians<f64> {
    #[inline]
    fn to_degrees(&self) -> Degrees<f64> {
        Degrees::new(self.get() / f64::consts::PI * 180.0)
    }
}

/// A position on the earth represented by latitude [rad], longitude [rad], and geoid height [m].
pub struct GeodeticCoord<E> {
    pub lat: Radians<f64>,
    pub lon: Radians<f64>,
    pub hgt: f64,
    _ellipsoid: PhantomData<E>,
}

impl<E> GeodeticCoord<E> {
    /// Construct a geodetic coord with latitude/longitude in radians.
    #[inline]
    pub fn new(lat: Radians<f64>, lon: Radians<f64>, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat,
            lon: lon,
            hgt: hgt,
            _ellipsoid: PhantomData,
        }
    }

    /// Construct a geodetic coord with latitude/longitude in degrees.
    #[inline]
    pub fn from_degrees(lat: Degrees<f64>, lon: Degrees<f64>, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat.to_radians(),
            lon: lon.to_radians(),
            hgt: hgt,
            _ellipsoid: PhantomData,
        }
    }
}

impl<E> PartialEq for GeodeticCoord<E> {
    #[inline]
    fn eq(&self, other: &GeodeticCoord<E>) -> bool {
        let &GeodeticCoord { lat: slat, lon: slon, hgt: shgt, .. } = self;
        let &GeodeticCoord { lat: olat, lon: olon, hgt: ohgt, .. } = other;
        (slat, slon, shgt) == (olat, olon, ohgt)
    }
}

impl<E> ApproxEq for GeodeticCoord<E> {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn default_max_relative() -> f64 {
        f64::default_epsilon()
    }

    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn relative_eq(&self, other: &GeodeticCoord<E>, epsilon: f64, max_relative: f64) -> bool {
        self.lat.get().relative_eq(&other.lat.get(), epsilon, max_relative) &&
            self.lon.get().relative_eq(&other.lon.get(), epsilon, max_relative) &&
            self.hgt.relative_eq(&other.hgt, epsilon, max_relative)
    }

    fn ulps_eq(&self, other: &GeodeticCoord<E>, epsilon: f64, max_ulps: u32) -> bool {
        self.lat.get().ulps_eq(&other.lat.get(), epsilon, max_ulps) &&
            self.lon.get().ulps_eq(&other.lon.get(), epsilon, max_ulps) &&
            self.hgt.ulps_eq(&other.hgt, epsilon, max_ulps)
    }
}

impl<E: Ellipsoid> Debug for GeodeticCoord<E> {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        f.debug_struct(&format!("GeodeticCoord<{}>", E::name()))
            .field("lat", &self.lat)
            .field("lon", &self.lon)
            .field("hgt", &self.hgt)
            .finish()
    }
}

impl<E> Clone for GeodeticCoord<E> {
    #[inline]
    fn clone(&self) -> Self {
        GeodeticCoord {
            lat: self.lat,
            lon: self.lon,
            hgt: self.hgt,
            _ellipsoid: PhantomData,
        }
    }
}

impl<E> Copy for GeodeticCoord<E> {
}

/// A position on the earth represented by x [m], y [m], and z [m] in geocentric (ECEF) coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeocentricCoord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl GeocentricCoord {
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        GeocentricCoord {
            x: x,
            y: y,
            z: z,
        }
    }
}

impl ApproxEq for GeocentricCoord {
    type Epsilon = f64;

    fn default_epsilon() -> f64 {
        f64::default_epsilon()
    }

    fn default_max_relative() -> f64 {
        f64::default_epsilon()
    }

    fn default_max_ulps() -> u32 {
        f64::default_max_ulps()
    }

    fn relative_eq(&self, other: &GeocentricCoord, epsilon: f64, max_relative: f64) -> bool {
        self.x.relative_eq(&other.x, epsilon, max_relative) &&
            self.y.relative_eq(&other.y, epsilon, max_relative) &&
            self.z.relative_eq(&other.z, epsilon, max_relative)
    }

    fn ulps_eq(&self, other: &GeocentricCoord, epsilon: f64, max_ulps: u32) -> bool {
        self.x.ulps_eq(&other.x, epsilon, max_ulps) &&
            self.y.ulps_eq(&other.y, epsilon, max_ulps) &&
            self.z.ulps_eq(&other.z, epsilon, max_ulps)
    }
}

impl<E: Ellipsoid> From<GeocentricCoord> for GeodeticCoord<E> {
    fn from(xyz: GeocentricCoord) -> GeodeticCoord<E> {
        if xyz ==
            (GeocentricCoord {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            })
        {
            GeodeticCoord::new(Radians::new(0.0), Radians::new(0.0), 0.0)
        } else {
            let GeocentricCoord { x, y, z } = xyz;

            let radius = E::radius();
            let polar_radius = E::polar_radius();

            // eccentricity
            let e = E::ecc();
            // second eccentricity ^2
            let e2 = E::second_ecc().powi(2);

            let t = f64::atan2(radius * z, polar_radius * f64::sqrt(x.powi(2) + y.powi(2)));

            let lat = f64::atan2(z + e2 * polar_radius * f64::sin(t).powi(3),
                                 f64::sqrt(x.powi(2) + y.powi(2)) -
                                 e.powi(2) * radius * f64::cos(t).powi(3));

            // Prime vertical's radius of carvative
            let pv_rc = radius / f64::sqrt(1.0 - e.powi(2) * f64::sin(lat).powi(2));

            let lon = f64::atan2(y, x);
            let hgt = f64::sqrt((x.powi(2) + y.powi(2))) / f64::cos(lat) - pv_rc;

            GeodeticCoord::new(Radians::new(lat), Radians::new(lon), hgt)
        }
    }
}

impl<E: Ellipsoid> From<GeodeticCoord<E>> for GeocentricCoord {
    fn from(blh: GeodeticCoord<E>) -> GeocentricCoord {
        let (lat, lon) = (blh.lat.get(), blh.lon.get());

        // Prime vertical's radius of carvative
        let pv_rc = E::radius() / f64::sqrt(1.0 - E::ecc().powi(2) * f64::sin(lat).powi(2));

        GeocentricCoord {
            x: (pv_rc + blh.hgt) * f64::cos(lat) * f64::cos(lon),
            y: (pv_rc + blh.hgt) * f64::cos(lat) * f64::sin(lon),
            z: (pv_rc * (1.0 - E::ecc().powi(2)) + blh.hgt) * f64::sin(lat),
        }
    }
}

#[test]
fn test_pos_conversion() {
    use ellipsoid::WGS84;

    let blh = GeodeticCoord::<WGS84>::from_degrees(Degrees::new(31.0), Degrees::new(131.0), 100.0);
    assert_relative_eq!(GeodeticCoord::from(GeocentricCoord::from(blh)), blh);

    let xyz = GeocentricCoord::new(-3957314.620, 3310254.137, 3737540.043);
    assert_relative_eq!(GeocentricCoord::from(GeodeticCoord::<WGS84>::from(xyz)), xyz,
        epsilon = 0.1);
}
