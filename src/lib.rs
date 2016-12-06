extern crate euclid;

use euclid::{Degrees, Radians};

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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodeticCoord {
    pub lat: Radians<f64>,
    pub lon: Radians<f64>,
    pub hgt: f64,
}

impl GeodeticCoord {
    /// Construct a geodetic coord with latitude/longitude in radians.
    #[inline]
    pub fn new(lat: Radians<f64>, lon: Radians<f64>, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat,
            lon: lon,
            hgt: hgt,
        }
    }

    /// Construct a geodetic coord with latitude/longitude in degrees.
    #[inline]
    pub fn from_degrees(lat: Degrees<f64>, lon: Degrees<f64>, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat.to_radians(),
            lon: lon.to_radians(),
            hgt: hgt,
        }
    }
}

/// A position on the earth represented by x [m], y [m], and z [m] in geocentric (ECEF) coordinate.
#[derive(Clone, Copy)]
pub struct GeocentricCoord<E> {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    ellipsoid: PhantomData<E>,
}

impl<E: Ellipsoid> Debug for GeocentricCoord<E> {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        f.debug_struct(&format!("GeocentricCoord<{}>", E::name()))
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}

// This can't be derived because `E` can be an non-`PartialEq` type.
impl<E> PartialEq for GeocentricCoord<E> {
    #[inline]
    fn eq(&self, other: &GeocentricCoord<E>) -> bool {
        let &GeocentricCoord { x: sx, y: sy, z: sz, .. } = self;
        let &GeocentricCoord { x: ox, y: oy, z: oz, .. } = other;
        (sx, sy, sz) == (ox, oy, oz)
    }
}

impl<E> GeocentricCoord<E> {
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        GeocentricCoord {
            x: x,
            y: y,
            z: z,
            ellipsoid: PhantomData,
        }
    }
}

impl<E: Ellipsoid> From<GeocentricCoord<E>> for GeodeticCoord {
    fn from(xyz: GeocentricCoord<E>) -> GeodeticCoord {
        if xyz ==
            (GeocentricCoord {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                ellipsoid: PhantomData,
            })
        {
            GeodeticCoord::new(Radians::new(0.0), Radians::new(0.0), 0.0)
        } else {
            let GeocentricCoord { x, y, z, .. } = xyz;

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

impl<E: Ellipsoid> From<GeodeticCoord> for GeocentricCoord<E> {
    fn from(blh: GeodeticCoord) -> GeocentricCoord<E> {
        let (lat, lon) = (blh.lat.get(), blh.lon.get());

        // Prime vertical's radius of carvative
        let pv_rc = E::radius() / f64::sqrt(1.0 - E::ecc().powi(2) * f64::sin(lat).powi(2));

        GeocentricCoord {
            x: (pv_rc + blh.hgt) * f64::cos(lat) * f64::cos(lon),
            y: (pv_rc + blh.hgt) * f64::cos(lat) * f64::sin(lon),
            z: (pv_rc * (1.0 - E::ecc().powi(2)) + blh.hgt) * f64::sin(lat),
            ellipsoid: PhantomData,
        }
    }
}

#[test]
fn test_pos_conversion() {
    use ellipsoid::{WGS84, GRS80};

    let blh = GeodeticCoord::from_degrees(Degrees::new(31.0), Degrees::new(131.0), 100.0);
    assert_eq!(GeodeticCoord::from(GeocentricCoord::<WGS84>::from(blh)), blh);

    let xyz = GeodeticCoord::from_degrees(Degrees::new(33.43), Degrees::new(135.7), 100.0);
    assert_eq!(GeodeticCoord::from(GeocentricCoord::<GRS80>::from(xyz)), xyz);
}
