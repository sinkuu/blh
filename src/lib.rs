pub mod ellipsoid;
use ellipsoid::Ellipsoid;

use std::f64;
use std::marker::PhantomData;

/// A position on the earth represented by latitude [rad], longitude [rad], and geoid height [m].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodeticCoord {
    pub lat: f64,
    pub lon: f64,
    pub hgt: f64,
}

impl GeodeticCoord {
    #[inline]
    pub fn new(lat: f64, lon: f64, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat,
            lon: lon,
            hgt: hgt,
        }
    }

    #[inline]
    pub fn from_deg(lat: f64, lon: f64, hgt: f64) -> Self {
        GeodeticCoord {
            lat: lat.to_radians(),
            lon: lon.to_radians(),
            hgt: hgt,
        }
    }
}

/// A position on the earth represented by x [m], y[m], and z [m] in ECEF coordinate.
#[derive(Debug, Clone, Copy)]
pub struct GeocentricCoord<E> {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    ellipsoid: PhantomData<E>,
}

// This can't be derived because `E` can be an non-`PartialEq` type.
impl<E> PartialEq for GeocentricCoord<E> {
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
            GeodeticCoord {
                lat: 0.0,
                lon: 0.0,
                hgt: 0.0,
            }
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

            GeodeticCoord::from_rad(lat, lon, hgt)
        }
    }
}

impl<E: Ellipsoid> From<GeodeticCoord> for GeocentricCoord<E> {
    fn from(blh: GeodeticCoord) -> GeocentricCoord<E> {
        // Prime vertical's radius of carvative
        let pv_rc = E::radius() / f64::sqrt(1.0 - E::ecc().powi(2) * f64::sin(blh.lat).powi(2));

        GeocentricCoord {
            x: (pv_rc + blh.hgt) * f64::cos(blh.lat) * f64::cos(blh.lon),
            y: (pv_rc + blh.hgt) * f64::cos(blh.lat) * f64::sin(blh.lon),
            z: (pv_rc * (1.0 - E::ecc().powi(2)) + blh.hgt) * f64::sin(blh.lat),
            ellipsoid: PhantomData,
        }
    }
}

#[test]
fn test_pos_conversion() {
    use ellipsoid::{WGS84, GRS80};

    let blh = GeodeticCoord::from_deg(31.0, 131.0, 100.0);
    assert_eq!(GeodeticCoord::from(GeocentricCoord::<WGS84>::from(blh)),
               blh);

    let xyz = GeodeticCoord::from_deg(33.43, 135.7, 100.0);
    assert_eq!(GeodeticCoord::from(GeocentricCoord::<GRS80>::from(xyz)),
               xyz);
}
