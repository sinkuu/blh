use std::f64;

pub trait Ellipsoid {
    /// Equatorial radius [m].
    fn radius() -> f64;

    /// Inverse flattening.
    fn flattening_inv() -> f64;

    /// Polar radius [m].
    #[inline]
    fn polar_radius() -> f64 {
        Self::radius() * (1.0 - Self::flattening())
    }

    /// Flattening.
    #[inline]
    fn flattening() -> f64 {
        1.0 / Self::flattening_inv()
    }

    /// Eccentricity.
    #[inline]
    fn ecc() -> f64 {
        f64::sqrt(Self::flattening() * (2.0 - Self::flattening()))
    }

    /// Second eccentricity.
    #[inline]
    fn second_ecc() -> f64 {
        f64::sqrt((Self::radius().powi(2) - Self::polar_radius().powi(2)) /
                  Self::polar_radius().powi(2))
    }

    /// Third eccentricity.
    #[inline]
    fn third_ecc() -> f64 {
        f64::sqrt((Self::radius().powi(2) - Self::polar_radius().powi(2)) /
                  (Self::radius().powi(2) + Self::radius().powi(2)))
    }

    /// Angular eccentricity.
    #[inline]
    fn angular_ecc() -> f64 {
        f64::acos(Self::polar_radius() / Self::radius())
    }
}

// Earth's ellipsoid constants based on WGS-84.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WGS84 {}

impl Ellipsoid for WGS84 {
    #[inline]
    fn radius() -> f64 {
        6_378_137.0
    }

    #[inline]
    fn flattening_inv() -> f64 {
        298.257_223_563
    }
}

// Earth's ellipsoid constants based on GRS80.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GRS80 {}

impl Ellipsoid for GRS80 {
    #[inline]
    fn radius() -> f64 {
        6_378_137.0
    }

    #[inline]
    fn flattening_inv() -> f64 {
        298.257_222_101
    }
}
