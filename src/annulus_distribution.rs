//! Implementation of a uniform distribuition of points on a two-dimensional
//! annulus.
use rand::distributions::Distribution;
use rand::Rng;
use std::f64::consts::PI;
pub use Point;

/// The uniform distribution of 2D points on an annulus `{x: r_1 <= |x| <= r_2}`.
pub struct AnnulusDist {
    r1_sq: f64,
    r2_sq: f64,
}

impl AnnulusDist {
    /// Construct a new `AnnulusDist` with the given inner and outer radius
    /// `r1`, `r2`. Panics if not `0 < r1 < r2`.
    pub fn new(r1: f64, r2: f64) -> AnnulusDist {
        assert!(0. < r1, "AnnulusDist::new called with `r1 <= 0`");
        assert!(r1 < r2, "AnnulusDist::new called with `r2 <= r1`");
        AnnulusDist {
            r1_sq: r1 * r1,
            r2_sq: r2 * r2,
        }
    }
}

impl Distribution<Point> for AnnulusDist {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Point {
        // For points to be uniformly distributed in the annulus, the area of the disk with radius
        // equal to the distance of the point from the origin is distributed uniformly between r₁²
        // and r₂².
        let r = (self.r1_sq + rng.gen::<f64>() * (self.r2_sq - self.r1_sq)).sqrt();
        // The angle is uniform between 0 and 2π.
        let (y, x) = (2. * PI * rng.gen::<f64>()).sin_cos();
        Point(r * x, r * y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_in_annulus() {
        let r1 = 12.;
        let r2 = 58.;

        assert!(
            AnnulusDist::new(r1, r2)
                .sample_iter(&mut ::rand::thread_rng())
                .take(1000)
                .all(|p| {
                    let d = p.dist(Point(0., 0.));
                    r1 <= d && d <= r2
                })
        );
    }
}
