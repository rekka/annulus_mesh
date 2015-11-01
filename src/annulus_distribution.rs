//! Implementation of a uniform distribuition of points on a two-dimensional
//! annulus.
use rand::{Closed01, Rng};
use rand::distributions::{Sample, IndependentSample};
pub use Point;
use std::f64::consts::PI;

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

    /// Iterator over independent samples of this distribution.
    pub fn ind_iter<'a, R: Rng>(&'a self, rng: &'a mut R) -> AnnulusDistIterator<'a, R> {
        AnnulusDistIterator {
            dist: self,
            rng: rng,
        }
    }
}

impl Sample<Point> for AnnulusDist {
    fn sample<R: Rng>(&mut self, rng: &mut R) -> Point {
        self.ind_sample(rng)
    }
}

impl IndependentSample<Point> for AnnulusDist {
    fn ind_sample<R: Rng>(&self, rng: &mut R) -> Point {
        let Closed01(t) = rng.gen::<Closed01<f64>>();
        let r = (self.r1_sq + t * (self.r2_sq - self.r1_sq)).sqrt();
        let (y, x) = (2. * PI * rng.gen::<f64>()).sin_cos();
        Point(r * x, r * y)
    }
}

struct AnnulusDistIterator<'a, R: 'a> {
    dist: &'a AnnulusDist,
    rng: &'a mut R,
}

impl<'a, R: 'a + Rng> Iterator for AnnulusDistIterator<'a, R> {
    type Item = Point;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.dist.ind_sample(self.rng))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_in_annulus() {
        let r1 = 12.;
        let r2 = 58.;

        assert!(AnnulusDist::new(r1, r2).ind_iter(&mut ::rand::thread_rng()).take(1000).all(|p| {
            let d = p.dist(Point(0., 0.));
            r1 <= d && d <= r2
        }));
    }
}
