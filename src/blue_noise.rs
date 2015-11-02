//! Generator of a Poisson-disc distribution in a plane (blue noise) using
//! Bridson's algorithm for Poisson-disc sampling.
//!
//! Bridson, Robert. "Fast Poisson disk sampling in arbitrary dimensions." ACM
//! SIGGRAPH. Vol. 2007. 2007.
use super::Point;
use rand::Rng;
use annulus_distribution::AnnulusDist;

/// Generate blue noise in a rectange given by the two corners in `bounding_box`. The points have
/// minimal distance `r`. Seed points can be given in `ps`.
///
/// # Panics
///
/// Panics if `r <= 0`. Might panic if the distance of seed points in `ps` is less than `r`.
///
pub fn generate_blue_noise<R: Rng>(rng: &mut R,
                                   ps: &mut Vec<Point>,
                                   r: f64,
                                   bounding_box: (Point, Point)) {
    generate_blue_noise_cull(rng, ps, r, bounding_box, |_| true);
}

/// Generate blue noise in a rectange given by the two corners in `bounding_box`. The points have
/// minimal distance `r`. Seed points can be given in `ps`. Points are culled if the provided
/// function `cull` return `false`.
///
/// # Panics
///
/// Panics if `r <= 0`. Might panic if the distance of seed points in `ps` is less than `r`.
///
pub fn generate_blue_noise_cull<R: Rng, F: FnMut(Point) -> bool>(rng: &mut R,
                                                                 ps: &mut Vec<Point>,
                                                                 r: f64,
                                                                 bounding_box: (Point, Point),
                                                                 mut cull: F) {

    assert!(r > 0.);

    let (p0, p1) = bounding_box;

    // initialize neighbor grid: the diagonal of the square cell should be
    // slightly less than r / sqrt(2) to prevent the danger of storing two particles
    // in the same cell.
    let mut ng = NeighborGrid::new(p0, p1, 0.99 * r / 2f64.sqrt());

    // generate a random point if no seed point provided
    if ps.is_empty() {
        let p = Point(ng.xmin + rng.gen_range(0., ng.width),
                      ng.ymin + rng.gen_range(0., ng.height));
        ps.push(p);
    }

    // fill the neighbor grid by the points provided
    for (i, &p) in ps.iter().enumerate() {
        ng.push(p, i as i32);
    }

    // all points are active initially
    let mut active = ps.clone();

    let annulus_dist = AnnulusDist::new(r, 2.0 * r);
    let n_gen = 30;

    while !active.is_empty() {
        // select a random active point
        let i = rng.gen_range(0, active.len());
        let c = active[i];

        // generate a fixed number of points in an (r, 2r)-annulus around c, within the
        // bounding
        // box, cull them if necessary, and select the first one that is further than r
        // away from
        // all previous points (if such a point exists)
        if let Some(p) = annulus_dist.ind_iter(rng)
                                     .map(|p| p + c)
                                     .take(n_gen)
                                     .filter(|&p| ng.is_within(p) && cull(p))
                                     .filter(|&p| {
                                         ng.neighbor_iter(p).all(|k| p.dist(ps[k as usize]) >= r)
                                     })
                                     .next() {
            // store the generated neighbor
            ng.push(p, ps.len() as i32);
            ps.push(p);
            active.push(p);
        } else {
            // no admissible neighbor generated, make point inactive
            active.swap_remove(i);
        }

    }
}

/// Data structure for efficiently checking neighborhood particles.
struct NeighborGrid {
    nx: usize,
    ny: usize,
    h: f64,
    xmin: f64,
    ymin: f64,
    width: f64,
    height: f64,
    grid: Vec<i32>,
}

struct NeighborIterator<'a> {
    ng: &'a NeighborGrid,
    i: usize,
    j: usize,
    ci: usize,
    cj: usize,
}

impl NeighborGrid {
    fn new(p0: Point, p1: Point, h: f64) -> NeighborGrid {
        let xmin = p0.0.min(p1.0);
        let ymin = p0.1.min(p1.1);
        let width = (p0.0 - p1.0).abs();
        let height = (p0.1 - p1.1).abs();
        let nx = (width / h).ceil() as usize + 1;
        let ny = (height / h).ceil() as usize + 1;

        let grid = vec![-1; nx * ny];

        NeighborGrid {
            nx: nx,
            ny: ny,
            h: h,
            xmin: xmin,
            ymin: ymin,
            width: width,
            height: height,
            grid: grid,
        }
    }

    /// Check if the point is withing the bounding box.
    fn is_within(&self, p: Point) -> bool {
        let Point(x, y) = p;
        self.xmin <= x && x <= self.xmin + self.width && self.ymin <= y &&
        y <= self.ymin + self.height
    }

    fn get_coord(&self, p: Point) -> (usize, usize) {
        let Point(x, y) = p;
        assert!(self.is_within(p));
        let kx = ((x - self.xmin) / self.h).floor() as usize;
        let ky = ((y - self.ymin) / self.h).floor() as usize;

        return (kx, ky);
    }

    /// Store a point in the grid.
    ///
    /// # Panics
    ///
    /// Panics if the grid cell is already occupied.
    fn push(&mut self, p: Point, i: i32) {
        let (kx, ky) = self.get_coord(p);
        let k = kx + self.nx * ky;

        assert_eq!(self.grid[k], -1);
        self.grid[k] = i;
    }

    /// Returns an iterator over a grid neighborhood of radius 2 (possible particle neighbors) and
    /// returns the indices in non-empty cells.
    fn neighbor_iter(&self, p: Point) -> NeighborIterator {
        use std::cmp::max;
        let (i, j) = self.get_coord(p);
        NeighborIterator {
            ng: self,
            i: i,
            j: j,
            ci: max(2, i) - 2,
            cj: max(2, j) - 2,
        }
    }
}

/// Iterator over a grid neighborhood of radius 2 that returns indices in non-empty cells.
impl<'a> NeighborIterator<'a> {
    fn next_inter(&mut self) -> Option<i32> {
        use std::cmp::{min, max};
        let ui = min(self.i + 2, self.ng.nx - 1);
        let uj = min(self.j + 2, self.ng.ny - 1);
        if self.cj > uj {
            return None;
        }
        let r = Some(self.ng.grid[self.ci + self.cj * self.ng.nx]);

        self.ci += 1;
        if self.ci > ui {
            self.cj += 1;
            self.ci = max(self.i, 2) - 2;
        }

        r
    }
}

impl<'a> Iterator for NeighborIterator<'a> {
    type Item = i32;

    fn next(&mut self) -> Option<i32> {
        while let Some(k) = self.next_inter() {
            if k >= 0 {
                return Some(k);
            }
        }
        None
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use super::NeighborGrid;
    use Point;

    #[test]
    fn test_neighbor_iterator() {
        let mut ng = NeighborGrid::new(Point(0., 0.), Point(1., 1.), 0.25);

        ng.push(Point(0., 0.), 0);
        ng.push(Point(1., 1.), 1);
        ng.push(Point(0.9, 0.99), 2);

        assert_eq!(ng.neighbor_iter(Point(0.3, 0.4)).collect::<Vec<_>>(),
                   vec![0, 2]);
    }

    #[test]
    fn test_new() {
        NeighborGrid::new(Point(0., 0.), Point(0., 0.), 0.25);
    }

    #[test]
    fn test_generate_blue_noise() {
        let mut rng = ::rand::thread_rng();
        let mut ps = vec![];
        generate_blue_noise(&mut rng, &mut ps, 0.125, (Point(-1., -2.), Point(3., 4.)));

        assert!(!ps.is_empty());
    }
}
