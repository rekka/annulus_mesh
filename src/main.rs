//! Uniform point distribution in an annulus using Bridson’s algorithm for
//! Poisson-disc sampling.
//!
//! Bridson, Robert. "Fast Poisson disk sampling in arbitrary dimensions." ACM
//! SIGGRAPH. Vol. 2007. 2007.
//!
//! See also [Visualizing Algorithms](http://bost.ocks.org/mike/algorithms/).
extern crate docopt;
extern crate gnuplot;
extern crate pcg_rand;
extern crate rand;
#[macro_use]
extern crate serde_derive;
extern crate delaunator;

mod annulus_distribution;
pub mod blue_noise;

use docopt::Docopt;
#[allow(unused_imports)]
use gnuplot::{
    AxesCommon, Caption, Color, Coordinate, DashType, Figure, Fix, PlotOption, TextColor,
};
use pcg_rand::seeds::PcgSeeder;
use pcg_rand::Pcg32;
use rand::SeedableRng;
use std::f64::consts::PI;
use std::ops::{Add, Div};

const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static USAGE: &'static str = r"
Simple annulus triangulation generator.

Usage: annulus_mesh [options]
       annulus_mesh (-h | --help | --version)

Options:
    -r INNER        Inner radius [default: 0.2]
    -d DIST         Point distance [default: 0.2]
    -l NUM          Laplace smoothing steps [default: 1]
    -h, --help      Print help
    --plot FILE     Output mesh plot
    --no-show       Do not plot anything
    --version       Print version
";

#[derive(Deserialize, Debug)]
struct Args {
    flag_d: f64,
    flag_r: f64,
    flag_l: u32,
    flag_no_show: bool,
    flag_plot: Option<String>,
}

/// A point in 2D.
#[derive(Copy, Clone, Debug)]
pub struct Point(f64, f64);

impl Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Point(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Div<f64> for Point {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Point(self.0 / rhs, self.1 / rhs)
    }
}

impl Point {
    /// Distance between two points.
    fn dist(self, y: Point) -> f64 {
        self.dist_sq(y).sqrt()
    }

    /// The square of the distance between two points.
    fn dist_sq(self, y: Point) -> f64 {
        (self.0 - y.0).powi(2) + (self.1 - y.1).powi(2)
    }
}

/// Compute the Delaunay triangulation of a set of points.
///
/// # Panics
///
/// Panics if a triangulation does not exist.
fn triangulation(ps: &[Point]) -> Vec<Vec<usize>> {
    let ps: Vec<_> = ps.iter().map(|&Point(x, y)| delaunator::Point { x, y }).collect();

    let tri = delaunator::triangulate(&ps).expect("No triangulation exists.");

    tri.triangles.chunks(3).map(|c| c.to_owned()).collect::<Vec<_>>()
}

fn main() {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.version(Some(VERSION.to_string())).deserialize())
        .unwrap_or_else(|e| e.exit());

    let seed: u64 = 1;
    let mut rng = Pcg32::from_seed(PcgSeeder::seed(seed));

    let mut ps = vec![];

    let r1 = args.flag_r;
    let r2 = 1.;

    let h = args.flag_d;
    let mean_dist_ratio = 0.7;

    // generate boundary points
    for &r in &[r1, r2] {
        let k = (mean_dist_ratio * 2. * PI * r / h).round() as i32;
        for i in 0..k {
            let (y, x) = (2. * PI * i as f64 / k as f64).sin_cos();
            ps.push(Point(r * x, r * y));
        }
    }

    let n_fixed = ps.len();

    let mut ps = blue_noise::generate_blue_noise_cull(
        &mut rng,
        ps,
        h,
        (Point(-r2, -r2), Point(r2, r2)),
        |p| {
            let d = p.dist(Point(0., 0.));
            r1 <= d && d <= r2
        },
    ).unwrap();

    println!("Generated {} nodes", ps.len());

    let indices = triangulation(&ps);

    // Cull triangles outside of the annulus by testing if the barycenter lies outside.
    let indices = indices.into_iter().filter(|t| (t.iter().map(|&i| ps[i]).fold(Point(0., 0.), Add::add) / t.len() as f64).dist(Point(0., 0.)) >= r1).collect::<Vec<_>>();

    println!("Delaunay triangulation has {} triangles", indices.len());

    let mut neighbors: Vec<Vec<usize>> = ps.iter().map(|_| vec![]).collect();
    for poly in &indices {
        for (&i, &j) in poly.iter().zip(poly.iter().cycle().skip(1)) {
            neighbors[i].push(j);
            neighbors[j].push(i);
        }
    }

    for list in &mut neighbors {
        list.sort();
        list.dedup();
    }

    // Laplace smoothing
    let n_laplace_step = args.flag_l;

    let laplace_eps = 1.0e-3;

    for _ in 0..n_laplace_step {
        let mut change = 0.;
        // Gauss-Seidel
        for i in n_fixed..ps.len() {
            let mut x = 0.;
            let mut y = 0.;
            for &j in &neighbors[i] {
                x += ps[j].0;
                y += ps[j].1;
            }
            let n = neighbors[i].len() as f64;
            let p = ps[i];
            let q = Point(x / n, y / n);
            ps[i] = q;
            change += p.dist_sq(q);
        }

        if change.sqrt() < laplace_eps {
            break;
        }
    }

    if !args.flag_no_show {
        let mut fg = Figure::new();
        if let Some(file) = args.flag_plot {
            fg.set_terminal("pngcairo size 720, 720", &file);
        }

        fg.clear_axes();
        {
            let axes = fg.axes2d();
            axes.set_aspect_ratio(Fix(1.));
            for poly in indices {
                axes.lines(
                    poly.iter()
                        .cycle()
                        .take(poly.len() + 1)
                        .map(|&i| ps[i as usize].0),
                    poly.iter()
                        .cycle()
                        .take(poly.len() + 1)
                        .map(|&i| ps[i as usize].1),
                    &[],
                );
            }
        }

        fg.show();
    }
}
