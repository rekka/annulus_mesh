//! Uniform point distribution in an annulus using Bridson’s algorithm for
//! Poisson-disc sampling
//! from [Visualizing Algorithms](http://bost.ocks.org/mike/algorithms/).
//!
//! Uses [qhull](http://qhull.org/) to generate a Delaunay triangulation.
extern crate docopt;
extern crate rand;
extern crate rustc_serialize;
extern crate gnuplot;
extern crate mersenne_twister;

mod annulus_distribution;

use std::f64::consts::PI;
#[allow(unused_imports)]
use gnuplot::{Figure, Caption, Color, Fix, AxesCommon, PlotOption, DashType, Coordinate, TextColor};
use rand::{Rng, SeedableRng};
use std::io::{Read, Write, Cursor, BufRead};
use std::process::{Command, Stdio};
use docopt::Docopt;
use mersenne_twister::MT19937_64;
use std::ops::Add;
use annulus_distribution::AnnulusDist;

const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static USAGE: &'static str = r"
Simple wave equation solver.
Usage: annulus_mesh [options]
       annulus_mesh (-h | --help | --version)

Options:
    -r INNER        Inner radius [default: 0.2]
    -d DIST         Point distance [default: 0.2]
    -h, --help      Print help
    --plot FILE     Output mesh plot
    --no-show       Do not plot anything
    --version       Print version
";

#[derive(RustcDecodable, Debug)]
struct Args {
    flag_d: f64,
    flag_r: f64,
    flag_no_show: bool,
    flag_plot: Option<String>,
}

#[derive(Copy, Clone, Debug)]
pub struct Point(f64, f64);

impl Add for Point {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Point(self.0 + rhs.0, self.1 + rhs.1)
    }
}

fn dist_sq(x: Point, y: Point) -> f64 {
    (x.0 - y.0).powi(2) + (x.1 - y.1).powi(2)
}

fn dist(x: Point, y: Point) -> f64 {
    ((x.0 - y.0).powi(2) + (x.1 - y.1).powi(2)).sqrt()
}

fn qhull_triangulation(ps: &[Point]) -> Vec<Vec<usize>> {
    let mut process = Command::new("qhull")
                          .stdin(Stdio::piped())
                          .stdout(Stdio::piped())
                          .arg("d")
                          .arg("i")
                          .spawn()
                          .unwrap_or_else(|e| panic!("couldn't spawn qhull: {:?}", e));

    if let Some(ref mut out) = process.stdin {
        writeln!(out, "2\n{}", ps.len()).unwrap();

        for &Point(x, y) in ps {
            writeln!(out, "{} {}", x, y).unwrap();
        }
    }

    // process.stdout.unwrap().read_to_string(&mut s);
    let output = process.wait_with_output().unwrap();

    let mut lines = Cursor::new(output.stdout).lines();

    let n = lines.next().unwrap().unwrap().parse::<usize>().unwrap();

    let mut res = vec![];

    for line in lines {
        let s = line.unwrap();
        let poly: Vec<usize> = s.split_whitespace().map(|e| e.parse::<usize>().unwrap()).collect();
        res.push(poly);
    }

    assert_eq!(res.len(), n);
    res
}

fn main() {
    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.version(Some(VERSION.to_string())).decode())
                         .unwrap_or_else(|e| e.exit());

    let seed: u64 = 0;
    let mut rng: MT19937_64 = SeedableRng::from_seed(seed);

    let mut ps = vec![];
    let mut active: Vec<Point> = vec![];

    let r1 = args.flag_r;
    let r2 = 1.;
    let n_gen = 30;

    let h = args.flag_d;
    let mean_dist_ratio = 0.7;

    let k1 = (mean_dist_ratio * 2. * PI * r1 / h).round() as i32;
    let k2 = (mean_dist_ratio * 2. * PI * r2 / h).round() as i32;

    // generate boundary points
    for &(k, r) in &[(k1, r1), (k2, r2)] {
        for i in 0..k {
            let (y, x) = (2. * PI * i as f64 / k as f64).sin_cos();
            ps.push(Point(r * x, r * y));
        }
    }

    // make boundary points active
    active.extend(&ps);

    let n_fixed = ps.len();

    let annulus_dist = AnnulusDist::new(h, 2.0 * h);

    while !active.is_empty() {
        // select a random active point
        let i = rng.gen_range(0, active.len());

        let c = active[i];

        if let Some(p) = annulus_dist.ind_iter(&mut rng)
                                     .map(|p| p + c)
                                     .filter(|&p| {
                                         let d = dist(p, Point(0., 0.));
                                         r1 <= d && d <= r2
                                     })
                                     .take(n_gen)
                                     .filter(|&p| ps.iter().all(|&q| dist(p, q) >= h))
                                     .next() {
            ps.push(p);
            active.push(p);
        } else {
            active.remove(i);
        }

    }

    let indices = qhull_triangulation(&ps);


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
    let n_laplace_step = 1;

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
            change += dist_sq(p, q);
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
                axes.lines(poly.iter().cycle().take(poly.len() + 1).map(|&i| ps[i as usize].0),
                           poly.iter().cycle().take(poly.len() + 1).map(|&i| ps[i as usize].1),
                           &[]);
            }
        }


        fg.show();
    }

}
