//! Uniform point distribution in an annulus using Bridson’s algorithm for Poisson-disc sampling
//! from [Visualizing Algorithms](http://bost.ocks.org/mike/algorithms/).
//!
//! Uses [qhull](http://qhull.org/) to generate a Delaunay triangulation.
extern crate docopt;
extern crate rand;
extern crate rustc_serialize;
extern crate gnuplot;
extern crate mersenne_twister;

use std::f64::consts::PI;
#[allow(unused_imports)]
use gnuplot::{Figure, Caption, Color, Fix, AxesCommon, PlotOption, DashType, Coordinate, TextColor};
use rand::{Rng, SeedableRng};
use rand::distributions::{IndependentSample, Range};
use std::io::{Read, Write, Cursor, BufRead};
use std::process::{Command, Stdio};
use docopt::Docopt;
use mersenne_twister::MT19937_64;

const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static USAGE: &'static str = "
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

fn dist(x: (f64, f64), y: (f64, f64)) -> f64 {
    ((x.0 - y.0).powi(2) + (x.1 - y.1).powi(2)).sqrt()
}

fn qhull_triangulation(ps: &[(f64, f64)]) -> Vec<Vec<usize>> {
    let mut process = Command::new("qhull")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .arg("d")
        .arg("i")
        .spawn()
        .unwrap_or_else(|e|
            panic!("couldn't spawn qhull: {:?}", e)
        );

    if let Some(ref mut out) = process.stdin {
        writeln!(out, "2\n{}", ps.len()).unwrap();

        for &(x,y) in ps {
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
    let mut active: Vec<(f64, f64)> = vec![];

    let r1 = args.flag_r;
    let r2 = 1.;
    let n_gen = 50;

    let h = args.flag_d;
    let mean_dist_ratio = 0.7;

    let k1 = (mean_dist_ratio * 2. * PI * r1 / h).round() as i32;
    let k2 = (mean_dist_ratio * 2. * PI * r2 / h).round() as i32;

    // generate boundary points
    for &(k, r) in &[(k1, r1), (k2, r2)] {
        for i in 0..k {
            ps.push((r * (2. * PI * i as f64 / k as f64).cos(),
                    r * (2. * PI * i as f64 / k as f64).sin()));
        }
    }

    // make boundary points active
    active.extend(&ps);

    let n_fixed = ps.len();

    let x_range = Range::new(h * h, 4. * h * h);
    let theta_range = Range::new(0., 2. * PI);

    while !active.is_empty() {
        // select a random active point
        let i = rng.gen_range(0, active.len());

        let c = active[i];
        let mut deactivate = true;

        for _ in 0..n_gen {
            let mut p;

            // generate a random point in an annulus [h, 2h] of the active point
            loop {
                let r = x_range.ind_sample(&mut rng).sqrt();
                let theta = theta_range.ind_sample(&mut rng);

                p = (r * theta.cos() + c.0, r * theta.sin() + c.1);
                let d = dist(p, (0., 0.));

                // make sure it lies in the domain
                if r1 <= d && d <= r2 {
                    break;
                }
            }

            if ps.iter().all(|&q| dist(p, q) >= h) {
                ps.push(p);
                active.push(p);
                deactivate = false;
                break;
            }
        }

        if deactivate {
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
            let (x0, y0) = ps[i];
            let (x1, y1) = (x / n, y / n);
            ps[i] = (x1, y1);
            change += (x0 - x1).powi(2) + (y0 - y1).powi(2);
        }

        if change.sqrt() < laplace_eps { break; }
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
                        poly.iter().cycle().take(poly.len() + 1).map(|&i| ps[i as usize].1), &[]);
            }
        }


        fg.show();
    }

}
