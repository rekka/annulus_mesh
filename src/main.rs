//! Uniform point distribution in an annulus using Bridson’s algorithm for Poisson-disc sampling
//! from [Visualizing Algorithms](http://bost.ocks.org/mike/algorithms/).
extern crate rand;
extern crate gnuplot;

use std::f64::consts::PI;
#[allow(unused_imports)]
use gnuplot::{Figure, Caption, Color, Fix, AxesCommon, PlotOption, DashType, Coordinate, TextColor};
use rand::Rng;
use rand::distributions::{IndependentSample, Range};

fn dist(x: (f64, f64), y: (f64, f64)) -> f64 {
    ((x.0 - y.0).powi(2) + (x.1 - y.1).powi(2)).sqrt()
}

fn main() {
    let mut rng = rand::thread_rng();

    let mut ps = vec![];
    let mut active: Vec<(f64, f64)> = vec![];

    let r1 = 0.1;
    let r2 = 1.;
    let n_gen = 50;

    let h = 0.025;
    let mean_dist_ratio = 0.7;

    let k1 = (mean_dist_ratio * 2. * PI * r1 / h).round() as i32;
    let k2 = (mean_dist_ratio * 2. * PI * r2 / h).round() as i32;

    for &(k, r) in &[(k1, r1), (k2, r2)] {
        for i in 0..k {
            ps.push((r * (2. * PI * i as f64 / k as f64).cos(),
                    r * (2. * PI * i as f64 / k as f64).sin()));
        }
    }

    active.extend(&ps);

    let x_range = Range::new(h * h, 4. * h * h);
    let theta_range = Range::new(0., 2. * PI);

    while !active.is_empty() {
        let i = rng.gen_range(0, active.len());

        let c = active[i];
        let mut deactivate = true;

        for _ in 0..n_gen {
            let mut p;

            loop {
                let r = x_range.ind_sample(&mut rng).sqrt();
                let theta = theta_range.ind_sample(&mut rng);

                p = (r * theta.cos() + c.0, r * theta.sin() + c.1);
                let d = dist(p, (0., 0.));

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


    let mut fg = Figure::new();
    fg.clear_axes();
    fg.axes2d()
        .set_aspect_ratio(Fix(1.))
        .points(ps.iter().map(|&x| x.0), ps.iter().map(|&x| x.1), &[]);


    fg.show();

}
