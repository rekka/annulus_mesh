# annulus_mesh

[![Build Status](https://travis-ci.org/rekka/annulus_mesh.svg)](https://travis-ci.org/rekka/annulus_mesh)

This [Rust] crate generates a triangulation of an annulus by first
generating nodes using [Bridson’s algorithm][BridsonVisualization] for
Poisson-disc sampling, and then creating a [Delaunay
triangulation][Delaunay].

![Example](d=0.1-r=0.2-l=0.svg)

## Usage

Install the latest stable version of [Rust] from the [website][Rust].

Build using cargo:
```bash
cargo build --release
```

Print command line options:
```
target/release/annulus_mesh -h
```

[Rust]: https://www.rust-lang.org/
[BridsonVisualization]: http://bl.ocks.org/mbostock/dbb02448b0f93e4c82c3
[Delaunay]: https://en.wikipedia.org/wiki/Delaunay_triangulation
