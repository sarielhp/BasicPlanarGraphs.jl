# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `is_planar`: the Left-Right planarity test of Brandes, imported from `Graphs.jl`
  and re-exported, so it is one function whichever package is loaded. For a
  `PlanarGraph` it asks about the graph, not the rotation. Note that the `Graphs.jl`
  implementation recurses once per DFS edge and is quadratic; see
  `is_planar_linear`.
- `is_planar_linear`: the same criterion, in this package's own iterative
  implementation. Both DFS passes run on an explicit stack and every stage is one
  pass over vertices and darts, so it is `O(nv + ne)` in time and space and
  survives graphs with long paths, where `is_planar` overflows the stack at about
  10^5 vertices and is quadratic.
- `planar_embedding`: computes a rotation system for a planar graph in linear
  time, or returns `nothing` when there is none. The result satisfies Euler's
  formula on every component, so it feeds `separator`, `triangulate` and
  `grid_drawing`; a `PlanarGraph` whose stored rotation is not an embedding gets
  a recomputed one, and a valid one is returned unchanged.
- `grid_drawing`, with the `GridDrawing` type and the queries `coordinates`,
  `width` and `height`: a straight-line crossing-free drawing with integer
  coordinates, computed in linear time by the Fraysseix–Pach–Pollack placement
  after augmenting the embedding to a triangulation. A connected graph with
  `n >= 3` fits in `2(n-2) × (n-2)`; disconnected graphs are drawn component by
  component and packed.
- Graph files: `write_planar_graph`, `read_planar_graph`, `write_graph` and
  `read_graph`, a sectioned plain-text format with an optional `vertices`
  directive, an `edges` section and a `rotation` section. Reading validates the
  syntax, rejects duplicate and self-loop edges, requires the two sections to
  describe the same graph, and reports errors as `source:line: message`.
- Documentation: new pages on planarity and embedding, grid drawing, and graph
  files; the API reference and the README cover the new names.
- `examples/test_planarity.jl`: an end-to-end run on a real road network from
  SNAP (TIGER/Line Pennsylvania, ~1.09M vertices). It downloads and caches the
  dataset, cleans the directed listing into a pure undirected graph, and for
  1000, 10,000, 100,000 induced vertices and the whole graph reports planarity,
  the computed embedding, the separator and its verification, with timings and
  allocated memory. `tools/gate` runs it offline on a small local graph, so a
  broken example fails the gate.

### Changed

- `Graphs` compatibility raised to `1.15`, the first release with
  `Graphs.is_planar`.

## [0.1.0] - 2026-10-24

First release.

### Added

- `PlanarGraph{T<:Integer}`: a `Graphs.jl` `SimpleGraph` plus a rotation system
  (counter-clockwise neighbour order per vertex). Both constructors validate the
  rotation against the edges, reject self-loops, parallel edges and non-integer
  entries, and infer the vertex type. `planar_graph` is an alias.
- Embedding queries: `neighbor_order`, `rotation`, `faces`, `nfaces`,
  `is_planar_embedding` (Euler's formula per component), and `dual`.
- `triangulate`: chord insertion inside faces, verified before it returns;
  raises `ArgumentError` for faces that are not simple cycles and for rotation
  systems that are not planar.
- `separator(g; balance = 2/3)`: the linear-time Lipton–Tarjan vertex separator.
  Returns a sorted partition `(A, B, C)` with no `A`–`B` edge,
  `max(|A|, |B|) <= balance * n`, and `|C| = O(sqrt(n))`. Handles disconnected
  and degenerate inputs, tries several BFS roots, and refines until balanced.
- Generators: `planar_path`, `planar_cycle`, `planar_grid` (optionally
  triangulated by diagonals), `planar_wheel`, and `random_triangulation`.
- Documentation in `docs/`: embedding, separator, generators and triangulation,
  API reference, and development notes.
- Tooling: `tools/build` (load the package), `tools/gate` (audit, static
  analysis, tests with Aqua, one build), and `tools/check` as its alias.
- CI: `Pkg.test()` on Linux (x64) for Julia 1.10 and the latest release, with
  coverage.

[Unreleased]: https://github.com/sarielhp/BasicPlanarGraphs.jl/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.1.0
