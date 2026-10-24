# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.2] - 2026-10-24

### Added

- `separator_partition(g, r)` and `SeparatorPartition`: the recursive separator
  partition (a global separator plus parts of at most `r` vertices with no edge
  between parts, `|S| = O(n/sqrt(r))`, `O(1 + n/r)` parts, `O(n (1 + log(n/r)))`
  time, `O(n)` space), with optional induced part graphs.

### Changed

- `separator` is now the exact Lipton–Tarjan algorithm: every component left by
  the separator has at most `2n/3` vertices and `|C| <= 2*sqrt(2n)`, for every
  input, in `O(n)` worst-case time. The heuristic cycle candidates and the
  refinement loop (which could be exponential in principle and gave no balance
  guarantee) are gone. Faces with four or more sides get virtual star vertices,
  so no triangulation of the input is needed and graphs with bridges and cut
  vertices are handled; the balanced cycle is found by an exact walk down the
  dual tree.
- A pruning pass moves redundant separator vertices into `A` or `B`: paths give
  one separator vertex, cycles two, and `k x k` grids `k - 1`.
- Matching the twin darts of a rotation system (`PlanarGraph` to `DartGraph`,
  `induced`) uses two counting sorts instead of a hash table: `O(n + m)` in the
  worst case.
- `separator` requires `balance` in `[2/3, 1)` (it used to accept `(0, 1)` and
  silently return an unbalanced result below `2/3`), and rejects a rotation
  system that is not a planar embedding with an `ArgumentError`.

## [0.2.1] - 2026-10-24

### Changed

- Release numbering: v0.2.1 supersedes the withdrawn v0.2.0 tag. No functional
  changes since v0.1.1.

## [0.1.1] - 2026-10-24

### Changed

- Internal ordering no longer uses comparison sorts. The Left-Right test orders
  each vertex's outgoing darts by nesting depth with two stable counting passes,
  the separator orders dual edges and components with counting sorts, and vertex
  subsets are produced by a mark-and-scan over `1:nv`. `is_planar_linear`,
  `planar_embedding`, `separator`, `cycle_sides` and rotation validation now meet
  their documented `O(nv + ne)` bound.

### Fixed

- `separator` returns both sides sorted for disconnected graphs, as documented,
  instead of leaving the smaller components in discovery order.

## [0.1.0] - 2026-10-24

First release.

### Added

- `PlanarGraph{T<:Integer}`: a `Graphs.jl` `SimpleGraph` plus a rotation system
  (counter-clockwise neighbour order per vertex). Both constructors validate the
  rotation against the edges, reject self-loops, parallel edges and non-integer
  entries, and infer the vertex type. `planar_graph` is an alias.
- Embedding queries: `neighbor_order`, `rotation`, `faces`, `nfaces`,
  `is_planar_embedding` (Euler's formula per component), and `dual`.
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
- `triangulate`: chord insertion inside faces, verified before it returns;
  raises `ArgumentError` for faces that are not simple cycles and for rotation
  systems that are not planar.
- `separator(g; balance = 2/3)`: the linear-time Lipton–Tarjan vertex separator.
  Returns a sorted partition `(A, B, C)` with no `A`–`B` edge,
  `max(|A|, |B|) <= balance * n`, and `|C| = O(sqrt(n))`. Handles disconnected
  and degenerate inputs, tries several BFS roots, and refines until balanced.
- `cycle_separator`: a balanced separator that is a **simple cycle of the input
  graph**, computed by the algorithm of Har-Peled and Nayyeri; its phases and the
  shortcut are linear in the size of the internal triangulation. The chords
  added by the triangulation are shortcut back to boundary
  arcs of the faces of the input — a step performed neither by the paper nor by
  OGDF's `SeparatorHarPeled`, which returns only a vertex set — so consecutive
  vertices of the result, including the last and the first, are adjacent in the
  graph that was passed in. Both sides hold at most `2/3` of the vertices, and
  the cycle is oriented with the smaller side on its left.
- `cycle_sides(g, cycle)`: the two sides of any simple cycle of a connected
  `PlanarGraph`. Because rotations are counter-clockwise, the vertices to the
  left of the traversal are the inside and those to the right the outside;
  reversing the cycle swaps the two. The cycle itself is in neither side, and
  both are sorted.
- `grid_drawing`, with the `GridDrawing` type and the queries `coordinates`,
  `width` and `height`: a straight-line crossing-free drawing with integer
  coordinates, computed in linear time by the Fraysseix–Pach–Pollack placement
  after augmenting the embedding to a triangulation. A connected graph with
  `n >= 3` fits in `2(n-2) × (n-2)`; disconnected graphs are drawn component by
  component and packed.
- Generators: `planar_path`, `planar_cycle`, `planar_grid` (optionally
  triangulated by diagonals), `planar_wheel`, and `random_triangulation`.
- Graph files: `write_planar_graph`, `read_planar_graph`, `write_graph` and
  `read_graph`, a sectioned plain-text format with an optional `vertices`
  directive, an `edges` section and a `rotation` section. Reading validates the
  syntax, rejects duplicate and self-loop edges, requires the two sections to
  describe the same graph, and reports errors as `source:line: message`.
- Documentation: new pages on planarity and embedding, grid drawing, the cycle
  separator, and graph files; the API reference and the README cover the new
  names. The README now states that the repository was written by an AI coding
  assistant and what that implies (tests are the authority, independent
  reimplementations of the cited algorithms, no warranty), and credits the
  papers and reference implementations consulted, with their licenses.
- `examples/test_planarity.jl`: an end-to-end run on a real road network from
  SNAP (TIGER/Line Pennsylvania, ~1.09M vertices). It downloads and caches the
  dataset, cleans the directed listing into a pure undirected graph, and for
  1000, 10,000, 100,000 induced vertices and the whole graph reports planarity,
  the computed embedding, the separator and its verification, with timings and
  allocated memory. `tools/gate` runs it offline on a small local graph, so a
  broken example fails the gate.
- Tooling: `tools/build` (load the package), `tools/gate` (audit, static
  analysis, tests with Aqua, one build), `tools/check` as its alias, and the
  running-time survey in `benchmark/` (audited and format-checked by the gate,
  never run by it).
- CI: `Pkg.test()` on Linux (x64) for Julia 1.10 and the latest release, with
  coverage.

### Fixed

- `triangulate` now triangulates each simple face with a single fan: it finds a
  vertex whose every fan chord is still free and inserts all `k - 3` chords at
  once, instead of splitting the face one chord at a time. This removes the
  repeated copying of the face boundary that made the original implementation
  quadratic (5.5 s and 18.7 GiB of allocation for a 50,000-vertex cycle). The
  chord insertion is 2-4x faster than the balanced split it replaces (57 ms to
  25 ms for that cycle), so `triangulate` takes about 0.09 s and
  `cycle_separator` drops from 5.6 s to about 0.11 s on the same input. A face
  whose boundary repeats a vertex is cut at the second appearance of that
  vertex (connecting the previous and next boundary vertices) before the fan is
  inserted. A fan chord or a repeat chord that already exists is never added:
  the next free apex or another chord is used instead, so no parallel edge can
  appear.

### Changed

- `Graphs` compatibility raised to `1.15`, the first release with
  `Graphs.is_planar`.

[Unreleased]: https://github.com/sarielhp/BasicPlanarGraphs.jl/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.2
[0.2.1]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.1
[0.1.1]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.1.0
