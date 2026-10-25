# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `separator_partition(g, r; weights = w)`: parts bounded by weight (see
  `docs/partition.md`). Without `weights` the result is unchanged.
- `add_edge_in_face!` joins two components when the corners share no face because the
  vertices are in different components.

### Fixed

- `cycle_sides` with a non-integer vector (for example `Float64`) threw a `MethodError`;
  it throws `ArgumentError` now.
- The `triangulate` test compared rotations up to reversal; it now requires the same
  direction, which `triangulate` preserves.

## [0.2.8] - 2026-10-25

### Added

- `PlanarGraph{T}` is now an `AbstractGraph{T}`: `connected_components`, `bfs_parents`,
  `dijkstra_shortest_paths` and the other read-only `Graphs.jl` algorithms accept it.
  `Graphs.rem_edge!` works. Functions that build a new graph of the same type
  (`induced_subgraph`, `complement`) are not supported.
- `delete_edge!`, `contract_edge!` and `add_edge_in_face!`: in-place edits that keep
  the embedding.
- `r_division(g, r)` and `RDivision`: regions of at most `r` vertices (`r >= 100`) with
  at most `12 sqrt(r)` boundary vertices each; the number of regions and the total
  boundary are measured, not proved (`docs/r_division.md`).

## [0.2.7] - 2026-10-25

### Added

- `separator(g; weights = w)`: balance by vertex weight (nonnegative integers). The
  weights of both sides and of every component of `g - C` are at most `balance` of
  the total; `|C| <= 2*sqrt(2n)` still counts vertices. Results without `weights`
  are unchanged (checked on 63 fixed cases). `separator_partition` stays unweighted.

## [0.2.6] - 2026-10-25

### Added

- `five_coloring(g)`: a proper vertex coloring of a planar `SimpleGraph` or
  `PlanarGraph` with the colors `1:5` (Matula–Shiloach–Tarjan), and
  `is_proper_coloring(g, colors)`. Near linear in practice; linear time is not
  proved. Non-planarity is reported only when the method gets stuck.

## [0.2.5] - 2026-10-25

### Fixed

- `cycle_separator` on a triangulation returned a cycle of length proportional to
  the BFS depth (`2L + 4` on a triangulated tube of length `L` with a large end)
  when every rung of the ladder was light, because the ring outside the last rung,
  a case of the paper, was missing. The `O(sqrt(n))` length guarantee now holds on
  that family (214 instead of 3 204 vertices at `L = 1600`); the tests now reach
  phases 2 and 3 with a small step.
- `separator` and `separator_partition` reject a rotation system that is not planar
  in **any** component (they checked only the component that was cut, and not at
  all for `r >= nv`).
- `cycle_separator` and `cycle_sides` take linear time when the walk returns to a
  high-degree vertex many times (a star, a windmill): they were quadratic.
- `triangulate` takes linear expected time when two adjacent hubs are opposite on
  many faces: it was quadratic. The edges are now kept in a hash set. The documented
  degree bound distinguishes simple faces (at most 5) from faces that repeat a vertex.
- Documentation: `triangulate` on disconnected input, the winding-value condition of
  `cycle_sides`, the cost of walks supplied by the caller.

## [0.2.4] - 2026-10-24

### Fixed

- `hash(::PlanarGraph)` agrees with `==` (equal embeddings with shifted or
  reversed rotations hashed differently).
- `dual` checks its preconditions and throws a clear `ArgumentError` for a dual with
  parallel edges or self-loops, a disconnected or non-planar input, and a face
  count that overflows the vertex type; its documentation no longer says that
  bridge-free is enough.
- `PlanarGraph(graph, rotation)` throws `ArgumentError` for non-integer entries;
  `add_dart!` validates its anchor before changing the graph;
  `random_triangulation` accepts any `Integer`.
- The graph readers refuse a vertex count or label above `max_vertices` (default
  10 000 000, a keyword) before allocating, honour the `vertices` line (a larger
  label or a second `vertices` line is an error), and report a rotation line
  without `:` and a vertex without a rotation line as `ArgumentError`.
- `faces`, `triangulate` and `planar_embedding` keep the vertex type of their input
  (`Int32` stays `Int32`), `PlanarGraph(graph, rotation)` copies the graph, and
  `planar_embedding` and `is_planar_linear` throw `ArgumentError` for directed graphs.
- The benchmark `phases` suite works again and the gate runs the benchmark once
  (`--quick`); `benchmark/README.md` has current numbers. The planarity example runs
  each distinct size once and exits non-zero when a planar result fails its check.
- Documentation: test count, the `Graphs.jl` functions a `PlanarGraph` accepts, the
  list of exports, the `write_*` return value, the packing of components in drawings.
- `grid_drawing` and `planar_embedding` are linear for graphs with many components
  (they were quadratic).

### Changed

- `triangulate` cuts each face by a ladder (vertex `i` joined to vertex `k+1-i`) and
  closes the quadrangles with a diagonal, so a vertex gains a constant number of
  edges per face; regions left by blocked rungs get a slanted ladder. A vertex gains at most 5
  chords per face. Still `O(n + m)`; no edge is added twice.

## [0.2.3] - 2026-10-24

### Changed

- `cycle_separator` returns a closed walk of the input and no longer requires a
  2-connected graph: bridges, cut vertices and large faces are accepted. Both
  sides are at most `2/3` of the vertices for every connected input, and the time
  is `O(n)` (a graph with one huge face used to take quadratic time in
  `triangulate`, and the chord shortcut could throw on valid input). On a
  triangulation the result is still the paper's simple `O(sqrt(n))` cycle, now
  checked and backed by an exact fallback; on other graphs the exact balanced
  walk of the virtual star triangulation is lifted by replacing the star edge
  with a face-boundary arc.
- `cycle_sides` accepts closed walks that repeat vertices and edges (by winding
  numbers) and rejects walks that cross themselves.
- The dual spanning tree of the paper's first phase is built without union-find.

- `triangulate` is `O(n + m)` in the worst case. It used to search each face for a
  fan apex whose chords were all free, which took quadratic time on a face whose
  vertices have chords elsewhere (a ladder, `K_{2,n}`), and it hashed adjacency.
  It now goes around each face clipping ears: the chord between the neighbours of
  a position is added unless they are adjacent, and then the next position is
  tried (among four consecutive vertices one of the two chords is always new, by
  planarity); repeated occurrences of a cut vertex are clipped the same way. After
  a clip the next position is skipped, so a face's chords are spread over its
  vertices (degree `O(log k)` per vertex for a face of length `k`) instead of
  fanning from one. No edge is added twice, and the same `ArgumentError`s are
  raised for a non-planar rotation system or a single-edge component.

### Removed

- The chord shortcut of `cycle_separator` and the internal
  `_triangulate_with_regions`.

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

[Unreleased]: https://github.com/sarielhp/BasicPlanarGraphs.jl/compare/v0.2.8...HEAD
[0.2.8]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.8
[0.2.7]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.7
[0.2.6]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.6
[0.2.5]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.5
[0.2.4]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.4
[0.2.3]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.3
[0.2.2]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.2
[0.2.1]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.2.1
[0.1.1]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/sarielhp/BasicPlanarGraphs.jl/releases/tag/v0.1.0
