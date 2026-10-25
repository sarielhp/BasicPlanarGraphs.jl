# BasicPlanarGraphs

[![CI](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Planar graphs with an abstract embedding (a rotation system): linear-time
planarity testing and embedding, the Lipton–Tarjan planar separator, straight-line
grid drawings, and a plain-text graph file format.

## This code was written by an AI

Everything in this repository — `src/`, `test/`, `docs/`, the example, the tooling
and this README — was written by a large language model through an AI coding
assistant, in sessions with the maintainer, who chose what to build, the
algorithms and the interfaces, read the results, asked for changes, and ran the
checks. It is a product of *vibe programming*: the code comes out of conversation
and iteration, not out of a human typing it.

Please weigh that when you use it:

- **The tests are the guarantee, not the prose.** Every claim in `docs/` should be
  backable by a test or by a run of the example; the suite (98,000+ tests,
  including an independent `Graphs.jl`-based validator and a brute-force optimal
  separator oracle) is what is actually checked, and `tools/gate` runs it. Where
  a document and a test disagree, the test wins.
- **It can be confidently wrong.** The AI has produced plausible-looking
  arguments that were wrong, and code that passed a review and failed on a
  corner case. The algorithms are reimplementations of published work and of the
  reference implementations credited below — not verbatim copies — but
  independent reimplementation is exactly the kind of task that hides bugs.
- **No warranty.** MIT says it plainly. Treat this as an unreviewed contribution
  until you have read the part you depend on.

The maintainer is the author of record and is responsible for what is here; the
AI is the one that typed it.

## Documentation

Full documentation lives in [`docs/`](docs/README.md):

| Page | Contents |
|---|---|
| [Embedding](docs/embedding.md) | the `PlanarGraph` type, rotation systems, faces, duality, Euler's formula, the dart core |
| [Planarity and embedding](docs/planarity.md) | `is_planar_linear`, `planar_embedding`, the Left-Right test, what an embedding guarantees |
| [Grid drawing](docs/drawing.md) | `grid_drawing`, `GridDrawing`, grid-size bounds, the Fraysseix–Pach–Pollack placement |
| [Separator](docs/separator.md) | guarantees, the Lipton–Tarjan stages, pruning, disconnected inputs, complexity |
| [Coloring](docs/coloring.md) | `five_coloring`, a proper 5-coloring of a planar graph, and `is_proper_coloring` |
| [Separator partition](docs/partition.md) | `separator_partition`, the recursive partition into parts of at most `r` vertices, guarantees and cost |
| [Cycle separator](docs/cycle_separator.md) | guarantees, the left/right convention, the Har-Peled–Nayyeri phases, lifting to closed walks |
| [Generators and triangulation](docs/generators.md) | `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation`, `triangulate` |
| [Graph files](docs/io.md) | the sectioned text format, the writers and readers, what is rejected |
| [API reference](docs/api.md) | every exported name with signatures, return types, and errors |
| [Development](docs/development.md) | build, the quality gate, the test suite, benchmarks, CI |

Release history is in [CHANGELOG.md](CHANGELOG.md).

## Installation

Requires Julia 1.10 or later. The package is not registered in the General
registry, so install it from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/sarielhp/BasicPlanarGraphs.jl")
```

## Usage

```julia
using Graphs, BasicPlanarGraphs

g = random_triangulation(1000)
A, B, C = separator(g)

isempty(intersect(A, B))            # true: the three parts are disjoint
maximum([length(A), length(B)])     # at most 2/3 n
length(C)                           # O(sqrt(n))
```

## Overview

`PlanarGraph{T}` is a `Graphs.jl` `SimpleGraph{T}` plus, for every vertex, the
counter-clockwise cyclic order of its neighbours:

```julia
triangle = PlanarGraph([[2, 3], [3, 1], [1, 2]])
nv(triangle), ne(triangle), nfaces(triangle)     # 3, 3, 2
faces(triangle)                                  # [[1, 2, 3], [1, 3, 2]]
dual(triangulate(planar_grid(3, 4)))             # 3-regular dual, Euler holds
```

The constructor rejects rotations that are not permutations of the neighbours, so
the stored data always describes a graph; `is_planar_embedding` then checks
Euler's formula on every component. Equality of embeddings is equality up to a
cyclic shift and reversal at each vertex. Details, including the dart-based core
underneath: [Embedding](docs/embedding.md).

`separator(g; balance = 2/3)` returns `(A, B, C)` with `A`, `B`, `C` a partition
of the vertices, no edge between `A` and `B`, `|A|, |B| <= balance * n`, and
`|C| <= 2*sqrt(2n)`. Both bounds hold for every input and the running time is
`O(n)` in the worst case; `balance` must lie in `[2/3, 1)`. The stages, the
guarantees and the proof sketch: [Separator](docs/separator.md).

`cycle_separator(g)` returns a balanced separator that is a **closed walk of `g`**
rather than a vertex set — consecutive vertices, including the last and the first,
are adjacent. `cycle_sides(g, C)` returns the vertices on the two sides;
rotations are counter-clockwise, so the left of the traversal is the inside, and
the walk is oriented with the smaller side inside. Both sides hold at most `2/3`
of the vertices and the running time is `O(n)`, for every connected graph. If
`g` is a triangulation the walk is a simple cycle with `O(sqrt(n))` vertices
(Har-Peled and Nayyeri); otherwise a fill edge of a triangulation becomes a path
along its face, so the walk may repeat vertices and edges and its length is not
bounded. Details: [Cycle separator](docs/cycle_separator.md).

```julia
g = planar_grid(20, 20)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)
maximum([length(inside), length(outside)]) <= 2 / 3 * nv(g)
```

Planarity is decided by the Left-Right test, `is_planar_linear`, in linear time
on an explicit stack (a road network with a million intersections is fine), and
`planar_embedding` turns the partition of the return edges into a rotation
system — or recomputes one for a stored rotation that is not an embedding:

```julia
is_planar_linear(complete_graph(5))                # false
e = planar_embedding(smallgraph(:dodecahedral))    # 20 vertices, 12 faces
```

`grid_drawing` draws an embedding in straight lines on an integer grid of area
`O(n) × O(n)` — at most `2(n-2)` wide and `n-2` high for a connected graph — and
`write_planar_graph`/`read_planar_graph` (or `write_graph`/`read_graph` for a
bare `SimpleGraph`) store it as text. Details:
[Planarity and embedding](docs/planarity.md), [Grid drawing](docs/drawing.md),
[Graph files](docs/io.md).

Generators, triangulation, and the full export list:
[Generators](docs/generators.md) and [API reference](docs/api.md).

[`examples/test_planarity.jl`](examples/test_planarity.jl) puts the pipeline
together on a real road network: it downloads the TIGER/Line road graph of
Pennsylvania, cleans it into a pure graph, and for 1000, 10,000, 100,000 induced
vertices and for all 1,088,092 of them reports the planarity verdict, the
computed embedding, and the separator with its verification, timing every stage.
On the whole graph that is 8.4 s for the planarity test, 8.5 s for the embedding
(454,218 faces), and 20.2 s for the separator, which returns `|C| = 2292` at a
balance of 49.9 %, with a peak of 4.0 GiB resident. It is described in
[Development](docs/development.md).

## Credits and sources

The algorithms are from the literature, and each Julia implementation was
written against a published reference implementation. Where a reference's
structure is mirrored — the iterative DFS stacks of the Left-Right test, the
shelling-order bookkeeping of the grid placement — the header of the
corresponding `src/` file names it.

| Algorithm | Source paper | Reference implementation consulted | Its license |
|---|---|---|---|
| Planarity testing and embedding (`src/planarity.jl`) | U. Brandes, *The Left-Right Planarity Test* (2009) | [NetworkX](https://networkx.org) `algorithms/planarity.py`: the iterative formulation of the two DFS passes, the nesting-depth rules and the conflict-stack operations | BSD-3-Clause |
| Straight-line grid drawing (`src/grid_drawing.jl`) | H. de Fraysseix, J. Pach, R. Pollack, *How to draw a planar graph on a grid* (Combinatorica 10, 1990) | [OGDF](https://ogdf.github.io) `src/ogdf/planarlayout/FPPLayout.cpp`: the shelling order over a candidate list, the level-by-level placement and the x-offset pass | GPL-2.0 or GPL-3.0 |
| Planar separator (`src/separator.jl`) | R. J. Lipton, R. E. Tarjan, *A Separator Theorem for Planar Graphs* (SIAM J. Appl. Math. 36(2), 1979) | — | — |
| Cycle separator (`src/cycle_separator.jl`, `src/cycle_walk.jl`) | S. Har-Peled, A. Nayyeri, *A Simple Algorithm for Computing a Cycle Separator* (Comput. Geom. 83, 2019) | [OGDF](https://ogdf.github.io) `src/ogdf/separator/SeparatorHarPeled.cpp`: the dual-tree fundamental cycle, the face levels, the nested level cycles and the ladder | GPL-2.0 or GPL-3.0 |
| `is_planar` (re-exported) | same criterion as above | [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) `src/planarity.jl` | BSD-2-Clause |

No source file was copied: the Julia code is an independent implementation, and
the licenses above describe the projects consulted, not the terms under which
this package is distributed (it is MIT, see [LICENSE](LICENSE)).

The example in [`examples/test_planarity.jl`](examples/test_planarity.jl) runs
on the SNAP road networks (J. Leskovec and A. Krevl, *SNAP Datasets*, Stanford
University), which are derived from the US Census Bureau's TIGER/Line files.

## Development

```bash
tools/build   # instantiate the environment and load the package
tools/gate    # audit, static analysis, examples, tests (with Aqua), build (~35 s)
```

Pushing does not block on CI: `tools/push` pushes and prints the run URL, and
`tools/check_ci` reports the verdict later (exit `0` green, `1` failed,
`2` pending, `3` unknown). The test suite, the optimal-separator oracle, the
style rules and the CI matrix are described in [Development](docs/development.md).

## License

MIT — see [LICENSE](LICENSE).
