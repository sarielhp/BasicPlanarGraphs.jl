# BasicPlanarGraphs documentation

Planar graphs with a combinatorial embedding, plus the linear-time Lipton–Tarjan
planar separator, Left-Right planarity testing and embedding, straight-line grid
drawings, and a plain-text graph file format.

This directory describes the library. Every statement here is about the code in
`src/`; the tests that pin the behaviour are in `test/`.

Like the rest of the repository, these pages were written by an AI coding
assistant working with the maintainer; where a page and a test disagree, the test
is right. See [the README](https://github.com/sarielhp/BasicPlanarGraphs.jl#this-code-was-written-by-an-ai)
for what that means in practice.

## Contents

| Page | What it covers |
|---|---|
| [Embedding](embedding.md) | The `PlanarGraph` type, rotation systems, faces, duality, Euler's formula, and the dart core underneath. |
| [Planarity and embedding](planarity.md) | `is_planar_linear`, `planar_embedding`, what an embedding guarantees, and the Left-Right test behind both. |
| [Grid drawing](drawing.md) | `grid_drawing`, the `GridDrawing` type, the grid-size bounds, and the Fraysseix–Pach–Pollack placement. |
| [Separator](separator.md) | What `separator` guarantees, the Lipton–Tarjan stages, refinement, disconnected inputs, complexity, and caveats. |
| [Generators and triangulation](generators.md) | `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation`, and `triangulate`. |
| [Graph files](io.md) | The sectioned text format, `write_planar_graph`, `read_planar_graph`, `write_graph`, `read_graph`, and what is rejected. |
| [API reference](api.md) | Every exported name, with signatures, return types, and errors. |
| [Development](development.md) | Building, the quality gate (including the examples), the test suite, the planarity example on real data, and CI. |

## Installation and quick start

Requires Julia 1.10 or later. The package is not registered, so install it from
the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/sarielhp/BasicPlanarGraphs.jl")
```

```julia
using Graphs, BasicPlanarGraphs

g = random_triangulation(1000)
A, B, C = separator(g)

isempty(intersect(A, B))            # true: the three parts are disjoint
maximum([length(A), length(B)])     # at most 2/3 n
length(C)                           # O(sqrt(n))
```

Planarity, embeddings, drawings and files:

```julia
is_planar_linear(complete_graph(5))                # false
e = planar_embedding(smallgraph(:dodecahedral))    # PlanarGraph, nfaces(e) == 12
d = grid_drawing(e)                                # GridDrawing, coordinates in O(n) x O(n)
write_planar_graph("dodecahedron.txt", e)          # text file, rotation included
read_planar_graph("dodecahedron.txt") == e         # true
```

## In one paragraph

A [`PlanarGraph`](embedding.md) is a `Graphs.jl` `SimpleGraph` together with a
rotation system: for each vertex, the counter-clockwise cyclic order of its
neighbours. The constructor validates
that the rotation is a permutation of the neighbours, so the stored data always
describes a graph, and it rejects parallel edges and self-loops. On top of that
representation the package offers face walks, the dual, Euler-formula
validation, triangulation, a handful of embedded generators, and
`separator`, which splits the vertex set into `A`, `B` and `C` with no edge
between `A` and `B`, both sides at most `2/3 n`, and `|C| = O(sqrt(n))`.

Three questions the embedding alone does not answer are handled too:
`is_planar_linear` decides planarity in linear time (iteratively, so a million-vertex
road network is fine), `planar_embedding` computes a
rotation system for a graph that has none (or repairs one that is not an
embedding), and `grid_drawing` turns an embedding into a straight-line drawing
with integer coordinates in an `O(n) × O(n)` grid. `write_planar_graph` and
`read_planar_graph` store a graph and its embedding as text, and `write_graph`
and `read_graph` do the same for a bare `SimpleGraph`.

The internal working representation is a dart graph (`DartGraph`, `src/embedded.jl`),
an embedded multigraph stored as a rotation system over oriented edges. It is
described in [Embedding](embedding.md#the-dart-core); it is not exported, and
the public type converts to it and back at the boundary.

## Scope and limitations

- Graphs are **simple**: no self-loops, no parallel edges. The dart layer could
  represent them, the public type rejects them.
- Vertices are the integers `1:n`; the element type is parametric
  (`PlanarGraph{Int32}` works) and is inferred from the rotation.
- The package draws but does not render: `grid_drawing` returns integer
  coordinates in an `O(n) × O(n)` grid, and turning them into lines is the
  caller's business. There is no layout for a fixed canvas, no edge routing, and
  no curved edges.
- The embedding of a `PlanarGraph` is not checked for planarity on construction;
  `is_planar_embedding` reports it and `planar_embedding` repairs it.
- `separator` takes a `PlanarGraph`. It uses only adjacency; the embedding
  influences BFS traversal order and therefore only tie-breaking.
