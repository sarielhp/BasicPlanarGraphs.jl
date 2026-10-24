# BasicPlanarGraphs documentation

Planar graphs with a combinatorial embedding, plus the linear-time Lipton–Tarjan
planar separator.

This directory describes the library. Every statement here is about the code in
`src/`; the tests that pin the behaviour are in `test/`.

## Contents

| Page | What it covers |
|---|---|
| [Embedding](embedding.md) | The `PlanarGraph` type, rotation systems, faces, duality, Euler's formula, and the dart core underneath. |
| [Separator](separator.md) | What `separator` guarantees, the Lipton–Tarjan stages, refinement, disconnected inputs, complexity, and caveats. |
| [Generators and triangulation](generators.md) | `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation`, and `triangulate`. |
| [API reference](api.md) | Every exported name, with signatures, return types, and errors. |
| [Development](development.md) | Building, the quality gate, tests, and CI. |

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

The internal working representation is a dart graph (`DartGraph`, `src/embedded.jl`),
an embedded multigraph stored as a rotation system over oriented edges. It is
described in [Embedding](embedding.md#the-dart-core); it is not exported, and
the public type converts to it and back at the boundary.

## Scope and limitations

- Graphs are **simple**: no self-loops, no parallel edges. The dart layer could
  represent them, the public type rejects them.
- Vertices are the integers `1:n`; the element type is parametric
  (`PlanarGraph{Int32}` works) and is inferred from the rotation.
- The package does no drawing or layout: an embedding is stored, not rendered.
- `separator` takes a `PlanarGraph`. It uses only adjacency; the embedding
  influences BFS traversal order and therefore only tie-breaking.
