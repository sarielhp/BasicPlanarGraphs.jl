# API reference

Every name exported by `BasicPlanarGraphs`, with its signature, what it returns,
and the errors it raises. Internal helpers (`DartGraph` and friends) are not
exported; they are described in [Embedding](embedding.md#the-dart-core).

```julia
using Graphs, BasicPlanarGraphs
```

## Types

### `PlanarGraph{T<:Integer}`

Undirected simple planar graph plus a combinatorial embedding.

| Constructor | Result |
|---|---|
| `PlanarGraph(graph::SimpleGraph{T}, rotation)` | graph plus embedding; the rotation is checked against the edges |
| `PlanarGraph(rotation)` | embedding plus the graph induced by it |
| `planar_graph(args...)` | alias for `PlanarGraph(args...)` |

`rotation[v]` lists the neighbours of `v` counter-clockwise. The element type is
inferred from the rotation (`Int` when it is empty).

Errors (`ArgumentError`): rotation length differs from `nv`; a rotation is not a
permutation of the neighbours of its vertex; a neighbour is out of range; a vertex
lists itself; a rotation entry is not an integer.

Equality is equality of graph and cyclic order (up to shift and reversal at each
vertex); `hash` agrees with it.

## Embedding queries

### `nv(g)`, `ne(g)`, `vertices(g)`, `edges(g)`, `neighbors(g, v)`, `degree(g, v)`, `has_edge(g, u, v)`

`Graphs.jl` interface, forwarded to the wrapped graph. Any `Graphs.jl` function
that only reads adjacency accepts a `PlanarGraph`.

### `neighbor_order(g, v) -> Vector{T}`

Neighbours of `v` in counter-clockwise order. This is the stored vector; mutating
it corrupts the embedding. Use `rotation(g)` for a copy.

### `rotation(g) -> Vector{Vector{T}}`

A fresh copy of the whole rotation system.

### `faces(g) -> Vector{Vector{T}}`

Boundary walks of the faces, one vector of vertex numbers per face. Vertices on
bridges or cut points repeat within a walk. An isolated vertex contributes no
walk, so a one-vertex graph returns an empty vector.

### `nfaces(g) -> Int`

`length(faces(g))`, without materialising the walks.

### `is_planar_embedding(g) -> Bool`

`true` when every connected component satisfies `n_c - m_c + f_c = 2` (a
component with no edges counts one face). This is exactly the planarity test for
a rotation system.

### `dual(g) -> PlanarGraph`

One vertex per face of `g`, one edge per primal edge; the dual is 3-regular when
`g` is a triangulation. Throws `ArgumentError` for a graph with a bridge, whose
dual edge would be a self-loop.

## Triangulation

### `triangulate(g) -> typeof(g)`

Adds chords inside faces until every face is a triangle, preserving vertices and
embedding. Raises `ArgumentError` when a face is not a simple cycle and cannot be
split, or when the rotation system is not planar. On a 2-connected input the
result is maximal planar (`3n - 6` edges).

## Separator

### `separator(g; balance = 2/3) -> (A, B, C)`

Vertex separator `C` with sides `A` and `B`:

- `A`, `B`, `C` partition `1:nv(g)` and are sorted;
- no edge joins `A` to `B`;
- `max(|A|, |B|) <= balance*nv(g)`;
- `|C| = O(sqrt(nv(g)))`, at most `2*sqrt(2 nv(g)) + 1` on the test families.

`balance` must be in `(0, 1)`; otherwise `ArgumentError`. Deterministic for a
given graph. See [Separator](separator.md).

## Generators

### `planar_path(n) -> PlanarGraph`

Path on `n >= 1` vertices. `ArgumentError` for `n < 1`.

### `planar_cycle(n) -> PlanarGraph`

Cycle on `n >= 3` vertices, counter-clockwise. `ArgumentError` for `n < 3`.

### `planar_grid(rows, cols; diagonals = false) -> PlanarGraph`

`rows x cols` grid, vertices numbered `(r-1)*cols + c`. With `diagonals = true`
every cell gets its top-left to bottom-right diagonal, so all inner faces are
triangles. `ArgumentError` unless both dimensions are positive.

### `planar_wheel(rim) -> PlanarGraph`

Wheel: hub `1`, rim `2:rim+1`. `ArgumentError` for `rim < 3`.

### `random_triangulation(n; rng = Random.default_rng()) -> PlanarGraph`

Random maximal planar graph on `n >= 3` vertices, built by inserting each new
vertex into a uniformly chosen face (a stacked triangulation).
`ArgumentError` for `n < 3`. Pass `rng` for reproducibility.

## Summary table

| Name | Kind | Returns |
|---|---|---|
| `PlanarGraph`, `planar_graph` | type / constructor | `PlanarGraph{T}` |
| `nv`, `ne`, `vertices`, `edges`, `neighbors`, `degree`, `has_edge` | `Graphs.jl` interface | graph data |
| `neighbor_order` | query | `Vector{T}` |
| `rotation` | query | `Vector{Vector{T}}` |
| `faces` | query | `Vector{Vector{T}}` |
| `nfaces` | query | `Int` |
| `is_planar_embedding` | query | `Bool` |
| `dual` | operation | `PlanarGraph` |
| `triangulate` | operation | `PlanarGraph` |
| `separator` | algorithm | `Tuple{Vector,Vector,Vector}` |
| `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation` | generators | `PlanarGraph` |
