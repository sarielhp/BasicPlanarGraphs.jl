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

## Planarity and embedding

### `is_planar(g) -> Bool`

Whether `g` is planar. This is `Graphs.is_planar`, imported and re-exported under
the same name (a program that loads both packages sees one function), so it
accepts any graph `Graphs.jl` understands. It is the Left-Right test of Brandes
and runs in `O(nv + ne)` time. For a `PlanarGraph` it asks about the graph and
ignores the stored rotation. See [Planarity and embedding](planarity.md).

### `planar_embedding(g) -> Union{PlanarGraph, Nothing}`

A planar embedding of `g` as a `PlanarGraph`, or `nothing` when `g` is not
planar. `O(nv + ne)` time. The rotation system satisfies Euler's formula on every
connected component, so the result is accepted by `separator`, `triangulate` and
`grid_drawing`. Vertex labels are preserved and the answer is deterministic.

For a `PlanarGraph` whose stored rotation is already a planar embedding the
argument itself is returned; otherwise the embedding is recomputed from the
underlying graph, which repairs a rotation system that came from a file.

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

## Drawing

### `grid_drawing(g) -> GridDrawing`

Straight-line crossing-free drawing of a planar graph on the integer grid, with
integer coordinates, a width of at most `2n - 4` and a height of at most `n - 2`
for a connected graph with `n >= 3`; disconnected graphs are packed component by
component. `O(nv + ne)` time. A `SimpleGraph` gets a computed embedding, a
`PlanarGraph` its stored one, recomputed when it is not an embedding; a
`ArgumentError` is thrown when the graph is not planar. Deterministic. See
[Grid drawing](drawing.md).

### `GridDrawing`

The result of `grid_drawing`: the coordinates and the bounding box.

| Query | Returns |
|---|---|
| `d[v]` | `(x, y)` of vertex `v`, a `Tuple{Int,Int}` |
| `coordinates(d)` | `Vector{Tuple{Int,Int}}` in vertex order |
| `width(d)`, `height(d)` | `Int`, the dimensions of the bounding box |
| `length(d)` | `Int`, the number of vertices |

`d[v]` and `coordinates(d)` throw the usual `BoundsError` for a vertex outside
`1:nv`. The fields `d.x`, `d.y`, `d.width` and `d.height` hold the same data; the
queries are the supported interface.

### `coordinates(d::GridDrawing) -> Vector{Tuple{Int,Int}}`

The positions of all vertices, indexed by vertex label.

### `width(d)`, `height(d) -> Int`

Dimensions of the bounding box: every coordinate satisfies
`0 <= x <= width(d)` and `0 <= y <= height(d)`.

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

## Graph files

The format is described in [Graph files](io.md).

### `write_planar_graph(io, g)`, `write_planar_graph(path, g)`

Write a `PlanarGraph` as a sectioned text file: `vertices`, `edges`, and the
rotation system, ended by `end`. Returns the `IO` when given one. The file reads
back with `read_planar_graph` into an equal graph.

### `write_graph(io, g)`, `write_graph(path, g)`

The same format for a `SimpleGraph`, without a rotation section.

### `read_planar_graph(io) -> PlanarGraph`, `read_planar_graph(path) -> PlanarGraph`

Read a graph and its embedding. The rotation goes through the `PlanarGraph`
constructor, and when the file has both sections they must describe the same
graph; otherwise `ArgumentError`. `ArgumentError` also when the file has no
`rotation` section. Planarity is not checked: an embedding that violates Euler's
formula is read as stored.

### `read_graph(io) -> SimpleGraph`, `read_graph(path) -> SimpleGraph`

The graph of a file, without its embedding. The edges are those of the `edges`
section, or, when the file has only a rotation section, the edges that rotation
induces.

Both readers take an optional second argument naming the source for error
messages (default `<io>`); errors are `ArgumentError` with the form
`source:line: message`.

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
| `is_planar` | decision | `Bool` |
| `planar_embedding` | algorithm | `PlanarGraph` or `nothing` |
| `grid_drawing`, `GridDrawing` | drawing | `GridDrawing` |
| `coordinates`, `width`, `height` | drawing query | positions and bounding box |
| `dual` | operation | `PlanarGraph` |
| `triangulate` | operation | `PlanarGraph` |
| `separator` | algorithm | `Tuple{Vector,Vector,Vector}` |
| `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation` | generators | `PlanarGraph` |
| `write_planar_graph`, `write_graph` | file output | `IO` or `nothing` |
| `read_planar_graph`, `read_graph` | file input | `PlanarGraph` / `SimpleGraph` |
