# Grid drawing

An embedding fixes the faces; a drawing fixes where the vertices go.
[`grid_drawing`](api.md) takes a planar graph and returns integer coordinates for
its vertices such that straight segments between adjacent vertices never cross.

```julia
using Graphs, BasicPlanarGraphs

d = grid_drawing(planar_wheel(8))
d                     # GridDrawing(9 vertices, 14x7 grid)
d[1]                  # (0, 0)
coordinates(d)[1:3]   # [(0, 0), (14, 0), (7, 7)]
```

| Call | Returns | Time |
|---|---|---|
| `grid_drawing(g::PlanarGraph)` | `GridDrawing` | `O(nv + ne)` |
| `grid_drawing(g::SimpleGraph)` | `GridDrawing` | `O(nv + ne)` |

`d[v]` is the position of vertex `v` as a tuple, `coordinates(d)` the whole
vector of positions in vertex order, `length(d)` the number of vertices, and
[`width`](api.md) and [`height`](api.md) the dimensions of the bounding box.

## What the drawing guarantees

For a graph with `n = nv(g)` vertices:

- every position is a pair of integers in `0:width(d) × 0:height(d)`;
- distinct vertices get distinct positions, and no vertex lies in the interior
  of an edge;
- two edges meet only at a shared endpoint;
- for a connected graph with `n >= 3`, `width(d) <= 2n - 4` and
  `height(d) <= n - 2`, so the drawing fits in a `2(n-2) × (n-2)` grid and the
  area is `O(n) × O(n)`;
- a triangulation uses exactly that box: `2(n-2)` wide and `n-2` high.

Disconnected graphs are drawn one component at a time, each component in its own
box, and the boxes are packed left to right with two free columns between them;
the width is then the sum over components and the height the tallest component,
so both stay `O(n)`. A component with one vertex is a point at height `0`, one
with two vertices is a horizontal segment of length `1`.

```julia
coordinates(grid_drawing(SimpleGraph(1)))    # [(0, 0)]
coordinates(grid_drawing(planar_path(4)))    # [(0, 0), (4, 0), (2, 2), (2, 1)]
grid_drawing(SimpleGraph(0))                 # GridDrawing(0 vertices, 0x0 grid)
```

## The algorithm

The placement is the one of de Fraysseix, Pach and Pollack, in the formulation of
OGDF's `FPPLayout`, which is what the implementation follows:

1. The embedding is augmented to a maximal planar graph by inserting chords
   inside faces, using the same `_add_chords!` as [`triangulate`](api.md). The
   added edges are not part of the drawing; they only make every face a triangle,
   which the placement needs.
2. A **shelling order** `v1, …, vn` is computed by peeling vertices off the outer
   face: at every step the outer face is a cycle, and vertices with no chord to
   another outer vertex can be removed. Vertex `v1` comes first, `v2` second, and
   each removed vertex gets the next number from the end.
3. The vertices are placed in order. `v1` goes to `(0, 0)`, `v2` to `(2n-4, 0)`,
   and each `vk` goes to the apex of a triangle over the edge of the current
   outer face that is associated with it. A level-by-level sweep then shifts
   every vertex to the right by the offset of its predecessor, which restores the
   `O(n)` width after the local placements.

Every step is a pass over vertices and darts, so the drawing is computed in linear
time; the augmentation adds only edges, never vertices.

## Which embedding

A `PlanarGraph` is drawn with its stored rotation system, so two embeddings of
the same graph generally give different drawings. When the stored rotation is not
a planar embedding, the embedding is recomputed with
[`planar_embedding`](planarity.md) and the graph is drawn with that one: the
drawing is always valid, and the input is left untouched.

```julia
k4 = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])  # not an embedding
is_planar_embedding(k4)      # false
grid_drawing(k4)             # GridDrawing(4 vertices, 4x2 grid)
is_planar_embedding(k4)      # still false: the argument was not modified
```

A `SimpleGraph` has no embedding, so one is computed first. Either way the graph
itself must be planar, and a non-planar graph raises `ArgumentError`:

```julia
grid_drawing(complete_graph(5))                 # ArgumentError
grid_drawing(complete_graph(4))                 # GridDrawing(4 vertices, 4x2 grid)
```

The drawing is deterministic: the same graph and rotation always give the same
coordinates, since every stage visits neighbours in the order the embedding
stores them.

## Rendering

The package produces coordinates only. Drawing the segments is up to the caller,
which keeps the dependency list short and lets the same coordinates go to
`Luxor`, `CairoMakie`, an SVG template, or a test:

```julia
g = planar_grid(3, 3; diagonals = true)
d = grid_drawing(g)
for e in edges(g)
    (x1, y1), (x2, y2) = d[src(e)], d[dst(e)]
    # draw the segment (x1, y1)-(x2, y2)
end
```

Y coordinates grow downwards in most plotting packages and upwards in this grid;
flip the sign of `y` when the orientation matters.
