# Planarity and embedding

Two questions, one algorithm. *Is this graph planar?* and *what is a rotation
system for it?* are both answered in linear time by the Left-Right planarity
test: [`is_planar_linear`](api.md) for the decision, [`planar_embedding`](api.md)
for the rotation system.

```julia
using Graphs, BasicPlanarGraphs

is_planar_linear(complete_graph(4))        # true
is_planar_linear(complete_graph(5))        # false
is_planar_linear(smallgraph(:petersen))    # false

g = planar_embedding(smallgraph(:dodecahedral))
nfaces(g)                           # 12
is_planar_embedding(g)              # true
```

| Call | Returns | Time |
|---|---|---|
| `is_planar_linear(g)` | `Bool` | `O(nv + ne)`, iterative |
| `is_planar(g)` | `Bool` (`Graphs.is_planar`) | quadratic, recursive |
| `planar_embedding(g::SimpleGraph)` | `PlanarGraph` or `nothing` | `O(nv + ne)`, iterative |
| `planar_embedding(g::PlanarGraph)` | `g` or a fresh `PlanarGraph`, or `nothing` | `O(nv + ne)` |

## The criterion

Both functions implement the
[Left-Right planarity criterion](https://en.wikipedia.org/wiki/Left-right_planarity_test)
of Brandes (2009). A DFS orients the edges; the lowpoints of the oriented edges
give each edge a nesting depth; a second DFS maintains a stack of conflict
intervals and decides, for every return edge, whether it is routed on the left or
on the right of the DFS tree. The graph is planar exactly when that partition can
always be built — and that same partition, read around every vertex, is a
rotation system, which is why the embedding costs no more than the decision.

## Two entry points, and why

`is_planar_linear` is this package's own implementation. Both DFS passes run on
an explicit stack, so a long path does not become a long Julia call stack, and
the graph is scanned a constant number of times, so `O(n + m)` is a real bound on
real instances.

`is_planar` is `Graphs.is_planar`, imported and re-exported so that a program
loading both packages sees a single function and code written against `Graphs.jl`
keeps working. It decides the same criterion, but its two DFS passes are
recursive and its nesting-depth sort rescans a dictionary once per vertex, which
makes it quadratic and — from roughly 10^5 vertices on a path — fatal:

```julia
is_planar_linear(path_graph(100_000))   # true, 0.3 s
is_planar(path_graph(100_000))          # StackOverflowError
```

The overflow is not a catchable error: Julia warns that the program state may be
corrupted, and it is. Use `is_planar_linear` at that size; keep `is_planar` for
interchange. The two agree on every graph small enough for both (the tests check
random graphs, the generators, and the classical obstructions).

Both answer about the graph, not about a stored rotation. For a `PlanarGraph`
`planar_embedding` returns a stored rotation that is already planar as it is and
recomputes a broken one, so a broken embedding still yields the correct verdict:

```julia
k4 = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])
is_planar_linear(k4)     # true: the graph is K4
is_planar_embedding(k4)  # false: that rotation gives only two faces
```

## Computing an embedding

`planar_embedding` runs the same test and turns the left/right partition into the
rotation system stored by a [`PlanarGraph`](embedding.md). The left and right
sides of each return edge, read in DFS order around every vertex, give the
counter-clockwise order of the neighbours.

- The result is a `PlanarGraph` whose rotation system satisfies Euler's formula
  on every connected component, verified by `is_planar_embedding`. It is
  therefore accepted by [`separator`](separator.md), [`triangulate`](api.md) and
  [`grid_drawing`](drawing.md).
- A non-planar graph gives `nothing`; nothing is thrown, so `nothing` can be
  tested directly.
- Vertices keep their numbers, and the result has the vertex type of the input
  (`PlanarGraph{Int32}` in, `PlanarGraph{Int32}` out).
- Disconnected graphs are handled one component at a time.
- The answer is deterministic: the same graph always gives the same rotation
  system, because the DFS visits neighbours in the order the graph stores them.

For a `PlanarGraph` the call is a check first: when the stored rotation is
already a planar embedding the argument is returned unchanged (not copied);
otherwise the embedding is recomputed from the graph, which is the way to repair
a rotation system that came from a file or from a previous experiment.

```julia
g = planar_grid(3, 3)
g === planar_embedding(g)        # true, nothing to do

k4 = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])
h = planar_embedding(k4)
h === k4                         # false: the rotation was rebuilt
is_planar_embedding(h)           # true
```

## Cost

`is_planar_linear` stops as soon as the conflict stack reports a contradiction.
`planar_embedding` needs the partition the test builds on the way, and then
materialises the rotations. Both are `O(nv + ne)` time and space: the state is a
handful of dictionaries keyed by oriented edges, at most `2·ne` entries each, and
every stage is one pass over vertices and darts. On a million-vertex road network
the embedding takes a few seconds and a few GiB, which is why the example program
in `examples/test_planarity.jl` prints timings and allocation for each size it
runs.

The construction follows the reference implementation in
[NetworkX](https://networkx.org/documentation/stable/reference/algorithms/planarity.html).

## What is not checked

The embedding is *combinatorial*: it fixes the cyclic order of the neighbours at
every vertex, and therefore the faces, but it says nothing about geometry. For
coordinates see [Grid drawing](drawing.md), which turns an embedding into a
straight-line drawing on an integer grid.
