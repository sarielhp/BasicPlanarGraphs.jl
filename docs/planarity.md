# Planarity and embedding

Two questions, one algorithm. *Is this graph planar?* and *what is a rotation
system for it?* are both answered in linear time by the Left-Right planarity
test, which is what [`is_planar`](api.md) (from `Graphs.jl`) and
[`planar_embedding`](api.md) use.

```julia
using Graphs, BasicPlanarGraphs

is_planar(complete_graph(4))        # true
is_planar(complete_graph(5))        # false
is_planar(smallgraph(:petersen))    # false

g = planar_embedding(smallgraph(:dodecahedral))
nfaces(g)                           # 12
is_planar_embedding(g)              # true
```

| Call | Returns | Time |
|---|---|---|
| `is_planar(g)` | `Bool` | `O(nv + ne)` |
| `is_planar(g::PlanarGraph)` | `Bool` for the underlying graph | `O(nv + ne)` |
| `planar_embedding(g::SimpleGraph)` | `PlanarGraph` or `nothing` | `O(nv + ne)` |
| `planar_embedding(g::PlanarGraph)` | `g` or a fresh `PlanarGraph`, or `nothing` | `O(nv + ne)` |

## Testing planarity

`is_planar` is `Graphs.is_planar`, imported and re-exported under the same name,
so a program that loads both packages sees a single function. It is the
[Left-Right planarity criterion](https://en.wikipedia.org/wiki/Left-right_planarity_test)
of Brandes (2009): an iterative DFS orients the edges, the lowpoints of the
oriented edges give a nesting depth, and a second DFS maintains a stack of
conflict intervals that decides, for every return edge, whether it can be routed
on the left or on the right of the DFS tree. The graph is planar exactly when
that partition can always be built.

The test reads adjacency only. For a `PlanarGraph` it answers about the graph
and ignores the stored rotation, so a broken embedding still yields the correct
verdict:

```julia
k4 = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])
is_planar(k4)            # true: the graph is K4
is_planar_embedding(k4)  # false: that rotation gives only two faces
```

## Computing an embedding

`planar_embedding` runs the same Left-Right test — its output *is* the partition
of the return edges into left and right — and turns the result into the rotation
system stored by a [`PlanarGraph`](embedding.md). The left and right sides of
each return edge, read in DFS order around every vertex, give the
counter-clockwise order of the neighbours.

- The result is a `PlanarGraph` whose rotation system satisfies Euler's formula
  on every connected component, verified by `is_planar_embedding`. It is
  therefore accepted by [`separator`](separator.md), [`triangulate`](api.md) and
  [`grid_drawing`](drawing.md).
- A non-planar graph gives `nothing`; nothing is thrown, so `nothing` can be
  tested directly.
- Vertex labels are preserved, and the element type follows the input.
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

## Why the same algorithm twice

`is_planar` answers the decision problem and stops as soon as the conflict stack
reports a contradiction; `planar_embedding` needs the partition that the test
builds on the way. Both are one pass over the vertices and one over the edges
per DFS, with dictionaries keyed by oriented edges, so both are `O(nv + ne)` in
time and space on planar input. Neither recurses: each DFS keeps an explicit
stack, so a graph with a long path does not exhaust the Julia stack.

The construction follows the reference implementation in
[NetworkX](https://networkx.org/documentation/stable/reference/algorithms/planarity.html);
the tests compare our verdict with `Graphs.is_planar` on random graphs, all small
graphs, and the classical obstructions.

## What is not checked

The embedding is *combinatorial*: it fixes the cyclic order of the neighbours at
every vertex, and therefore the faces, but it says nothing about geometry. For
coordinates see [Grid drawing](drawing.md), which turns an embedding into a
straight-line drawing on an integer grid.
