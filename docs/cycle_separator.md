# Cycle separator: Har-Peled–Nayyeri

```julia
cycle_separator(g) -> Vector{Int}
cycle_sides(g, cycle) -> (inside, outside)
```

A balanced **cycle** separator of an embedded planar graph in `O(n)` time. On a
triangulation it follows S. Har-Peled and A. Nayyeri, *A Simple Algorithm for
Computing a Cycle Separator*, Computational Geometry: Theory and Applications 83
(2019), 101–110 ([arXiv:1709.08122](https://arxiv.org/abs/1709.08122)); on any
other connected graph it lifts a cycle of a triangulation to a closed walk.

## Guarantees

`C = cycle_separator(g)` is a vector of vertices of the connected graph `g` such
that

1. `C` is a **closed walk of `g`**: at least three vertices, and consecutive
   entries — including the last and the first — are adjacent in `g`. It may repeat
   vertices and edges;
2. the vertices off the walk split into the two sides of the walk:
   `cycle_sides(g, C)` returns `(inside, outside)` with
   `sort(vcat(inside, outside, unique(C))) == collect(1:nv(g))`, the sets
   disjoint, and no edge between `inside` and `outside`;
3. the two sides are **balanced**: `max(|inside|, |outside|) <= 2/3 * nv(g)`. This
   always holds; it never depends on the heuristics below;
4. if every face of `g` is a triangle, `C` is a **simple cycle** with `O(sqrt(n))`
   vertices. For other graphs its length is not bounded.

The result is deterministic, and the running time is `O(nv)` in the worst case.

## Orientation

Rotations are counter-clockwise, so a walk read in the order returned has its
**left** side as the inside. `cycle_sides` reports the vertices to the left of the
traversal as `inside` and those to the right as `outside`; reversing the vector
swaps them. `cycle_separator` orients its result so that `inside` is the smaller
side.

For a walk that repeats vertices and edges, left and right come from winding
numbers: crossing an edge from its right to its left changes a face's number by
the times the walk uses the edge forward minus backward, so an edge used both ways
(a spike of the walk) can be crossed freely. The faces get two adjacent values,
the larger on the left; a walk with three values crosses itself and is rejected
with an `ArgumentError`, as is a vertex sequence that is not a closed walk of `g`.

```julia
g = planar_grid(20, 20)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)

isempty(intersect(inside, outside))                  # true
maximum([length(inside), length(outside)]) <= 2/3 * nv(g)
cycle_sides(g, reverse(C)) == (outside, inside)      # true
```

## Preconditions

`g` must be connected and have at least three vertices, and its rotation system
must be a planar embedding; otherwise an `ArgumentError`. Bridges, cut vertices
and large faces are fine.

## The algorithm

**Triangulated input.** All faces are triangles, and the three phases are the
paper's:

1. **A fundamental cycle.** Breadth-first search from a root gives a tree; the
   dual of the non-tree edges is a spanning tree of the faces (no union-find is
   needed, since the cotree of a spanning tree is a tree); an edge whose removal
   splits it most evenly gives a balanced fundamental cycle.
2. **Level cycles.** With the two tree paths to the endpoints of that edge forced,
   every face has a level, and the boundaries between consecutive levels are
   disjoint nested cycles.
3. **A ladder.** Every `Δ`-th level cycle, `Δ = ceil(sqrt(n/2))`, gives a small
   ladder; the first rung with `1/3..2/3` of the faces on both sides is returned,
   else the pieces of the heavy ring.

The result of phases 2 and 3 is **checked**: if it is not a walk with both sides
at most `2n/3`, the exact balanced walk below replaces it.

**Any other connected input.** The exact balanced walk comes from the virtual
star triangulation of the [separator](separator.md): every face with four or more
sides gets a star vertex joined to all its corners (one edge per corner, so a
vertex that appears several times on a face boundary gets several parallel star
edges, and graphs with bridges and cut vertices need no special case). The BFS
tree of `g` rooted at a pseudo-center is the spanning tree, with the stars as
leaves. The same walk down the dual tree that finds the Lipton–Tarjan cycle finds
a dual edge with at most `2/3` of the vertices strictly on each side, in time
linear in the number of triangles, and its fundamental cycle is a tree path plus
either a real edge or a star edge.

**Lifting.** A star edge of face `f` closes the tree path through `f`, and it is
replaced by the shorter boundary arc of `f` between the same two corners. This
only adds vertices to the walk: every vertex off the walk keeps the side it had,
so both sides stay at most `2/3`. The walk repeats a vertex when the tree path
touches the face boundary, and an edge when it runs along it.

## Complexity

Every step is a linear pass over the vertices, darts and faces (or triangles of
the virtual triangulation): the BFS, the dual tree, the balanced-edge walk
(at most four visits per triangle, see [Separator](separator.md)), the lifting
and the final check of the sides. `triangulate` is not called, so the
cost of a graph with a huge face is no longer dominated by it. Indicative timings
on one development machine: `0.02 s` for a `100 x 100` grid, `0.3 s` for a
`400 x 400` grid, and `0.05 s` for a ladder with `80 000` vertices.

## Usage notes

- The input must be a `PlanarGraph`; the rotation system is used for the
  embedding. A bare `SimpleGraph` has no embedding to work with (compute one with
  `planar_embedding`).
- For a vertex separator rather than a cycle, use [`separator`](separator.md): it
  has a size bound for every input, which a cycle cannot have.
- `cycle_sides` is a separate query: it works on any connected `PlanarGraph` and
  any closed walk you already have, not only on the result of `cycle_separator`.
