# Cycle separator: Har-Peled–Nayyeri

```julia
cycle_separator(g) -> Vector{Int}
cycle_sides(g, cycle) -> (inside, outside)
```

A balanced **cycle** separator of an embedded planar graph, in time linear in the
size of the internal triangulation, following S. Har-Peled and A. Nayyeri, *A
Simple Algorithm for Computing a Cycle Separator*, Computational Geometry:
Theory and Applications 83 (2019), 101–110
([arXiv:1709.08122](https://arxiv.org/abs/1709.08122)).

## Guarantees

`C = cycle_separator(g)` is a vector of vertices of `g` such that

1. `C` is a **simple cycle of `g`**: at least three vertices, no repeats, and
   consecutive entries — including the last and the first — are adjacent in `g`;
2. the vertices outside `C` split into the two sides of the cycle:
   `cycle_sides(g, C)` returns `(inside, outside)` with
   `sort(vcat(inside, outside, C)) == collect(1:nv(g))`, the three sets pairwise
   disjoint, and no edge between `inside` and `outside`;
3. the two sides are **balanced**: `max(|inside|, |outside|) <= 2/3 * nv(g)`;
4. `C` is **short in the triangulation used internally** (`O(sqrt(n))`), and on
   graphs whose faces are already triangles that bound applies to `C` itself.

The result is deterministic for a given graph, and repeated calls return the
same vector.

## Orientation

Rotations in this package are counter-clockwise, so a cycle read in the order
returned has its **left** side as the inside: `cycle_sides` reports the vertices
to the left of the traversal as `inside` and those to the right as `outside`.
Reversing the vector swaps the two sides. `cycle_separator` orients its result so
that `inside` is the smaller side, which fixes an otherwise arbitrary choice.

`cycle_sides` raises an `ArgumentError` when `cycle` has fewer than three
entries, repeats a vertex, contains a value outside `1:nv(g)`, or has two
consecutive entries that are not adjacent (including the wrap-around pair). It
also requires a connected graph.

```julia
g = planar_grid(20, 20)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)

sort(vcat(inside, outside, C)) == collect(1:nv(g))   # true
isempty(intersect(inside, outside))                  # true
maximum([length(inside), length(outside)]) <= 2/3 * nv(g)
cycle_sides(g, reverse(C)) == (outside, inside)      # true
```

## Preconditions

`g` must be connected, 2-connected, and have at least three vertices:
2-connectivity is what makes every face a simple cycle, which the internal
triangulation needs. A connected plane graph is 2-connected exactly when every
face boundary is a simple cycle, so the check is one pass over the darts; a graph
with a bridge or a cut vertex is rejected with an `ArgumentError` rather than
silently producing a non-cycle. A disconnected or too small graph is rejected
the same way.

## The algorithm

The three phases are the paper's, run on a triangulation `h` of `g` built
internally with [`triangulate`](api.md#triangulate-g-typeof-g):

1. **A fundamental cycle.** Breadth-first search from a root gives a tree. The
   dual of the non-tree edges is a spanning tree of the faces because the
   embedding is planar; an edge whose removal splits that tree most evenly gives
   a balanced fundamental cycle `S`. `S` is returned directly when it is already
   short.
2. **Level cycles.** With the two tree paths to the endpoints of the chosen edge
   forced, every face has a level (the largest level of its vertices). The
   boundary between faces of level `i` and `i+1` is a cycle `C_i`; these cycles
   are disjoint and nested.
3. **A ladder.** Taking every `Δ`-th level cycle, `Δ = ceil(sqrt(n/2))`, gives a
   small ladder. The number of faces inside each rung is tracked by one dual BFS
   with a potential; the first rung that leaves `1/3..2/3` of the faces on both
   sides is returned. If no rung qualifies, the two pieces of the *heavy* ring —
   the one that straddles the middle — are tried, then the inner region plus the
   first piece. The first simple cycle wins, and the fundamental cycle is the
   fallback.

The ladder step is an internal knob (`_cycle_separator(dg; delta = ...)`) used by
the tests to exercise that branch; the public interface always uses
`ceil(sqrt(n/2))`.

## Shortcutting the chords

The paper works on a triangulated graph, and so does OGDF's `SeparatorHarPeled`,
which returns only the *vertex set* of the cycle. The vector returned here,
however, must be a cycle of the caller's graph, which may have arbitrarily large
faces. The implementation therefore records, for every face of `h`, the face of
`g` that contains it, and **shortcuts** the chords of the triangulation cycle:

- a maximal run of fill edges is replaced by the boundary arc of the face of `g`
  that contains it, chosen away from the separator's interior;
- the global version marks each face of `g` as inside, outside, or crossed and
  traces the boundary of the inside faces (and of the complement), keeping the
  valid simple cycle with the best balance;
- if the face set is degenerate — every face of `g` crossed, as happens for a
  cycle graph — a second pass replaces one run at a time by the boundary arc of
  its face that contains no other vertex of the cycle. This is the path taken for
  `planar_cycle(n)`, whose two faces are both crossed.

Every candidate is checked to be a simple cycle of `g` before it is returned; if
the shortcut fails for both the ladder cycle and `S`, an error is raised rather
than returning something that is not a cycle.

This shortcut is not in the paper and not in OGDF; it is what makes the returned
vector a cycle of the input graph. The construction of the regions, the level
cycles and the ladder follows the paper and the OGDF implementation.

## Complexity

Every separator phase and the shortcut are linear passes over the vertices,
darts and faces of the internal triangulation, so the running time is
`O(nv + ne)` with the constant of a few dozen passes for an input that is
already triangulated. The preprocessing is [`triangulate`](api.md#triangulate-g-typeof-g);
its cost is linear on graphs whose faces are small (triangulations, grids) and is
worst for a single huge face, where the triangulation itself is the bottleneck.
Indicative timings on one development machine: `0.02 s` for a `100 x 100` grid,
`0.4 s` for a `300 x 300` grid, and `0.3 s` for a random triangulation with
`50 000` vertices.

## Usage notes

- The input must be a `PlanarGraph`; the rotation system is used for the
  embedding. A bare `SimpleGraph` has no embedding to work with (compute one with
  `planar_embedding`).
- For a vertex separator rather than a cycle, use [`separator`](separator.md); the cycle
  is often much smaller than the vertex set it separates, but the two guarantees
  differ.
- `cycle_sides` is a separate query: it works on any connected `PlanarGraph` and
  any simple cycle you already have, not only on the result of
  `cycle_separator`.
