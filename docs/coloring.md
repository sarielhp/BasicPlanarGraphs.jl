# Coloring

```julia
five_coloring(g) -> Vector{Int}
is_proper_coloring(g, colors) -> Bool
```

`five_coloring` returns `c` with `c[u] != c[v]` for every edge `uv` and
`1 <= c[v] <= 5`. The input is a `SimpleGraph` or a `PlanarGraph` (only the edges
are used); it may be disconnected and have isolated vertices.

## Method

Matula, Shiloach and Tarjan (1980). A planar graph has a vertex of degree at most
five. Remove such a vertex `v` and remember its neighbours. If `v` has five
neighbours, two of them are not adjacent (five pairwise adjacent neighbours and `v`
would make a `K5`); identify them, which keeps the graph planar. Undo the steps in
reverse order: an identified pair gets one color, and `v` then sees at most four
colors, so a fifth is free.

## Guarantees

- A returned coloring is always proper and uses at most five colors, also for a
  non-planar input that the method happens to finish.
- `ArgumentError` when no vertex of degree at most five remains or five neighbours
  of a removed vertex are pairwise adjacent. This proves the graph is not planar.
  The converse does not hold: for example `K3,3` is colored. Use
  [`is_planar_linear`](api.md) when planarity must be decided.
- Directed graphs and other `AbstractGraph` types give `ArgumentError`.

## Cost

Adjacency sets are hash sets, so an adjacency test is `O(1)`, and an identification
moves the smaller set into the larger. On random triangulations the elimination
time roughly doubles when `n` doubles (0.03 s at 50,000 vertices, 0.84 s at
800,000). Linear time is **not** proved; the paper reaches it with linked
adjacency lists and a bound on the degrees of the identified pair. The Left-Right
planarity test, if you run it, costs about ten times the coloring.
