# Separator partition

```julia
separator_partition(g, r; subgraphs = false) -> SeparatorPartition
```

Recursive separator partition of an embedded planar graph: one global separator
`S` and parts `P₁, …, Pₖ` of at most `r` vertices each, with no edge between two
different parts. The specification is
`~/papers/26/ai/triangulations/lib/planar_partition/specs.md`.

## Result

A `SeparatorPartition` with `separator`, `parts` and `part_of` (original vertex
numbers; `part_of[v] == 0` exactly on the separator). With `subgraphs = true` it
also holds `local_index` and `graphs`: the induced, inherited-embedding graph of
each part, vertex `j` of `graphs[i]` being `parts[i][j]`.

## Guarantees

1. `S, P₁, …, Pₖ` partition the vertices and every part is nonempty;
2. every part has at most `r` vertices;
3. no edge joins two different parts;
4. `|S| <= 16 n / sqrt(r)`;
5. `k <= 2n/r + 1`;
6. `part_of` is the direct lookup above;
7. the result is deterministic, and `g` is not modified.

The graph may be disconnected, with isolated vertices. `r >= nv(g)` gives an
empty separator and one part; the empty graph gives empty output; `r < 1` is an
`ArgumentError`, as is a rotation system that is not a planar embedding. A part
may contain several components of `g - S`.

## Weights

`separator_partition(g, r; weights = w)` takes one nonnegative integer per vertex, as
[`separator`](separator.md) does, and bounds a part by its total weight:

- every part weighs at most `r`, and `k <= 2W/r + 1` for the total weight `W`;
- a vertex heavier than `r` fits in no part and is in the separator;
- guarantees 1, 3 and 6 are unchanged, and every group is cut by the weighted separator
  until it weighs at most `r`;
- the separator bound of guarantee 4 is **not** claimed: the weight can sit on few
  vertices, and the separator bounds the number of vertices it removes, not their weight.

Without `weights` the result is identical to before (checked on fixed cases).

## Algorithm and cost

A group of `q > r` vertices is split with [`separator`](separator.md): its
separator goes into `S` and its two sides, each at most `2q/3` (the primitive's
guarantee, for disconnected groups too), become the child groups. Groups of at
most `r` vertices are terminal. Terminal groups are packed next-fit into parts of
at most `r` vertices.

- **Separation.** Distinct sides have no connecting edge, so distinct terminal
  groups have none, and packing keeps that.
- **Separator size.** A group of size `q` adds at most `2sqrt(2q)` vertices,
  charged to its `q` vertices. Along any vertex's path the sizes shrink by `2/3`
  per level and all exceed `r`, so its total charge is at most
  `2sqrt(2)/sqrt(r) * 1/(1 - sqrt(2/3)) < 16/sqrt(r)`.
- **Part count.** Two consecutive packed parts hold more than `r` vertices
  together, because the first terminal group of the later part did not fit in the
  earlier one; pairing them gives `k <= 2n/r + 1`.
- **Time.** There are `O(1 + log(n/r))` levels; groups of a level are disjoint, so
  their induced graphs, separator calls and packing total `O(n)` per level. The
  conversion from `PlanarGraph`, the induced graphs (`induced` matches twin darts
  with counting sorts) and the final packing are `O(n)`. Total
  `O(n (1 + log(n/r)))`.
- **Space.** A processed group is dropped once its children exist, and live groups
  are disjoint, so the live size is `O(n)`; the output is `O(n)`.

```julia
g = planar_grid(30, 30)
p = separator_partition(g, 50)
all(length(P) <= 50 for P in p.parts)
```
