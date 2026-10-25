# r-division

```julia
r_division(g, r) -> RDivision
```

A division of an embedded planar graph into regions of at most `r` vertices that
overlap only in a small *boundary*. It is the stronger relative of
[`separator_partition`](partition.md): there the parts share nothing and the
separator is one global set; here every region also has a bounded boundary.

## Result

A `RDivision` with

- `regions[i]`: the vertices of region `i`, boundary vertices included;
- `boundary[i]`: the vertices of region `i` that lie in another region too;
- `region_of[v]`: the region of `v` if it lies in only one, `0` for a boundary vertex;
- `boundary_vertices`: all vertices in more than one region.

## Guarantees

For every input and every `r >= 100`:

1. every vertex is in some region and every edge has both ends in one region;
2. a vertex is in more than one region exactly when it is a boundary vertex;
3. every region has at most `r` vertices;
4. every region has at most `12 sqrt(r)` boundary vertices;
5. the result is deterministic and the input is not modified.

The number of regions `k` and the total boundary are expected to be `O(n/r)` and
`O(n/sqrt(r))`. That is argued below and measured, not proved:

| graph | n | r | regions (n/r) | boundary vertices (n/sqrt r) | largest boundary (limit) |
|---|---|---|---|---|---|
| grid | 90,000 | 200 | 709 (450) | 12,848 (6,364) | 66 (170) |
| grid | 90,000 | 1000 | 128 (90) | 5,348 (2,846) | 122 (379) |
| random triangulation | 20,000 | 200 | 149 (100) | 204 (1,414) | 11 (170) |
| triangulated grid | 14,400 | 300 | 67 (48) | 1,959 (831) | 88 (208) |
| wheel | 20,001 | 200 | 139 (100) | 140 (1,414) | 3 (170) |

## Method

Two phases of recursive separators, both through the machinery of
[`separator`](separator.md), with the weights added in 0.2.7.

1. A piece with more than `r` vertices is cut by a balanced separator. Each child
   holds its side and the separator vertices that have a neighbour there; a
   separator vertex that lands in both children is *shared*. Children have at most
   `2q/3 + 2 sqrt(2q) < q` vertices for `q > 72`, so this stops.
2. A piece with at most `r` vertices but more than `12 sqrt(r)` shared vertices
   is cut by a separator with weight 1 on the shared vertices. Both children keep at
   most `2/3` of the shared vertices; the new separator adds at most `2 sqrt(2r)`.
   The count `b` becomes at most `2b/3 + 2 sqrt(2r) < b` whenever `b > 6 sqrt(2) sqrt(r)`,
   which is below the threshold `12 sqrt(r)`, so every piece ends below it.

The final pieces are packed next-fit into regions: a piece joins the current region
while the size stays at most `r` and the sum of the shared counts stays at most
`12 sqrt(r)`. The sum bounds the real boundary of the union, so property 4 holds.

## Why the counts are expected to be small

Along any vertex's path in the recursion the sizes shrink by a constant factor and
stay above `r` in phase 1, so a separator of size `O(sqrt q)` charged to the `q` vertices
of its piece costs `O(1/sqrt r)` per vertex in all: `O(n/sqrt(r))` boundary vertices
are created. Nested phase-1 pieces differ by a constant factor, so pieces of similar
size are disjoint and the pieces of size in `[x, 1.2x)` number at most `n'/x`: the
phase-1 tree has `O(n/r)` pieces. The step that is not proved is the number of phase-2
cuts, and the bound on the number of regions after packing rests on it (two
consecutive regions together exceed `r` vertices or `12 sqrt(r)` shared vertices).

`r >= 100` is required because a separator of `2 sqrt(2q)` vertices must be smaller than
a third of `q` for a piece to shrink; an `ArgumentError` is thrown otherwise.

## Cost

Each phase-1 level is linear in the total size of its pieces, and there are
`O(log(n/r))` levels, so `O(n (1 + log(n/r)))`; a phase-2 cut costs time linear in a piece
of at most `r` vertices. The measured time is 0.9 s for the 90,000-vertex grid with
`r = 200`. No region is built as a `PlanarGraph`: the induced graphs of the regions can be
derived from `regions` and the input.
