# Separator: Lipton–Tarjan

```julia
separator(g; balance = 2/3, weights = nothing) -> (A, B, C)
```

Computes a balanced vertex separator of a planar graph in linear time, following
R. J. Lipton and R. E. Tarjan, *A Separator Theorem for Planar Graphs*, SIAM J.
Appl. Math. 36(2), 1979. The implementation follows the proof of the theorem, so
the guarantees below are worst case, for every input, and need no retries.

## Guarantees

For `n = nv(g)` the three returned vectors:

1. **partition** the vertices: `sort(vcat(A, B, C)) == collect(1:n)`, pairwise
   disjoint, each sorted ascending;
2. **separate** `A` from `B`: no edge has one end in `A` and the other in `B`
   (`C` may have edges to either);
3. are **balanced**: `max(|A|, |B|) <= balance * n`, always;
4. have a **small separator**: `|C| <= 2*sqrt(2n)`, always.

`balance` must lie in `[2/3, 1)`; anything else throws an `ArgumentError`. `2/3`
is what the theorem proves, so a smaller value cannot be promised. A larger
value does not change the algorithm, it only lets the final pruning step move
more separator vertices into `A` and `B`.

## Weights

`weights` is a vector with one nonnegative integer per vertex. Balance is then by
weight: with total weight `W`, the weights of `A` and `B` are at most `balance * W`,
and so is the weight of every component of `g - C`. The bound `|C| <= 2*sqrt(2n)`
still counts vertices, not weight. Without `weights` every vertex weighs 1 and the
result is identical to before.

Lipton and Tarjan's proof carries over with three changes: the median level `l1` is
where the prefix sum of the weights reaches `W/2` (the sizes `|L(l)|` and `k` that
choose `l0` and `l2` still count vertices); the band is balanced by weight, and so
is the search for the dual edge (the proof that a balanced edge exists only
compares sums, so it is unchanged); and the packing and pruning steps add weights.
A vertex heavier than `W/3`, or than a third of the weight of the component being
cut, is returned alone as the separator: the rest of that component then weighs
less than two thirds of it.

Weights must be integers so that every comparison is exact and the incremental
sums cannot drift; scale real weights first. The total may not exceed
`typemax(Int) ÷ 3`. A wrong length, a negative or non-integer entry or a total
above the limit is an `ArgumentError`. `separator_partition` takes the same
`weights`; see [Separator partition](partition.md#weights).

The running time is `O(n)` in the worst case: no hashing, no comparison sorts,
no data-dependent restarts (every loop is a pass over vertices, darts, triangles or
levels; see [Complexity](#complexity)).
The result is deterministic.

The input must be a planar embedding. A rotation system of higher genus is
detected (Euler's formula fails, or the dual of the non-tree edges is not a
tree) and reported as an `ArgumentError`.

## The algorithm

### 1. Level structure

Breadth-first search from a root gives levels `L(0), L(1), ...`:

- `l1` is the first level whose prefix reaches `n/2` vertices, and `k` is the
  number of vertices on levels `0..l1`;
- `l0 <= l1` is the largest level with `|L(l0)| + 2(l1 - l0) <= 2*sqrt(k)`;
- `l2 >= l1 + 1` is the smallest level with
  `|L(l2)| + 2(l2 - l1 - 1) <= 2*sqrt(n - k)`.

Both levels exist (this is the counting argument of the theorem), and the tests
are done in integer arithmetic. The vertices below `l0` are fewer than `n/2`, and
so are those above `l2`; the **band** is `L(l0+1) .. L(l2-1)`.

### 2. The base separator

`C = L(l0) ∪ L(l2)`. If the band has at most `2n/3` vertices, every component of
`G - C` has at most `2n/3` vertices and the separator is complete. Otherwise a
cycle through the band is added.

### 3. The balanced cycle (Lemma 2)

Lemma 2 needs a triangulation. Instead of adding edges, every face with four or
more sides gets a virtual **star** vertex joined to all its corners. This works
for any connected plane graph, including bridges and cut vertices (a repeated
corner just gives parallel star edges). The spanning tree is the BFS tree of the
input, with each star hanging off it as a leaf, so fundamental cycles keep the
BFS depth bound: the band vertices on one are at most `2(l2 - l0 - 1)`, and
`|C| <= [|L(l0)| + 2(l1 - l0)] + [|L(l2)| + 2(l2 - l1 - 1)] <= 2sqrt(k) + 2sqrt(n-k)
<= 2*sqrt(2n)`.

The non-tree edges form a spanning tree of the dual, whose nodes are the
triangles (at most three neighbours each). A dual edge splits the band vertices
into those strictly on one side (all their triangles there), those strictly on
the other, and those on the cycle. An edge is *balanced* when each strict side
has at most `2/3` of the band vertices. The existence argument: orient every
unbalanced edge toward its heavy side; walking those orientations ends at a
triangle with at most three branches, and a vertex is strictly inside at most
`k - 1` of its `k` sides (every band vertex has at least two triangles, so it
has one in some branch), so the `k` sides cannot all be heavy: heavy means
strictly more than `2/3` of the band, and `k - 1 <= 2k/3` for `k <= 3`. Hence a balanced edge exists, and the walk finds it.

The walk starts at a root triangle and moves to the child whose side is heavy.
Per-vertex counters of how many triangles lie in the current region give the
strict counts exactly: the largest child is measured by *removing* the current
triangle and the small subtrees from the region, the small children from
scratch. A step costs the size of its small subtrees, and the region shrinks by
them, so the total over the walk is `O(#triangles)`.

The cycle's band vertices (the BFS tree path between the endpoints of the
balanced edge, restricted to the band) are added to `C`. The Jordan curve of the
cycle separates the strictly-inside band vertices from the strictly-outside
ones, and the band touches the rest of the graph only through `L(l0)` and
`L(l2)`, both in `C`.

### 4. Packing

Every component of `G - C` has at most `2n/3` vertices, and together they have
at most `n`. A component with at least `n/3` vertices goes alone into `A` and the
rest into `B`; otherwise components fill `A` until it passes `n/3` (so it stays
below `2n/3`) and the rest go to `B`. Both bins are at most `2n/3`.

### 5. Pruning

A vertex of `C` with no neighbour in `B` can join `A` (and symmetrically), as
long as the side stays within `balance * n`. One pass over `C` reaches a fixed
point, because moving a vertex into a side can only block other vertices from the
opposite side, never free them. This removes redundant vertices: a path gives
one separator vertex, a cycle two, and a `k x k` grid `k - 1`.

### 6. Choosing the root

The size of the level structure depends on the BFS root, so four roots are tried
(the original root, a farthest vertex, the midpoint of a BFS path, and the
farthest vertex from that midpoint) and the partition with the smallest
separator after pruning is kept. Every root satisfies the guarantees; the choice
is a constant factor in time, not in the bound.

### 7. Disconnected inputs

Components are labelled first. If the largest component has at most `2n/3`
vertices no separator is needed, and packing alone splits the graph with
`C = ∅`. Otherwise the separator of the largest component (as above, with `n`
replaced by its size) is the separator of the graph: its pieces have at most
`2/3` of that component, and everything else has fewer than `n/3`, so packing
works as in step 4. `n <= 2` is handled directly (pruning leaves one separator
vertex for a single edge).

## Complexity

Every stage is a pass over vertices, darts or triangles of the virtual
triangulation (fewer than `2m` of them):

- BFS and levels: `O(n + m)` per BFS, at most four roots to try and three more BFS
  runs to choose them;
- faces and triangles: `O(n + m)`; the dual tree has no union-find, the cotree
  of a spanning tree already is a tree;
- the balanced-edge walk: `O(#triangles)`, as above;
- packing and pruning: `O(n + m)`; conversion from `PlanarGraph` matches twin
  darts with two counting sorts, so it is `O(n + m)` in the worst case, not in
  expectation.

Indicative timings on one development machine: about 0.3 s for a `500 x 500`
grid and 1.4 s for `10^6` vertices; random triangulations with a random vertex
numbering are slower per vertex because of cache misses, not because of the
algorithm.

## Usage notes

- The input must be a `PlanarGraph`, not a bare `SimpleGraph`.
- Vertices in `A`, `B` and `C` are sorted.
- If you need components rather than a partition, use the graph interface:
  `connected_components(SimpleGraph(g))` from `Graphs.jl`.

```julia
g = planar_grid(20, 20)
A, B, C = separator(g)

isempty(intersect(A, B))                     # partition and separation
maximum([length(A), length(B)]) <= 2/3 * nv(g)
length(C) <= 2 * sqrt(2 * nv(g))
```
