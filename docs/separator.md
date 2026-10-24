# Separator: Lipton–Tarjan

```julia
separator(g; balance = 2/3) -> (A, B, C)
```

Computes a balanced vertex separator of a planar graph in linear time, following
R. J. Lipton and R. E. Tarjan, *A Separator Theorem for Planar Graphs*, SIAM J.
Appl. Math. 36(2), 1979.

## Guarantees

For `n = nv(g)` the three returned vectors:

1. **partition** the vertices: `sort(vcat(A, B, C)) == collect(1:n)`, pairwise
   disjoint, each sorted ascending;
2. **separate** `A` from `B`: no edge has one end in `A` and the other in `B`
   (`C` may have edges to either);
3. are **balanced**: `max(|A|, |B|) <= balance * n`;
4. have a **small separator**: `|C| = O(sqrt(n))`. On the graphs in the test
   suite, `|C| <= 2*sqrt(2n) + 1`.

`balance` must lie in `(0, 1)`; anything else throws an `ArgumentError`.

The result is deterministic: the same graph gives the same partition on every
call. The algorithm uses only adjacency, so the rotation system affects BFS
visit order and therefore tie-breaking, not correctness.

Two structural notes:

- Validity is **structural, not asserted**. Every returned partition comes from
  `_split_and_pack`, which takes the connected components of `G \ C` and packs
  them into two bins of nearly equal size. Components are separated by
  construction, so there can be no `A`–`B` edge whatever `C` contains.
- The work goes into *balance* and *size*: the algorithm chooses `C`, then keeps
  growing it until the packing fits.

## The algorithm

### 1. Level structure

Breadth-first search from a root gives levels `L(0), L(1), ...`. The separator
combines two levels with a middle band:

- `l1` is the first level whose prefix reaches `n/2` vertices;
- `l0 <= l1` is the largest level with `|L(l0)| + 2(l1 - l0) <= 2*sqrt(k)`,
  where `k` is the number of vertices on levels `0..l1`;
- `l2 >= l1 + 1` is the smallest level with
  `|L(l2)| + 2(l2 - l1 - 1) <= 2*sqrt(n - k)`.

This is Theorem 4's choice: the two levels are small enough, and everything
below `l0` or above `l2` is already a modest fraction of the graph.

### 2. The base separator

`C = L(l0) ∪ L(l2)`; the band is `L(l0+1) .. L(l2-1)`. If the band alone has at
most `balance * n` vertices, the base separator is already enough and the packing
decides the split.

### 3. Cutting the band

Otherwise the band is too big, and a cycle through the band is added to the
separator (Lemmas 2 and 3):

- the band is **triangulated** in place, keeping the rotation system;
- the darts that belong to the BFS tree of the whole search are marked; the
  remaining darts form a spanning tree of the dual;
- each face is **weighted** by the band vertices incident to it, so the weight of
  the faces inside a cycle bounds the band vertices strictly inside it;
- a post-order pass computes the subtree weight of every dual edge, and the dual
  edges whose removal splits the dual most evenly are ranked;
- the fundamental cycle of a non-tree edge splits the band into two parts; its
  primal path is a candidate separator. When the path runs through the shrunk
  levels (already part of `C`), only the two tree paths are kept.

At most six candidates are examined (`MAX_CYCLE_CANDIDATES`): the first one that
makes both sides fit wins, otherwise the best one by largest side is kept.

A band that cannot be triangulated (it is not 2-connected — a face that is not a
simple cycle) contributes no cycle candidates; the level separator has to do on
its own, and refinement takes over.

### 4. Refinement

After packing, if a side is still larger than `balance * n` the separator is
grown: the largest remaining component is separated recursively by the same
procedure and the two separators are merged. Components shrink by a constant
factor per pass, so the added vertices sum to `O(sqrt(n))`. The loop is bounded
by `MAX_REFINEMENTS = 32`, and it stops early when a pass adds nothing.

### 5. Choosing the root

The quality of the level structure depends strongly on the BFS root. Four roots
are tried — the original root, a farthest vertex, the midpoint of a BFS path (the
centroid on tree-like graphs) and the farthest vertex from that midpoint — and
the best result is kept. Candidates are ranked by

1. balanced before unbalanced,
2. smaller `|C|`,
3. smaller `max(|A|, |B|)`.

That is a constant four passes, not a factor of `n`.

### 6. Disconnected inputs

Disconnected graphs follow the end of the proof of Theorem 4. With components
sorted by size and `n` the total:

| Largest component | Result |
|---|---|
| larger than `balance * n` | separate that component recursively; its two sides join the corresponding sides and the rest of the graph goes to the larger one |
| at least `n/3` | that component is one side, the other components the other side, `C = ∅` |
| smaller | components are packed greedily into two sides, in decreasing size, `C = ∅` |

Degenerate sizes are handled directly: `n == 0` returns three empty vectors, and
`n <= 2` puts every vertex into `C`.

## Complexity

Each stage is one pass over vertices, darts or faces:

- BFS and level bookkeeping: `O(n + m)` per root, four roots;
- band triangulation and dual spanning tree: `O(n + m)` for the band;
- packing and refinement: `O(n + m)` per pass, at most 32 passes, and passes stop
  as soon as the partition is balanced.

So the running time is `O(n + m)` with a small constant, not a function of the
number of passes in practice. Indicative timings on one development machine:
about 0.2 s for a `320 x 320` grid (`n = 102400`) and for a random triangulation
with `n = 32000`.

## Usage notes

- The input must be a `PlanarGraph`, not a bare `SimpleGraph`. The embedding is
  used for adjacency only, so a graph built from a planar layout works even if
  you never needed the rotation system otherwise.
- A looser `balance` (for example `0.75`) generally means a smaller `C`; a
  tighter one may force refinement and a larger separator, and extremely small
  values can make the guarantee impossible to meet.
- Vertices in `A`, `B` and `C` are sorted; use `findall`/indexing rather than
  assuming anything about the order beyond that.
- If you need components rather than a partition, use the graph interface:
  `connected_components(SimpleGraph(g))` from `Graphs.jl`.

```julia
g = planar_grid(20, 20)
A, B, C = separator(g)

isempty(intersect(A, B))                     # partition and separation
maximum([length(A), length(B)]) <= 2/3 * nv(g)
length(C) <= 2 * sqrt(2 * nv(g)) + 1
```
