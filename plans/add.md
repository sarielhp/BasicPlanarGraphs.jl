# Possible additions

Status: five-coloring (v0.2.6), weighted `separator` (v0.2.7; not
`separator_partition`), a true r-division, embedding-preserving edits and
`PlanarGraph <: AbstractGraph` are done; the rest are proposals. Written 2026-10-24 after the v0.2.4 release,
while the independent review of the separator, cycle separator and triangulation
code was still open; revisit after its findings are fixed.

The maintainer's pick so far: **5-coloring** (below). The rest are candidates.

## 1. Five-coloring (maintainer's request) — done in v0.2.6

Deviations from the plan below: no planarity pre-check (it cost ten times the
coloring); `ArgumentError` only when the method gets stuck; hash-set adjacency with
small-to-large merging, near linear in practice but not proved linear; no Kempe-chain
oracle. The paper's linked-list scheme is the way to a proof.

`five_coloring(g) -> Vector{Int}`: a proper vertex coloring of a planar graph with
colors `1:5`, in `O(n)` time. (Six colors is the trivial degeneracy bound; four is
the Four Color Theorem and has no simple algorithm, so five is the right target.)

Input: a `SimpleGraph` or a `PlanarGraph`. The algorithm does not need the
embedding, only planarity; a `PlanarGraph` input can use its rotation to choose the
pair to merge, a `SimpleGraph` is checked with `is_planar_linear` (or the algorithm
fails when it meets a minimum degree above 5, which proves non-planarity, and
throws `ArgumentError`).

### Algorithm (Matula–Shiloach–Tarjan 1980; Chiba–Nishizeki–Saito 1981)

1. Every planar graph has a vertex of degree at most 5. Keep buckets of vertices by
   degree (degrees `0..5`) and repeatedly remove a vertex `v` of degree at most 5,
   pushing it on a stack together with its neighbours at that moment.
2. Degree at most 4: nothing special; `v` will see at most 4 colors.
3. Degree exactly 5: its neighbours `a_1..a_5` cannot be pairwise adjacent (that
   would be a `K5` with `v`), so two of them, `a`, `b`, are non-adjacent. Merge
   `a` and `b` into one vertex (edge contraction through `v`, or identification);
   the graph stays planar, so the process continues. Record the merge.
4. Unwind the stack: color the merged vertex, give both `a` and `b` its color, and
   `v` then sees at most 4 colors among its neighbours (`a`, `b` count once), so a
   fifth is free.
5. Linear time needs care: after merges a vertex can reach large degree, so keep
   merged adjacency lists as linked lists with a merge-the-shorter rule and choose
   the pair to merge among neighbours of small degree (the published argument shows
   a pair with degrees at most 11 exists), so each merge costs `O(1)` amortized.
   The simpler Kempe-chain version is `O(n^2)` and is fine as a test oracle.

### Work items

- `src/coloring.jl`, `export five_coloring`, docstring, `docs/coloring.md`, row in
  `docs/api.md` and the README table, `CHANGELOG.md` entry.
- Tests (`test/coloring_test.jl`): proper coloring checked independently (every edge
  has different colors, colors in `1:5`) on `K4`, wheels (odd and even rim),
  grids, trees, random triangulations, thinned triangulations, disconnected graphs
  with isolated vertices, `n = 0, 1, 2`; non-planar input (`K5`, `K3,3`, Petersen
  subdivisions) gives `ArgumentError`; a Kempe-chain oracle on small graphs; a
  scaling test (many components and a long path).
- Decide whether to also expose `six_coloring` (degeneracy order, 10 lines) and
  `is_proper_coloring(g, colors)` as a public validator. Recommended: yes to the
  validator, it makes the tests and users' checks trivial.

## 2. Weighted separators — `separator` done, `separator_partition` not

`separator(g; weights)` is implemented (integer weights, see `docs/separator.md`).
A weighted `separator_partition` is left out: a vertex heavier than `r` fits no part,
so its size bounds would need a different statement.

Optional vertex weights for `separator` and `separator_partition`, so that
"balanced" means by weight. Lipton–Tarjan handles weights directly; the exact
balanced-walk search and the bin packing count vertices in a few places only.

## 3. A true r-division

`separator_partition` bounds the total separator and the part sizes; an r-division
also bounds each part's boundary by `O(sqrt r)`. Read the spec
(`/home/sariel/papers/26/ai/triangulations/lib/planar_partition/specs.md`) first:
add only if it needs the per-part boundary bound.

## 4. Embedding-preserving edits

Public `delete_edge`, `contract_edge`, `add_edge_in_face` on `PlanarGraph`, on top
of the dart core, for incremental algorithms.

## 5. `PlanarGraph <: AbstractGraph`

Today only seven `Graphs.jl` functions are forwarded; `connected_components` and
`bfs_parents` throw `MethodError`. Subtyping needs `is_directed`, `edgetype`,
`inneighbors`, `outneighbors` and friends, and then every `Graphs.jl` algorithm works.

## 6. Smaller items

- SVG writer for `grid_drawing`.
- `verify_embedding(g)` that returns the reason an embedding is invalid, not a `Bool`.
- A uniformly random planar graph generator (`random_triangulation` is a stacked
  triangulation, not uniform).
- Biconnected components.
