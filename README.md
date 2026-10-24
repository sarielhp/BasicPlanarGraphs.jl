# BasicPlanarGraphs

[![CI](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Planar graphs with an abstract embedding (a rotation system), and a linear-time
planar separator based on Lipton and Tarjan, *A Separator Theorem for Planar
Graphs* (SIAM J. Appl. Math. 36(2), 1979).

## Installation

Requires Julia 1.10 or later. The package is not registered yet, so install it
from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/sarielhp/BasicPlanarGraphs.jl")
```

## Usage

```julia
using Graphs, BasicPlanarGraphs

g = random_triangulation(1000)
A, B, C = separator(g)

isempty(intersect(A, B))            # true: the three parts are disjoint
maximum([length(A), length(B)])     # at most 2/3 n
length(C)                           # O(sqrt(n))
```

## The type

`PlanarGraph{T}` is a `Graphs.jl` `SimpleGraph{T}` plus, for every vertex, the
counter-clockwise cyclic order of its neighbours:

```julia
triangle = PlanarGraph([[2, 3], [3, 1], [1, 2]])
nv(triangle), ne(triangle), nfaces(triangle)     # 3, 3, 2
faces(triangle)                                  # both walks are 1-2-3
dual(triangulate(planar_grid(3, 4)))             # 3-regular dual, Euler holds
```

The constructor rejects rotations that are not permutations of the neighbours,
so the stored data always describes a graph: every rotation entry is an edge and
every edge appears twice. Two embeddings are equal when their rotations agree up
to a cyclic shift and reversal at each vertex, that is when they describe the
same embedding.

`is_planar_embedding` checks Euler's formula `n_c - m_c + f_c = 2` on every
connected component; the constructor does not, which is what makes it possible
to test that a non-planar rotation system is rejected.

Additional operations: `neighbor_order`, `rotation`, `faces`, `nfaces`, `dual`,
`triangulate`, and the generators `planar_path`, `planar_cycle`, `planar_grid`,
`planar_wheel`, `random_triangulation`.

## The separator

`separator(g; balance = 2/3)` returns `(A, B, C)` with

- `A`, `B`, `C` a partition of the vertices,
- no edge between `A` and `B`,
- `|A|, |B| <= balance * n`,
- `|C| = O(sqrt(n))`; on the test families `|C| <= 2 sqrt(2 n)`.

Every returned partition is derived from the connected components of `G \ C`, so
validity is structural; the work goes into balance and size.

How it works, following the paper:

1. Breadth-first search on the input graph gives the levels `L(0), L(1), ...`.
   A middle level `l1` and two levels `l0 < l1 < l2` are chosen so that the
   vertex separators `L(l0)` and `L(l2)` are small (Theorem 4).
2. If the middle band (levels `l0+1 .. l2-1`) contains more than `balance * n`
   vertices, it is triangulated and its fundamental cycles are ranked by the
   number of vertices they enclose, using the spanning tree of the dual made of
   the non-tree edges (Lemma 2 and 3). The best cycles become extra separator
   candidates.
3. The best candidate is then refined: if a component is still too large, the
   same procedure is applied inside that component and the two separators are
   merged. Components shrink geometrically, so the added vertices sum to
   `O(sqrt(n))`.
4. Several BFS roots (the ends and the middle of a BFS pseudo-diameter, which is
   the centroid of a tree) are tried and the best result is kept; the level
   structure of a sparse graph depends strongly on the root.
5. Disconnected inputs follow the end of the proof of Theorem 4: a component
   larger than `2/3 n` is handled recursively, a component of at least `n/3` is
   one side, and otherwise the components are packed greedily.

Running time is linear in `n + m`; the constant is small (about 0.2 s for a
`320 x 320` grid, `n = 102400`, and for a random triangulation with `n = 32000`).

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
julia-audit --strict --include-tests
julia-static-analysis
```

These commands and the audit tooling are described in [AGENTS.md](AGENTS.md).
CI runs `Pkg.test()` on Linux x64 for Julia 1.10 and the latest release.

```bash
tools/build            # instantiate the environment and load the package
tools/gate             # audit, static analysis, tests (with Aqua), build (~35 s)
tools/gate --full      # .. plus a standalone Aqua run
```

Pushing does not block on CI. `tools/push` pushes and prints the run URL;
`tools/check_ci` reports the verdict later (exit `0` green, `1` failed,
`2` still running or not registered yet, `3` unknown):

```bash
tools/push            # push, do not wait
tools/check_ci --jobs # later: what happened, per job
```

The test suite covers the embedding layer (rotations, faces, dual, Euler),
triangulation, the generators, and the separator: validity and balance on
structured, random, sparse, and disconnected graphs, the size bound, a
brute-force optimal separator as an oracle on small graphs, and a generous
linear-scaling timing test. `Aqua.test_all` runs as part of `Pkg.test()`.

## License

MIT — see [LICENSE](LICENSE).
