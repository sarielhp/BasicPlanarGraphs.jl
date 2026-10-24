# BasicPlanarGraphs

[![CI](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sarielhp/BasicPlanarGraphs.jl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Planar graphs with an abstract embedding (a rotation system), and a linear-time
planar separator based on Lipton and Tarjan, *A Separator Theorem for Planar
Graphs* (SIAM J. Appl. Math. 36(2), 1979).

## Documentation

Full documentation lives in [`docs/`](docs/README.md):

| Page | Contents |
|---|---|
| [Embedding](docs/embedding.md) | the `PlanarGraph` type, rotation systems, faces, duality, Euler's formula, the dart core |
| [Separator](docs/separator.md) | guarantees, the Lipton–Tarjan stages, refinement, disconnected inputs, complexity |
| [Generators and triangulation](docs/generators.md) | `planar_path`, `planar_cycle`, `planar_grid`, `planar_wheel`, `random_triangulation`, `triangulate` |
| [API reference](docs/api.md) | every exported name with signatures, return types, and errors |
| [Development](docs/development.md) | build, the quality gate, the test suite, CI |

Release history is in [CHANGELOG.md](CHANGELOG.md).

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

## Overview

`PlanarGraph{T}` is a `Graphs.jl` `SimpleGraph{T}` plus, for every vertex, the
counter-clockwise cyclic order of its neighbours:

```julia
triangle = PlanarGraph([[2, 3], [3, 1], [1, 2]])
nv(triangle), ne(triangle), nfaces(triangle)     # 3, 3, 2
faces(triangle)                                  # [[1, 2, 3], [1, 3, 2]]
dual(triangulate(planar_grid(3, 4)))             # 3-regular dual, Euler holds
```

The constructor rejects rotations that are not permutations of the neighbours, so
the stored data always describes a graph; `is_planar_embedding` then checks
Euler's formula on every component. Equality of embeddings is equality up to a
cyclic shift and reversal at each vertex. Details, including the dart-based core
underneath: [Embedding](docs/embedding.md).

`separator(g; balance = 2/3)` returns `(A, B, C)` with `A`, `B`, `C` a partition
of the vertices, no edge between `A` and `B`, `|A|, |B| <= balance * n`, and
`|C| = O(sqrt(n))` — at most `2*sqrt(2n) + 1` on the test families. Running time is
linear in `n + m`, with a small constant. The stages, the guarantees and the
caveats: [Separator](docs/separator.md).

Generators, triangulation, and the full export list:
[Generators](docs/generators.md) and [API reference](docs/api.md).

## Development

```bash
tools/build   # instantiate the environment and load the package
tools/gate    # audit, static analysis, tests (with Aqua), build (~35 s)
```

Pushing does not block on CI: `tools/push` pushes and prints the run URL, and
`tools/check_ci` reports the verdict later (exit `0` green, `1` failed,
`2` pending, `3` unknown). The test suite, the optimal-separator oracle, the
style rules and the CI matrix are described in [Development](docs/development.md).

## License

MIT — see [LICENSE](LICENSE).
