# Development

## Build, gate, and tests

The project follows the Julia standards in `~/prog/standards/julia` and the
review-cycle contract in `~/prog/standards/workflow`; the settings live in
`tools/review_cycle.json`.

```bash
tools/build            # instantiate the environment and load the package
tools/gate             # audit, static analysis, tests (with Aqua), build (~35 s)
tools/gate --full      # .. plus a standalone Aqua run
tools/gate --no-build  # .. without the build step
```

`tools/gate` runs, in order:

| Step | Command | What it checks |
|---|---|---|
| audit | `julia-audit --strict --include-tests` | nesting ≤ 4, cognitive complexity ≤ 15, method tiers (80 / 120 lines), file size, argument count, idiom and trap rules |
| static analysis | `julia-static-analysis` | every file parses, JuliaFormatter would change nothing, JET type inference, ExplicitImports |
| tests | `Pkg.test()` | the suite below, including `Aqua.test_all` |
| build | `tools/build` | the package still loads |

Or run the suite directly:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Malformed input is rejected with `ArgumentError`, never by silently returning a
wrong object; if you add a path that cannot honour its contract, throw.

## Test suite

`test/runtests.jl` only includes the other files; each runs standalone and seeds
its generators.

| File | Covers |
|---|---|
| `embedding_test.jl` | construction and validation, cyclic-order equality, faces and `nfaces`, bridges and cut vertices, `dual`, Euler's formula, `Int32` element types |
| `triangulate_test.jl` | chord insertion, maximal-planarity counts, rejection of non-2-connected and non-planar inputs |
| `generators_test.jl` | vertex and edge counts, degrees, face counts, determinism under a seed |
| `separator_test.jl` | validity, balance and size on structured, random, sparse and disconnected graphs, `balance` handling, determinism |
| `planarity_test.jl` | `is_planar` and `planar_embedding` against `Graphs.is_planar` on random graphs and the classical obstructions, Euler's formula on the result, recomputation of broken rotations, edge and vertex preservation, linear scaling |
| `drawing_test.jl` | integer coordinates inside the bounding box, distinct positions, no vertex on an edge, no crossing edges, the `2n-4 × (n-2)` bounds, components, trivial graphs, rejection of non-planar input |
| `graph_io_test.jl` | round trips for embedded and bare graphs, the file format, isolated vertices and empty graphs, every rejected input |
| `scaling_test.jl` | quadrupling `n` does not blow up the running time (generous factor for loaded machines) |

Two properties make the separator tests meaningful rather than circular:

- **Independent validation.** `check_separator` re-derives the partition check
  with `Graphs.jl` primitives only (`sort(vcat(A, B, C))`, edge endpoints by
  side), so it shares no code with the algorithm.
- **An optimal oracle.** `minimum_balanced_separator` enumerates every vertex
  subset for `n <= 16` and returns the smallest balanced separator, so the test
  can assert that the algorithm is never better than optimal — the direction that
  catches an invalid `C`.

The scaling test compares `separator` on a `60 x 60` and a `120 x 120` grid and
on random triangulations of 2000 and 8000 vertices, allowing a factor of 30 for a
4× size increase.

## CI

`.github/workflows/ci.yml` runs `Pkg.test()` through
`julia-actions/julia-runtest` on Linux (x64) for Julia 1.10 and the latest
release, with coverage from the latest-release job.

CI is fire-and-forget: pushing does not wait for the run.

```bash
tools/push             # push and print the run URL
tools/check_ci         # later: verdict for HEAD (0 green, 1 failed, 2 pending/unknown run, 3 unknown)
tools/check_ci --jobs  # .. with a per-job breakdown
```

## Style

- JuliaFormatter's default style (`.JuliaFormatter.toml`); formatting-only
  changes go in their own commit.
- Explicit imports (`using Graphs: ...`), typed struct fields, concrete
  containers, `!` on mutating functions.
- Library code throws; only scripts exit.
- Numerical behaviour is pinned by tests: a refactor that changes a result is a
  separate, deliberate change, not a cleanup.
