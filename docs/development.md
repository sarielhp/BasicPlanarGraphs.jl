# Development

## Build, gate, and tests

The project follows the Julia standards in `~/prog/standards/julia` and the
review-cycle contract in `~/prog/standards/workflow`; the settings live in
`tools/review_cycle.json`.

```bash
tools/build            # instantiate the environment and load the package
tools/gate             # audit, static analysis, examples, tests (with Aqua), build (~45 s)
tools/gate --full      # .. plus a standalone Aqua run
tools/gate --no-build  # .. without the build step
```

`tools/gate` runs, in order:

| Step | Command | What it checks |
|---|---|---|
| audit | `julia-audit --strict --include-tests src test examples benchmark` | nesting ≤ 4, cognitive complexity ≤ 15, method tiers (80 / 120 lines), file size, argument count, idiom and trap rules |
| static analysis | `julia-static-analysis src test examples benchmark` | every file parses, JuliaFormatter would change nothing, JET type inference, ExplicitImports |
| examples | `examples/test_planarity.jl --graph … --sizes 10,all` | the planarity example runs end to end on a 13-vertex wheel written by the gate: read, planarity test, embedding, output, separator, verification |
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

`examples/` is part of the checked code: the audit and the format check cover it,
and `tools/gate` runs `examples/test_planarity.jl` on a 13-vertex wheel that the
gate writes itself, so a broken example fails the gate without any network access.
`benchmark/` is audited and format-checked the same way, but the gate never runs
it, because its timings are machine-dependent.

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

## The planarity example

`examples/test_planarity.jl` is the end-to-end demo on real data: it downloads a
SNAP road network (TIGER/Line Pennsylvania, ~1.09M vertices), cleans the directed
listing into a pure undirected graph, then for 1000, 10,000, 100,000 induced
vertices and the whole graph tests planarity, computes and writes an embedding,
runs the separator and verifies it, printing the time and memory of every stage.
Generated files are cached under `/tmp/<user>/BasicPlanarGraphs/`.

```bash
julia --project=. examples/test_planarity.jl                 # the four sizes
julia --project=. examples/test_planarity.jl --sizes 1000,10000
julia --project=. examples/test_planarity.jl --dataset roadNet-TX
julia --project=. examples/test_planarity.jl --graph my.graph.txt   # offline
```

It never runs in CI (network download, several GiB at the full size); the gate
only runs the small `--graph` smoke test.

On an otherwise idle machine (58 GB RAM) the four sizes give:

| vertices | n | m | planar | check | embed | write | separator | \|C\| | balance |
|---|---|---|---|---|---|---|---|---|---|
| first 1000 | 1000 | 1276 | true | 0.00 s | 0.04 s | 0.00 s | 0.00 s | 0 | 50.6 % |
| first 10000 | 10000 | 14209 | true | 0.04 s | 0.05 s | 0.01 s | 0.11 s | 197 | 56.3 % |
| first 100000 | 100000 | 143993 | true | 0.53 s | 0.59 s | 0.13 s | 2.20 s | 2726 | 57.6 % |
| all | 1088092 | 1541898 | true | 8.42 s | 8.47 s | 1.16 s | 20.22 s | 2292 | 49.9 % |

Every separator is verified by the program itself (the three parts partition
`1:n`, no edge joins `A` and `B`, and each side is within `2/3 n`), and every
written embedding satisfies Euler's formula on every component. Peak resident
memory for the whole run is about 4 GiB.

## Benchmarks

`benchmark/run.jl` compares the public entry points that the documentation calls
linear on the same graph, and reports `separator`'s scaling separately. It is a
diagnostic, not a gate: the timings are wall-clock and machine-dependent, so
`tools/gate` and CI never run it. Methodology, the reference numbers and the
interpretation live in [`benchmark/README.md`](../benchmark/README.md).

```bash
julia --project=. benchmark/run.jl --quick          # smoke test
julia --project=. benchmark/run.jl                  # full survey
julia --project=. benchmark/run.jl --suite=phases   # internals of separator
```

## Releases

A release is a Git tag plus a GitHub Release object whose notes are the
matching `CHANGELOG.md` section. The package is not registered in the General
registry and there is no plan to register it for the foreseeable future;
installation stays `Pkg.add(url = "https://github.com/sarielhp/BasicPlanarGraphs.jl")`.

1. Set `version` in `Project.toml` and add the matching
   `## [X.Y.Z] - YYYY-MM-DD` section to `CHANGELOG.md`.
2. Run `tools/gate`.
3. Commit, then create an annotated tag: `git tag -a vX.Y.Z -m "vX.Y.Z"`.
4. `tools/push` pushes the branch; push the tag explicitly with
   `git push origin vX.Y.Z`.
5. Create the GitHub Release from the changelog section, for example
   `gh release create vX.Y.Z --verify-tag --title "vX.Y.Z" --notes-file NOTES.md`
   with `NOTES.md` holding that section's text.

CI never creates a tag, a GitHub release, or a registry entry.

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
