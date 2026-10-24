## Julia Tooling and Authorization

Follow the Julia [GUIDELINES](~/prog/standards/julia/GUIDELINES.md) and
[RATIONALE](~/prog/standards/julia/RATIONALE.md). Code is formatted with
JuliaFormatter's default style (`.JuliaFormatter.toml`); formatting changes go in their
own commit.

Use the exact command forms when applicable:

| Command | Use |
|---|---|
| `rtk proxy julia-audit --strict` | Depth, complexity, size, idiom and trap audit |
| `rtk proxy julia-static-analysis [paths]` | Syntax, formatting; JET and imports for a project |
| `rtk proxy julia-static-analysis --aqua` | Aqua package checks (also run by the test suite) |
| `rtk proxy julia --project=. -e 'using Pkg; Pkg.test()'` | Run the test suite |

Run routine audits and tests for authorized work. Existing command approval allows
execution; it does not authorize commits, releases, registration, or changes outside
the task. Reuse authorization already given in the session. Keep sandboxing and manual
approval review enabled. Bootstrap the analysis packages separately, when authorized, with
`julia-install-tools`; the gates never install anything. Never enlarge an allow marker
(`# julia-audit: allow-<rule>`) to hide a new finding.

Before any refactor of numerical code, record the outputs on fixed inputs and require
bit-identical results afterwards; a result change is its own deliberate commit.

## Package notes

- Public repository: <https://github.com/sarielhp/BasicPlanarGraphs.jl>, MIT licensed,
  default branch `main`. CI (`.github/workflows/ci.yml`) runs `Pkg.test()` via
  `julia-actions/julia-runtest` on Julia 1.10 and the latest release, Linux x64, macOS
  aarch64 and Windows x64; build and test locally before pushing. `Manifest.toml` is
  git-ignored (it is a package, not an application), so a fresh clone resolves its own.
- `src/embedded.jl` holds the dart-based rotation-system core (`DartGraph`); the public
  `PlanarGraph` in `src/planar_graph.jl` wraps a `Graphs.jl` `SimpleGraph` plus rotations.
- The separator (`src/separator.jl`) follows Lipton–Tarjan: BFS levels on the input graph,
  middle-band triangulation for cycle candidates, and a refinement loop that separates the
  largest remaining component. Every result comes from `_split_and_pack`, so the returned
  partition is always valid; only balance and separator size need testing.
- Tests use `Graphs.jl` primitives as an independent check and a brute-force optimal
  separator (n <= 16 only) as an oracle.
