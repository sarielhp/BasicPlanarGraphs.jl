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

## CI is fire-and-forget

Do not block on GitHub Actions. Pushing or dispatching a run means issuing the remote
command and returning; the verdict is a separate, later question.

- Never `gh run watch`, never poll a run in a loop, and never sleep or wait inside a task
  for CI to finish. Report that a run was triggered and move on.
- The one place that inspects results is `tools/check_ci` (Ruby, one shot). It asks about
  HEAD by default (`--commit REV` for another commit, `--any` for the newest run on the
  branch), and its exit status carries the verdict: `0` success, `1` failed/cancelled,
  `2` still running or no run registered yet, `3` the check could not be performed.
  `--jobs` adds a per-job verdict; `--wait` is an explicit opt-in and is not used in
  normal flow.
- `tools/push` pushes the current branch and prints the run URL without waiting.
- Any other tool in `tools/` that touches CI follows the same rule: trigger remotely,
  print the URL, leave the verdict to `tools/check_ci`.
- Surface a CI failure only when the task at hand depends on it, when the user asks for
  the verdict, or when a run was already known to be red.

## Repository lifecycle

The project follows the review-cycle contract in
`~/prog/standards/workflow/README.md`; the settings live in `tools/review_cycle.json`
(target `.`, gate and build, `push: false`). Every entry point is Ruby and self-documenting
(`--help`); none installs anything or touches Git except where stated.

| Command | Behavior |
|---|---|
| `tools/build` | Instantiate the environment and load the package. No tests, no linters, no installs. |
| `tools/gate` | Fast gate: `julia-audit --strict --include-tests src test examples`, `julia-static-analysis src test examples`, the `examples/test_planarity.jl` smoke run (offline, on a 13-vertex wheel), `Pkg.test()` (with Aqua), then one `tools/build`. `--full` adds a standalone Aqua run; `--no-build` skips the build. |
| `tools/check` | Compatibility alias for `tools/gate` (a symlink). |
| `tools/push` | `git push` and return; reports the run URL, never waits. |
| `tools/check_ci` | Verdict of the latest Actions run: exit `0` green, `1` failed, `2` pending, `3` unknown. |
| `tools/bump` | Cut a release: minor by default, changelog, gate, commit, tag, push, GitHub Release. |

`tools/gate` is the entry point for "is this good"; `tools/build` for "does it load". A
gate failure stops at the failing step and exits `1`; a missing tool exits `3`.

## Autonomous review cycle (`review-cycle`)

The engine and review-audit live in ~/prog/26/ruby/review_cycle and run from PATH. The
repository contract is ~/prog/standards/workflow/README.md; project choices belong in
tools/review_cycle.json. Read project/language guidelines, agent instructions, and tooling
settings before remediation.

| Command | Effect |
|---|---|
| `rtk proxy review-cycle --doctor` | Check host and configured sandbox prerequisites. |
| `rtk proxy review-cycle --setup` | Requested onboarding: configure local bws profiles/tools and project review files. |
| `rtk proxy review-cycle --no-push` | Requested audit/remediation with local commits, gate/build, and differential verification. |
| `rtk proxy review-cycle --loop --no-push` | Bounded loop across lenses; default maximum seven cycles. |
| `rtk proxy review-cycle --push` | Authorized upstream publication with remote SHA verification. |
| `rtk proxy review-cycle --status` | Inspect rotation and history. |
| `rtk proxy review-cycle --dry-run --no-push` | Audit/report without remediation or commits; writes reports/logs and runs pre-flight checks. |

Commits remain local by default; inspect existing explicit push:true settings and use
--no-push when needed. Review commits do not contain lifecycle structured verification or
create release tags. Run `tools/gate` before authorized publication. Install only through
the explicitly authorized project installation workflow.

Command approval permits execution, not unrelated work. Reuse session authorization. Keep
sandboxing and manual approval review enabled; never enlarge writable roots, modify
approval rules, install dependencies, or publish merely because these instructions are
present.

## Examples

`examples/` is held to the same standards as `src` and `test`: it is audited and
format-checked, and `tools/gate` executes `examples/test_planarity.jl` on a small
local graph, so a broken example fails the gate. The example never runs in CI
against the real dataset (1M vertices, network download); run it by hand:

```bash
julia --project=. examples/test_planarity.jl            # 1000, 10000, 100000, all
julia --project=. examples/test_planarity.jl --help
```

## Releases

- Cut a release with `tools/bump` (Ruby). It bumps `Project.toml`, moves the
  `[Unreleased]` changelog entries under `## [X.Y.Z] - YYYY-MM-DD`, runs
  `tools/gate`, commits `release: vX.Y.Z`, creates the annotated tag, pushes the
  branch and tag, and creates the GitHub Release from those notes.
- The default increment is **minor** (`0.1.0 -> 0.2.0`) while the package is
  pre-1.0. Use `--patch` for fixes only; `--major` (1.0.0) is an explicit
  decision. `--dry-run` prints the plan; `--no-push` stops after the local
  commit and tag.
- **Do not register this package in the General registry**, and do not open a
  registry PR, for the foreseeable future. Installation stays
  `Pkg.add(url = "https://github.com/sarielhp/BasicPlanarGraphs.jl")`.
- CI never creates tags or GitHub releases; the maintainer uses `tools/bump`.

## Package notes

- Public repository: <https://github.com/sarielhp/BasicPlanarGraphs.jl>, MIT licensed,
  default branch `main`. CI (`.github/workflows/ci.yml`) runs `Pkg.test()` via
  `julia-actions/julia-runtest` on Linux x64 only, for Julia 1.10 and the latest release;
  build and test locally before pushing. `Manifest.toml` is git-ignored (it is a package,
  not an application), so a fresh clone resolves its own.
- CI must not generate or release binaries. The workflow only builds the package
  in-place to run the tests: no compiled/binary artifacts, no `Pkg.build` outputs or
  uploaded build products, no GitHub releases, tags, or uploads to a package registry.
  This package is pure Julia and has no native code to compile.
- `src/embedded.jl` holds the dart-based rotation-system core (`DartGraph`); the public
  `PlanarGraph` in `src/planar_graph.jl` wraps a `Graphs.jl` `SimpleGraph` plus rotations.
- The separator (`src/separator.jl`) follows Lipton–Tarjan: BFS levels on the input graph,
  middle-band triangulation for cycle candidates, and a refinement loop that separates the
  largest remaining component. Every result comes from `_split_and_pack`, so the returned
  partition is always valid; only balance and separator size need testing.
- Tests use `Graphs.jl` primitives as an independent check and a brute-force optimal
  separator (n <= 16 only) as an oracle.
