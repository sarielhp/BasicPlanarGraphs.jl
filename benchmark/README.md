# Benchmarks

`run.jl` is a running-time survey of the public entry points that the
documentation describes as linear. It answers one question: on the same input,
which of them is disproportionately expensive?

```bash
julia --project=. benchmark/run.jl                  # every suite, default sizes
julia --project=. benchmark/run.jl --quick          # small graphs, smoke test
julia --project=. benchmark/run.jl --suite=operations --reps=5
julia --project=. benchmark/run.jl --suite=separator
julia --project=. benchmark/run.jl --suite=phases   # internals, diagnostic
julia --project=. benchmark/run.jl --help
```

Every call is warmed once and the best of `--reps` wall-clock runs is reported,
so compilation and GC pauses do not dominate. The numbers are machine-dependent
and are **not** part of `tools/gate` or CI; the only timing assertion in the
suite is the generous linear-scaling check in `test/scaling_test.jl`.

## Reference run

AMD Ryzen 7 5700G (16 CPUs), Linux, Julia 1.13.2, `--reps=3`, `--seed=3`, at
version 0.2.3 plus the fixes of the review branch (`triangulate` by ladder).

### Operations on the same graph

| operation | random triangulation, n=20 000, m=59 994 | grid 200×200, n=40 000, m=79 600 | wheel, n=50 001, m=100 000 |
|---|---:|---:|---:|
| `faces` | 18.2 ms | 7.5 ms | 9.6 ms |
| `nfaces` | 12.9 ms | 7.0 ms | 9.4 ms |
| `is_planar_embedding` | 12.2 ms | 9.6 ms | 11.5 ms |
| `dual` | 27.7 ms | n/a (parallel edges) | 27.5 ms |
| `separator` | 57.5 ms | 23.3 ms | 46.1 ms |
| `planar_embedding(::PlanarGraph)` | 16.6 ms | 9.7 ms | 11.3 ms |
| `triangulate` | 41.1 ms | 61.2 ms | 74.8 ms |
| `grid_drawing` | 55.9 ms | 63.7 ms | 85.1 ms |
| `cycle_separator` | 37.9 ms | 38.7 ms | 31.1 ms |
| `planar_embedding(::SimpleGraph)` | 146.9 ms | 150.2 ms | 216.8 ms |
| `is_planar_linear(::SimpleGraph)` | 132.0 ms | 148.9 ms | 198.9 ms |

`dual` is `n/a` for a grid because a corner cell shares two edges with the outer
face, so the dual has parallel edges, which a simple `PlanarGraph` cannot store.

### `separator` scaling

| n | random triangulation | grid | wheel |
|---:|---:|---:|---:|
| 5 000 | 9.9 ms | 2.2 ms | 3.2 ms |
| 10 000 | 20.6 ms | 5.1 ms | 6.6 ms |
| 20 000 | 45.7 ms | 9.6 ms | 13.8 ms |
| 40 000 | 141.4 ms | 20.5 ms | 34.0 ms |

### Where `separator` spends it

`--suite=phases`, random triangulation, n=10 000 (`--quick`): the virtual
triangulation takes 2.9 ms, the choice of roots 1.6 ms, the level separator of
each of the four roots about 3 ms, and the whole call 29 ms.

## Reading the numbers

- Every operation costs between about 90 and 2 500 ns per edge. The slowest is
  `planar_embedding(::SimpleGraph)` (and so `is_planar_linear`): it runs the
  Left-Right test and builds the rotation system, while `is_planar_embedding` only
  validates one that exists, so a factor of about 10 between them is expected.
- `separator` on a random triangulation takes about 2 µs per vertex at 5 000 to
  20 000 vertices and 3.5 µs at 40 000; the algorithm is linear, and the rise is
  most likely cache behaviour, which this survey did not investigate. On a grid
  or a wheel the cost per vertex is flat within 30 %.
- The timings are machine-dependent and vary from run to run; the only timing
  assertion in the test suite is the generous scaling check in
  `test/scaling_test.jl`.
