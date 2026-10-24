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

AMD Ryzen 7 5700G (16 CPUs), Linux, Julia 1.13.2, `--reps=3`, `--seed=3`.

### Operations on the same graph

| operation | random triangulation, n=20 000, m=59 994 | grid 200×200, n=40 000, m=79 600 | wheel, n=50 001, m=100 000 |
|---|---:|---:|---:|
| `faces` | 19.2 ms | 11.8 ms | 23.3 ms |
| `nfaces` | 13.9 ms | 11.3 ms | 16.4 ms |
| `is_planar_embedding` | 15.3 ms | 14.7 ms | 17.7 ms |
| `dual` | 22.4 ms | n/a (multigraph) | 32.7 ms |
| `separator` | **281.3 ms** | 47.9 ms | 179.9 ms |
| `planar_embedding(::PlanarGraph)` | 14.4 ms | 16.7 ms | 17.8 ms |
| `triangulate` | 41.2 ms | 56.1 ms | 76.4 ms |
| `grid_drawing` | 52.5 ms | 65.7 ms | 91.5 ms |
| `cycle_separator` | 76.9 ms | 94.7 ms | 95.6 ms |
| `planar_embedding(::SimpleGraph)` | 150.9 ms | 169.5 ms | 239.4 ms |
| `is_planar_linear(::SimpleGraph)` | 143.8 ms | 164.8 ms | 238.9 ms |

`dual` is `n/a` for a grid because a corner cell shares two edges with the outer
face, so the dual has parallel edges, which a simple `PlanarGraph` cannot store.

### `separator` scaling

| n | random triangulation | grid | wheel |
|---:|---:|---:|---:|
| 5 000 | 55.2 ms | 4.0 ms | 11.5 ms |
| 10 000 | 132.7 ms | 7.8 ms | 28.3 ms |
| 20 000 | 260.7 ms | 20.3 ms | 57.2 ms |
| 40 000 | 585.8 ms | 38.3 ms | 123.8 ms |

### Where `separator` spends it

`--suite=phases`, random triangulation, n=40 000:

```
root=1      |band|=39664  triangulate=failed  tri=33.9 ms  band_cycles=56.5 ms
root=28305  |band|=39491  triangulate=failed  tri=27.4 ms  band_cycles=61.7 ms
root=2      |band|=39959  triangulate=ok      tri=45.0 ms  band_cycles=94.7 ms
root=35649  |band|=39810  triangulate=failed  tri=20.2 ms  band_cycles=56.8 ms
```

## Reading the numbers

- **`separator` is the outlier**, by roughly 5× over the next slowest operation
  on the same graph. The cause is input shape, not a superlinear step: a random
  triangulation has `O(log n)` BFS diameter, so Lipton–Tarjan's middle band is
  almost the whole graph, while on a grid it is a thin strip. `_best_connected`
  repeats the band triangulation for up to four roots, and three of the four
  bands here are not 2-connected, so the work is discarded after `triangulate`
  throws. When a band does work, up to six candidate cycles are each evaluated by
  a full `_split_and_pack` pass.
- The timings are sensitive to the seed even at fixed `n`: the same 20 000-vertex
  triangulation with another seed took 740 ms instead of 281 ms, because a
  different root makes the band triangulation succeed and the candidate cycles
  get evaluated. Scaling stays linear in `n` in every observed case.
- One stage of `separator` is formally `O(n log n)`: `_balanced_darts`
  (`src/separator.jl`, line 394) sorts all dual-tree edges to find the most
  balanced ones.
  The file header claims `O(n + m)`, so this sort is the one place where the
  claim is stronger than the implementation.
- `planar_embedding(::SimpleGraph)` is about 10× `is_planar_embedding` per edge.
  That is expected in kind: the former runs the Left-Right test and builds a
  rotation system, the latter only validates one that already exists.
- Everything else sits between about 150 and 1 300 ns per edge, which is the
  range these operations are expected to occupy.
