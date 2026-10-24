# Generators and triangulation

Small embedded constructors, mainly for examples, benchmarks and tests, plus the
operation that turns any embedding into a triangulation. All of them return a
`PlanarGraph` whose rotation system is a valid planar embedding for the intended
drawing.

## `planar_path(n)`

Path `1 - 2 - ... - n`, embedded in the plane (both sides of the path are the two
faces). `n >= 1`, otherwise `ArgumentError`. For an interior vertex the rotation
is `[v-1, v+1]`; `nfaces(planar_path(n)) == 1` for `n >= 2`, and a single vertex
has no faces:

```julia
nv(planar_path(1)), nfaces(planar_path(1))   # (1, 0)
ne(planar_path(5))                           # 4
```

## `planar_cycle(n)`

Simple cycle on `n >= 3` vertices, counter-clockwise, with an inner and an outer
face (`nfaces == 2`).

```julia
faces(planar_cycle(4))    # [[1, 4, 3, 2], [1, 2, 3, 4]]
```

## `planar_grid(rows, cols; diagonals = false)`

The `rows x cols` grid, vertices numbered `(r-1)*cols + c`, neighbours listed
counter-clockwise (right, up, left, down). `rows, cols >= 1`, otherwise
`ArgumentError`. The outer face contains the whole boundary.

With `diagonals = true` each cell also gets its top-left to bottom-right
diagonal, so every inner face is a triangle:

```julia
ne(planar_grid(2, 2))                      # 4    (a square)
ne(planar_grid(2, 2; diagonals = true))    # 5
nfaces(planar_grid(2, 2; diagonals = true))# 3
```

A single row is a path, a single column is a path, and a `1 x 1` grid is one
vertex.

## `planar_wheel(rim)`

Wheel with `rim >= 3` rim vertices and a hub. The hub is vertex `1`, the rim is
`2:rim+1` in cyclic order, and every rim vertex lists `[1, prev, next]`. The hub
has degree `rim`; the graph has `rim + 1` vertices and `2*rim` edges.

```julia
nv(planar_wheel(4)), ne(planar_wheel(4)), nfaces(planar_wheel(4))   # (5, 8, 5)
degree(planar_wheel(4), 1)                                          # 4
```

## `random_triangulation(n; rng = Random.default_rng())`

A random **maximal planar** graph on `n >= 3` vertices, built incrementally:

1. start from the embedded triangle `1-2-3`;
2. for `v = 4:n`, pick a face uniformly at random (through `rng`) and insert `v`
   into it, replacing one triangular face by three.

The result is always a triangulation: `ne == 3n - 6` and `nfaces == 2n - 4` for
`n >= 3`, and `is_planar_embedding` holds. It is a uniform random *stacked*
triangulation, not a uniform random planar graph — the construction only reaches
graphs that can be built by this face-splitting process.

```julia
g = random_triangulation(8; rng = MersenneTwister(1))
ne(g), nfaces(g)          # (18, 12)
```

Pass an explicit `rng` for reproducibility; the default is
`Random.default_rng()`.

## Triangulation

```julia
triangulate(g) -> PlanarGraph   # or DartGraph for a DartGraph
```

Adds edges inside faces until every face is a triangle, preserving the vertex set
and the embedding. Chords are inserted inside the faces of the **given** rotation
system, so the result is verified before it is returned:

- if a face is not a simple cycle and no chord can split it, or
- if the rotation system does not describe a planar embedding in the first
  place,

an `ArgumentError` is raised instead of returning a graph that only looks like a
triangulation. For a 2-connected planar graph the result is maximal planar.

```julia
ne(planar_cycle(4))              # 4
ne(triangulate(planar_cycle(4))) # 6   both diagonals, maximal planar
nfaces(triangulate(planar_cycle(4)))# 4

d = dual(triangulate(planar_grid(3, 4)))
all(==(3), degree.(Ref(d), 1:nv(d)))   # true: 3-regular dual, Euler holds
```

The output type follows the input: `triangulate(::PlanarGraph)` returns a
`PlanarGraph`, `triangulate(::DartGraph)` a `DartGraph`. The separator relies on
this: the middle band must be triangulated before its fundamental cycles are
enumerated, and a band that cannot be triangulated is detected here rather than
by an implicit assumption later.

## Errors at a glance

| Call | Condition | Error |
|---|---|---|
| `planar_path(0)` | `n < 1` | `ArgumentError` |
| `planar_cycle(2)` | `n < 3` | `ArgumentError` |
| `planar_grid(0, 3)` | a dimension `< 1` | `ArgumentError` |
| `planar_wheel(2)` | `rim < 3` | `ArgumentError` |
| `random_triangulation(2)` | `n < 3` | `ArgumentError` |
| `triangulate(g)` | face is not a simple cycle, or the embedding is not planar | `ArgumentError` |
