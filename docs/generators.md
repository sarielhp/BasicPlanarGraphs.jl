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
and the embedding, in `O(nv + ne)` time. No edge is added twice, so a simple graph
stays simple, and the result is maximal planar (`3n - 6` edges) for a connected
graph with at least three vertices. Bridges and cut vertices are handled: a face
boundary that repeats a vertex is cut at every repeat first.

An `ArgumentError` is raised instead of returning a graph that only looks like a
triangulation if the rotation system does not describe a planar embedding, or if
a component is a single edge (no chord can make its face a triangle).

A boundary that repeats a vertex (a bridge or a cut vertex) first loses its
repeated occurrences: such an occurrence can always be clipped, because the
vertex appears twice, so its two neighbours lie in different blocks around it
and are not adjacent. After a clip the next occurrence is skipped, so the
chords do not chain through one vertex.

The remaining simple cycle `p_1 … p_k` is cut by a **ladder**: `p_i` is joined to
`p_{k+1-i}`. The rungs are nested, so they cut the face into quadrangles, and
each quadrangle `a, b, c, d` takes the diagonal `ac` or `bd`: two edges outside
the face joining `a, c` and `b, d` would cross, so one of them is new. A vertex
gets at most one rung and two diagonals, a constant degree, instead of a fan's
`k - 3` at one vertex.

A rung whose ends are already adjacent outside the face is skipped, which merges
its quadrangles into a larger region. All rungs of that region are edges outside
the face, and each rung of a **slanted ladder** (started one position later, so
the vertices are paired at a slant) crosses one of them in the sense of
alternating ends on the boundary. By planarity it is therefore not blocked, so
the slanted ladder cuts the region into quadrangles and a final triangle or
quadrangle too. A vertex gains at most 5 chords per face (3 unless it ends a
rung of the first ladder), and the bound is attained. The only edges that can
block a chord are non-crossing edges of the disk outside the face, which is all
the argument uses. If the input is not a planar embedding, a quadrangle may have
both diagonals blocked and an `ArgumentError` is raised.

The adjacency test scans the smaller of the two degrees. A test that succeeds is
charged to the chord it adds, one that fails to the edge it finds, and the sum
over the edges of the smaller degree is `O(n)` in a planar graph (Chiba and
Nishizeki). So the whole triangulation is `O(nv + ne)`, with no hashing and no
search for a free apex.

```julia
ne(planar_cycle(4))              # 4
ne(triangulate(planar_cycle(4))) # 6   both diagonals, maximal planar
nfaces(triangulate(planar_cycle(4)))# 4

d = dual(triangulate(planar_grid(3, 4)))
all(==(3), degree.(Ref(d), 1:nv(d)))   # true: 3-regular dual, Euler holds
```

The output type follows the input: `triangulate(::PlanarGraph)` returns a
`PlanarGraph`, `triangulate(::DartGraph)` a `DartGraph`.

## Errors at a glance

| Call | Condition | Error |
|---|---|---|
| `planar_path(0)` | `n < 1` | `ArgumentError` |
| `planar_cycle(2)` | `n < 3` | `ArgumentError` |
| `planar_grid(0, 3)` | a dimension `< 1` | `ArgumentError` |
| `planar_wheel(2)` | `rim < 3` | `ArgumentError` |
| `random_triangulation(2)` | `n < 3` | `ArgumentError` |
| `triangulate(g)` | the embedding is not planar, or a component is a single edge | `ArgumentError` |
