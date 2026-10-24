# Graph files

Planar graphs are written and read as sectioned plain text: one line per record,
one optional section for the embedding. A file is readable by a person, stable
under version control, and reads back into the same graph.

```julia
using Graphs, BasicPlanarGraphs

g = planar_grid(3, 3; diagonals = true)
write_planar_graph("grid.txt", g)

h = read_planar_graph("grid.txt")   # PlanarGraph, same rotation
s = read_graph("grid.txt")          # SimpleGraph, the same edges
```

Both functions also take an `IO` as their first argument, so they work with
`IOBuffer`, pipes, and [`Base.open`](@ref).

| Call | Writes | Reads |
|---|---|---|
| `write_planar_graph(io_or_path, g::PlanarGraph)` | vertices, edges, rotation | — |
| `write_graph(io_or_path, g::SimpleGraph)` | vertices, edges | — |
| — | | `read_planar_graph(io_or_path) -> PlanarGraph` |
| — | | `read_graph(io_or_path) -> SimpleGraph` |

`read_graph` on a file that has a rotation section ignores the rotation when the
file also has an `edges` section; a file with *only* a rotation section describes
its graph through that rotation, and the returned graph is the one the rotation
induces. `read_planar_graph` requires a rotation section.

## The format

```
# BasicPlanarGraphs 1
vertices 4
edges
1 2
1 4
2 3
3 4
rotation
1: 2 4
2: 3 1
3: 4 2
4: 1 3
end
```

- The first line is a comment that names the format and its version; it is the
  only comment the writers emit. Comments start with `#` and run to the end of
  the line, blank lines are ignored, and both are accepted anywhere.
- `vertices <n>` fixes the vertex count, which is what makes isolated vertices
  and the empty graph representable. It is optional — without it the count is the
  largest label in the file — and it may appear anywhere.
- `edges` opens a section of edge lines. Each line is `u v`, the two labels in
  either order; an edge may be listed once per file.
- `rotation` opens a section of rotation lines. Each line is `v: w1 w2 …`, the
  neighbours of `v` counter-clockwise — the same order as
  [`rotation`](api.md) and [`faces`](embedding.md) use. A line `v:` gives `v` an
  empty rotation. Lines may appear in any order, and every vertex with a
  rotation needs exactly one line.
- A section lasts until the next `edges` or `rotation` line, or until `end`,
  which closes the file; a file that stops after its last record is also
  accepted.

Vertex labels are `1:n`; they must be positive integers. The two writers emit
`vertices`, `edges`, and (for a `PlanarGraph`) `rotation`, then `end`.

## What is rejected

Reading validates the syntax and the consistency of the file, and raises
`ArgumentError` naming the source and the line:

| Input | Error |
|---|---|
| `1 2` with no section open | `unexpected line ... (missing a section keyword?)` |
| `edges` then `1` or `1 x` | ``an edge is `u v` `` / `` `x` is not an integer `` |
| label `0`, or a self-loop `1 1` | `vertex labels start at 1` / `self-loops are not supported` |
| the same edge twice | `duplicate edge 1-2` |
| `vertices many`, `vertices -1` | `` `many` is not an integer `` / `vertex count must be non-negative` |
| two rotation lines for one vertex | `duplicate rotation for vertex 2` |
| a rotation that lists the vertex itself | `vertex 2 lists itself as a neighbour` |
| fewer rotation lines than vertices | `rotation has 3 lines for 4 vertices` |
| `edges` and `rotation` differing | ``the `edges` and `rotation` sections describe different graphs`` |
| `read_planar_graph` without a `rotation` section | ``no `rotation` section; use read_graph for an unembedded graph`` |

The rotation system goes through the [`PlanarGraph`](api.md)
constructor, so a rotation that is not a permutation of the neighbours of its
vertex is an error. The format **does not** check planarity: an embedding that
violates Euler's formula is stored as written, and
[`is_planar_embedding`](api.md) reports it. Use
[`planar_embedding`](planarity.md) to recompute one.

## Round trip

Writing a `PlanarGraph` and reading it back returns an equal graph, rotation
included; writing a `SimpleGraph` and reading it back returns the same graph
without an embedding. The rotation is written in the direction the graph stores
it, so the file is a faithful picture:

```julia
io = IOBuffer()
write_planar_graph(io, PlanarGraph([[2, 3], [3, 1], [1, 2]]))
print(String(take!(io)))
```

```
# BasicPlanarGraphs 1
vertices 3
edges
1 2
1 3
2 3
rotation
1: 2 3
2: 3 1
3: 1 2
end
```

Edges are written in the order `Graphs.jl` enumerates them, so two runs on the
same graph produce byte-identical files.
