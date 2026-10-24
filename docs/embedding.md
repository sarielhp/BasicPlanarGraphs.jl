# Embedding: the rotation system

A `PlanarGraph` is an undirected simple graph plus a combinatorial embedding.
Nothing is drawn: the embedding is the counter-clockwise cyclic order of the
neighbours around every vertex, called a **rotation system**. Two rotation
systems describe the same embedding when they differ only by a cyclic shift or a
reversal at each vertex.

## The type

```julia
struct PlanarGraph{T<:Integer}
    graph::SimpleGraph{T}
    rotation::Vector{Vector{T}}   # rotation[v] = neighbours of v, counter-clockwise
end
```

`T` is the vertex type. It is inferred from the rotation, so `PlanarGraph{Int32}`
works, and the `Graphs.jl` interface is forwarded to the wrapped graph.

## Constructors

```julia
PlanarGraph(graph, rotation)   # graph::SimpleGraph, rotation of matching size
PlanarGraph(rotation)          # build the graph from the rotation
planar_graph(args...)          # alias for PlanarGraph(args...)
```

`PlanarGraph(rotation)` derives the edges from the rotation, so the two agree by
construction. It validates as it builds:

| Input | Error |
|---|---|
| a neighbour `w` with `w < 1` or `w > n` | `ArgumentError("neighbour $w of vertex $v is out of range")` |
| `v` listing itself | `ArgumentError("vertex $v lists itself as a neighbour")` |
| a non-integer entry | `ArgumentError("rotation entries must be integers, ...")` |
| the same neighbour twice | `ArgumentError("rotation of vertex $v is not a permutation of its neighbours")` |

`PlanarGraph(graph, rotation)` performs the mirror check against an existing
graph: the rotation must have one entry per vertex, and `rotation[v]` must be a
permutation of `neighbors(graph, v)`. A rotation that is missing an edge,
mentions a non-neighbour, or repeats an entry is rejected with an
`ArgumentError` naming the vertex. The constructor therefore *always* yields a
consistent graph-plus-embedding pair; it does **not** check planarity, which is
what [`is_planar_embedding`](#eulers-formula-and-planarity) is for.

```julia
t = PlanarGraph([[2, 3], [3, 1], [1, 2]])   # counter-clockwise triangle
nv(t), ne(t), nfaces(t)                     # (3, 3, 2)
faces(t)                                    # [[1, 2, 3], [1, 3, 2]]
```

## Equality and hashing

`==` requires the same underlying graph and the same cyclic order at every
vertex, up to a starting point and direction. So the same triangle written with
the rotations reversed is equal, while the same graph with a different embedding
is not. `hash` is consistent with this, so `PlanarGraph`s can be dictionary keys.

## Queries

| Function | Result |
|---|---|
| `neighbor_order(g, v)` | the stored neighbour vector of `v`, counter-clockwise. It is the internal vector: do not mutate it. |
| `rotation(g)` | a fresh copy of the whole rotation system. |
| `faces(g)` | boundary walks, one `Vector{T}` per face. |
| `nfaces(g)` | `length(faces(g))`. |
| `is_planar_embedding(g)` | Euler's formula on every component. |
| `dual(g)` | the dual embedding. |
| `triangulate(g)` | maximal planar supergraph; see [Generators and triangulation](generators.md#triangulation). |

Because `PlanarGraph` forwards `nv`, `ne`, `vertices`, `edges`, `neighbors`,
`degree` and `has_edge` to its graph, everything that works on a `SimpleGraph`
works here — `Graphs.jl` algorithms that only read adjacency accept a
`PlanarGraph` directly.

## Faces

A face is traced by walking darts: for a dart `d`, the next dart of the face on
its left is `prev_rot[twin[d]]` in the dart core (§ below), which in the public
type is the same rule expressed through the stored rotations.

`faces(g)` returns the vertices along each walk, starting at an arbitrary vertex.
Properties worth knowing:

- The number of walks is the number of faces of the embedding.
- A vertex on a **bridge** or the articulation point of a **cut** appears more
  than once in a walk, because the walk traverses it twice.
- An **isolated vertex** has no darts and therefore contributes no walk: a
  one-vertex graph has `nfaces == 0` and `faces(g) == []`.
- For a connected graph, `ne - nv + nfaces == 2`; this is the invariant that
  `is_planar_embedding` checks component by component.

```julia
c = planar_cycle(4)
faces(c)          # [[1, 4, 3, 2], [1, 2, 3, 4]]  (inner and outer walk)
nfaces(c)         # 2
```

## Euler's formula and planarity

`is_planar_embedding(g)` returns `true` when every connected component satisfies

```
n_c - m_c + f_c = 2
```

with `f_c` counted per component; a component with no edges is a lone vertex and
counts as one face. A rotation system satisfying this on every component is a
planar embedding, so the check is exact rather than heuristic.

The constructor does not call it. That separation is deliberate: it is what makes
it possible to build a rotation system that describes a non-planar graph and test
that it is rejected, and what lets `dual` and `triangulate` report a bad input
instead of silently producing nonsense.

## Duality

`dual(g)` has one vertex per face walk of `g` and one edge per primal edge: each
dart and its twin lie on the two faces that share the edge, so the dual of a
triangulation is 3-regular. The rotation of the dual is built from the primal
face orbits, so the dual is itself a `PlanarGraph` and the operation can be
iterated.

The primal graph must be **bridge-free**. A bridge has the same face on both
sides, so its dual edge would be a self-loop, which `PlanarGraph` does not
represent; `dual(planar_path(3))` therefore throws an `ArgumentError` rather than
returning a broken object.

```julia
d = dual(triangulate(planar_grid(2, 2)))
nv(d), ne(d), all(==(3), degree.(Ref(d), 1:nv(d)))    # (4, 6, true)
```

## The dart core

`src/embedded.jl` holds `DartGraph`, an embedded **multigraph** used by the
algorithms. A *dart* is an oriented edge; every edge contributes a dart `d` and
its twin, pointing the other way. Darts are numbered `1:2m`.

```julia
mutable struct DartGraph
    n::Int
    deg::Vector{Int}       # degree of each vertex
    first::Vector{Int}     # one dart leaving each vertex
    tail::Vector{Int}      # tail[d]
    head::Vector{Int}      # head[d]
    twin::Vector{Int}      # twin[d]
    next_rot::Vector{Int}  # cyclic doubly linked neighbour list, counter-clockwise
    prev_rot::Vector{Int}
end
```

The rotation is the pair of linked lists: around a vertex `v`, following
`next_rot` from `first[v]` visits every dart leaving `v` in counter-clockwise
order and comes back. The face on the left of a dart is then traced by

```
face_successor(d) = prev_rot[twin[d]]
```

and `face_orbits` collects those cycles, returning the orbits and a
`face_of` array that maps each dart to its face number.

The core is deliberately more permissive than the public type: it allows
parallel edges, and `add_dart!` inserts a dart into a rotation without requiring
the twin, so structures can be built incrementally (`_match_pending_twins!`
pairs them up afterwards and rejects parallel edges and self-loops).

Supporting operations, all linear in the size of the graph:

| Function | Purpose |
|---|---|
| `add_dart!(g, u, v, anchor=0)` | insert `u -> v` after `anchor` in the rotation of `u` |
| `add_edge!(g, u, v, anchor_u=0, anchor_v=0)` | add both darts of `{u, v}`, linked as twins |
| `darts_of(g, v)`, `rotation_lists(g)` | in-order darts, and rotations as vertex lists |
| `face_orbits(g)` | face walks and `face_of` |
| `euler_holds(g)` | Euler's formula per component |
| `count_components(g)`, `component_labels(g)` | connectivity |
| `components_without(g, removed)` | components after deleting a vertex mask |
| `bfs_tree(g, roots; keep)` | `parent`, `parent_dart`, `depth`, BFS `order` |
| `level_structure(g, roots; keep)` | just what the separator needs from a BFS |
| `restrict(g, keep)`, `induced(g, keep)` | induced subgraphs, keeping or renumbering vertices |
| `copy(g)` | deep copy |

`PlanarGraph(g::DartGraph)` and `DartGraph(g::PlanarGraph)` convert between the
public type and the core; the conversion is linear and validates what the public
type guarantees (no parallel edges).

## Complexity

Every embedding query is linear in `n + m`:

| Operation | Time | Space |
|---|---|---|
| construction and validation | `O(n + m)` | `O(n + m)` |
| `faces`, `nfaces` | `O(n + m)` | `O(n + m)` |
| `is_planar_embedding` | `O(n + m)` | `O(n + m)` |
| `dual` | `O(n + m)` | `O(n + m)` |
| `triangulate` | each face is examined once; inserting a chord copies the two new orbits | `O(n + m)` |

Validation uses a sort per vertex to compare against the neighbour list, so
construction is `O(sum_v d_v log d_v)` in the worst case; for the bounded-degree
graphs this package targets that is linear.
