# Public planar graph type: a `Graphs.jl` graph plus a rotation system.

"""
    PlanarGraph(graph, rotation)
    PlanarGraph(rotation)

An undirected planar graph together with a combinatorial embedding, given by
the counter-clockwise cyclic order of the neighbours of every vertex
(`rotation[v]`).

`graph` must contain exactly the edges that occur in `rotation`; the
constructor validates this. Use [`is_planar_embedding`](@ref) to query the
embedding check explicitly.

A `PlanarGraph` is a `Graphs.jl` `AbstractGraph` (undirected, read-only through the
generic interface), so the algorithms of `Graphs.jl` accept it: `connected_components`,
`bfs_parents`, `dijkstra_shortest_paths`, `greedy_color`, `adjacency_matrix`, and so
on. Functions that build a new graph of the same type (`induced_subgraph`, `complement`)
are not supported; use `g.graph` for those. The mutating `Graphs.rem_edge!` is, and
keeps the embedding; see [`delete_edge!`](@ref).

```julia
g = PlanarGraph([[2, 3], [3, 1], [1, 2]])   # triangle, counter-clockwise
nv(g), ne(g), nfaces(g)
```

Two embeddings are equal when they have the same underlying graph and their
rotations agree up to a cyclic shift and reversal at every vertex, that is
when they describe the same embedding of the same graph.
"""
struct PlanarGraph{T<:Integer} <: AbstractGraph{T}
    graph::SimpleGraph{T}
    rotation::Vector{Vector{T}}

    function PlanarGraph{T}(
        graph::SimpleGraph{T},
        rotation::Vector{Vector{T}},
    ) where {T<:Integer}
        _check_rotation(graph, rotation)
        return new{T}(graph, rotation)
    end
end

# The graph is copied, so changing the argument later cannot break the invariant
# that `graph` and `rotation` describe the same edges.
function PlanarGraph(graph::SimpleGraph{T}, rotation::AbstractVector) where {T<:Integer}
    return PlanarGraph{T}(copy(graph), _convert_rotation(rotation, T))
end

function PlanarGraph(rotation::AbstractVector)
    T = _rotation_eltype(rotation)
    n = length(rotation)
    graph = SimpleGraph{T}(n)
    for v = 1:n
        for w in rotation[v]
            (1 <= w <= n) ||
                throw(ArgumentError("neighbour $w of vertex $v is out of range"))
            w == v && throw(ArgumentError("vertex $v lists itself as a neighbour"))
            Graphs.add_edge!(graph, T(min(v, w)), T(max(v, w)))
        end
    end
    return PlanarGraph(graph, rotation)
end

"""    planar_graph(args...) = PlanarGraph(args...)"""
planar_graph(args...) = PlanarGraph(args...)

# Element type for the rotation system: the type of the entries, promoted over
# all of them, and `Int` when there are no entries at all.
function _rotation_eltype(rotation)
    T = Union{}
    for r in rotation, x in r
        x isa Integer ||
            throw(ArgumentError("rotation entries must be integers, got $(typeof(x))"))
        T = T === Union{} ? typeof(x) : promote_type(T, typeof(x))
    end
    return T === Union{} ? Int : T
end

function _convert_rotation(rotation, ::Type{T}) where {T<:Integer}
    all(r -> r isa AbstractVector, rotation) ||
        throw(ArgumentError("rotation must be a vector of vectors"))
    _rotation_eltype(rotation)
    return [T[_to_vertex(T, x) for x in r] for r in rotation]
end

function _to_vertex(::Type{T}, x::Integer) where {T<:Integer}
    typemin(T) <= x <= typemax(T) ||
        throw(ArgumentError("rotation entry $x does not fit the vertex type $T"))
    return T(x)
end

function _check_rotation(graph::SimpleGraph, rot::Vector{Vector{T}}) where {T}
    n = nv(graph)
    length(rot) == n ||
        throw(ArgumentError("rotation has $(length(rot)) entries for $n vertices"))
    # A permutation check with stamps, not two sorts: mark the neighbours of
    # `v`, consume each rotation entry, and leave `-v` behind so a repeat fails.
    mark = zeros(Int, n)
    for v = 1:n
        rot_v = rot[v]
        expected = degree(graph, v)
        length(rot_v) == expected || throw(
            ArgumentError(
                "rotation of vertex $v has $(length(rot_v)) entries, expected $expected",
            ),
        )
        for w in neighbors(graph, v)
            mark[w] = v
        end
        for w in rot_v
            (1 <= w <= n) || throw(
                ArgumentError(
                    "rotation of vertex $v is not a permutation of its neighbours",
                ),
            )
            mark[w] == v || throw(
                ArgumentError(
                    "rotation of vertex $v is not a permutation of its neighbours",
                ),
            )
            mark[w] = -v
        end
    end
    return nothing
end

# --- Graphs.jl interface -------------------------------------------------

Graphs.nv(g::PlanarGraph) = nv(g.graph)
Graphs.ne(g::PlanarGraph) = ne(g.graph)
Graphs.vertices(g::PlanarGraph) = vertices(g.graph)
Graphs.edges(g::PlanarGraph) = edges(g.graph)
Graphs.neighbors(g::PlanarGraph, v::Integer) = neighbors(g.graph, v)
Graphs.degree(g::PlanarGraph, v::Integer) = degree(g.graph, v)
Graphs.has_edge(g::PlanarGraph, u::Integer, v::Integer) = has_edge(g.graph, u, v)
Graphs.has_vertex(g::PlanarGraph, v::Integer) = Graphs.has_vertex(g.graph, v)
Graphs.is_directed(::Type{<:PlanarGraph}) = false
Graphs.edgetype(::PlanarGraph{T}) where {T} = Graphs.SimpleEdge{T}
Graphs.outneighbors(g::PlanarGraph, v::Integer) = neighbors(g.graph, v)
Graphs.inneighbors(g::PlanarGraph, v::Integer) = neighbors(g.graph, v)

"""Neighbours of `v` in counter-clockwise order (the embedding's rotation)."""
neighbor_order(g::PlanarGraph, v::Integer) = g.rotation[v]

"""Copy of the rotation system."""
rotation(g::PlanarGraph) = [copy(r) for r in g.rotation]

# --- equality ------------------------------------------------------------

function Base.:(==)(a::PlanarGraph, b::PlanarGraph)
    a.graph == b.graph || return false
    length(a.rotation) == length(b.rotation) || return false
    return all(v -> same_cyclic_order(a.rotation[v], b.rotation[v]), eachindex(a.rotation))
end

# Equal embeddings have equal graphs, whatever shift or direction their rotations
# have, so the hash looks at the graph only.
Base.hash(g::PlanarGraph, h::UInt) = hash(g.graph, hash(:PlanarGraph, h))

"""
    same_cyclic_order(x, y) -> Bool

`y` is `x` read from some starting point, possibly in the opposite direction.
This is the equality of cyclic orders that the rotation system stores.
"""
function same_cyclic_order(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || return false
    isempty(x) && return true
    backwards = reverse(y)
    return any(
        p -> _shift_matches(x, y, p) || _shift_matches(x, backwards, p),
        eachindex(x),
    )
end

# `y` is `x` read from position `p` onwards, wrapping around.
function _shift_matches(x::AbstractVector, y::AbstractVector, p::Integer)
    m = length(x)
    return all(i -> x[mod1(p + i - 1, m)] == y[i], eachindex(y))
end

# --- embedding queries ---------------------------------------------------

"""
    faces(g) -> Vector{Vector{T}}

Boundary walks of the faces of the embedding. Vertices on a bridge or a cut
vertex occur more than once in a walk; the number of walks is the number of
faces.
"""
function faces(g::PlanarGraph{T}) where {T}
    dg = DartGraph(g)
    orbits, _ = face_orbits(dg)
    return [T[tail(dg, d) for d in orbit] for orbit in orbits]
end

"""
    nfaces(g) -> Int

Number of face walks of the embedding, that is `length(faces(g))`. A component
without edges (a single isolated vertex) has no darts and therefore contributes
no walk.
"""
nfaces(g::PlanarGraph) = length(face_orbits(DartGraph(g))[1])

"""
    is_planar_embedding(g) -> Bool

`true` when every connected component satisfies Euler's formula
`n_c - m_c + f_c = 2`. The rotation system is assumed to be well formed (the
constructor guarantees that).
"""
is_planar_embedding(g::PlanarGraph) = euler_holds(DartGraph(g))

"""
    dual(g) -> PlanarGraph

Dual of the embedding: one vertex per face, one edge per primal edge.

A `PlanarGraph` is simple, so the dual must be too. An `ArgumentError` is thrown
when the graph is not connected or has no edge, when the embedding is not planar,
when the dual would have a self-loop (a bridge), when it would have parallel edges
(two faces sharing more than one edge, as around a vertex of degree two or at a
2-edge cut), or when the number of faces does not fit the vertex type. A
triangulation with at least four vertices always qualifies.
"""
function dual(g::PlanarGraph{T}) where {T}
    dg = DartGraph(g)
    orbits, face_of = face_orbits(dg)
    _check_dualizable(dg, length(orbits), T)
    mark = zeros(Int, length(orbits))
    rot = Vector{Vector{T}}()
    for (f, orbit) in enumerate(orbits)
        around = T[]
        for d in orbit
            h = face_of[dg.twin[d]]
            h == f && throw(ArgumentError("the dual has a self-loop (a bridge)"))
            mark[h] == f && throw(ArgumentError("the dual has parallel edges"))
            mark[h] = f
            push!(around, T(h))
        end
        push!(rot, around)
    end
    return PlanarGraph(rot)
end

function _check_dualizable(dg::DartGraph, nfaces::Int, ::Type{T}) where {T}
    ndarts(dg) > 0 || throw(ArgumentError("the dual needs a graph with an edge"))
    euler_holds(dg) ||
        throw(ArgumentError("the dual needs a planar embedding (Euler's formula fails)"))
    dg.n - nedges(dg) + nfaces == 2 ||
        throw(ArgumentError("the dual needs a connected graph"))
    nfaces <= typemax(T) ||
        throw(ArgumentError("$nfaces faces do not fit the vertex type $T"))
    return nothing
end

# --- conversions ---------------------------------------------------------

function DartGraph(g::PlanarGraph)
    dg = DartGraph(nv(g))
    for v = 1:dg.n
        anchor = 0
        for w in g.rotation[v]
            anchor = add_dart!(dg, v, Int(w), anchor)
        end
    end
    _match_pending_twins!(dg)
    return dg
end

# The embedding of `g` with the vertex type `T`; `g` itself when it has it already.
_retype(::Type{T}, g::PlanarGraph{T}) where {T<:Integer} = g

function _retype(::Type{T}, g::PlanarGraph) where {T<:Integer}
    graph = SimpleGraph{T}(nv(g))
    for e in edges(g.graph)
        Graphs.add_edge!(graph, T(src(e)), T(dst(e)))
    end
    return PlanarGraph(graph, g.rotation)
end

function PlanarGraph(dg::DartGraph)
    rot = rotation_lists(dg)
    for v = 1:dg.n
        length(unique(rot[v])) == length(rot[v]) || throw(
            ArgumentError(
                "embedded graph has parallel edges; it is not a simple PlanarGraph",
            ),
        )
    end
    return PlanarGraph(rot)
end
