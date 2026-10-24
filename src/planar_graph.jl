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

```julia
g = PlanarGraph([[2, 3], [3, 1], [1, 2]])   # triangle, counter-clockwise
nv(g), ne(g), nfaces(g)
```

Two embeddings are equal when they have the same underlying graph and their
rotations agree up to a cyclic shift and reversal at every vertex, that is
when they describe the same embedding of the same graph.
"""
struct PlanarGraph{T<:Integer}
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

function PlanarGraph(graph::SimpleGraph{T}, rotation::AbstractVector) where {T<:Integer}
    return PlanarGraph{T}(graph, _convert_rotation(rotation, T))
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
    length(rotation) >= 0 || throw(ArgumentError("rotation must be a vector of vectors"))
    return [T[x for x in r] for r in rotation]
end

function _check_rotation(graph::SimpleGraph, rot::Vector{Vector{T}}) where {T}
    n = nv(graph)
    length(rot) == n ||
        throw(ArgumentError("rotation has $(length(rot)) entries for $n vertices"))
    for v = 1:n
        expected = sort!(collect(neighbors(graph, v)))
        got = sort!(copy(rot[v]))
        length(got) == length(expected) || throw(
            ArgumentError(
                "rotation of vertex $v has $(length(got)) entries, expected $(length(expected))",
            ),
        )
        got == expected || throw(
            ArgumentError("rotation of vertex $v is not a permutation of its neighbours"),
        )
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

Base.hash(g::PlanarGraph, h::UInt) = hash(g.rotation, hash(:PlanarGraph, h))

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
function faces(g::PlanarGraph)
    dg = DartGraph(g)
    orbits, _ = face_orbits(dg)
    return [[tail(dg, d) for d in orbit] for orbit in orbits]
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
Requires a bridge-free primal graph (otherwise the dual has self-loops).
"""
function dual(g::PlanarGraph{T}) where {T}
    dg = DartGraph(g)
    orbits, face_of = face_orbits(dg)
    rot = Vector{Vector{T}}()
    for orbit in orbits
        push!(rot, T[face_of[dg.twin[d]] for d in orbit])
    end
    return PlanarGraph(rot)
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
