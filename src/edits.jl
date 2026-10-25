# Edits of a PlanarGraph that keep the embedding.
#
# Every function changes `g` in place and returns it. Vertex numbers never change,
# so a contracted vertex is left isolated instead of being renumbered away. The
# rotation lists are plain vectors, so a change costs `O(degree)` per vertex
# touched (and finding the entry of one vertex in the list of another is a scan).

function _require_edge(g::PlanarGraph, u::Integer, v::Integer)
    has_edge(g.graph, u, v) || throw(ArgumentError("$u-$v is not an edge of the graph"))
    return nothing
end

# Remove the first entry `x` of `r`.
function _remove_entry!(r::Vector, x)
    deleteat!(r, findfirst(==(x), r)::Int)
    return r
end

"""
    delete_edge!(g, u, v) -> g

Remove the edge `uv`, keeping the cyclic order of the remaining edges at every
vertex, so the result is again a planar embedding of the smaller graph. Throws an
`ArgumentError` when `uv` is not an edge. `Graphs.rem_edge!(g, u, v)` is the same
operation returning a `Bool`. The cost is `O(deg(u) + deg(v))`.
"""
function delete_edge!(g::PlanarGraph, u::Integer, v::Integer)
    _require_edge(g, u, v)
    _remove_entry!(g.rotation[u], v)
    _remove_entry!(g.rotation[v], u)
    Graphs.rem_edge!(g.graph, u, v)
    return g
end

function Graphs.rem_edge!(g::PlanarGraph, u::Integer, v::Integer)
    has_edge(g.graph, u, v) || return false
    delete_edge!(g, u, v)
    return true
end

"""
    contract_edge!(g, u, v) -> g

Contract the edge `uv` into `u`: the neighbours of `v` become neighbours of `u`, in
the position `v` held in the embedding, and `v` is left isolated (it keeps its
number). A neighbour of both `u` and `v` would get two parallel edges to `u`; one
of them is dropped, so the graph stays simple and the embedding stays planar (the
triangle `uvw` of a triangulation disappears: three edges fewer in all). Throws an
`ArgumentError` when `uv` is not an edge. The cost is `O(deg(u) + deg(v))` plus a
scan of the rotation of every neighbour of `v`.

```julia
g = planar_wheel(6)
contract_edge!(g, 1, 2)
is_planar_embedding(g)      # true
```
"""
function contract_edge!(g::PlanarGraph{T}, u::Integer, v::Integer) where {T}
    _require_edge(g, u, v)
    around_u = _after(g.rotation[u], v)
    around_v = _after(g.rotation[v], u)
    shared = Set{T}(w for w in around_v if has_edge(g.graph, u, w))
    for w in around_v
        _hand_over!(g, w, T(u), T(v), w in shared)
    end
    g.rotation[u] = vcat(around_u, T[w for w in around_v if !(w in shared)])
    g.rotation[v] = T[]
    Graphs.rem_edge!(g.graph, u, v)
    return g
end

# The entries of `r` that follow `x` cyclically, without `x`.
function _after(r::Vector, x)
    i = findfirst(==(x), r)::Int
    return vcat(r[(i+1):end], r[1:(i-1)])
end

# Neighbour `w` of `v` now sees `u` in place of `v`; when it saw `u` already, the
# entry of `v` is dropped.
function _hand_over!(g::PlanarGraph{T}, w::T, u::T, v::T, shared::Bool) where {T}
    w == u && return nothing
    r = g.rotation[w]
    shared ? _remove_entry!(r, v) : (r[findfirst(==(v), r)::Int] = u)
    Graphs.rem_edge!(g.graph, v, w)
    Graphs.add_edge!(g.graph, u, w)
    return nothing
end

"""
    add_edge_in_face!(g, u, v; after_u = nothing, after_v = nothing) -> g

Add the edge `uv` inside a face. At `u` the new edge goes right after the neighbour
`after_u` in the counter-clockwise rotation, that is into the corner between
`after_u` and the next neighbour; `after_v` does the same at `v`. A vertex of degree
0 needs no corner and one of degree 1 has only one; a larger degree requires it.
The two corners must lie on one face, which is checked by walking that face (the
cost is the sum of the degrees along it). An endpoint of degree 0 can go into any
face. Edges between two components are not supported: build such a graph with
`PlanarGraph(rotation)`.

Throws an `ArgumentError` for `u == v`, an existing edge, an unknown vertex, a
corner that is missing, wrong or on another face.

```julia
g = planar_cycle(6)
add_edge_in_face!(g, 1, 4; after_u = 2, after_v = 5)    # a chord of the inner face
```
"""
function add_edge_in_face!(
    g::PlanarGraph{T},
    u::Integer,
    v::Integer;
    after_u = nothing,
    after_v = nothing,
) where {T}
    _check_new_edge(g, u, v)
    pu = _corner(g, u, after_u)
    pv = _corner(g, v, after_v)
    _share_face(g, (u, pu), (v, pv)) ||
        throw(ArgumentError("the corners of $u and $v are not on one face"))
    insert!(g.rotation[u], pu + 1, T(v))
    insert!(g.rotation[v], pv + 1, T(u))
    Graphs.add_edge!(g.graph, u, v)
    return g
end

function _check_new_edge(g::PlanarGraph, u::Integer, v::Integer)
    (Graphs.has_vertex(g.graph, u) && Graphs.has_vertex(g.graph, v)) ||
        throw(ArgumentError("$u or $v is not a vertex of the graph"))
    u == v && throw(ArgumentError("a self-loop cannot be added"))
    has_edge(g.graph, u, v) && throw(ArgumentError("$u-$v is already an edge"))
    return nothing
end

# Position in the rotation of `u` after which a new edge goes: 0 when `u` has no
# edge.
function _corner(g::PlanarGraph, u::Integer, after)
    r = g.rotation[u]
    isempty(r) && return after === nothing ? 0 : throw(_bad_corner(u, after))
    after === nothing && length(r) == 1 && return 1
    after === nothing &&
        throw(ArgumentError("vertex $u has $(length(r)) edges: give after_u"))
    i = findfirst(==(after), r)
    i === nothing && throw(_bad_corner(u, after))
    return i
end

_bad_corner(u, after) = ArgumentError("$after is not a neighbour of $u")

# The face of the corner after position `p` at `u` contains the dart from `u` to
# `r[p]`; the face to the left of a dart `x -> y` goes on along `y -> pred(x)`.
function _share_face(g::PlanarGraph, (u, pu), (v, pv))
    (pu == 0 || pv == 0) && return true
    first = (u, g.rotation[u][pu])
    target = (v, g.rotation[v][pv])
    dart = first
    while true
        dart == target && return true
        dart = _next_on_face(g, dart)
        dart == first && return false
    end
end

function _next_on_face(g::PlanarGraph, (x, y))
    r = g.rotation[y]
    i = findfirst(==(x), r)::Int
    return (y, r[mod1(i - 1, length(r))])
end
