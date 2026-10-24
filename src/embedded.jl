# Internal dart-based representation of an embedded multigraph.
#
# A dart is an oriented edge. Every edge contributes two darts, `d` and its
# twin `twin[d]`, which point in opposite directions. Around every vertex the
# darts leaving it are kept in a cyclic doubly linked list (`next_rot`,
# `prev_rot`), which is the rotation system of the embedding.
#
# The rotation order is read counter-clockwise. The face to the left of a dart
# is therefore traced by
#
#     face_successor(d) = prev_rot[twin[d]]
#
# See `PlanarGraph` for the public, `Graphs.jl`-backed type.

"""
    DartGraph(n)

Mutable embedded multigraph with `n` vertices, stored as a rotation system
over darts. Vertices are numbered `1:n`; darts are numbered `1:2m`.

This is the working representation of the algorithms in this package. It
accepts parallel edges, which the public [`PlanarGraph`](@ref) does not.
"""
mutable struct DartGraph
    n::Int
    deg::Vector{Int}
    first::Vector{Int}
    tail::Vector{Int}
    head::Vector{Int}
    twin::Vector{Int}
    next_rot::Vector{Int}
    prev_rot::Vector{Int}
end

function DartGraph(n::Integer)
    n >= 0 || throw(ArgumentError("number of vertices must be non-negative, got $n"))
    return DartGraph(
        Int(n),
        zeros(Int, n),
        zeros(Int, n),
        Int[],
        Int[],
        Int[],
        Int[],
        Int[],
    )
end

Base.copy(g::DartGraph) = DartGraph(
    g.n,
    copy(g.deg),
    copy(g.first),
    copy(g.tail),
    copy(g.head),
    copy(g.twin),
    copy(g.next_rot),
    copy(g.prev_rot),
)

ndarts(g::DartGraph) = length(g.head)

nedges(g::DartGraph) = length(g.head) >>> 1

tail(g::DartGraph, d::Integer) = g.tail[d]

function _insert_after!(g::DartGraph, d::Integer, anchor::Integer, v::Integer)
    if anchor == 0
        g.first[v] == 0 || throw(ArgumentError("anchor required for non-empty rotation"))
        g.first[v] = d
        g.next_rot[d] = d
        g.prev_rot[d] = d
        return g
    end
    g.tail[anchor] == v || throw(ArgumentError("dart $anchor does not leave vertex $v"))
    nxt = g.next_rot[anchor]
    g.next_rot[anchor] = d
    g.prev_rot[d] = anchor
    g.next_rot[d] = nxt
    g.prev_rot[nxt] = d
    return g
end

"""
    add_dart!(g, u, v, anchor=0) -> Int

Insert the dart `u -> v` into the rotation of `u`, immediately after `anchor`
(or as the first dart when `anchor == 0`). The twin is added separately.
"""
function add_dart!(g::DartGraph, u::Integer, v::Integer, anchor::Integer = 0)
    (1 <= u <= g.n && 1 <= v <= g.n) ||
        throw(ArgumentError("vertices out of range: $u, $v"))
    u == v && throw(ArgumentError("self-loops are not supported"))
    d = length(g.head) + 1
    push!(g.tail, Int(u))
    push!(g.head, Int(v))
    push!(g.twin, 0)
    push!(g.next_rot, d)
    push!(g.prev_rot, d)
    _insert_after!(g, d, anchor, Int(u))
    g.deg[u] += 1
    return d
end

"""
    add_edge!(g, u, v, anchor_u=0, anchor_v=0) -> Int

Add an edge `{u, v}` and return its first dart `u -> v`.
"""
function add_edge!(
    g::DartGraph,
    u::Integer,
    v::Integer,
    anchor_u::Integer = 0,
    anchor_v::Integer = 0,
)
    d = add_dart!(g, u, v, anchor_u)
    e = add_dart!(g, v, u, anchor_v)
    g.twin[d] = e
    g.twin[e] = d
    return d
end

"""Rotations of `g` as vertex neighbours, in counter-clockwise order."""
function rotation_lists(g::DartGraph)
    rots = [Int[] for _ = 1:g.n]
    for v = 1:g.n
        d = g.first[v]
        d == 0 && continue
        start = d
        while true
            push!(rots[v], g.head[d])
            d = g.next_rot[d]
            d == start && break
        end
    end
    return rots
end

"""Darts leaving `v`, in counter-clockwise order."""
function darts_of(g::DartGraph, v::Integer)
    d = g.first[v]
    d == 0 && return Int[]
    out = Int[d]
    x = g.next_rot[d]
    while x != d
        push!(out, x)
        x = g.next_rot[x]
    end
    return out
end

face_successor(g::DartGraph, d::Integer) = g.prev_rot[g.twin[d]]

"""
    face_orbits(g) -> (orbits, face_of_dart)

Decompose the darts of `g` into face boundary walks.
"""
function face_orbits(g::DartGraph)
    face_of = zeros(Int, ndarts(g))
    orbits = Vector{Vector{Int}}()
    for d = 1:ndarts(g)
        face_of[d] == 0 || continue
        orbit = Int[]
        x = d
        while face_of[x] == 0
            face_of[x] = length(orbits) + 1
            push!(orbit, x)
            x = face_successor(g, x)
        end
        push!(orbits, orbit)
    end
    return orbits, face_of
end

"""Number of connected components of `g` (isolated vertices counted)."""
function count_components(g::DartGraph)
    seen = falses(g.n)
    total = 0
    for s = 1:g.n
        seen[s] && continue
        total += 1
        stack = Int[s]
        seen[s] = true
        while !isempty(stack)
            v = pop!(stack)
            for d in darts_of(g, v)
                w = g.head[d]
                if !seen[w]
                    seen[w] = true
                    push!(stack, w)
                end
            end
        end
    end
    return total
end

"""
    bfs_tree(g, roots; keep=nothing) -> (parent, parent_dart, depth, order)

Breadth-first traversal seeded with `roots`. Vertices outside `keep` are not
visited. `parent[v] == 0` for the roots and for unreached vertices. `order`
lists the visited vertices in BFS order.
"""
function bfs_tree(g::DartGraph, roots; keep = nothing)
    parent = zeros(Int, g.n)
    parent_dart = zeros(Int, g.n)
    depth = fill(-1, g.n)
    order = Int[]
    queue = Int[]
    for r in roots
        depth[r] == -1 || continue
        depth[r] = 0
        push!(queue, r)
        push!(order, r)
    end
    head = 1
    while head <= length(queue)
        v = queue[head]
        head += 1
        for d in darts_of(g, v)
            w = g.head[d]
            depth[w] == -1 || continue
            keep === nothing || keep[w] || continue
            depth[w] = depth[v] + 1
            parent[w] = v
            parent_dart[w] = d
            push!(queue, w)
            push!(order, w)
        end
    end
    return parent, parent_dart, depth, order
end

"""
    LevelStructure

What the level-based separator needs from a breadth-first search: `parent[v]`
(0 at the roots and at unreached vertices), `depth[v]` (the BFS level) and
`order` (the reached vertices in BFS order).
"""
struct LevelStructure
    parent::Vector{Int}
    depth::Vector{Int}
    order::Vector{Int}
end

"""
    level_structure(g, roots; keep = nothing) -> LevelStructure

Breadth-first search seeded with `roots`, discarding everything a
[`LevelStructure`](@ref) does not need.
"""
function level_structure(g::DartGraph, roots; keep = nothing)
    parent, _, depth, order = bfs_tree(g, roots; keep)
    return LevelStructure(parent, depth, order)
end

"""
    restrict(g, keep) -> DartGraph

Induced subgraph on the vertices selected by `keep`, with vertex numbering
unchanged and the rotation system inherited from `g`.
"""
function restrict(g::DartGraph, keep::AbstractVector{Bool})
    h = DartGraph(g.n)
    for v = 1:g.n
        keep[v] || continue
        anchor = 0
        for d in darts_of(g, v)
            keep[g.head[d]] || continue
            anchor = add_dart!(h, v, g.head[d], anchor)
        end
    end
    _match_pending_twins!(h)
    return h
end

"""
    induced(g, keep) -> (DartGraph, Vector{Int})

Subgraph induced by `keep`, renumbered to `1:count(keep)`. The second value
maps each new vertex back to its index in `g`.
"""
function induced(g::DartGraph, keep::AbstractVector{Bool})
    ids = Int[v for v = 1:g.n if keep[v]]
    index = zeros(Int, g.n)
    for (i, v) in enumerate(ids)
        index[v] = i
    end
    h = DartGraph(length(ids))
    for (i, v) in enumerate(ids)
        anchor = 0
        for d in darts_of(g, v)
            w = g.head[d]
            keep[w] || continue
            anchor = add_dart!(h, i, index[w], anchor)
        end
    end
    _match_pending_twins!(h)
    return h, ids
end

# Pair up darts of `h` that were created without twins. Rotations built by
# `DartGraph`, `restrict` and `induced` come from simple graphs, so each
# ordered pair of ends occurs once.
function _match_pending_twins!(h::DartGraph)
    base = h.n + 1
    pending = Dict{Int,Int}()
    sizehint!(pending, ndarts(h) >>> 1)
    for d = 1:ndarts(h)
        key = tail(h, d) * base + h.head[d]
        haskey(pending, key) && throw(ArgumentError("parallel edges are not supported"))
        pending[key] = d
    end
    for d = 1:ndarts(h)
        h.twin[d] == 0 || continue
        e = get(pending, h.head[d] * base + tail(h, d), 0)
        e == 0 && throw(ArgumentError("inconsistent rotation system"))
        e == d && throw(ArgumentError("self-loops are not supported"))
        h.twin[d] = e
        h.twin[e] = d
    end
    return h
end

"""
    euler_holds(g) -> Bool

`true` when every connected component satisfies Euler's formula
`n_c - m_c + f_c = 2`, which is what makes a rotation system a planar
embedding. A component without edges is a lone vertex: it has one face and no
face walk.
"""
function euler_holds(g::DartGraph)
    label, ncomps = component_labels(g)
    vertex_count = zeros(Int, ncomps)
    edge_count = zeros(Int, ncomps)
    for v = 1:g.n
        vertex_count[label[v]] += 1
    end
    for d = 1:ndarts(g)
        d < g.twin[d] || continue
        edge_count[label[tail(g, d)]] += 1
    end
    face_count = zeros(Int, ncomps)
    orbits, _ = face_orbits(g)
    for orbit in orbits
        face_count[label[tail(g, orbit[1])]] += 1
    end
    return all(c -> _euler_component(vertex_count, edge_count, face_count, c), 1:ncomps)
end

_euler_component(n, m, f, c::Integer) = n[c] - m[c] + (m[c] == 0 ? 1 : f[c]) == 2

"""Map each vertex of `g` to its component id."""
function component_labels(g::DartGraph)
    label = zeros(Int, g.n)
    current = 0
    for s = 1:g.n
        label[s] == 0 || continue
        current += 1
        stack = Int[s]
        label[s] = current
        while !isempty(stack)
            v = pop!(stack)
            for d in darts_of(g, v)
                w = g.head[d]
                if label[w] == 0
                    label[w] = current
                    push!(stack, w)
                end
            end
        end
    end
    return label, current
end

"""Connected components of `g` after deleting `removed` (isolated vertices kept)."""
function components_without(g::DartGraph, removed::AbstractVector{Bool})
    seen = falses(g.n)
    comps = Vector{Vector{Int}}()
    for s = 1:g.n
        (removed[s] || seen[s]) && continue
        comp = Int[]
        stack = Int[s]
        seen[s] = true
        while !isempty(stack)
            v = pop!(stack)
            push!(comp, v)
            for d in darts_of(g, v)
                w = g.head[d]
                (seen[w] || removed[w]) && continue
                seen[w] = true
                push!(stack, w)
            end
        end
        push!(comps, comp)
    end
    return comps
end
