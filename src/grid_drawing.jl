# Straight-line drawings on an integer grid in O(n) × O(n) area.
#
# The placement is the one of de Fraysseix, Pach and Pollack, in the formulation
# used by OGDF (`src/ogdf/planarlayout/FPPLayout.cpp`, GPL-2.0 or GPL-3.0): the
# graph is augmented
# to a maximal planar graph by inserting chords into the faces of its embedding,
# a shelling order of that triangulation is computed by peeling vertices off the
# outer face, and the vertices are then placed level by level, each new vertex
# at the apex of a triangle over an edge of the current outer face. Every
# coordinate is an integer, the width is at most `2(n-2)` and the height at most
# `n-2`, so the total area is quadratic in the number of vertices.
#
# Disconnected graphs are drawn component by component and the components are
# packed side by side, which keeps the linear bound on both dimensions.

"""
    GridDrawing

Straight-line drawing of a planar graph on the integer grid. `x[v]` and `y[v]`
are the coordinates of vertex `v`, and `width`/`height` are the dimensions of
the bounding box, so that every vertex lies in `0:width × 0:height`.

Use [`coordinates`](@ref) for the labelled positions and `d[v]` for a single
vertex.
"""
struct GridDrawing
    x::Vector{Int}
    y::Vector{Int}
    width::Int
    height::Int
end

"""Number of vertices of a [`GridDrawing`](@ref)."""
Base.length(d::GridDrawing) = length(d.x)

"""Position of vertex `v` in a [`GridDrawing`](@ref)."""
Base.getindex(d::GridDrawing, v::Integer) = (d.x[v], d.y[v])

"""Width of the bounding box of a [`GridDrawing`](@ref)."""
width(d::GridDrawing) = d.width

"""Height of the bounding box of a [`GridDrawing`](@ref)."""
height(d::GridDrawing) = d.height

"""    coordinates(d::GridDrawing) -> Vector{Tuple{Int,Int}}"""
coordinates(d::GridDrawing) = [(d.x[v], d.y[v]) for v in eachindex(d.x)]

function Base.show(io::IO, d::GridDrawing)
    print(io, "GridDrawing(", length(d), " vertices, ", d.width, "x", d.height, " grid)")
    return nothing
end

"""
    grid_drawing(g::PlanarGraph) -> GridDrawing
    grid_drawing(g::SimpleGraph) -> GridDrawing

Straight-line crossing-free drawing of the planar graph `g` on the integer grid,
with integer coordinates in `O(n) × O(n)`.

The embedding of `g` is used as given; it is augmented to a maximal planar graph
by inserting chords inside faces (the added edges are not drawn, they only guide
the placement), and the resulting triangulation is drawn by the
Fraysseix–Pach–Pollack algorithm. Disconnected graphs are drawn component by
component and packed side by side.

For a `SimpleGraph` the embedding is computed first with
[`planar_embedding`](@ref), and for a `PlanarGraph` whose stored rotation is not
a planar embedding it is recomputed; an `ArgumentError` is thrown only when the
graph itself is not planar.

```julia
g = planar_grid(3, 3)
d = grid_drawing(g)
width(d), height(d)
d[1]
```
"""
function grid_drawing(g::PlanarGraph)
    return _grid_drawing(DartGraph(_drawable_embedding(g)))
end

function grid_drawing(g::SimpleGraph)
    return _grid_drawing(DartGraph(_drawable_embedding(g)))
end

# The embedding to draw, repaired if the stored one is not an embedding.
function _drawable_embedding(g)
    embedded = planar_embedding(g)
    embedded === nothing && throw(ArgumentError("the graph is not planar"))
    return embedded
end

function _grid_drawing(dg::DartGraph)
    n = dg.n
    x = zeros(Int, n)
    y = zeros(Int, n)
    n == 0 && return GridDrawing(x, y, 0, 0)

    label, ncomps = component_labels(dg)
    offset = 0
    for c = 1:ncomps
        keep = [label[v] == c for v = 1:n]
        sub, ids = induced(dg, keep)
        if length(ids) <= 2
            offset = _place_trivial_component!(x, y, ids, offset)
            continue
        end
        placed = _triangulate_component(sub)
        cx, cy = _fpp_coordinates(placed)
        minx, miny = minimum(cx), minimum(cy)
        for (i, v) in enumerate(ids)
            x[v] = cx[i] - minx + offset
            y[v] = cy[i] - miny
        end
        offset += maximum(cx) - minx + 2
    end

    return GridDrawing(x, y, maximum(x), maximum(y))
end

# One vertex, or two vertices joined by an edge: a single column, or two.
function _place_trivial_component!(x, y, ids::Vector{Int}, offset::Int)
    for (i, v) in enumerate(ids)
        x[v] = offset + i - 1
        y[v] = 0
    end
    return offset + length(ids) + 1
end

# Chords are inserted until every face is a triangle, which is what the
# shelling order and the placement need.
function _triangulate_component(h::DartGraph)
    placed = _add_chords!(copy(h))
    orbits, _ = face_orbits(placed)
    all(orbit -> length(orbit) == 3, orbits) || throw(
        ArgumentError("cannot draw on a grid: the embedding could not be triangulated"),
    )
    return placed
end

# --- Fraysseix–Pach–Pollack ----------------------------------------------

"""
    _fpp_coordinates(h::DartGraph) -> (Vector{Int}, Vector{Int})

Coordinates of a triangulated connected embedded graph with at least three
vertices, given as the OGDF implementation computes them.
"""
function _fpp_coordinates(h::DartGraph)
    n = h.n
    orbits, _ = face_orbits(h)
    e12 = orbits[1][1]
    e2n = face_successor(h, e12)
    en1 = face_successor(h, e2n)
    v1, v2, vn = tail(h, e12), tail(h, e2n), tail(h, en1)
    length(unique((v1, v2, vn))) == 3 ||
        throw(ArgumentError("cannot draw on a grid: the outer face is not a triangle"))

    sh = _shelling_order(h, e12, e2n, en1)
    return _fpp_place(h, sh.num, sh.e_wp, sh.e_wq)
end

# A linked list of vertices that may be peeled off next. The shelling order
# needs deletion of an arbitrary vertex and popping from the front, both in
# constant time.
mutable struct _NodeList
    next::Vector{Int}
    prev::Vector{Int}
    present::Vector{Bool}
    head::Int
    tail::Int
end

_NodeList(n::Integer) = _NodeList(zeros(Int, n), zeros(Int, n), falses(n), 0, 0)

function _push_back!(list::_NodeList, v::Int)
    list.prev[v] = list.tail
    list.next[v] = 0
    list.tail == 0 ? (list.head = v) : (list.next[list.tail] = v)
    list.tail = v
    list.present[v] = true
    return list
end

function _delete!(list::_NodeList, v::Int)
    list.present[v] || throw(ArgumentError("vertex $v is not in the candidate list"))
    p, q = list.prev[v], list.next[v]
    p == 0 ? (list.head = q) : (list.next[p] = q)
    q == 0 ? (list.tail = p) : (list.prev[q] = p)
    list.present[v] = false
    return list
end

function _pop_front!(list::_NodeList)
    v = list.head
    v == 0 && throw(ArgumentError("the candidate list is empty"))
    _delete!(list, v)
    return v
end


# The state of the peeling: the order as far as it is known, for every vertex
# the number of chords that still keep it from being removed, whether it is on
# the current outer face, the two outer-face darts that will place it, and the
# queue of vertices that may be removed next.
mutable struct _Shelling
    num::Vector{Int}
    num_diag::Vector{Int}
    outer::BitVector
    e_wp::Vector{Int}
    e_wq::Vector{Int}
    possible::_NodeList
    v1::Int
    v2::Int
end

function _Shelling(n::Integer, v1::Int, v2::Int)
    return _Shelling(
        zeros(Int, n),
        zeros(Int, n),
        falses(n),
        zeros(Int, n),
        zeros(Int, n),
        _NodeList(n),
        v1,
        v2,
    )
end
# Peel the outer face down to the initial triangle, numbering the vertices in
# reverse removal order; `num` is the shelling order, and `e_wp`/`e_wq` are the
# darts at every vertex pointing to its predecessor and successor on the outer
# face at the time it is removed.
function _shelling_order(h::DartGraph, e12::Int, e2n::Int, en1::Int)
    n = h.n
    v1, v2, vn = tail(h, e12), tail(h, e2n), tail(h, en1)
    sh = _Shelling(n, v1, v2)

    sh.num[v1] = 1
    sh.num[v2] = 2
    sh.outer[v1] = true
    sh.outer[v2] = true
    sh.outer[vn] = true
    _push_back!(sh.possible, vn)

    sh.e_wq[v1] = h.twin[en1]
    sh.e_wp[v2] = e2n
    sh.e_wq[vn] = h.twin[e2n]
    sh.e_wp[vn] = en1

    for k = n:-1:3
        vk = _pop_front!(sh.possible)
        sh.num[vk] = k
        wq = h.head[sh.e_wq[vk]]
        wp = h.head[sh.e_wp[vk]]
        sh.outer[vk] = false

        _drop_chord!(h, sh, wp, wq)
        sh.e_wq[wp] = h.next_rot[sh.e_wq[wp]]
        sh.e_wp[wq] = h.prev_rot[sh.e_wp[wq]]
        _expose_outer_face!(h, sh, wp, wq)
    end
    return sh
end

# `wp` and `wq` lose a chord when they are still joined by one.
function _drop_chord!(h::DartGraph, sh::_Shelling, wp::Int, wq::Int)
    h.head[h.next_rot[sh.e_wq[wp]]] == wq || return nothing
    _release_candidate!(sh, wp)
    _release_candidate!(sh, wq)
    return nothing
end

function _release_candidate!(sh::_Shelling, v::Int)
    sh.num_diag[v] -= 1
    sh.num_diag[v] == 0 && _push_back!(sh.possible, v)
    return nothing
end

# The vertices between `wp` and `wq` join the outer face; the chords they bring
# are the constraints that keep them from being peeled off.
function _expose_outer_face!(h::DartGraph, sh::_Shelling, wp::Int, wq::Int)
    e = sh.e_wq[wp]
    u = h.head[e]
    while u != wq
        e = _expose_outer_vertex!(h, sh, u, e)
        u = h.head[e]
    end
    return nothing
end

function _expose_outer_vertex!(h::DartGraph, sh::_Shelling, u::Int, e::Int)
    sh.outer[u] = true
    sh.e_wp[u] = h.twin[e]
    sh.e_wq[u] = h.next_rot[h.next_rot[sh.e_wp[u]]]
    e = sh.e_wq[u]
    _count_new_chords!(h, sh, u)
    sh.num_diag[u] == 0 && _push_back!(sh.possible, u)
    return e
end

function _count_new_chords!(h::DartGraph, sh::_Shelling, u::Int)
    e2 = h.prev_rot[sh.e_wp[u]]
    while e2 != sh.e_wq[u]
        w = h.head[e2]
        sh.outer[w] && _count_chord!(sh, u, w)
        e2 = h.prev_rot[e2]
    end
    return nothing
end

# A chord to an outer vertex counts at both ends, except along the root edge.
function _count_chord!(sh::_Shelling, u::Int, w::Int)
    sh.num_diag[u] += 1
    (w == sh.v1 || w == sh.v2) && return nothing
    sh.num_diag[w] += 1
    sh.num_diag[w] == 1 && _delete!(sh.possible, w)
    return nothing
end

# Place the vertices in shelling order: every vertex goes above the edge of the
# current outer face that it replaces, on the first odd level that keeps all
# distances even, and the vertices in between are shifted apart.
function _fpp_place(h::DartGraph, num::Vector{Int}, e_wp::Vector{Int}, e_wq::Vector{Int})
    n = h.n
    x = zeros(Int, n)
    y = zeros(Int, n)
    x_rel = zeros(Int, n)
    upper = zeros(Int, n)
    nxt = zeros(Int, n)
    order = zeros(Int, n)
    for v = 1:n
        num[v] == 0 &&
            throw(ArgumentError("cannot draw on a grid: the shelling order is incomplete"))
        order[num[v]] = v
    end

    v1, v2 = order[1], order[2]
    x_rel[v1] = 0
    x_rel[v2] = 0
    y[v1] = 0
    y[v2] = 0
    nxt[v1] = v2
    nxt[v2] = 0

    for k = 3:n
        vk = order[k]
        wp = h.head[e_wp[vk]]
        wq = h.head[e_wq[vk]]

        xq = 2
        w = wp
        while true
            w = nxt[w]
            xq += x_rel[w]
            w == wq && break
        end

        x_rel[vk] = div(xq + y[wq] - y[wp], 2)
        y[vk] = div(xq + y[wq] + y[wp], 2)
        x_rel[wq] = xq - x_rel[vk]

        dx = 1
        w = nxt[wp]
        while w != wq
            dx += x_rel[w]
            x[w] = dx - x_rel[vk]
            upper[w] = vk
            w = nxt[w]
        end

        nxt[wp] = vk
        nxt[vk] = wq
    end

    vn = order[n]
    x[vn] = x_rel[vn]
    x[v2] = x[vn] + x_rel[v2]
    x[v1] = 0
    for k = (n-1):-1:3
        x[order[k]] += x[upper[order[k]]]
    end
    return x, y
end
