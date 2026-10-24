# Triangulation of an embedded graph by inserting chords inside faces.
#
# Two linear passes over each face, both going around its boundary and clipping
# ears: the chord `ac` between the neighbours of a boundary position `b` is
# added and `b` leaves the boundary.
#
#   * A boundary that repeats a vertex (a bridge or a cut vertex) first loses its
#     repeated occurrences. Such an occurrence can always be clipped: the vertex
#     appears twice, so `a` and `c` lie in different blocks around it and are
#     not adjacent. After a clip the next occurrence is skipped, so the chords
#     are spread over the vertices instead of chaining through one of them.
#   * The remaining simple cycle is clipped at every position whose neighbours
#     `a` and `c` are not adjacent yet; otherwise the next position is tried.
#     Among four consecutive vertices `a, b, c, d` the chord `ac` or `bd` is
#     always new, since two edges outside the face joining `a, c` and `b, d`
#     would cross, so a failure is followed by a success. After a clip the next
#     vertex is skipped, so every pass removes a constant fraction of the face.
#
# The adjacency test scans the smaller of the two degrees. A test that succeeds
# is charged to the chord it adds, one that fails to the edge it finds, and the
# sum over the edges of the smaller degree is `O(n)` in a planar graph (Chiba and
# Nishizeki), so the whole triangulation is `O(n + m)`, with no hashing, and never
# adds an edge that is already there.

"""
    triangulate(g) -> same type as g

Add edges to `g` until every face is a triangle, preserving the embedding and
the vertex set, in `O(nv + ne)` time. No edge is added twice, so a simple graph
stays simple, and the new edges are spread over the vertices of a face instead of
fanning from one of them. The result is a maximal planar graph when `g` is
connected with at least three vertices.

A face boundary that repeats a vertex (a bridge or a cut vertex) is handled. An
`ArgumentError` is raised when the rotation system does not describe a planar
embedding, or when a component has fewer than three vertices but an edge (no
chord can make its face a triangle). For a [`PlanarGraph`](@ref) the result is
returned as a `PlanarGraph`.
"""
function triangulate(g::DartGraph)
    h = _add_chords!(copy(g))
    _verify_triangulation(h)
    return h
end

triangulate(g::PlanarGraph) = PlanarGraph(triangulate(DartGraph(g)))

function _add_chords!(h::DartGraph)
    euler_holds(h) || throw(
        ArgumentError(
            "cannot triangulate: the rotation system does not describe a planar embedding",
        ),
    )
    orbits, _ = face_orbits(h)
    occurrences = zeros(Int, h.n)
    for orbit in orbits
        length(orbit) <= 3 && continue
        orbit = _peel_repeats!(h, orbit, occurrences)
        length(orbit) <= 3 && continue
        _clip_ears!(h, orbit)
    end
    return h
end

# --- faces that repeat a vertex ------------------------------------------

# Clip the occurrences of repeated vertices off the boundary walk `orbit`, going
# around the list of those occurrences and skipping the one next to a clip.
# `occurrences[v]` counts the occurrences of `v` that are left; it is all zero
# again on return. Return the remaining boundary walk, which is a simple cycle.
function _peel_repeats!(h::DartGraph, orbit::Vector{Int}, occurrences::Vector{Int})
    verts = Int[tail(h, d) for d in orbit]
    for v in verts
        occurrences[v] += 1
    end
    candidates = findall(t -> occurrences[verts[t]] >= 2, eachindex(verts))
    result =
        isempty(candidates) ? orbit :
        _clip_repeats!(h, orbit, verts, occurrences, candidates)
    for v in verts
        occurrences[v] = 0
    end
    return result
end

function _clip_repeats!(
    h::DartGraph,
    orbit::Vector{Int},
    verts::Vector{Int},
    occurrences::Vector{Int},
    candidates::Vector{Int},
)
    k = length(orbit)
    prevp = [mod1(t - 1, k) for t = 1:k]
    nextp = [mod1(t + 1, k) for t = 1:k]
    dartp = copy(orbit)
    removed = falses(k)
    ring = _Ring(k, candidates)
    cur = candidates[1]
    while ring.size > 0
        t = cur
        cur = _ring_next(ring, t)
        _ring_remove!(ring, t)
        occurrences[verts[t]] >= 2 || continue
        b = prevp[t]
        c = nextp[t]
        verts[b] == verts[c] && throw(_not_triangulable())
        dartp[b] = _add_chord!(h, verts[b], verts[c], dartp[b], dartp[c])
        nextp[b] = c
        prevp[c] = b
        removed[t] = true
        occurrences[verts[t]] -= 1
        cur == c && ring.size > 0 && (cur = _ring_next(ring, c))
    end
    return _remaining_walk(dartp, nextp, removed)
end

# A circular doubly linked list over some of the positions `1..k`.
mutable struct _Ring
    next::Vector{Int}
    prev::Vector{Int}
    size::Int
end

function _Ring(k::Int, members::Vector{Int})
    ring = _Ring(zeros(Int, k), zeros(Int, k), length(members))
    for (i, t) in enumerate(members)
        ring.next[t] = members[mod1(i + 1, length(members))]
        ring.prev[t] = members[mod1(i - 1, length(members))]
    end
    return ring
end

_ring_next(ring::_Ring, t::Int) = ring.next[t]

function _ring_remove!(ring::_Ring, t::Int)
    ring.next[ring.prev[t]] = ring.next[t]
    ring.prev[ring.next[t]] = ring.prev[t]
    ring.size -= 1
    return nothing
end

function _remaining_walk(dartp::Vector{Int}, nextp::Vector{Int}, removed::BitVector)
    start = findfirst(!, removed)
    out = Int[]
    t = start
    while true
        push!(out, dartp[t])
        t = nextp[t]
        t == start && break
    end
    return out
end

_not_triangulable() =
    ArgumentError("cannot triangulate: a face boundary cannot be made a simple cycle")

# Insert the chord `u v` inside a face: the new dart at `u` goes after `du`, the
# dart leaving `u` along the face, and the one at `v` after `dv`. Returns the
# dart `u -> v`.
function _add_chord!(h::DartGraph, u::Int, v::Int, du::Int, dv::Int)
    forward = add_dart!(h, u, v, du)
    backward = add_dart!(h, v, u, dv)
    h.twin[forward] = backward
    h.twin[backward] = forward
    return forward
end

# --- simple faces --------------------------------------------------------

# Clip the simple face with the darts `orbit` down to a triangle, going around
# it. A test that failed is remembered as the pair of neighbours it was made
# with, so it is not repeated until a neighbour changes. Each pass removes at
# least a third of the vertices, so there are at most about `3k` visits; the
# bound only guards against an input that is not a planar embedding.
function _clip_ears!(h::DartGraph, orbit::Vector{Int})
    k = length(orbit)
    verts = Int[tail(h, d) for d in orbit]
    dart = copy(orbit)
    prevp = [mod1(t - 1, k) for t = 1:k]
    nextp = [mod1(t + 1, k) for t = 1:k]
    failed = zeros(Int, k, 2)
    remaining = k
    cur = 1
    for _ = 1:(20k+20)
        remaining > 3 || return nothing
        a = prevp[cur]
        c = nextp[cur]
        if _fails(h, failed, cur, verts[a], verts[c])
            cur = c
        else
            dart[a] = _add_chord!(h, verts[a], verts[c], dart[a], dart[c])
            nextp[a] = c
            prevp[c] = a
            remaining -= 1
            cur = nextp[c]
        end
    end
    throw(_not_triangulable())
end

# True when the chord between the vertices `u` and `w` cannot be added, because
# they are adjacent; the answer is kept for `b` while its neighbours stay.
function _fails(h::DartGraph, failed::Matrix{Int}, b::Int, u::Int, w::Int)
    (failed[b, 1] == u && failed[b, 2] == w) && return true
    _adjacent(h, u, w) || return false
    failed[b, 1] = u
    failed[b, 2] = w
    return true
end

# Adjacent vertices, found by scanning the one of smaller degree.
function _adjacent(h::DartGraph, u::Int, w::Int)
    h.deg[u] <= h.deg[w] || return _adjacent(h, w, u)
    for d in incident_darts(h, u)
        h.head[d] == w && return true
    end
    return false
end

# --- verification --------------------------------------------------------

# The result must be a triangulation of a valid embedding.
function _verify_triangulation(h::DartGraph)
    orbits, _ = face_orbits(h)
    all(orbit -> length(orbit) == 3, orbits) ||
        throw(ArgumentError("cannot triangulate: some face is not a simple cycle"))
    euler_holds(h) ||
        throw(ArgumentError("cannot triangulate: the rotation system is not planar"))
    return nothing
end
