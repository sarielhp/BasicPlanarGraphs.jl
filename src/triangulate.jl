# Triangulation of an embedded graph by inserting chords inside faces.
#
# A face boundary that repeats a vertex (a bridge or a cut vertex) first loses its
# repeated occurrences: such an occurrence can always be clipped, since the vertex
# appears twice, so its neighbours `a` and `c` lie in different blocks around it
# and are not adjacent. After a clip the next occurrence is skipped, so the chords
# are spread over the vertices instead of chaining through one of them.
#
# The remaining simple cycle `p_1 … p_k` is cut by a ladder: `p_i` is joined to
# `p_{k+1-i}`. The rungs are nested, so they cut the face into quadrangles, and
# each quadrangle `a, b, c, d` takes the diagonal `ac` or `bd`: both outside the
# face would cross, so one is always new. Every vertex gets at most one rung and
# two diagonals, a constant degree. A rung whose ends are already adjacent
# outside the face is skipped, which merges its quadrangles into a larger
# region. All rungs of such a region are edges outside the face, and each rung of
# a ladder started one position later (a slanted ladder) crosses one of them, so
# none is blocked: the slanted ladder cuts the region into quadrangles, and a
# vertex gains at most 5 chords from a face whose boundary is a simple cycle (3
# unless it ends a first-ladder rung); the peeling of repeated vertices before it
# can give a vertex a few more (9 was the most seen), and every bound is a constant.
# On an input that is not a planar embedding a quadrangle can have both diagonals
# blocked, and the triangulation is rejected.
#
#
# The edges are kept in a hash set, so an adjacency test is `O(1)` in expectation
# and no edge is added twice. A scan of the smaller degree would be quadratic: two
# adjacent hubs joined by many paths are opposite on every face, so the same edge
# would be found, and charged, once per face. With the set, every test is `O(1)`
# and there are `O(n + m)` of them, so the whole triangulation is `O(n + m)`
# expected time.
"""
    triangulate(g) -> same type as g

Add edges to `g` until every face is a triangle, preserving the embedding and
the vertex set, in `O(nv + ne)` expected time. No edge is added twice, so a simple graph
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

function triangulate(g::PlanarGraph{T}) where {T}
    return _retype(T, PlanarGraph(triangulate(DartGraph(g))))
end

function _add_chords!(h::DartGraph)
    euler_holds(h) || throw(
        ArgumentError(
            "cannot triangulate: the rotation system does not describe a planar embedding",
        ),
    )
    orbits, _ = face_orbits(h)
    occurrences = zeros(Int, h.n)
    b = _Builder(h)
    for orbit in orbits
        length(orbit) <= 3 && continue
        orbit = _peel_repeats!(b, orbit, occurrences)
        length(orbit) <= 3 && continue
        _ladder!(b, orbit)
    end
    return h
end

# The embedded graph under construction and the set of its edges, with every
# chord added so far, so that no chord repeats an edge.
struct _Builder
    graph::DartGraph
    edges::Set{Int}
end

function _Builder(h::DartGraph)
    edges = Set{Int}()
    sizehint!(edges, 3 * h.n)
    for d = 1:ndarts(h)
        h.tail[d] < h.head[d] && push!(edges, _edge_key(h.n, h.tail[d], h.head[d]))
    end
    return _Builder(h, edges)
end

_edge_key(n::Int, u::Int, w::Int) = min(u, w) * (n + 1) + max(u, w)

# --- faces that repeat a vertex ------------------------------------------

# Clip the occurrences of repeated vertices off the boundary walk `orbit`, going
# around the list of those occurrences and skipping the one next to a clip.
# `occurrences[v]` counts the occurrences of `v` that are left; it is all zero
# again on return. Return the remaining boundary walk, which is a simple cycle.
function _peel_repeats!(b::_Builder, orbit::Vector{Int}, occurrences::Vector{Int})
    verts = Int[tail(b.graph, d) for d in orbit]
    for v in verts
        occurrences[v] += 1
    end
    candidates = findall(t -> occurrences[verts[t]] >= 2, eachindex(verts))
    result =
        isempty(candidates) ? orbit :
        _clip_repeats!(b, orbit, verts, occurrences, candidates)
    for v in verts
        occurrences[v] = 0
    end
    return result
end

function _clip_repeats!(
    builder::_Builder,
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
        dartp[b] = _add_chord!(builder, verts[b], verts[c], dartp[b], dartp[c])
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
function _add_chord!(b::_Builder, u::Int, v::Int, du::Int, dv::Int)
    h = b.graph
    forward = add_dart!(h, u, v, du)
    backward = add_dart!(h, v, u, dv)
    h.twin[forward] = backward
    h.twin[backward] = forward
    push!(b.edges, _edge_key(h.n, u, v))
    return forward
end

# --- simple faces --------------------------------------------------------

# A rung of the ladder: the chord between the positions `s` and `t`, with its two
# darts. The base rung is the boundary edge between the positions `1` and `k`.
struct _Rung
    s::Int
    t::Int
    forward::Int
    backward::Int
end

# Triangulate the simple face with the darts `orbit` (a dart leaves each boundary
# vertex along the face) with a ladder, then close the regions it leaves. A region
# of at least five vertices has all its rungs blocked by edges outside the face; a
# ladder started one position later pairs the vertices at a slant, which those
# edges do not block, so a slanted ladder leaves only quadrangles and a triangle.
function _ladder!(b::_Builder, orbit::Vector{Int}, slanted::Bool = false)
    for region in _regions(orbit, _add_rungs!(b, orbit))
        _close_region!(b, region, slanted)
    end
    return nothing
end

# Join the positions `i + 1` and `k - i` for `i = 1, 2, …` whenever the two
# vertices are not adjacent yet. The base rung comes first.
function _add_rungs!(b::_Builder, orbit::Vector{Int})
    k = length(orbit)
    verts = Int[tail(b.graph, d) for d in orbit]
    rungs = [_Rung(1, k, 0, orbit[k])]
    for i = 1:((k-3)÷2)
        s = i + 1
        t = k - i
        _adjacent(b, verts[s], verts[t]) && continue
        forward = _add_chord!(b, verts[s], verts[t], orbit[s], orbit[t])
        push!(rungs, _Rung(s, t, forward, b.graph.twin[forward]))
    end
    return rungs
end

# The boundary darts of the region between each two consecutive rungs, and of the
# one inside the last rung. A rung's forward dart bounds the region outside it, its
# backward dart the region inside.
function _regions(orbit::Vector{Int}, rungs::Vector{_Rung})
    regions = Vector{Vector{Int}}()
    for j = 2:length(rungs)
        a, b = rungs[j-1], rungs[j]
        push!(regions, vcat(orbit[b.t:(a.t-1)], a.backward, orbit[a.s:(b.s-1)], b.forward))
    end
    a = rungs[end]
    push!(regions, vcat(orbit[a.s:(a.t-1)], a.backward))
    return regions
end

# A triangle is done, a quadrangle takes the diagonal that is new, and a larger
# region, left by skipped rungs, gets a slanted ladder.
function _close_region!(b::_Builder, region::Vector{Int}, slanted::Bool)
    k = length(region)
    k <= 3 && return nothing
    if k > 4
        slanted && throw(_not_triangulable())
        return _ladder!(b, vcat(region[2:end], region[1]), true)
    end
    verts = Int[tail(b.graph, d) for d in region]
    for (s, t) in ((1, 3), (2, 4))
        _adjacent(b, verts[s], verts[t]) && continue
        _add_chord!(b, verts[s], verts[t], region[s], region[t])
        return nothing
    end
    throw(_not_triangulable())
end


# Adjacent vertices, by a lookup in the set of edges. The same pair is asked about
# once per face it is opposite on, so a scan of the smaller degree would charge
# one edge again and again (two adjacent hubs joined by many paths); the set
# makes every test `O(1)` in expectation.
_adjacent(b::_Builder, u::Int, w::Int) = _edge_key(b.graph.n, u, w) in b.edges

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
