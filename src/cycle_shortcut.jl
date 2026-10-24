# Project a cycle of the internal triangulation onto a cycle of the input graph.
#
# The Har-Peled / Nayyeri pipeline in cycle_separator.jl works on a
# triangulation. Its cycle can use chords that do not exist in the input, so
# the helpers here replace runs of fill edges by boundary arcs of the faces of
# the input, verify that the result is a simple cycle of the input, and orient
# it. Everything is linear in the size of the triangulation.

# --- cycles of the input graph -------------------------------------------

# The directed darts of `cyc`, or `nothing` when two consecutive vertices are
# not adjacent.
function _cycle_darts(g::DartGraph, cyc::AbstractVector{Int})
    k = length(cyc)
    k >= 3 || return nothing
    darts = Vector{Int}(undef, k)
    for i = 1:k
        a = cyc[i]
        b = cyc[mod1(i+1, k)]
        d = 0
        for x in darts_of(g, a)
            if g.head[x] == b
                d = x
                break
            end
        end
        d == 0 && return nothing
        darts[i] = d
    end
    return darts
end

function _cycle_edge_mask(g::DartGraph, cyc::AbstractVector{Int})
    darts = _cycle_darts(g, cyc)
    darts === nothing && return nothing
    mask = falses(ndarts(g))
    for d in darts
        mask[d] = true
        mask[g.twin[d]] = true
    end
    return mask
end

function _is_simple_cycle(g::DartGraph, cyc::AbstractVector{Int})
    length(cyc) >= 3 || return false
    allunique(cyc) || return false
    return _cycle_darts(g, cyc) !== nothing
end

function _check_cycle(g::DartGraph, cyc::Vector{Int})
    k = length(cyc)
    k >= 3 || throw(ArgumentError("a cycle needs at least 3 vertices, got $k"))
    all(v -> 1 <= v <= g.n, cyc) ||
        throw(ArgumentError("cycle vertices must lie in 1:$(g.n)"))
    allunique(cyc) || throw(ArgumentError("the cycle repeats a vertex"))
    _cycle_darts(g, cyc) === nothing && throw(
        ArgumentError(
            "consecutive cycle vertices must be adjacent, including the last " *
            "and the first",
        ),
    )
    return nothing
end

# The two sides of `cyc`, with the left side first. `nothing` when `cyc` is not
# a cycle of `g`.
function _cycle_sides(g::DartGraph, orbits, face_of::Vector{Int}, cyc::Vector{Int})
    darts = _cycle_darts(g, cyc)
    darts === nothing && return nothing
    blocked = falses(ndarts(g))
    seen = falses(g.n)
    for v in cyc
        seen[v] = true
    end
    for d in darts
        blocked[d] = true
        blocked[g.twin[d]] = true
    end
    left = _face_set(g, orbits, face_of, face_of[darts[1]], blocked)
    inside = Int[]
    for (f, orbit) in enumerate(orbits)
        left[f] || continue
        for d in orbit
            v = g.tail[d]
            seen[v] && continue
            seen[v] = true
            push!(inside, v)
        end
    end
    outside = Int[v for v = 1:g.n if !seen[v]]
    return _sorted_ids(inside, g.n), outside
end

# Replace the chords of a cycle of the triangulation by boundary paths of the
# input graph. `region[d]` is the face of `g` containing the face left of `d`.
function _shortcut_cycle(g::DartGraph, h::DartGraph, region::Vector{Int}, cyc::Vector{Int})
    length(cyc) >= 3 || return nothing
    darts = _cycle_darts(h, cyc)
    darts === nothing && return nothing

    h_orbits, h_face_of = face_orbits(h)
    inside = _cycle_inside_faces(h, h_orbits, h_face_of, darts)
    g_orbits, g_face_of = face_orbits(g)
    nf = length(g_orbits)
    inside_count, total_count = _region_counts(h_orbits, region, inside, nf)

    best = Int[]
    best_key = nothing
    for keep in _side_face_sets(inside_count, total_count)
        for candidate in _boundary_cycles(g, g_orbits, g_face_of, keep)
            key = _evaluate_cycle(g, g_orbits, g_face_of, candidate)
            key === nothing && continue
            if best_key === nothing || key < best_key
                best = candidate
                best_key = key
            end
        end
    end
    return isempty(best) ? nothing : best
end

# Faces of `h` on the left of the cycle whose darts are given.
function _cycle_inside_faces(h::DartGraph, orbits, face_of::Vector{Int}, darts::Vector{Int})
    blocked = falses(ndarts(h))
    for d in darts
        blocked[d] = true
        blocked[h.twin[d]] = true
    end
    return _face_set(h, orbits, face_of, face_of[darts[1]], blocked)
end

# For each face of `g`, how many faces of `h` it contains and how many of those
# are inside the cycle.
function _region_counts(orbits, region::Vector{Int}, inside, nf::Int)
    inside_count = zeros(Int, nf)
    total_count = zeros(Int, nf)
    for (f, orbit) in enumerate(orbits)
        r = region[orbit[1]]
        total_count[r] += 1
        inside[f] && (inside_count[r] += 1)
    end
    return inside_count, total_count
end

# The two candidate face sets of `g`: the faces reached by the separator and
# their complement.
function _side_face_sets(inside_count::Vector{Int}, total_count::Vector{Int})
    nf = length(inside_count)
    reached = BitVector(inside_count[r] > 0 for r = 1:nf)
    outside = BitVector(inside_count[r] < total_count[r] for r = 1:nf)
    return reached, outside
end

# Rank a candidate cycle of `g` (or return `nothing` when it is not one),
# orienting it so that its smaller side is on the left.
function _evaluate_cycle(g::DartGraph, orbits, face_of::Vector{Int}, candidate::Vector{Int})
    (length(candidate) >= 3 && allunique(candidate)) || return nothing
    sides = _cycle_sides(g, orbits, face_of, candidate)
    sides === nothing && return nothing
    inside, outside = sides
    if length(inside) > length(outside)
        reverse!(candidate)
        inside, outside = outside, inside
    end
    large = max(length(inside), length(outside))
    return (large > 2 * g.n / 3 ? 1 : 0, large, length(candidate))
end

# Shortcut one chord run at a time, keeping the original edges of the cycle.
# The triangulation is a copy of `g` with darts appended, so a dart with index
# at most `ndarts(g)` is an edge of `g`.
function _local_shortcut(
    g::DartGraph,
    h::DartGraph,
    region::Vector{Int},
    cyc::Vector{Int},
    g_orbits,
)
    darts = _cycle_darts(h, cyc)
    darts === nothing && return nothing
    k = length(darts)
    original = BitVector(d <= ndarts(g) for d in darts)
    any(original) || return nothing
    start = findfirst(original)
    start === nothing && return nothing
    darts = circshift(darts, 1 - start)
    original = circshift(original, 1 - start)
    verts = Int[h.tail[d] for d in darts]
    on_cycle = falses(h.n)
    for v in verts
        on_cycle[v] = true
    end

    out = Int[]
    i = 1
    while i <= k
        if original[i]
            push!(out, verts[i])
            i += 1
            continue
        end
        r = region[darts[i]]
        j = i
        while j < k && !original[j+1] && region[darts[j+1]] == r
            j += 1
        end
        a = verts[i]
        b = j < k ? verts[j+1] : verts[1]
        arc = _clean_arc(g, g_orbits[r], a, b, on_cycle)
        arc === nothing && return nothing
        append!(out, arc[1:(end-1)])
        i = j + 1
    end
    return _is_simple_cycle(g, out) ? out : nothing
end

# Boundary path of the face with the given orbit from `a` to `b` that contains
# no vertex of the cycle except its endpoints. When both paths qualify, take
# the shorter one.
function _clean_arc(g::DartGraph, orbit, a::Int, b::Int, on_cycle)
    bv = Int[g.tail[d] for d in orbit]
    length(bv) >= 3 || return nothing
    pa = findfirst(==(a), bv)
    pb = findfirst(==(b), bv)
    (pa === nothing || pb === nothing) && return nothing

    forward = Int[]
    p = pa
    while true
        push!(forward, bv[p])
        p == pb && break
        p = mod1(p + 1, length(bv))
    end
    backward = Int[]
    p = pa
    while true
        push!(backward, bv[p])
        p == pb && break
        p = mod1(p - 1, length(bv))
    end

    clean(arc) = all(v -> v == arc[1] || v == arc[end] || !on_cycle[v], arc)
    forward_clean = clean(forward)
    backward_clean = clean(backward)
    if !forward_clean
        return backward_clean ? backward : nothing
    end
    backward_clean || return forward
    return length(forward) <= length(backward) ? forward : backward
end

# Orient `cyc` so that cycle_sides reports the smaller side as the inside.
function _orient_cycle(g::DartGraph, orbits, face_of::Vector{Int}, cyc::Vector{Int})
    sides = _cycle_sides(g, orbits, face_of, cyc)
    sides === nothing && return cyc
    inside, outside = sides
    length(inside) > length(outside) && reverse!(cyc)
    return cyc
end

# Boundary cycles of a set of faces of `g`.
function _boundary_cycles(g::DartGraph, orbits, face_of::Vector{Int}, keep)
    bnd = _directed_boundary(g, orbits, face_of, keep)
    cycles = Vector{Vector{Int}}()
    visited = falses(ndarts(g))
    is_boundary = d -> bnd[d]
    for start = 1:ndarts(g)
        bnd[start] || continue
        visited[start] && continue
        cyc = Int[]
        d = start
        closed = false
        for _ = 1:ndarts(g)
            visited[d] = true
            push!(cyc, g.tail[d])
            d = _next_boundary_dart(g, d, is_boundary)
            d == 0 && break
            if d == start
                closed = true
                break
            end
            visited[d] && break
        end
        closed && push!(cycles, cyc)
    end
    return cycles
end
