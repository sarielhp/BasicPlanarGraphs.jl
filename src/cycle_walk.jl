# Closed walks as cycle separators.
#
# A cycle separator of a triangulation turns into a closed walk of any connected
# plane graph: a fill edge (a, b) inside a face is replaced by a path along the
# boundary of that face. The walk only gains vertices, so the vertices off the
# walk keep the side they had; the walk may repeat vertices and edges.
#
# Here the triangulation is the virtual one of separator.jl: every face with
# four or more sides gets a star vertex, and the star edge of the balanced dual
# edge closes the BFS tree path through the face. Everything is linear.

# --- walks of the input graph --------------------------------------------

# The directed darts of the closed walk `cyc`, or `nothing` when two
# consecutive vertices are not adjacent.
function _cycle_darts(g::DartGraph, cyc::AbstractVector{Int})
    k = length(cyc)
    k >= 3 || return nothing
    darts = Vector{Int}(undef, k)
    for i = 1:k
        a = cyc[i]
        b = cyc[mod1(i+1, k)]
        d = 0
        for x in incident_darts(g, a)
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
    _cycle_darts(g, cyc) === nothing && throw(
        ArgumentError(
            "consecutive cycle vertices must be adjacent, including the last " *
            "and the first",
        ),
    )
    return nothing
end

# --- the two sides of a closed walk --------------------------------------

# Winding numbers of the faces: crossing an edge from its right to its left
# changes the number by the times the walk uses the edge forward minus backward,
# so an edge used both ways (a spike of the walk) is free to cross. The face
# right of the first dart gets 0. A walk that keeps to one side of itself has
# two adjacent values, the larger one on its left. A vertex off the walk touches
# faces of one winding number only.
function _winding(g::DartGraph, orbits, face_of::Vector{Int}, darts::Vector{Int})
    weight = zeros(Int, ndarts(g))
    for d in darts
        weight[d] += 1
        weight[g.twin[d]] -= 1
    end
    wind = fill(typemin(Int), length(orbits))
    seed = face_of[g.twin[darts[1]]]
    wind[seed] = 0
    queue = Int[seed]
    for f in queue
        for d in orbits[f]
            h = face_of[g.twin[d]]
            expected = wind[f] - weight[d]
            if wind[h] == typemin(Int)
                wind[h] = expected
                push!(queue, h)
            elseif wind[h] != expected
                return nothing
            end
        end
    end
    return wind
end

# The two sides of the closed walk `cyc`, left side first. `nothing` when `cyc`
# is not a closed walk of `g` that keeps to one side of itself.
function _cycle_sides(g::DartGraph, orbits, face_of::Vector{Int}, cyc::Vector{Int})
    darts = _cycle_darts(g, cyc)
    darts === nothing && return nothing
    wind = _winding(g, orbits, face_of, darts)
    wind === nothing && return nothing
    lo, hi = extrema(wind)
    hi - lo <= 1 || return nothing
    seen = falses(g.n)
    for v in cyc
        seen[v] = true
    end
    inside = Int[]
    for (f, orbit) in enumerate(orbits)
        (hi > lo && wind[f] == hi) || continue
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

# Orient `cyc` so that `cycle_sides` reports the smaller side as the inside.
function _orient_cycle(g::DartGraph, orbits, face_of::Vector{Int}, cyc::Vector{Int})
    sides = _cycle_sides(g, orbits, face_of, cyc)
    sides === nothing && return cyc
    length(sides[1]) > length(sides[2]) && reverse!(cyc)
    return cyc
end

# True when `cyc` is a closed walk whose two sides each hold at most `2n/3`.
function _is_balanced_walk(g::DartGraph, orbits, face_of::Vector{Int}, cyc::Vector{Int})
    sides = _cycle_sides(g, orbits, face_of, cyc)
    sides === nothing && return false
    return 3 * max(length(sides[1]), length(sides[2])) <= 2 * g.n
end

# --- the balanced walk ---------------------------------------------------

# Lemma 2 on the virtual triangulation with every vertex of cost one: the
# fundamental cycle of the balanced dual edge has at most `2n/3` vertices
# strictly on each side. The BFS tree is rooted at a pseudo-center so that its
# paths are short.
function _balanced_walk(dg::DartGraph)
    tri = _triangles(dg)
    root = _pseudo_center(dg, _bfs(dg, 1).order[end])
    tree = _bfs(dg, root)
    dual = _dual_tree(dg, tri, tree.parent_dart)
    parent, size = _root_dual(dual, tri.ntri)
    search = _Search(
        tri,
        dual,
        trues(dg.n),
        dg.n,
        parent,
        size,
        Int[],
        _Region(dg.n),
        _Region(dg.n),
    )
    e = _search_edge(search)
    walk = _tree_walk(tree, dual.ep[e], dual.eq[e])
    dual.ef[e] == 0 && return walk
    return _close_through_face(walk, tri.orbits[dual.ef[e]], dg, dual.ej[e])
end

# The BFS tree path from `p` to `q`.
function _tree_walk(tree, p::Int, q::Int)
    left = Int[]
    right = Int[]
    a, b = p, q
    while a != b
        if tree.depth[a] >= tree.depth[b]
            push!(left, a)
            a = tree.parent[a]
        else
            push!(right, b)
            b = tree.parent[b]
        end
    end
    push!(left, a)
    return append!(left, reverse!(right))
end

# `walk` runs from corner 1 to corner `j` of the face with the given darts. The
# star edge of the balanced dual edge is replaced by the shorter boundary arc
# from corner `j` back to corner 1, unless that leaves fewer than three vertices.
function _close_through_face(walk::Vector{Int}, orbit::Vector{Int}, dg::DartGraph, j::Int)
    k = length(orbit)
    corner(i) = dg.tail[orbit[i]]
    forward = Int[corner(i) for i = (j+1):k]
    backward = Int[corner(i) for i = (j-1):-1:2]
    short, long =
        length(forward) <= length(backward) ? (forward, backward) : (backward, forward)
    return append!(walk, length(walk) + length(short) >= 3 ? short : long)
end
