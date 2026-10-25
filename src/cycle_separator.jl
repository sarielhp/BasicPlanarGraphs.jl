# Har-Peled / Nayyeri cycle separator.
#
# Reference: S. Har-Peled and A. Nayyeri, "A Simple Algorithm for Computing a
# Cycle Separator", Computational Geometry: Theory and Applications 83 (2019),
# 101-110, arXiv:1709.08122. The three phases follow the paper:
#
#   1. a (possibly long) cycle separator `S` from a balanced edge of the dual
#      spanning tree (Section 3.1, after Lipton and Tarjan);
#   2. nested level cycles `C_0, C_1, ...` around the BFS root (Section 3.2);
#   3. a short cycle from a small ladder of those cycles (Section 3.3).
#
# The paper works on a triangulation, and phases 2 and 3 give its `O(sqrt(n))`
# bound there. For any other connected graph the cycle comes from the virtual
# star triangulation of separator.jl and is lifted to a closed walk of the
# input (cycle_walk.jl); it is balanced, but its length is not bounded. On a
# triangulation the result of phases 2 and 3 is checked and the exact balanced walk
# is the fallback, so balance never depends on their heuristics; the exact walk
# itself is balanced by Lemma 2 (and its lifting keeps that), and is not checked.
# Phase 3 includes the ring outside the last rung, which the length bound needs.

# --- public interface ----------------------------------------------------

"""
    cycle_separator(g) -> Vector{Int}

A balanced cycle separator of the connected embedded planar graph `g`, in
`O(nv(g))` time.

The result is a closed walk of `g`: at least three vertices, and consecutive
entries, including the last and the first, are adjacent in `g`. Removing its
vertices leaves two sides, and use [`cycle_sides`](@ref) to obtain them: no edge
joins them, and each has at most `2/3` of the vertices. The walk is oriented so
that the side `cycle_sides` calls the inside is the smaller one.

If every face of `g` is a triangle the walk is a simple cycle with `O(sqrt(n))`
vertices (S. Har-Peled and A. Nayyeri, 2019). For other graphs the walk may
repeat vertices and edges (a fill edge of a triangulation becomes a path along
the face boundary) and its length is not bounded; both sides are still at most
`2/3` of the vertices.

```julia
g = planar_grid(20, 20)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)
isempty(intersect(inside, outside))              # true
max(length(inside), length(outside)) <= 2 / 3 * nv(g)
```

The graph must be connected and have at least three vertices.
"""
function cycle_separator(g::PlanarGraph)
    dg = DartGraph(g)
    dg.n >= 3 ||
        throw(ArgumentError("a cycle separator needs at least 3 vertices, got $(dg.n)"))
    count_components(dg) == 1 ||
        throw(ArgumentError("cycle_separator requires a connected graph"))
    return _cycle_separator(dg)
end

"""
    cycle_sides(g, cycle) -> (inside, outside)

Split the vertices of the connected graph `g` by the closed walk `cycle` and
return the two sides. The walk is read in the order given: because rotations
are counter-clockwise, the vertices locally to the **left** of the traversal
are `inside` and those to the right are `outside`. Reversing the walk swaps the
two results; neither side is privileged beyond this convention.

`cycle` must be a sequence of at least three vertices of `g` whose consecutive
entries, including the last and the first, are adjacent, and which may repeat
vertices and edges (as the result of [`cycle_separator`](@ref) may); otherwise,
or when the faces get more than two winding values (the walk crosses itself,
goes around twice, or touches itself with opposite orientations), an
`ArgumentError` is raised. Vertices of the
walk itself are in neither result, and both vectors are sorted ascending.

```julia
g = planar_wheel(10)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)
sort(vcat(inside, outside, C)) == collect(1:nv(g))   # true
```
"""
function cycle_sides(g::PlanarGraph, cycle::AbstractVector{<:Integer})
    dg = DartGraph(g)
    count_components(dg) == 1 ||
        throw(ArgumentError("cycle_sides requires a connected graph"))
    cyc = Int[Int(v) for v in cycle]
    _check_cycle(dg, cyc)
    orbits, face_of = face_orbits(dg)
    sides = _cycle_sides(dg, orbits, face_of, cyc)
    sides === nothing && throw(
        ArgumentError(
            "the given vertices do not form a closed walk that keeps to one side of itself",
        ),
    )
    return sides
end

# --- the pipeline ---------------------------------------------------------

function _cycle_separator(g::DartGraph; delta::Union{Nothing,Int} = nothing)
    orbits, face_of = face_orbits(g)
    if all(orbit -> length(orbit) == 3, orbits)
        step = delta === nothing ? max(1, ceil(Int, sqrt(g.n / 2))) : delta
        step >= 1 || throw(ArgumentError("delta must be positive, got $step"))
        result = _cycle_in_triangulation(g, step)
        for candidate in (result.cycle, result.fundamental)
            _is_balanced_walk(g, orbits, face_of, candidate) &&
                return _orient_cycle(g, orbits, face_of, candidate)
        end
    end
    return _orient_cycle(g, orbits, face_of, _balanced_walk(g))
end

# --- phase 1: a possibly long cycle separator ----------------------------

# BFS tree of the triangulation, a balanced edge of the dual of the non-tree
# edges, and the root of the fundamental cycle through that edge.
function _fundamental_cycle(h::DartGraph)
    parent, parent_dart, depth, _ = bfs_tree(h, [1])
    tree = _tree_edges(h, parent_dart)
    face_parent, face_dart, order = _dual_spanning_tree(h, tree)
    best = _balanced_dual_dart(face_parent, face_dart, order)
    u = h.tail[best]
    v = h.head[best]
    return (parent = parent, depth = depth, u = u, v = v, root = _lca(parent, depth, u, v))
end

# Darts of the BFS tree, in both directions.
function _tree_edges(h::DartGraph, parent_dart::Vector{Int})
    tree = falses(ndarts(h))
    for v = 1:h.n
        d = parent_dart[v]
        d == 0 && continue
        tree[d] = true
        tree[h.twin[d]] = true
    end
    return tree
end

# The dual of the non-tree edges is a spanning tree when the embedding is
# planar, so no union-find is needed. Returns its parent pointers and a BFS
# order over the faces.
function _dual_spanning_tree(h::DartGraph, tree::BitVector)
    orbits, face_of = face_orbits(h)
    nf = length(orbits)
    adj = [Tuple{Int,Int}[] for _ = 1:nf]
    edges = 0
    for d = 1:ndarts(h)
        (d < h.twin[d] && !tree[d]) || continue
        f1 = face_of[d]
        f2 = face_of[h.twin[d]]
        push!(adj[f1], (f2, d))
        push!(adj[f2], (f1, d))
        edges += 1
    end
    edges == nf - 1 || throw(_not_planar_error())

    face_parent = zeros(Int, nf)
    face_dart = zeros(Int, nf)
    order = Int[1]
    seen = falses(nf)
    seen[1] = true
    for f in order
        for (g, d) in adj[f]
            seen[g] && continue
            seen[g] = true
            face_parent[g] = f
            face_dart[g] = d
            push!(order, g)
        end
    end
    length(order) == nf || throw(_not_planar_error())
    return face_parent, face_dart, order
end

# Dart of the dual edge whose removal splits the faces most evenly.
function _balanced_dual_dart(
    face_parent::Vector{Int},
    face_dart::Vector{Int},
    order::Vector{Int},
)
    nf = length(face_parent)
    subtree = ones(Int, nf)
    for i = length(order):-1:2
        f = order[i]
        subtree[face_parent[f]] += subtree[f]
    end
    best = 0
    best_key = typemax(Int)
    for i = 2:length(order)
        f = order[i]
        key = max(subtree[f], nf - subtree[f])
        if key < best_key
            best_key = key
            best = face_dart[f]
        end
    end
    best == 0 && throw(_not_planar_error())
    return best
end

_not_planar_error() = ArgumentError(
    "cannot compute a separator: the rotation system does not describe " *
    "a planar embedding",
)

function _lca(parent::Vector{Int}, depth::Vector{Int}, u::Int, v::Int)
    a, b = u, v
    while depth[a] > depth[b]
        a = parent[a]
    end
    while depth[b] > depth[a]
        b = parent[b]
    end
    while a != b
        a = parent[a]
        b = parent[b]
    end
    return a
end

# Root-to-`u` and root-to-`v` tree paths, both including the root.
function _tree_paths(parent::Vector{Int}, root::Int, u::Int, v::Int)
    pu = Int[u]
    x = u
    while x != root
        x = parent[x]
        push!(pu, x)
    end
    reverse!(pu)
    pv = Int[v]
    x = v
    while x != root
        x = parent[x]
        push!(pv, x)
    end
    reverse!(pv)
    return pu, pv
end

# BFS tree from `root` with the two given paths forced: this is the tree of
# Lemma 3.1, whose levels phase 2 uses.
function _forced_tree(h::DartGraph, root::Int, pu::Vector{Int}, pv::Vector{Int})
    parent, _, depth, _ = bfs_tree(h, [root])
    for path in (pu, pv), i = 2:length(path)
        parent[path[i]] = path[i-1]
    end
    return parent, depth
end

# --- phase 2: the nested level cycles ------------------------------------

# Level of a face: the largest level of its three vertices.
function _face_levels(h::DartGraph, orbits, depth::Vector{Int})
    return [maximum(depth[h.tail[d]] for d in orbit) for orbit in orbits]
end

# The faces of the triangulation and their levels, bundled for the phases that
# need all three.
struct _FaceLevels
    orbits::Vector{Vector{Int}}
    face_of::Vector{Int}
    flevel::Vector{Int}
end

_is_level_boundary(h::DartGraph, faces::_FaceLevels, i::Int, d::Int) =
    faces.flevel[faces.face_of[d]] == i && faces.flevel[faces.face_of[h.twin[d]]] == i + 1

# The cycle `C_i` through `w`, traced with the faces of level at most `i` on
# the left. Returns an empty vector when `w` is not on that boundary.
function _trace_level_cycle(h::DartGraph, faces::_FaceLevels, i::Int, w::Int)
    start = 0
    for d in darts_of(h, w)
        if _is_level_boundary(h, faces, i, d)
            start = d
            break
        end
    end
    start == 0 && return Int[]
    is_boundary = d -> _is_level_boundary(h, faces, i, d)
    return _trace_ordered_cycle(h, start, is_boundary)
end

# --- phase 3: a short cycle from a ladder --------------------------------

function _cycle_in_triangulation(h::DartGraph, step::Int)
    first = _fundamental_cycle(h)
    pu, pv = _tree_paths(first.parent, first.root, first.u, first.v)
    parent, depth = _forced_tree(h, first.root, pu, pv)
    hmax = max(depth[first.u], depth[first.v])
    fundamental = vcat(pu, reverse(pv)[1:(end-1)])
    length(fundamental) <= 2 * step &&
        return (cycle = fundamental, fundamental = fundamental)

    orbits, face_of = face_orbits(h)
    faces = _FaceLevels(orbits, face_of, _face_levels(h, orbits, depth))
    psi = depth[first.u] >= depth[first.v] ? first.u : first.v
    psi_path = Int[psi]
    x = psi
    while x != first.root
        x = parent[x]
        push!(psi_path, x)
    end
    reverse!(psi_path)

    cycles = Vector{Vector{Int}}(undef, hmax)
    for i = 0:(hmax-1)
        cycles[i+1] = i == 0 ? Int[] : _trace_level_cycle(h, faces, i, psi_path[i+1])
    end

    sums = zeros(Int, step - 1)
    for i0 = 1:(step-1)
        total = 0
        for alpha = i0:step:(hmax-1)
            total += length(cycles[alpha+1])
        end
        sums[i0] = total
    end
    isempty(sums) && return (cycle = fundamental, fundamental = fundamental)
    i0 = argmin(sums)
    ladder = Int[i0+j*step for j = 0:((hmax-1-i0)÷step)]
    ladder_cycles = [cycles[alpha+1] for alpha in ladder]
    any(isempty, ladder_cycles) && return (cycle = fundamental, fundamental = fundamental)

    rungs = (levels = ladder, cycles = ladder_cycles)
    cycle = _ladder_cycle(h, faces, rungs, fundamental, first.root)
    _is_simple_cycle(h, cycle) || return (cycle = fundamental, fundamental = fundamental)
    return (cycle = cycle, fundamental = fundamental)
end

# The short cycle of Section 3.3: a balanced ladder cycle, one of the two
# pieces of the heavy ring, or the inner region together with one piece.
function _ladder_cycle(
    h::DartGraph,
    faces::_FaceLevels,
    rungs,
    fundamental::Vector{Int},
    root::Int,
)
    ladder_cycles = rungs.cycles
    k = length(rungs.levels)
    index = _ring_index(h, faces, rungs, root)
    index === nothing && return Int[]
    nf = length(faces.orbits)

    histogram = zeros(Int, k + 1)
    for r in index
        histogram[r+1] += 1
    end
    inside = zeros(Int, k + 1)
    for j = 1:k
        inside[j+1] = inside[j] + histogram[j]
    end

    for j = 1:k
        if nf / 3 <= inside[j+1] <= 2nf / 3
            return ladder_cycles[j]
        end
    end

    # The heavy ring lies between two rungs, inside the first rung (heavy = 1), or,
    # when every rung is light, outside the last one (heavy = k + 1): the case of
    # the trivial outermost cycle. Its two pieces are cut by the fundamental
    # cycle, which crosses the ring, so the outside is empty and a piece with a
    # third of the faces exists.
    heavy = findfirst(j -> inside[j] < nf / 3 && inside[j+1] > 2nf / 3, 1:k)
    if heavy === nothing
        (k >= 1 && inside[k+1] < nf / 3) || return Int[]
        heavy = k + 1
    end
    parts =
        (inner = BitVector(index .< heavy - 1), annulus = BitVector(index .== heavy - 1))
    seed = ladder_cycles[min(heavy, k)][1]
    return _heavy_ring_cycle(h, faces, fundamental, seed, parts)
end

# Number of the ladder ring that contains each face: 0 inside the innermost
# ring, `k` outside the outermost one. `nothing` if the ladder is not nested as
# the construction needs it to be.
function _ring_index(h::DartGraph, faces::_FaceLevels, rungs, root::Int)
    ladder = rungs.levels
    ladder_cycles = rungs.cycles
    nf = length(faces.orbits)
    delta = zeros(Int8, ndarts(h))
    for (j, alpha) in enumerate(ladder)
        darts = _cycle_darts(h, ladder_cycles[j])
        darts === nothing && return nothing
        for d in darts
            if faces.flevel[faces.face_of[d]] == alpha
                delta[d] = 1
                delta[h.twin[d]] = -1
            else
                delta[d] = -1
                delta[h.twin[d]] = 1
            end
        end
    end

    index = zeros(Int, nf)
    visited = falses(nf)
    seed = faces.face_of[darts_of(h, root)[1]]
    visited[seed] = true
    queue = Int[seed]
    head = 1
    while head <= length(queue)
        f = queue[head]
        head += 1
        for d in faces.orbits[f]
            g = faces.face_of[h.twin[d]]
            visited[g] && continue
            visited[g] = true
            index[g] = index[f] + delta[d]
            push!(queue, g)
        end
    end
    all(r -> 0 <= r <= length(ladder), index) || return nothing
    return index
end

function _heavy_ring_cycle(
    h::DartGraph,
    faces::_FaceLevels,
    fundamental::Vector{Int},
    seed_vertex::Int,
    parts,
)
    edge_mask = _cycle_edge_mask(h, fundamental)
    edge_mask === nothing && return Int[]
    first_dart = findfirst(edge_mask)
    first_dart === nothing && return Int[]
    in_fundamental =
        _face_set(h, faces.orbits, faces.face_of, faces.face_of[first_dart], edge_mask)
    first = parts.annulus .& in_fundamental
    second = parts.annulus .& .!in_fundamental
    nf = length(faces.orbits)
    if count(first) >= nf / 3
        return _cycle_of_face_set(h, faces.orbits, faces.face_of, first, seed_vertex)
    end
    if count(second) >= nf / 3
        return _cycle_of_face_set(h, faces.orbits, faces.face_of, second, seed_vertex)
    end
    return _cycle_of_face_set(
        h,
        faces.orbits,
        faces.face_of,
        parts.inner .| first,
        seed_vertex,
    )
end

# --- face sets and their boundaries --------------------------------------

# Faces reachable from `seed` without crossing a dart marked in `blocked`.
function _face_set(h::DartGraph, orbits, face_of::Vector{Int}, seed::Int, blocked)
    seen = falses(length(orbits))
    seen[seed] = true
    queue = Int[seed]
    head = 1
    while head <= length(queue)
        f = queue[head]
        head += 1
        for d in orbits[f]
            blocked[d] && continue
            g = face_of[h.twin[d]]
            seen[g] && continue
            seen[g] = true
            push!(queue, g)
        end
    end
    return seen
end

# Darts of edges with a face of `face_set` on the left and another face on the
# right, that is, the boundary of the region oriented with the region left.
function _directed_boundary(h::DartGraph, orbits, face_of::Vector{Int}, face_set)
    bnd = falses(ndarts(h))
    for (f, orbit) in enumerate(orbits)
        face_set[f] || continue
        for d in orbit
            face_set[face_of[h.twin[d]]] || (bnd[d] = true)
        end
    end
    return bnd
end

# Cycle that bounds the region made of `face_set`, starting at a dart at
# `seed_vertex`. Empty when that vertex is not on the boundary.
function _cycle_of_face_set(h::DartGraph, orbits, face_of, face_set, seed_vertex)
    bnd = _directed_boundary(h, orbits, face_of, face_set)
    start = 0
    for d in darts_of(h, seed_vertex)
        bnd[d] && (start = d; break)
    end
    start == 0 && return Int[]
    return _trace_boundary_cycle(h, bnd, start)
end

# --- tracing cycles along a boundary -------------------------------------

# Follow the boundary at each vertex: rotate clockwise from the opposite dart
# until a boundary dart is found.
function _next_boundary_dart(h::DartGraph, d::Int, is_boundary)
    v = h.head[d]
    x = h.prev_rot[h.twin[d]]
    for _ = 1:h.deg[v]
        is_boundary(x) && return x
        x = h.prev_rot[x]
    end
    return 0
end

function _trace_ordered_cycle(h::DartGraph, start::Int, is_boundary)
    cyc = Int[]
    d = start
    for _ = 1:ndarts(h)
        push!(cyc, h.tail[d])
        d = _next_boundary_dart(h, d, is_boundary)
        d == 0 && return cyc
        d == start && return cyc
    end
    return cyc
end

function _trace_boundary_cycle(h::DartGraph, bnd::BitVector, start::Int)
    is_boundary = d -> bnd[d]
    return _trace_ordered_cycle(h, start, is_boundary)
end
