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
# The paper works on a triangulated graph. This implementation triangulates the
# input internally and then *shortcuts* the chords of the resulting cycle: every
# maximal run of fill edges is replaced by the part of the boundary of the face
# of the input graph that contains it. The shortcut is neither in the paper nor
# in OGDF's `SeparatorHarPeled` (which returns only a vertex set); it is what
# makes the returned vector a cycle of the caller's graph.

# --- public interface ----------------------------------------------------

"""
    cycle_separator(g) -> Vector{Int}

A balanced cycle separator of the connected, 2-connected embedded planar graph
`g`, computed by the algorithm of S. Har-Peled and A. Nayyeri (*A Simple
Algorithm for Computing a Cycle Separator*, 2019). The separator phases and the
shortcut are linear in the size of the internal triangulation.

The result is a simple cycle of `g`, listed in traversal order: consecutive
vertices, including the last and the first, are adjacent in `g`, and no vertex
is repeated. Removing the cycle splits `g` into the vertices on its two sides;
use [`cycle_sides`](@ref) to obtain them. The cycle is oriented so that the
left side, which [`cycle_sides`](@ref) calls the inside, is the smaller of the
two sides.

Both sides have at most `2/3` of the vertices, and the cycle has `O(sqrt(n))`
vertices in the triangulation used internally. On graphs whose faces are already
triangles no shortcut is needed and that bound applies to the returned cycle.

```julia
g = planar_grid(20, 20)
C = cycle_separator(g)
inside, outside = cycle_sides(g, C)
isempty(intersect(inside, outside))              # true
max(length(inside), length(outside)) <= 2 / 3 * nv(g)
```

The graph must be connected and 2-connected (so every face is a simple cycle,
which is what [`triangulate`](@ref) needs), with at least three vertices.
"""
function cycle_separator(g::PlanarGraph)
    dg = DartGraph(g)
    dg.n >= 3 ||
        throw(ArgumentError("a cycle separator needs at least 3 vertices, got $(dg.n)"))
    count_components(dg) == 1 ||
        throw(ArgumentError("cycle_separator requires a connected graph"))
    _faces_are_simple_cycles(dg) || throw(
        ArgumentError(
            "cycle_separator requires a 2-connected graph: a face boundary " *
            "repeats a vertex, so the graph has a bridge or a cut vertex",
        ),
    )
    return _cycle_separator(dg)
end

# A connected plane graph is 2-connected exactly when every face boundary is a
# simple cycle; a bridge or a cut vertex makes some boundary walk repeat a
# vertex. `triangulate` does not reject every such graph, so check it here.
function _faces_are_simple_cycles(g::DartGraph)
    orbits, _ = face_orbits(g)
    stamp = zeros(Int, g.n)
    for (f, orbit) in enumerate(orbits)
        for d in orbit
            v = g.tail[d]
            stamp[v] == f && return false
            stamp[v] = f
        end
    end
    return true
end

"""
    cycle_sides(g, cycle) -> (inside, outside)

Split the vertices of the connected graph `g` by the simple cycle `cycle` and
return the two sides. The cycle is read in the order given: because rotations
are counter-clockwise, the vertices locally to the **left** of the traversal
are `inside` and those to the right are `outside`. Reversing the cycle swaps the
two results; neither side is privileged beyond this convention.

`cycle` must be a sequence of at least three distinct vertices of `g` whose
consecutive entries, including the last and the first, are adjacent; otherwise
an `ArgumentError` is raised. Vertices of the cycle itself are in neither
result, and both vectors are sorted ascending.

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
    sides === nothing && throw(ArgumentError("the given vertices do not form a cycle"))
    return sides
end

# --- the pipeline ---------------------------------------------------------

function _cycle_separator(g::DartGraph; delta::Union{Nothing,Int} = nothing)
    h, region = _triangulate_with_regions(g)
    step = delta === nothing ? max(1, ceil(Int, sqrt(h.n / 2))) : delta
    step >= 1 || throw(ArgumentError("delta must be positive, got $step"))
    result = _cycle_in_triangulation(h, step)
    g_orbits, g_face_of = face_orbits(g)
    for candidate in (result.cycle, result.fundamental)
        cyc = _shortcut_cycle(g, h, region, candidate)
        cyc === nothing || return _orient_cycle(g, g_orbits, g_face_of, cyc)
    end
    for candidate in (result.cycle, result.fundamental)
        cyc = _local_shortcut(g, h, region, candidate, g_orbits)
        cyc === nothing || return _orient_cycle(g, g_orbits, g_face_of, cyc)
    end
    error(
        "cycle separator: the cycle of the triangulation could not be reduced " *
        "to a cycle of the input graph",
    )
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
# planar. Returns its parent pointers and a BFS order over the faces.
function _dual_spanning_tree(h::DartGraph, tree::BitVector)
    orbits, face_of = face_orbits(h)
    nf = length(orbits)
    dsu = collect(1:nf)
    adj = [Tuple{Int,Int}[] for _ = 1:nf]
    for d = 1:ndarts(h)
        d < h.twin[d] || continue
        tree[d] && continue
        f1 = face_of[d]
        f2 = face_of[h.twin[d]]
        f1 == f2 && continue
        a = _find_root!(dsu, f1)
        b = _find_root!(dsu, f2)
        a == b && continue
        dsu[a] = b
        push!(adj[f1], (f2, d))
        push!(adj[f2], (f1, d))
    end

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

function _find_root!(dsu::Vector{Int}, x::Int)
    while dsu[x] != x
        dsu[x] = dsu[dsu[x]]
        x = dsu[x]
    end
    return x
end

_not_planar_error() = ArgumentError(
    "cannot compute a cycle separator: the rotation system does not describe " *
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

    heavy = findfirst(j -> inside[j] < nf / 3 && inside[j+1] > 2nf / 3, 1:k)
    heavy === nothing && return Int[]
    parts =
        (inner = BitVector(index .< heavy - 1), annulus = BitVector(index .== heavy - 1))
    return _heavy_ring_cycle(h, faces, fundamental, ladder_cycles[heavy][1], parts)
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
