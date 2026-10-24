# Lipton-Tarjan planar separator.
#
# Reference: R. J. Lipton and R. E. Tarjan, "A Separator Theorem for Planar
# Graphs", SIAM J. Appl. Math. 36(2), 1979. Lemma 2 (fundamental cycles of a
# bounded-radius tree), Lemma 3 (two levels plus the middle band) and Theorem 4
# (the algorithm).
#
# Every stage is a single pass over vertices, darts or faces, so the running
# time is O(n + m).

const DEFAULT_BALANCE = 2 / 3

# Candidate fundamental cycles to examine around the middle band. The first
# balanced one wins; this is a small constant, not a factor of n.
const MAX_CYCLE_CANDIDATES = 6

"""
    separator(g; balance = 2/3) -> (A, B, C)

Balanced vertex separator of the embedded planar graph `g`, computed by the
linear-time Lipton-Tarjan algorithm.

The three returned vectors partition `vertices(g)`: no edge joins `A` to `B`,
`|A|` and `|B|` are at most `balance * nv(g)`, and `|C| = O(sqrt(nv(g)))`
(at most `2*sqrt(2*nv(g))` for well-behaved inputs).

```julia
g = planar_grid(20, 20)
A, B, C = separator(g)
isempty(intersect(A, B))            # true
length(C) <= 2 * sqrt(2 * nv(g))    # true
```
"""
function separator(g::PlanarGraph; balance::Real = DEFAULT_BALANCE)
    0 < balance < 1 || throw(ArgumentError("balance must lie in (0, 1), got $balance"))
    n = nv(g)
    n == 0 && return (Int[], Int[], Int[])
    dg = DartGraph(g)
    label, ncomps = component_labels(dg)
    ncomps <= 1 && return _best_connected(dg, balance)
    return _separator_disconnected(dg, label, ncomps, balance)
end

# --- connected case ------------------------------------------------------

# Bounds on the refinement loop: each pass separates the largest remaining
# component, so the sizes shrink fast and a handful of passes always suffice.
const MAX_REFINEMENTS = 32

# The level structure, and with it the quality of the separator, depends
# strongly on the BFS root, so a few roots are tried and the best result is
# kept. The candidates are the endpoints and the midpoint of a BFS
# pseudo-diameter: on tree-like graphs that is the centroid.
function _best_connected(dg::DartGraph, balance::Real)
    best = _separator_connected(dg, 1, balance, 0)
    best_key = _candidate_key(best, dg.n, balance)
    for root in _root_candidates(dg, 1)[2:end]
        candidate = _separator_connected(dg, root, balance, 0)
        key = _candidate_key(candidate, dg.n, balance)
        if key < best_key
            best = candidate
            best_key = key
        end
    end
    return best
end

# Balanced results first, then the smallest separator, then the best split.
function _candidate_key(candidate, n::Integer, balance::Real)
    a, b, c = candidate
    unbalanced = _is_balanced(candidate, n, balance) ? 0 : 1
    return (unbalanced, length(c), max(length(a), length(b)))
end

function _root_candidates(dg::DartGraph, root::Integer)
    far, = _farthest(dg, root)
    mid = _midpoint(dg, far)
    far2, = _farthest(dg, mid)
    return unique(Int[root, far, mid, far2])
end

# A farthest vertex from `source`, with the BFS tree that reaches it.
function _farthest(dg::DartGraph, source::Integer)
    parent, _, depth, order = bfs_tree(dg, [source])
    best = source
    for v in order
        depth[v] > depth[best] && (best = v)
    end
    return best, parent, depth
end

# Vertex half way along the BFS path from the farthest vertex back to `source`.
function _midpoint(dg::DartGraph, source::Integer)
    far, parent, depth = _farthest(dg, source)
    v = far
    for _ = 1:(depth[far]÷2)
        v = parent[v]
    end
    return v
end

function _separator_connected(dg::DartGraph, root::Integer, balance::Real)
    return _separator_connected(dg, root, balance, 0)
end

function _separator_connected(dg::DartGraph, root::Integer, balance::Real, depth::Integer)
    n = dg.n
    n == 0 && return (Int[], Int[], Int[])
    n <= 2 && return (Int[], Int[], collect(1:n))
    removed = _initial_separator(dg, root, balance)
    return _refine(dg, removed, balance, depth)
end

# Separator of Theorem 4 before refinement: the two levels, plus a balanced
# cycle of the middle band when the band alone is larger than `balance * n`.
function _initial_separator(dg::DartGraph, root::Integer, balance::Real)
    n = dg.n
    levels = level_structure(dg, [root])
    l0, l2 = _level_plan(levels.depth, n)
    _, mid, _, base = _level_parts(levels.depth, l0, l2)
    length(mid) <= balance * n && return base
    best = base
    best_score = _largest_side(dg, base)
    for cycle in _band_cycles(dg, levels, l0, l2)
        candidate = _merge_vertices(base, cycle)
        score = _largest_side(dg, candidate)
        score <= balance * n && return candidate
        if score < best_score
            best = candidate
            best_score = score
        end
    end
    return best
end

function _largest_side(dg::DartGraph, removed::Vector{Int})
    a, b, _ = _split_and_pack(dg, removed)
    return max(length(a), length(b))
end

# Grow the separator until both sides fit into `balance * n`: the largest
# remaining component is split by the same algorithm. The components shrink by
# a constant factor each pass, so the added separators sum to O(sqrt(n)).
function _refine(dg::DartGraph, removed::Vector{Int}, balance::Real, depth::Integer)
    n = dg.n
    candidate = _split_and_pack(dg, removed)
    _is_balanced(candidate, n, balance) && return candidate
    (depth >= MAX_REFINEMENTS || length(removed) >= n) && return candidate
    comps = components_without(dg, _mask(removed, n))
    isempty(comps) && return candidate
    biggest = comps[argmax(length.(comps))]
    length(biggest) < 3 && return candidate
    sub, ids = induced(dg, _mask(biggest, n))
    _, _, extra = _separator_connected(sub, 1, balance, depth + 1)
    merged = _merge_vertices(removed, Int[ids[v] for v in extra])
    length(merged) == length(removed) && return candidate
    return _refine(dg, merged, balance, depth + 1)
end

function _mask(vertices::Vector{Int}, n::Integer)
    mask = falses(n)
    for v in vertices
        mask[v] = true
    end
    return mask
end

# Levels l0 and l2 of Theorem 4.
function _level_plan(level::Vector{Int}, n::Integer)
    counts = _level_counts(level)
    l1 = _middle_level(counts, n)
    k = _cumulative(counts, l1)
    l0 = _lower_level(counts, l1, k)
    l2 = _upper_level(counts, l1, n - k)
    return l0, l2
end

"""
Split the vertices by level: below `l0`, strictly between `l0` and `l2`, above
`l2`, and the two boundary levels themselves (the separator).
"""
function _level_parts(level::Vector{Int}, l0::Int, l2::Int)
    low = Int[]
    mid = Int[]
    high = Int[]
    base = Int[]
    for v in eachindex(level)
        lv = level[v]
        if lv == l0 || lv == l2
            push!(base, v)
        elseif lv < l0
            push!(low, v)
        elseif lv > l2
            push!(high, v)
        else
            push!(mid, v)
        end
    end
    return low, mid, high, base
end

# --- level bookkeeping ---------------------------------------------------

# counts[l + 1] = number of vertices on BFS level l.
function _level_counts(level::Vector{Int})
    counts = zeros(Int, maximum(level) + 1)
    for v in eachindex(level)
        counts[level[v]+1] += 1
    end
    return counts
end

# Number of vertices on levels 0 through l (levels outside 0:r are empty).
function _cumulative(counts::Vector{Int}, l::Int)
    l < 0 && return 0
    l >= length(counts) && return sum(counts)
    return sum(counts[1:(l+1)])
end

function _level_size(counts::Vector{Int}, l::Int)
    return 0 <= l < length(counts) ? counts[l+1] : 0
end

# Level whose prefix first reaches half of the vertices.
function _middle_level(counts::Vector{Int}, n::Integer)
    acc = 0
    for l = 0:(length(counts)-1)
        acc += counts[l+1]
        acc >= n / 2 && return l
    end
    return length(counts) - 1
end

# Level l0 <= l1 with |L(l0)| + 2(l1 - l0) <= 2 sqrt(k), as in Theorem 4.
function _lower_level(counts::Vector{Int}, l1::Int, k::Integer)
    limit = 2 * sqrt(k)
    l = l1
    while l >= 0 && _level_size(counts, l) + 2 * (l1 - l) > limit
        l -= 1
    end
    return l
end

# Level l2 >= l1 + 1 with |L(l2)| + 2(l2 - l1 - 1) <= 2 sqrt(n - k).
function _upper_level(counts::Vector{Int}, l1::Int, rest::Integer)
    limit = 2 * sqrt(rest)
    last = length(counts) - 1
    l = l1 + 1
    while l <= last && _level_size(counts, l) + 2 * (l - l1 - 1) > limit
        l += 1
    end
    return l
end

# --- middle band: fundamental cycles -------------------------------------

"""
Candidate cycles that split the middle band. The band (levels `l0+1 .. l2-1`)
is triangulated in place and cut by the fundamental cycle of the non-tree edge
whose dual subtree is most balanced, as in Lemma 2.
"""
function _band_cycles(dg::DartGraph, levels::LevelStructure, l0::Int, l2::Int)
    n = dg.n
    keep = BitVector(l0 < levels.depth[v] < l2 for v = 1:n)
    any(keep) || return Vector{Vector{Int}}()
    h = _band_triangulation(restrict(dg, keep))
    h === nothing && return Vector{Vector{Int}}()
    tree = _tree_dart_mask(h, levels, keep)
    _, face_of = face_orbits(h)
    weights = _face_weights(h, face_of, levels, l0, l2)
    sum(weights) == 0 && return Vector{Vector{Int}}()
    adj = _dual_spanning_tree(h, face_of, tree)
    roots = _forest_roots(levels, l0)
    return _candidate_cycles(h, adj, weights, levels, roots)
end

# Triangulating the band needs it to be 2-connected; when it is not, the cycles
# of the band cannot be used and the level separator has to do on its own.
function _band_triangulation(band::DartGraph)
    try
        return triangulate(band)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
end

# Mark the darts of `h` that belong to the BFS tree of the whole graph.
function _tree_dart_mask(h::DartGraph, levels::LevelStructure, keep::BitVector)
    mask = falses(ndarts(h))
    for v in eachindex(levels.parent)
        keep[v] || continue
        p = levels.parent[v]
        (p == 0 || !keep[p]) && continue
        for d in darts_of(h, v)
            h.head[d] == p || continue
            mask[d] = true
            mask[h.twin[d]] = true
            break
        end
    end
    return mask
end

# Each band vertex contributes one unit to one incident face, so the weight of
# the faces inside a cycle bounds the band vertices strictly inside it.
function _face_weights(
    h::DartGraph,
    face_of::Vector{Int},
    levels::LevelStructure,
    l0::Int,
    l2::Int,
)
    weights = zeros(Int, maximum(face_of; init = 0))
    for v = 1:h.n
        d = h.first[v]
        d == 0 && continue
        (l0 < levels.depth[v] < l2) || continue
        weights[face_of[d]] += 1
    end
    return weights
end

# Spanning tree of the dual made of the primal non-tree darts.
function _dual_spanning_tree(h::DartGraph, face_of::Vector{Int}, tree::BitVector)
    nf = maximum(face_of; init = 0)
    dsu = collect(1:nf)
    adj = [Tuple{Int,Int}[] for _ = 1:nf]
    for d = 1:ndarts(h)
        (d < h.twin[d] && !tree[d]) || continue
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
    return adj
end

function _find_root!(dsu::Vector{Int}, x::Int)
    while dsu[x] != x
        dsu[x] = dsu[dsu[x]]
        x = dsu[x]
    end
    return x
end

# Root of the BFS subtree a band vertex belongs to (0 outside the band).
function _forest_roots(levels::LevelStructure, l0::Int)
    roots = zeros(Int, length(levels.parent))
    for v in levels.order
        levels.depth[v] > l0 || continue
        p = levels.parent[v]
        roots[v] = (p == 0 || levels.depth[p] <= l0) ? v : roots[p]
    end
    return roots
end

# Primal darts whose dual edge splits the dual tree most evenly.
function _balanced_darts(adj::Vector{Vector{Tuple{Int,Int}}}, weights::Vector{Int})
    nf = length(adj)
    nf <= 1 && return Int[]
    parent = zeros(Int, nf)
    parent_dart = zeros(Int, nf)
    order = Int[1]
    seen = falses(nf)
    seen[1] = true
    for f in order
        for (g, d) in adj[f]
            seen[g] && continue
            seen[g] = true
            parent[g] = f
            parent_dart[g] = d
            push!(order, g)
        end
    end
    length(order) == nf || return Int[]
    subtree = copy(weights)
    for i = length(order):-1:2
        f = order[i]
        subtree[parent[f]] += subtree[f]
    end
    total = sum(weights)
    ranked = Tuple{Int,Int}[]
    for i = 2:length(order)
        f = order[i]
        s = subtree[f]
        push!(ranked, (max(s, total - s), parent_dart[f]))
    end
    sort!(ranked; by = first)
    return Int[d for (_, d) in ranked]
end

function _candidate_cycles(
    h::DartGraph,
    adj,
    weights::Vector{Int},
    levels::LevelStructure,
    roots::Vector{Int},
)
    candidates = Vector{Vector{Int}}()
    for d in _balanced_darts(adj, weights)
        x = tail(h, d)
        z = h.head[d]
        roots[x] == 0 && continue
        if roots[x] == roots[z]
            push!(candidates, _tree_path(levels.parent, levels.depth, x, z))
        else
            # The fundamental cycle passes through the shrunk levels, which are
            # part of the separator already: keep only the two tree paths.
            push!(
                candidates,
                vcat(
                    _up_path(levels.parent, x, roots[x]),
                    _up_path(levels.parent, z, roots[z]),
                ),
            )
        end
        length(candidates) >= MAX_CYCLE_CANDIDATES && break
    end
    return candidates
end

# Path from `v` up to its tree root `r`, both included.
function _up_path(parent::Vector{Int}, v::Int, r::Int)
    path = Int[v]
    while path[end] != r
        push!(path, parent[path[end]])
    end
    return path
end

# Tree path from x to z (both in the same tree), given parent pointers and
# depths.
function _tree_path(parent::Vector{Int}, depth::Vector{Int}, x::Int, z::Int)
    left = Int[]
    right = Int[]
    a, b = x, z
    while depth[a] > depth[b]
        push!(left, a)
        a = parent[a]
    end
    while depth[b] > depth[a]
        push!(right, b)
        b = parent[b]
    end
    while a != b
        push!(left, a)
        push!(right, b)
        a = parent[a]
        b = parent[b]
    end
    push!(left, a)
    return vcat(left, reverse(right))
end

# --- final partition -----------------------------------------------------

function _merge_vertices(base::Vector{Int}, extra::Vector{Int})
    isempty(extra) && return sort(copy(base))
    merged = vcat(base, extra)
    sort!(merged)
    unique!(merged)
    return merged
end

# Partition `vertices(g) \\ separator` into connected components, then put the
# components into two bins of nearly equal size.
function _split_and_pack(dg::DartGraph, separator::Vector{Int})
    comps = components_without(dg, _mask(separator, dg.n))
    sort!(comps; by = length, rev = true)
    a = Int[]
    b = Int[]
    size_a = 0
    size_b = 0
    for comp in comps
        if size_a <= size_b
            append!(a, comp)
            size_a += length(comp)
        else
            append!(b, comp)
            size_b += length(comp)
        end
    end
    return sort(a), sort(b), sort(separator)
end

_is_balanced(candidate, n::Integer, balance::Real) =
    max(length(candidate[1]), length(candidate[2])) <= balance * n

# --- disconnected case (Theorem 4, end of the proof) ---------------------

function _separator_disconnected(
    dg::DartGraph,
    label::Vector{Int},
    ncomps::Int,
    balance::Real,
)
    n = dg.n
    comps = [Int[] for _ = 1:ncomps]
    for v = 1:n
        push!(comps[label[v]], v)
    end
    sort!(comps; by = length, rev = true)
    biggest = comps[1]
    length(biggest) > balance * n && return _recurse_into_component(dg, comps, balance)
    length(biggest) >= n / 3 && return (sort(biggest), _flatten(comps[2:end]), Int[])
    return _pack_components(comps, n)
end

function _recurse_into_component(dg::DartGraph, comps, balance::Real)
    sub, ids = induced(dg, _mask(comps[1], dg.n))
    a, b, c = _best_connected(sub, balance)
    part_a = sort(Int[ids[v] for v in a])
    part_b = sort(Int[ids[v] for v in b])
    rest = _flatten(comps[2:end])
    if length(part_a) >= length(part_b)
        return part_a, sort(vcat(part_b, rest)), sort(Int[ids[v] for v in c])
    end
    return part_b, sort(vcat(part_a, rest)), sort(Int[ids[v] for v in c])
end

function _pack_components(comps, n::Integer)
    threshold = n / 3
    cut = length(comps)
    cumulative = 0
    for (i, comp) in enumerate(comps)
        cumulative += length(comp)
        if cumulative > threshold
            cut = i
            break
        end
    end
    return sort(_flatten(comps[1:cut])), sort(_flatten(comps[(cut+1):end])), Int[]
end

function _flatten(comps)
    out = Int[]
    for comp in comps
        append!(out, comp)
    end
    return out
end
