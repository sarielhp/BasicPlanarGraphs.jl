# Lipton-Tarjan planar separator.
#
# Reference: R. J. Lipton and R. E. Tarjan, "A Separator Theorem for Planar
# Graphs", SIAM J. Appl. Math. 36(2), 1979. Lemma 2 (a balanced fundamental
# cycle), Lemma 3 (two levels plus the middle band) and Theorem 4 (the
# algorithm).
#
# The implementation follows the proof, so both guarantees are worst case and
# deterministic, with no retries and no refinement loop:
#
#   * every component left by the separator has at most `2n/3` vertices, and
#   * the separator before pruning has at most `2*sqrt(2n)` vertices.
#
# Running time is `O(n)`: every stage is a pass over vertices, darts or the
# triangles of the virtual triangulation, and the search for the balanced cycle
# touches each triangle a constant number of times (see `_search_edge`).
#
# Lemma 2 needs a triangulation. Instead of adding edges to the graph, every
# face with four or more sides gets a virtual *star* vertex joined to all its
# corners. This works for any connected plane graph, including graphs with
# bridges and cut vertices, and it keeps the BFS tree of the input graph as the
# spanning tree: each star hangs off the tree as a leaf. The triangles of the
# result are the nodes of the dual tree that the search walks.

const DEFAULT_BALANCE = 2 / 3

"""
    separator(g; balance = 2/3) -> (A, B, C)

Balanced vertex separator of the embedded planar graph `g`, computed by the
Lipton-Tarjan algorithm in `O(nv(g))` time.

The three returned vectors partition `vertices(g)`: no edge joins `A` to `B`,
`|A|` and `|B|` are at most `balance * nv(g)`, and `|C| <= 2*sqrt(2*nv(g))`.
All three hold for every input, not only for well-behaved ones. `balance` must
lie in `[2/3, 1)`: `2/3` is what the theorem guarantees, and a larger value only
lets the final pruning step move more separator vertices into `A` and `B`.

```julia
g = planar_grid(20, 20)
A, B, C = separator(g)
isempty(intersect(A, B))            # true
length(C) <= 2 * sqrt(2 * nv(g))    # true
```
"""
function separator(g::PlanarGraph; balance::Real = DEFAULT_BALANCE)
    2 / 3 <= balance < 1 ||
        throw(ArgumentError("balance must lie in [2/3, 1), got $balance"))
    n = nv(g)
    n == 0 && return (Int[], Int[], Int[])
    cap = max(floor(Int, balance * n), (2 * n) ÷ 3)
    return _unpack(_best_partition(_planar_dart_graph(g), cap))
end

# The darts of `g`, after a check that every component is planarly embedded. The
# check has to see the whole graph: the search itself cuts only the largest
# component, so a rotation system of higher genus elsewhere would go unnoticed.
function _planar_dart_graph(g::PlanarGraph)
    dg = DartGraph(g)
    euler_holds(dg) || throw(_not_planar_error())
    return dg
end

# --- candidates and partitions -------------------------------------------

# Every root gives a valid separator; the smallest result after packing and
# pruning wins. The number of roots is a constant, so the work stays linear.
function _best_partition(dg::DartGraph, cap::Int)
    label, ncomps = component_labels(dg)
    candidates =
        ncomps == 1 ? _connected_separators(dg) : _big_component_separators(dg, label)
    best = Int8[]
    best_size = typemax(Int)
    for separator in candidates
        side = _partition(dg, separator, cap)
        size = count(iszero, side)
        if size < best_size
            best = side
            best_size = size
        end
    end
    return best
end

# Only a component with more than `2n/3` vertices has to be cut; the others
# already fit, and the packing step distributes them.
function _big_component_separators(dg::DartGraph, label::Vector{Int})
    sizes = zeros(Int, maximum(label))
    for l in label
        sizes[l] += 1
    end
    big = argmax(sizes)
    3 * sizes[big] <= 2 * dg.n && return [Int[]]
    sub, ids = induced(dg, label .== big)
    return [Int[ids[v] for v in sep] for sep in _connected_separators(sub)]
end

function _unpack(side::Vector{Int8})
    return findall(==(1), side), findall(==(2), side), findall(iszero, side)
end

# --- connected case ------------------------------------------------------

# The level structure, and with it the size of the separator, depends on the
# BFS root, so a few roots are tried: the endpoints and the midpoint of a BFS
# pseudo-diameter, which is the centroid on tree-like graphs.
function _connected_separators(dg::DartGraph)
    dg.n <= 2 && return [collect(1:dg.n)]
    tri = _triangles(dg)
    return [_level_separator(dg, tri, root) for root in _root_candidates(dg)]
end

function _root_candidates(dg::DartGraph)
    far = _bfs(dg, 1).order[end]
    mid = _pseudo_center(dg, far)
    return unique(Int[1, far, mid, _bfs(dg, mid).order[end]])
end

# Midpoint of a BFS path that starts at `far`, a farthest vertex: the centroid on
# tree-like graphs.
function _pseudo_center(dg::DartGraph, far::Int)
    sweep = _bfs(dg, far)
    end_vertex = sweep.order[end]
    return _ancestor(sweep.parent, end_vertex, sweep.depth[end_vertex] ÷ 2)
end

function _ancestor(parent::Vector{Int}, v::Int, steps::Int)
    for _ = 1:steps
        v = parent[v]
    end
    return v
end

# Breadth-first search from `root`; the graph is connected, so `order` lists
# every vertex.
function _bfs(dg::DartGraph, root::Int)
    parent = zeros(Int, dg.n)
    parent_dart = zeros(Int, dg.n)
    depth = fill(-1, dg.n)
    depth[root] = 0
    order = Int[root]
    sizehint!(order, dg.n)
    for v in order
        for d in incident_darts(dg, v)
            w = dg.head[d]
            depth[w] == -1 || continue
            depth[w] = depth[v] + 1
            parent[w] = v
            parent_dart[w] = d
            push!(order, w)
        end
    end
    return (; parent, parent_dart, depth, order)
end

# Theorem 4: the vertices of levels `l0` and `l2`, plus the band vertices of a
# balanced cycle when the band between them is larger than `2n/3`.
function _level_separator(dg::DartGraph, tri, root::Int)
    n = dg.n
    tree = _bfs(dg, root)
    l0, l2 = _level_plan(_level_counts(tree.depth), n)
    base = Int[v for v = 1:n if tree.depth[v] == l0 || tree.depth[v] == l2]
    band = BitVector(l0 < tree.depth[v] < l2 for v = 1:n)
    width = count(band)
    3 * width <= 2 * n && return base
    return vcat(base, _band_cycle(dg, tri, tree, band, width))
end

# --- level bookkeeping ---------------------------------------------------

# counts[l + 1] = number of vertices on BFS level l.
function _level_counts(depth::Vector{Int})
    counts = zeros(Int, maximum(depth) + 1)
    for lv in depth
        counts[lv+1] += 1
    end
    return counts
end

# Levels l0 and l2 of Theorem 4. Level `l1` is where the prefix sum first
# reaches `n/2`, so the levels below `l0` and above `l2` each hold at most half
# of the vertices.
function _level_plan(counts::Vector{Int}, n::Int)
    l1 = _middle_level(counts, n)
    k = sum(view(counts, 1:(l1+1)))
    return _lower_level(counts, l1, k), _upper_level(counts, l1, n - k)
end

function _level_size(counts::Vector{Int}, l::Int)
    return 0 <= l < length(counts) ? counts[l+1] : 0
end

function _middle_level(counts::Vector{Int}, n::Int)
    acc = 0
    for l = 0:(length(counts)-1)
        acc += counts[l+1]
        2 * acc >= n && return l
    end
    return length(counts) - 1
end

# Level l0 <= l1 with |L(l0)| + 2(l1 - l0) <= 2 sqrt(k). Squared, so exact.
function _lower_level(counts::Vector{Int}, l1::Int, k::Int)
    l = l1
    while l >= 0 && (_level_size(counts, l) + 2 * (l1 - l))^2 > 4 * k
        l -= 1
    end
    return l
end

# Level l2 >= l1 + 1 with |L(l2)| + 2(l2 - l1 - 1) <= 2 sqrt(n - k).
function _upper_level(counts::Vector{Int}, l1::Int, rest::Int)
    last = length(counts) - 1
    l = l1 + 1
    while l <= last && (_level_size(counts, l) + 2 * (l - l1 - 1))^2 > 4 * rest
        l += 1
    end
    return l
end

# --- the virtual triangulation -------------------------------------------

# Triangles of the graph in which every face with four or more sides is
# replaced by a fan around a star vertex. A face with `k >= 4` darts
# `d_1 .. d_k` (corners `a_j = tail(d_j)`) becomes the triangles
# `(s, a_j, a_{j+1})`, numbered `tbase + j - 1`; a triangular face stays one
# triangle. Stars have no cost, so a triangle lists only its real vertices
# (`vc == 0` for the triangles of a star).
struct _Triangles
    ntri::Int
    orbits::Vector{Vector{Int}}
    tbase::Vector{Int}
    dart_tri::Vector{Int}   # the triangle on the left of each dart
    deg::Vector{Int}        # number of triangles at each real vertex
    va::Vector{Int}
    vb::Vector{Int}
    vc::Vector{Int}
end

function _triangles(dg::DartGraph)
    orbits, _ = face_orbits(dg)
    dg.n - nedges(dg) + length(orbits) == 2 || throw(_not_planar_error())
    tbase = Vector{Int}(undef, length(orbits))
    ntri = 0
    for (f, orbit) in enumerate(orbits)
        tbase[f] = ntri + 1
        ntri += length(orbit) == 3 ? 1 : length(orbit)
    end
    tri = _Triangles(
        ntri,
        orbits,
        tbase,
        zeros(Int, ndarts(dg)),
        zeros(Int, dg.n),
        zeros(Int, ntri),
        zeros(Int, ntri),
        zeros(Int, ntri),
    )
    for f in eachindex(orbits)
        _fill_face!(tri, dg, f)
    end
    return tri
end

function _fill_face!(tri::_Triangles, dg::DartGraph, f::Int)
    orbit = tri.orbits[f]
    k = length(orbit)
    k >= 3 || throw(ArgumentError("a face with fewer than three sides"))
    corner(j) = dg.tail[orbit[mod1(j, k)]]
    for j = 1:k
        tri.deg[corner(j)] += k == 3 ? 1 : 2
        tri.dart_tri[orbit[j]] = k == 3 ? tri.tbase[f] : tri.tbase[f] + j - 1
    end
    for j = 1:(k==3 ? 1 : k)
        t = tri.tbase[f] + j - 1
        tri.va[t] = corner(j)
        tri.vb[t] = corner(j + 1)
        k == 3 && (tri.vc[t] = corner(3))
    end
    return nothing
end

# --- the dual tree -------------------------------------------------------

# The edges of the virtual triangulation that are not in the spanning tree are
# the edges of a spanning tree of the dual (the triangles). A triangle has three
# sides, so it has at most three dual neighbours; `ep[e], eq[e]` are the real
# vertices whose tree path closes the fundamental cycle of dual edge `e`.
struct _Dual
    nbr::Vector{Int}
    eid::Vector{Int}
    cnt::Vector{Int}
    ep::Vector{Int}
    eq::Vector{Int}
    ef::Vector{Int}   # face of a star edge, 0 for an edge of the graph
    ej::Vector{Int}   # corner of that face the star edge ends at
end

# `ends = (p, q, f, j)`: the real vertices of the cycle of the new dual edge and
# the face and corner of its star edge, `(0, 0)` for an edge of the graph.
function _link!(dual::_Dual, t1::Int, t2::Int, ends::NTuple{4,Int})
    push!(dual.ep, ends[1])
    push!(dual.eq, ends[2])
    push!(dual.ef, ends[3])
    push!(dual.ej, ends[4])
    for (a, b) in ((t1, t2), (t2, t1))
        slot = 3 * (a - 1) + (dual.cnt[a] += 1)
        dual.nbr[slot] = b
        dual.eid[slot] = length(dual.ep)
    end
    return nothing
end

function _dual_tree(dg::DartGraph, tri::_Triangles, parent_dart::Vector{Int})
    nt = tri.ntri
    dual =
        _Dual(zeros(Int, 3nt), zeros(Int, 3nt), zeros(Int, nt), Int[], Int[], Int[], Int[])
    in_tree = falses(ndarts(dg))
    for d in parent_dart
        d == 0 && continue
        in_tree[d] = true
        in_tree[dg.twin[d]] = true
    end
    for d = 1:ndarts(dg)
        e = dg.twin[d]
        (d < e && !in_tree[d]) || continue
        _link!(dual, tri.dart_tri[d], tri.dart_tri[e], (dg.tail[d], dg.head[d], 0, 0))
    end
    _link_stars!(dual, dg, tri)
    length(dual.ep) == nt - 1 || throw(_not_planar_error())
    return dual
end

# The first star edge of every star is its tree edge; the others join the
# consecutive triangles of the fan, and their cycle closes through the tree
# path from the first corner to the corner of that star edge.
function _link_stars!(dual::_Dual, dg::DartGraph, tri::_Triangles)
    for (f, orbit) in enumerate(tri.orbits)
        length(orbit) == 3 && continue
        first_corner = dg.tail[orbit[1]]
        for j = 2:length(orbit)
            t = tri.tbase[f] + j - 1
            _link!(dual, t - 1, t, (first_corner, dg.tail[orbit[j]], f, j))
        end
    end
    return nothing
end

# --- counting the vertices on each side of a dual edge -------------------

# The triangles of a region of the dual tree, as seen by the band vertices: how
# many triangles of each vertex lie in the region. A vertex is strictly inside
# when all its triangles do (`full`), and outside when none does. Stamps make a
# reset free.
mutable struct _Region
    count::Vector{Int}
    stamp::Vector{Int}
    epoch::Int
    full::Int
    touched::Int
end

_Region(n::Int) = _Region(zeros(Int, n), zeros(Int, n), 1, 0, 0)

function _reset!(r::_Region)
    r.epoch += 1
    r.full = 0
    r.touched = 0
    return r
end

function _enter!(r::_Region, v::Int, deg::Vector{Int})
    if r.stamp[v] != r.epoch
        r.stamp[v] = r.epoch
        r.count[v] = 0
    end
    c = r.count[v] + 1
    r.count[v] = c
    c == 1 && (r.touched += 1)
    c == deg[v] && (r.full += 1)
    return nothing
end

function _leave!(r::_Region, v::Int, deg::Vector{Int})
    c = r.count[v]
    c == deg[v] && (r.full -= 1)
    r.count[v] = c - 1
    c == 1 && (r.touched -= 1)
    return nothing
end

mutable struct _Search
    tri::_Triangles
    dual::_Dual
    band::BitVector
    width::Int
    parent::Vector{Int}
    size::Vector{Int}
    stack::Vector{Int}
    region::_Region
    scratch::_Region
end

function _apply_triangle!(op, r::_Region, s::_Search, t::Int)
    for v in (s.tri.va[t], s.tri.vb[t], s.tri.vc[t])
        (v != 0 && s.band[v]) && op(r, v, s.tri.deg)
    end
    return nothing
end

# Apply `op` to the triangles of the subtree of dual node `t`.
function _each_below!(op, r::_Region, s::_Search, t::Int)
    empty!(s.stack)
    push!(s.stack, t)
    while !isempty(s.stack)
        u = pop!(s.stack)
        _apply_triangle!(op, r, s, u)
        for i = 1:s.dual.cnt[u]
            w = s.dual.nbr[3*(u-1)+i]
            w == s.parent[u] || push!(s.stack, w)
        end
    end
    return nothing
end

# Parents and subtree sizes of the dual tree rooted at triangle 1.
function _root_dual(dual::_Dual, nt::Int)
    parent = zeros(Int, nt)
    order = Int[1]
    sizehint!(order, nt)
    for u in order
        for i = 1:dual.cnt[u]
            w = dual.nbr[3*(u-1)+i]
            w == parent[u] && continue
            parent[w] = u
            push!(order, w)
        end
    end
    length(order) == nt || throw(_not_planar_error())
    size = ones(Int, nt)
    for i = nt:-1:2
        size[parent[order[i]]] += size[order[i]]
    end
    return parent, size
end

# --- the balanced cycle (Lemma 2) ----------------------------------------

# Band vertices on the fundamental cycle of the balanced dual edge, that is on
# the BFS tree path between the endpoints `p` and `q`.
function _band_cycle(dg::DartGraph, tri::_Triangles, tree, band::BitVector, width::Int)
    dual = _dual_tree(dg, tri, tree.parent_dart)
    parent, size = _root_dual(dual, tri.ntri)
    search =
        _Search(tri, dual, band, width, parent, size, Int[], _Region(dg.n), _Region(dg.n))
    e = _search_edge(search)
    return _band_path(tree, band, dual.ep[e], dual.eq[e])
end

function _band_path(tree, band::BitVector, p::Int, q::Int)
    path = Int[]
    a, b = p, q
    while a != b
        if tree.depth[a] >= tree.depth[b]
            band[a] && push!(path, a)
            a = tree.parent[a]
        else
            band[b] && push!(path, b)
            b = tree.parent[b]
        end
    end
    band[a] && push!(path, a)
    return path
end

# A dual edge is balanced when at most `2/3` of the band vertices lie strictly
# on either side. Exactly one side can hold more than `2/3`, so every bad edge
# points at its heavy side. Walk down the dual tree from the root along the
# edges whose heavy side is below; at the end of the walk every incident bad
# edge points at the current triangle, and then the three (at most) sides
# cannot all be heavy, because a vertex is strictly inside at most `k - 1` of
# `k` sides. So some child edge is balanced.
#
# Each step measures the children: the small ones from scratch, and the largest
# one by removing the current triangle and the small subtrees from the region
# it already holds. The work of a step is bounded by the size of the small
# subtrees, which shrink the region by a constant factor or are removed from it,
# so the whole walk is linear.
function _search_edge(s::_Search)
    _each_below!(_enter!, s.region, s, 1)
    cur = 1
    while true
        e, next = _descend!(s, cur)
        e != 0 && return e
        next != 0 || error("no balanced dual edge: the embedding is not planar")
        cur = next
    end
end

function _kids(s::_Search, cur::Int)
    kids = Tuple{Int,Int}[]
    for i = 1:s.dual.cnt[cur]
        w = s.dual.nbr[3*(cur-1)+i]
        w == s.parent[cur] || push!(kids, (w, s.dual.eid[3*(cur-1)+i]))
    end
    return kids
end

_balanced(s::_Search, r::_Region) =
    3 * r.full <= 2 * s.width && 3 * (s.width - r.touched) <= 2 * s.width

_heavy(s::_Search, r::_Region) = 3 * r.full > 2 * s.width

# Returns `(edge, 0)` when a child edge of `cur` is balanced, else `(0, next)`
# with the region of `s` rebuilt for the child that holds the heavy side.
function _descend!(s::_Search, cur::Int)
    kids = _kids(s, cur)
    isempty(kids) && error("no balanced dual edge: the embedding is not planar")
    big = kids[argmax([s.size[c] for (c, _) in kids])]
    next = 0
    for (c, e) in kids
        c == big[1] && continue
        _reset!(s.scratch)
        _each_below!(_enter!, s.scratch, s, c)
        _balanced(s, s.scratch) && return (e, 0)
        _heavy(s, s.scratch) && (next = c)
    end
    _peel!(s, cur, kids, big[1])
    _balanced(s, s.region) && return (big[2], 0)
    _heavy(s, s.region) && return (0, big[1])
    next == 0 && return (0, 0)
    _reset!(s.region)
    _each_below!(_enter!, s.region, s, next)
    return (0, next)
end

# Leave the region of `cur` with the subtree of its child `big` only.
function _peel!(s::_Search, cur::Int, kids, big::Int)
    _apply_triangle!(_leave!, s.region, s, cur)
    for (c, _) in kids
        c == big || _each_below!(_leave!, s.region, s, c)
    end
    return nothing
end

# --- packing the components and pruning the separator ---------------------

# `side[v]`: 0 for the separator, 1 for `A`, 2 for `B`.
function _partition(dg::DartGraph, separator::Vector{Int}, cap::Int)
    comp, sizes = _components_without(dg, separator)
    bins = _bins(sizes, dg.n)
    side = Int8[c == 0 ? Int8(0) : bins[c] for c in comp]
    _prune!(dg, side, cap)
    max(count(==(1), side), count(==(2), side)) <= cap ||
        error("internal error: the separator leaves a side above $cap vertices")
    return side
end

# Component labels of the graph without `separator` (0 on the separator) and
# the component sizes.
function _components_without(dg::DartGraph, separator::Vector{Int})
    comp = zeros(Int, dg.n)
    for v in separator
        comp[v] = -1
    end
    sizes = Int[]
    stack = Int[]
    for s = 1:dg.n
        comp[s] == 0 || continue
        push!(sizes, _flood!(comp, stack, dg, s, length(sizes) + 1))
    end
    for v in separator
        comp[v] = 0
    end
    return comp, sizes
end

function _flood!(comp::Vector{Int}, stack::Vector{Int}, dg::DartGraph, s::Int, id::Int)
    comp[s] = id
    push!(stack, s)
    size = 0
    while !isempty(stack)
        v = pop!(stack)
        size += 1
        for d in incident_darts(dg, v)
            w = dg.head[d]
            comp[w] == 0 || continue
            comp[w] = id
            push!(stack, w)
        end
    end
    return size
end

# Two bins, each of at most `2n/3` vertices, for components that are all at most
# `2n/3`: a component of at least `n/3` goes alone into the first bin and the
# rest into the second; otherwise components fill the first bin until it passes
# `n/3`, which leaves it below `2n/3`, and the second bin has the rest.
function _bins(sizes::Vector{Int}, n::Int)
    bins = fill(Int8(2), length(sizes))
    isempty(sizes) && return bins
    big = argmax(sizes)
    if 3 * sizes[big] >= n
        bins[big] = 1
        return bins
    end
    cum = 0
    for c in eachindex(sizes)
        bins[c] = 1
        cum += sizes[c]
        3 * cum > n && break
    end
    return bins
end

# A separator vertex with no neighbour in one side can join the other side,
# while that side stays within `cap`. One pass reaches a fixed point, because a
# vertex that joins a side can only block, never free, the others.
function _prune!(dg::DartGraph, side::Vector{Int8}, cap::Int)
    sa = count(==(1), side)
    sb = count(==(2), side)
    for v = 1:dg.n
        side[v] == 0 || continue
        has_a, has_b = _neighbor_sides(dg, side, v)
        can_a = !has_b && sa < cap
        can_b = !has_a && sb < cap
        if can_a && (!can_b || sa <= sb)
            side[v] = 1
            sa += 1
        elseif can_b
            side[v] = 2
            sb += 1
        end
    end
    return side
end

function _neighbor_sides(dg::DartGraph, side::Vector{Int8}, v::Int)
    has_a = false
    has_b = false
    for d in incident_darts(dg, v)
        s = side[dg.head[d]]
        has_a |= s == 1
        has_b |= s == 2
    end
    return has_a, has_b
end
