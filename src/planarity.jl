# Linear-time planarity testing and embedding.
#
# The algorithm is the Left-Right planarity criterion of Ulrik Brandes, *The
# Left-Right Planarity Test* (2009): a DFS orients the graph, the nesting graph
# of its lowpoints is built, a stack of conflict intervals is refined in a
# second DFS, and the resulting left/right partition of the back edges yields
# the rotation system. Every stage is a pass over vertices and edges, so the
# running time is O(n + m).
#
# The translation follows the reference implementation in NetworkX
# (`networkx/algorithms/planarity.py`, BSD-3-Clause), which is the published
# algorithm with the iterative formulations that avoid deep recursion.
#
# Edges are keyed by oriented pairs `(tail, head)` in `Dict`s; the graph has at
# most `3n - 6` edges when it is planar, so the dictionaries stay linear.

# --- conflict intervals --------------------------------------------------

mutable struct _Interval
    low::Union{Nothing,Tuple{Int,Int}}
    high::Union{Nothing,Tuple{Int,Int}}
end

_Interval() = _Interval(nothing, nothing)

_interval_empty(i::_Interval) = i.low === nothing && i.high === nothing

_copy_interval(i::_Interval) = _Interval(i.low, i.high)

# A return edge in an interval constrains all edges inside it to the same side.
function _conflicting(i::_Interval, b::Tuple{Int,Int}, lowpt)
    _interval_empty(i) && return false
    high = i.high::Tuple{Int,Int}
    return lowpt[high] > lowpt[b]
end

mutable struct _ConflictPair
    left::_Interval
    right::_Interval
end

_ConflictPair() = _ConflictPair(_Interval(), _Interval())

function _swap_pair!(p::_ConflictPair)
    p.left, p.right = p.right, p.left
    return p
end

function _lowest(p::_ConflictPair, lowpt)
    _interval_empty(p.left) && return lowpt[p.right.low::Tuple{Int,Int}]
    _interval_empty(p.right) && return lowpt[p.left.low::Tuple{Int,Int}]
    return min(lowpt[p.left.low::Tuple{Int,Int}], lowpt[p.right.low::Tuple{Int,Int}])
end

_stack_top(stack) = isempty(stack) ? nothing : stack[end]

# The lowest lowpoint on the conflict stack, or `nothing` when it is empty.
function _stack_lowest(st)
    isempty(st.stack) && return nothing
    return _lowest(last(st.stack), st.lowpt)
end

# Whether the return edge `d` ends at `u`; `nothing` points nowhere.
_points_to(d::Nothing, u::Int) = false
_points_to(d::Tuple{Int,Int}, u::Int) = d[2] == u

# Where the constraint that starts at `d` continues; a missing entry is the end.
_reference(st, d::Nothing) = nothing
_reference(st, d::Tuple{Int,Int}) = get(st.ref, d, nothing)

# --- embedding under construction ---------------------------------------

# Half-edge lists in both directions around every vertex: `cw[v][w]` is the
# neighbour after `w` clockwise, `ccw[v][w]` the one before it. `leftmost[v]` is
# the neighbour that `add_half_edge_first` marks, which is the insertion point
# for the next tree edge.
mutable struct _HalfEdges
    cw::Vector{Dict{Int,Int}}
    ccw::Vector{Dict{Int,Int}}
    leftmost::Vector{Int}
end

_HalfEdges(n::Integer) =
    _HalfEdges([Dict{Int,Int}() for _ = 1:n], [Dict{Int,Int}() for _ = 1:n], zeros(Int, n))

function _add_half_edge!(
    h::_HalfEdges,
    v::Int,
    w::Int;
    cw_ref::Union{Nothing,Int} = nothing,
    ccw_ref::Union{Nothing,Int} = nothing,
)
    if isempty(h.cw[v])
        h.cw[v][w] = w
        h.ccw[v][w] = w
        h.leftmost[v] = w
        return h
    end

    if cw_ref !== nothing
        ref_ccw = h.ccw[v][cw_ref]
        h.cw[v][w] = cw_ref
        h.ccw[v][w] = ref_ccw
        h.cw[v][ref_ccw] = w
        h.ccw[v][cw_ref] = w
        return h
    end

    ccw_ref === nothing && error("_add_half_edge!: a reference neighbour is required")
    ref_cw = h.cw[v][ccw_ref]
    h.cw[v][w] = ref_cw
    h.ccw[v][w] = ccw_ref
    h.ccw[v][ref_cw] = w
    h.cw[v][ccw_ref] = w
    return h
end

function _add_half_edge_first!(h::_HalfEdges, v::Int, w::Int)
    ref = isempty(h.cw[v]) ? nothing : h.leftmost[v]
    _add_half_edge!(h, v, w; cw_ref = ref)
    h.leftmost[v] = w
    return h
end

# Neighbours of `v` in clockwise order, read off the half-edge lists.
function _rotation_of(h::_HalfEdges, v::Integer)
    isempty(h.cw[v]) && return Int[]
    start = h.leftmost[v]
    order = Int[start]
    x = h.cw[v][start]
    while x != start
        push!(order, x)
        x = h.cw[v][x]
    end
    return order
end

# --- algorithm state -----------------------------------------------------

mutable struct _LRState
    adj::Vector{Vector{Int}}
    height::Vector{Int}
    parent_edge::Vector{Union{Nothing,Tuple{Int,Int}}}
    out_adjs::Vector{Vector{Int}}
    oriented::Set{Tuple{Int,Int}}
    lowpt::Dict{Tuple{Int,Int},Int}
    lowpt2::Dict{Tuple{Int,Int},Int}
    nesting_depth::Dict{Tuple{Int,Int},Int}
    ordered_adjs::Vector{Vector{Int}}
    ref::Dict{Tuple{Int,Int},Tuple{Int,Int}}
    side::Dict{Tuple{Int,Int},Int}
    stack::Vector{_ConflictPair}
    stack_bottom::Dict{Tuple{Int,Int},Union{Nothing,_ConflictPair}}
    lowpt_edge::Dict{Tuple{Int,Int},Tuple{Int,Int}}
    left_ref::Vector{Int}
    right_ref::Vector{Int}
    halves::_HalfEdges
end

function _LRState(adj::Vector{Vector{Int}})
    n = length(adj)
    return _LRState(
        adj,
        fill(-1, n),
        fill(nothing, n),
        [Int[] for _ = 1:n],
        Set{Tuple{Int,Int}}(),
        Dict{Tuple{Int,Int},Int}(),
        Dict{Tuple{Int,Int},Int}(),
        Dict{Tuple{Int,Int},Int}(),
        [Int[] for _ = 1:n],
        Dict{Tuple{Int,Int},Tuple{Int,Int}}(),
        Dict{Tuple{Int,Int},Int}(),
        _ConflictPair[],
        Dict{Tuple{Int,Int},Union{Nothing,_ConflictPair}}(),
        Dict{Tuple{Int,Int},Tuple{Int,Int}}(),
        zeros(Int, n),
        zeros(Int, n),
        _HalfEdges(n),
    )
end

_has_orientation(st::_LRState, v::Int, w::Int) =
    (v, w) in st.oriented || (w, v) in st.oriented

# --- phase 1: orientation and lowpoints ----------------------------------

function _update_lowpoints!(
    st::_LRState,
    vw::Tuple{Int,Int},
    e::Union{Nothing,Tuple{Int,Int}},
)
    e === nothing && return nothing
    low, second = st.lowpt[vw], st.lowpt2[vw]
    if low < st.lowpt[e]
        st.lowpt2[e] = min(st.lowpt[e], second)
        st.lowpt[e] = low
    elseif low > st.lowpt[e]
        st.lowpt2[e] = min(st.lowpt2[e], low)
    else
        st.lowpt2[e] = min(st.lowpt2[e], second)
    end
    return nothing
end

# What a DFS edge turned out to be: already oriented (ignore it), a back edge to
# an earlier vertex, or a tree edge, which the caller must descend into.
const _ORIENTED = 0
const _BACK_EDGE = 1
const _TREE_EDGE = 2

function _dfs_orientation!(st::_LRState, root::Int)
    ind = zeros(Int, length(st.adj))
    skip_init = Set{Tuple{Int,Int}}()
    dfs_stack = Int[root]

    while !isempty(dfs_stack)
        v = pop!(dfs_stack)
        e = st.parent_edge[v]
        while ind[v] < length(st.adj[v])
            w = st.adj[v][ind[v]+1]
            vw = (v, w)
            kind = vw in skip_init ? _BACK_EDGE : _orient_edge!(st, v, w, vw)
            if kind == _TREE_EDGE
                push!(dfs_stack, v)
                push!(dfs_stack, w)
                push!(skip_init, vw)
                break
            end
            kind == _BACK_EDGE && _nest_and_lowpoint!(st, v, vw, e)
            ind[v] += 1
        end
    end
    return st
end

# Orient `v -> w` and set its lowpoints; a new vertex is a tree edge.
function _orient_edge!(st::_LRState, v::Int, w::Int, vw::Tuple{Int,Int})
    _has_orientation(st, v, w) && return _ORIENTED
    push!(st.oriented, vw)
    push!(st.out_adjs[v], w)
    st.lowpt[vw] = st.height[v]
    st.lowpt2[vw] = st.height[v]
    if st.height[w] == -1
        st.parent_edge[w] = vw
        st.height[w] = st.height[v] + 1
        return _TREE_EDGE
    end
    st.lowpt[vw] = st.height[w]
    return _BACK_EDGE
end

# The nesting depth of `vw`, and what it does to the lowpoints of its parent.
function _nest_and_lowpoint!(st::_LRState, v::Int, vw::Tuple{Int,Int}, e)
    chordal = st.lowpt2[vw] < st.height[v] ? 1 : 0
    st.nesting_depth[vw] = 2 * st.lowpt[vw] + chordal
    _update_lowpoints!(st, vw, e)
    return nothing
end

# --- phase 2: the LR test ------------------------------------------------

# `ref[x] = nothing` is the reference implementation's way of forgetting the
# constraint that pointed from `x` onwards; a missing entry means the same.
_set_ref!(st, low::Tuple{Int,Int}, high::Nothing) = (delete!(st.ref, low); st)
_set_ref!(st, low::Tuple{Int,Int}, high::Tuple{Int,Int}) = (st.ref[low] = high; st)

# Point the constraint that starts where `interval` starts, if it starts at all.
_chain_ref!(st, interval::_Interval, high) = _chain_ref!(st, interval.low, high)
_chain_ref!(st, low::Nothing, high) = st
_chain_ref!(st, low::Tuple{Int,Int}, high) = _set_ref!(st, low, high)

function _add_constraints!(st::_LRState, ei::Tuple{Int,Int}, parent::Tuple{Int,Int})
    pair = _ConflictPair()
    _merge_return_edges!(st, pair, ei, parent) || return false
    _merge_conflict_edges!(st, pair, ei) || return false
    _interval_empty(pair.left) && _interval_empty(pair.right) || push!(st.stack, pair)
    return true
end

# Merge the intervals above the stack bottom of `ei` into `pair.right`.
function _merge_return_edges!(st::_LRState, pair::_ConflictPair, ei, parent)
    while true
        isempty(st.stack) && return false
        q = pop!(st.stack)
        _interval_empty(q.left) || _swap_pair!(q)
        _interval_empty(q.left) || return false
        _merge_below!(st, pair, q, parent)
        _stack_top(st.stack) == st.stack_bottom[ei] && break
    end
    return true
end

# One interval below the lowpoint of the parent edge: merge it, or align it.
function _merge_below!(st::_LRState, pair::_ConflictPair, q::_ConflictPair, parent)
    if st.lowpt[q.right.low::Tuple{Int,Int}] > st.lowpt[parent]
        if _interval_empty(pair.right)
            pair.right = _copy_interval(q.right)
        else
            _chain_ref!(st, pair.right, q.right.high)
        end
        pair.right.low = q.right.low
    else
        st.ref[q.right.low::Tuple{Int,Int}] = st.lowpt_edge[parent]
    end
    return nothing
end

# Merge the pairs that conflict with `ei` into `pair.left`.
function _merge_conflict_edges!(st::_LRState, pair::_ConflictPair, ei)
    while _conflicts_with_top(st, ei)
        q = pop!(st.stack)
        _conflicting(q.right, ei, st.lowpt) && _swap_pair!(q)
        _conflicting(q.right, ei, st.lowpt) && return false
        _chain_ref!(st, pair.right, q.right.high)
        q.right.low !== nothing && (pair.right.low = q.right.low)
        _prepend_left!(st, pair, q)
    end
    return true
end

function _conflicts_with_top(st::_LRState, ei::Tuple{Int,Int})
    isempty(st.stack) && return false
    return _pair_conflicts(st, last(st.stack), ei)
end

function _pair_conflicts(st::_LRState, pair::_ConflictPair, ei::Tuple{Int,Int})
    return _conflicting(pair.left, ei, st.lowpt) || _conflicting(pair.right, ei, st.lowpt)
end

function _prepend_left!(st::_LRState, pair::_ConflictPair, q::_ConflictPair)
    if _interval_empty(pair.left)
        pair.left = _copy_interval(q.left)
    else
        _chain_ref!(st, pair.left, q.left.high)
    end
    pair.left.low = q.left.low
    return nothing
end

function _remove_back_edges!(st::_LRState, e::Tuple{Int,Int})
    u = e[1]
    while _stack_lowest(st) == st.height[u]
        q = pop!(st.stack)
        q.left.low !== nothing && (st.side[q.left.low::Tuple{Int,Int}] = -1)
    end

    if !isempty(st.stack)
        q = _trim_interval!(st, pop!(st.stack), u)
        push!(st.stack, q)
    end

    st.lowpt[e] < st.height[u] || return nothing
    top = _stack_top(st.stack)
    top === nothing && return nothing
    hl = top.left.high
    hr = top.right.high
    if hl !== nothing && (hr === nothing || st.lowpt[hl] > st.lowpt[hr])
        st.ref[e] = hl
    else
        _set_ref!(st, e, hr)
    end
    return nothing
end

function _trim_interval!(st::_LRState, pair::_ConflictPair, u::Int)
    _trim_side!(st, pair.left, pair.right, u)
    _trim_side!(st, pair.right, pair.left, u)
    return pair
end

# Drop the entries that return to `u`; when a side empties, the remaining side
# takes over the reference and the emptied low edge gets side -1.
function _trim_side!(st::_LRState, side::_Interval, other::_Interval, u::Int)
    while _points_to(side.high, u)
        side.high = _reference(st, side.high)
    end
    if side.high === nothing && side.low !== nothing
        _chain_ref!(st, side, other.low)
        st.side[side.low::Tuple{Int,Int}] = -1
        side.low = nothing
    end
    return nothing
end

function _dfs_testing!(st::_LRState, root::Int)
    ind = zeros(Int, length(st.adj))
    skip_init = Set{Tuple{Int,Int}}()
    dfs_stack = Int[root]

    while !isempty(dfs_stack)
        v = pop!(dfs_stack)
        e = st.parent_edge[v]
        skip_final = false
        while ind[v] < length(st.ordered_adjs[v])
            w = st.ordered_adjs[v][ind[v]+1]
            ei = (v, w)
            descend = !(ei in skip_init) && _start_testing_edge!(st, w, ei)
            if descend
                push!(dfs_stack, v)
                push!(dfs_stack, w)
                push!(skip_init, ei)
                skip_final = true
                break
            end
            _integrate_edge!(st, v, w, ei, e) || return false
            ind[v] += 1
        end
        skip_final || e === nothing || _remove_back_edges!(st, e)
    end
    return true
end

# Record where the conflict stack stands for `ei`. A tree edge is descended into
# and its constraint work is left for the second visit; a back edge gets a pair
# of its own.
function _start_testing_edge!(st::_LRState, w::Int, ei::Tuple{Int,Int})
    st.stack_bottom[ei] = _stack_top(st.stack)
    ei == st.parent_edge[w] && return true
    st.lowpt_edge[ei] = ei
    push!(st.stack, _ConflictPair(_Interval(), _Interval(ei, ei)))
    return false
end

# Account for the return edges of `ei`, which ends at `w`.

function _integrate_edge!(st::_LRState, v::Int, w::Int, ei::Tuple{Int,Int}, e)
    st.lowpt[ei] < st.height[v] || return true
    if w == first(st.ordered_adjs[v])
        e === nothing || (st.lowpt_edge[e] = st.lowpt_edge[ei])
        return true
    end
    e === nothing && return true
    return _add_constraints!(st, ei, e)
end

# --- phase 3: sides and the embedding ------------------------------------

function _sign!(st::_LRState, e::Tuple{Int,Int})
    old_ref = Dict{Tuple{Int,Int},Union{Nothing,Tuple{Int,Int}}}()
    dfs_stack = Tuple{Int,Int}[e]

    while !isempty(dfs_stack)
        current = pop!(dfs_stack)
        reference = get(st.ref, current, nothing)
        if reference !== nothing
            push!(dfs_stack, current)
            push!(dfs_stack, reference)
            old_ref[current] = reference
            delete!(st.ref, current)
            continue
        end
        previous = get(old_ref, current, nothing)
        factor = previous === nothing ? 1 : get(st.side, previous, 1)
        st.side[current] = get(st.side, current, 1) * factor
    end
    return get(st.side, e, 1)
end

# Every oriented edge `v -> w` contributes the half-edge `v -> w` to the
# rotation of `v`, in nesting-depth order; the embedding pass below then inserts
# the half-edge `w -> v` of every edge into the rotation of `w`.
function _initialize_rotation!(halves::_HalfEdges, v::Int, ordered::Vector{Int})
    isempty(ordered) && return halves
    previous = 0
    for w in ordered
        _add_half_edge!(halves, v, w; ccw_ref = previous == 0 ? nothing : previous)
        previous = w
    end
    return halves
end

function _dfs_embedding!(st::_LRState, root::Int)
    ind = zeros(Int, length(st.adj))
    dfs_stack = Int[root]
    halves = st.halves

    while !isempty(dfs_stack)
        v = pop!(dfs_stack)
        while ind[v] < length(st.ordered_adjs[v])
            w = st.ordered_adjs[v][ind[v]+1]
            ind[v] += 1
            ei = (v, w)

            if ei == st.parent_edge[w]
                _add_half_edge_first!(halves, w, v)
                st.left_ref[v] = w
                st.right_ref[v] = w
                push!(dfs_stack, v)
                push!(dfs_stack, w)
                break
            end

            if get(st.side, ei, 1) == 1
                _add_half_edge!(halves, w, v; ccw_ref = st.right_ref[w])
            else
                _add_half_edge!(halves, w, v; cw_ref = st.left_ref[w])
                st.left_ref[w] = v
            end
        end
    end
    return st
end

# --- driver --------------------------------------------------------------

function _rotation_system(adj::Vector{Vector{Int}})
    n = length(adj)
    m = sum(length, adj; init = 0) ÷ 2
    n > 2 && m > 3 * n - 6 && return nothing

    st = _LRState(adj)
    roots = Int[]
    for v = 1:n
        st.height[v] == -1 || continue
        st.height[v] = 0
        push!(roots, v)
        _dfs_orientation!(st, v)
    end

    for v = 1:n
        st.ordered_adjs[v] =
            sort(st.out_adjs[v]; by = w -> st.nesting_depth[(v, w)], alg = MergeSort)
    end
    for r in roots
        _dfs_testing!(st, r) || return nothing
    end

    for vw in st.oriented
        st.nesting_depth[vw] = _sign!(st, vw) * st.nesting_depth[vw]
    end
    for v = 1:n
        st.ordered_adjs[v] =
            sort(st.out_adjs[v]; by = w -> st.nesting_depth[(v, w)], alg = MergeSort)
        _initialize_rotation!(st.halves, v, st.ordered_adjs[v])
    end
    for r in roots
        _dfs_embedding!(st, r)
    end

    return [_rotation_of(st.halves, v) for v = 1:n]
end

_adjacency(graph::SimpleGraph) = [Int.(collect(neighbors(graph, v))) for v = 1:nv(graph)]

"""
    is_planar(g::PlanarGraph) -> Bool

Whether the underlying graph of `g` is planar. This is `Graphs.is_planar`, the
Left-Right planarity test; the stored rotation system plays no part, so a
`PlanarGraph` with a broken embedding still reports whether its graph is planar.
Use [`is_planar_embedding`](@ref) for the embedding question, and
[`planar_embedding`](@ref) to compute one.

!!! warning
    `Graphs.is_planar` recurses once per DFS edge and sorts the nesting depths
    by scanning them once per vertex, so its stack depth is the depth of the DFS
    forest and its running time is quadratic. It overflows the Julia stack on
    graphs with long paths — a `path_graph(100_000)`, or a road network — and
    the session is no longer trustworthy after that. Use
    [`is_planar_linear`](@ref) for anything larger than a few thousand vertices.
"""
is_planar(g::PlanarGraph) = is_planar(g.graph)

"""
    is_planar_linear(g) -> Bool

Whether `g` is planar, decided in `O(nv(g) + ne(g))` time by this package's own
iterative Left-Right test; the same criterion as `Graphs.is_planar`, without its
two scalability limits.

The two DFS passes run on an explicit stack, so the depth of the DFS forest does
not become the depth of the Julia call stack, and the graph is scanned a constant
number of times, so the running time is linear. Both matter in practice:
`is_planar` (which is `Graphs.is_planar`) throws a `StackOverflowError` from
about 10^5 vertices on a path, and is quadratic well before that, whereas this
function handles a road network with a million intersections.

Accepts a `SimpleGraph` and a `PlanarGraph`; for the latter a valid stored
rotation already proves planarity, so the answer is `true` without running the
test, and an invalid one is recomputed from the underlying graph.

```julia
julia> is_planar_linear(path_graph(100_000))
true
```

Use [`planar_embedding`](@ref) when the rotation system is wanted as well.
"""
is_planar_linear(g::Union{PlanarGraph,SimpleGraph}) = planar_embedding(g) !== nothing

"""
    planar_embedding(g::SimpleGraph) -> Union{PlanarGraph, Nothing}
    planar_embedding(g::PlanarGraph) -> Union{PlanarGraph, Nothing}

Compute a planar embedding of `g`, in `O(n + m)` time, or `nothing` when the
graph is not planar. The result is a [`PlanarGraph`](@ref) whose rotation system
satisfies Euler's formula on every component, so it can be passed to
[`separator`](@ref), [`triangulate`](@ref) and [`grid_drawing`](@ref).

For a `PlanarGraph` whose stored rotation is not planar, the embedding is
recomputed from the underlying graph. Vertices keep their labels.
"""
function planar_embedding(g::SimpleGraph)
    rotation = _rotation_system(_adjacency(g))
    rotation === nothing && return nothing
    return PlanarGraph(rotation)
end

function planar_embedding(g::PlanarGraph)
    is_planar_embedding(g) && return g
    return planar_embedding(g.graph)
end
