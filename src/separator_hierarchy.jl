# Persistent vertex-disjoint separator hierarchy.
#
# Every node owns a balanced Lipton-Tarjan separator of its induced subproblem and
# has the two sides as children. Separators are removed, not copied, so the owned
# lists partition the vertices. Nodes are numbered in depth-first preorder, which
# makes every subtree an interval of node ids.

"""
    HierarchyNode

One node of a [`SeparatorHierarchy`](@ref): `owned` (a view into one packed vertex
array), `children` (zero, one or two node ids, in traversal order), `last` (the last
node id of the subtree, inclusive) and `size` (the number of original vertices of the
subproblem, owned ones included).
"""
struct HierarchyNode
    owned::SubArray{Int,1,Vector{Int},Tuple{UnitRange{Int}},true}
    children::Vector{Int}
    last::Int
    size::Int
end

"""
    SeparatorHierarchy

Result of [`separator_hierarchy`](@ref): `root` (`1`, or `0` for the empty graph),
`owner[v]` (the node that owns original vertex `v`) and `nodes`, indexed by preorder
id.
"""
struct SeparatorHierarchy
    root::Int
    owner::Vector{Int}
    nodes::Vector{HierarchyNode}
end

"""
    separator_hierarchy(g::PlanarGraph; validate = false) -> SeparatorHierarchy
    separator_hierarchy(rotations; validate = false) -> SeparatorHierarchy

Recursive balanced-separator hierarchy of the embedded planar graph `g`, or of the
simple planar graph whose vertex `v` has the neighbours `rotations[v]` in
counter-clockwise cyclic order (the starting position is arbitrary). The graph may be
disconnected. Neither the graph nor the rotations are modified, and the result is
deterministic.

For a node `t` let `S_t = h.nodes[t].owned`, `U_t` the union of the owned lists in the
subtree of `t`, and `q_t = |U_t| = h.nodes[t].size`. Then:

1. the lists `S_t` partition `1:n`, and `h.owner[v] == t` exactly for `v ∈ S_t`;
2. `U_root = 1:n`; for `n = 0` there are no nodes and `root == 0`;
3. `S_t`, `U_a`, `U_b` partition `U_t` for the children `a, b` (zero, one or two,
   all nonempty);
4. no edge of `g` joins two different child subproblems, so the owners of the ends of
   an edge are equal or ancestor and descendant;
5. `3 * q_c <= 2 * q_t` for every child `c`;
6. `|S_t| <= 2 * sqrt(2 * q_t)`, the bound of [`separator`](@ref);
7. a singleton subproblem is owned by a leaf; a node without children owns its whole
   subproblem, which obeys 6;
8. `S_t` may be empty, and then `t` has two children;
9. nodes are numbered in depth-first preorder, the subtree of `t` being `t:last`.

There are at most `2n - 1` nodes and the depth is `O(log(n + 1))`. Along a root-to-leaf
path `∑ |S_t| <= 2 * sqrt(2n) / (1 - sqrt(2/3))`.

# Cost
Time `O(n * (1 + log(n + 1)))` and space `O(n)`, deterministic. The conversion to darts
is linear (counting sorts, no hashing). A node with `q` vertices costs `O(q)`: the
separator, the induced subgraphs of its sides (renumbered by maps of local size) and
the original-id lists. The subproblems of one depth are disjoint and the depth is
`O(log(n + 1))`. A subproblem is dropped as soon as its children exist, so the live
graphs are disjoint and total `O(n)`.

The embedding is a precondition. With `validate = true` Euler's formula is checked on
every component, in `O(n)` additional time; a violation throws an `ArgumentError`.
Index errors, loops and parallel edges always throw an `ArgumentError`.

See [`child_containing`](@ref) for routing.

```julia
h = separator_hierarchy(planar_grid(10, 10))
h.nodes[h.root].size                # 100
child_containing(h, h.root, h.owner[57])
```
"""
function separator_hierarchy(g::PlanarGraph; validate::Bool = false)
    return _hierarchy(DartGraph(g), validate)
end

function separator_hierarchy(
    rotations::AbstractVector{<:AbstractVector{<:Integer}};
    validate::Bool = false,
)
    return _hierarchy(_rotation_dart_graph(rotations), validate)
end

function _rotation_dart_graph(rotations)
    n = length(rotations)
    dg = DartGraph(n)
    for v = 1:n
        anchor = 0
        for w in rotations[v]
            (1 <= w <= n) ||
                throw(ArgumentError("neighbour $w of vertex $v is out of range"))
            w == v && throw(ArgumentError("vertex $v lists itself as a neighbour"))
            anchor = add_dart!(dg, v, Int(w), anchor)
        end
    end
    _match_pending_twins!(dg)
    return dg
end

# Preorder construction. The recursion depth is `O(log n)`, since every child has at
# most two thirds of the vertices of its parent.
mutable struct _HierarchyBuilder
    packed::Vector{Int}
    owner::Vector{Int}
    start::Vector{Int}
    count::Vector{Int}
    children::Vector{Vector{Int}}
    last::Vector{Int}
    size::Vector{Int}
end

function _hierarchy(dg::DartGraph, validate::Bool)
    validate && !euler_holds(dg) && throw(_not_planar_error())
    n = dg.n
    n == 0 && return SeparatorHierarchy(0, Int[], HierarchyNode[])
    b = _HierarchyBuilder(
        sizehint!(Int[], n),
        zeros(Int, n),
        Int[],
        Int[],
        Vector{Int}[],
        Int[],
        Int[],
    )
    _build!(b, dg, collect(1:n))
    nodes = [
        HierarchyNode(
            view(b.packed, b.start[t]:(b.start[t]+b.count[t]-1)),
            b.children[t],
            b.last[t],
            b.size[t],
        ) for t in eachindex(b.start)
    ]
    return SeparatorHierarchy(1, b.owner, nodes)
end

# Create the node of the subproblem `dg` (original ids `ids`), then its children.
# `dg` and `ids` are not referenced after the children are built from them.
function _build!(b::_HierarchyBuilder, dg::DartGraph, ids::Vector{Int})
    t = length(b.start) + 1
    q = dg.n
    side = q == 1 ? Int8[0] : _best_partition(dg, (2 * q) ÷ 3)
    push!(b.start, length(b.packed) + 1)
    for v in findall(iszero, side)
        push!(b.packed, ids[v])
        b.owner[ids[v]] = t
    end
    push!(b.count, length(b.packed) - b.start[t] + 1)
    push!(b.children, Int[])
    push!(b.last, t)
    push!(b.size, q)
    for s in (Int8(1), Int8(2))
        mask = side .== s
        any(mask) || continue
        sub, local_ids = induced(dg, mask)
        push!(b.children[t], length(b.start) + 1)
        _build!(b, sub, Int[ids[v] for v in local_ids])
    end
    b.last[t] = length(b.start)
    return nothing
end

"""
    child_containing(h, t, label) -> Int

The child of node `t` whose subtree contains the node `label`, or `0` when
`label == t`. Worst-case `O(1)` time, without allocation. Throws an `ArgumentError`
when `label` is not in the subtree of `t`; an invalid `t` raises a bounds error. To
route an original vertex `v` pass `h.owner[v]`.
"""
function child_containing(h::SeparatorHierarchy, t::Integer, label::Integer)
    node = h.nodes[t]
    t <= label <= node.last ||
        throw(ArgumentError("node $label is not in the subtree of node $t"))
    label == t && return 0
    for c in node.children
        label <= h.nodes[c].last && return c
    end
    throw(ArgumentError("node $label is not in the subtree of node $t"))
end
