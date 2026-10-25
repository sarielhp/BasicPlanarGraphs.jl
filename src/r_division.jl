# r-division of an embedded planar graph.
#
# Regions of at most `r` vertices, each with a boundary of `O(sqrt(r))` vertices
# shared with other regions, and `O(n/r)` regions (Frederickson 1987). Two phases
# of recursive separators, both built on `_best_partition`:
#
#   1. a piece with more than `r` vertices is cut by a balanced separator into two
#      pieces, each holding the separator vertices it needs; those vertices are
#      shared, that is they occur in both pieces;
#   2. a piece of at most `r` vertices with more than `12 sqrt(r)` shared vertices
#      is cut by a separator that balances the shared vertices (a weight of 1 on
#      them), so both children keep at most two thirds of them plus the new
#      separator, which has at most `2 sqrt(2r)` vertices. The count falls
#      whenever it is above `6 sqrt(2) sqrt(r) < 12 sqrt(r)`.
#
# The pieces left are packed next-fit into regions that respect both limits.

const MIN_R_DIVISION = 100
const BOUNDARY_SQUARED = 144

"""
    RDivision

Result of [`r_division`](@ref); all vertex numbers are those of the input graph.

- `regions`: the vertex list of every region, boundary vertices included;
- `boundary`: the boundary vertices of every region, a subset of `regions[i]`;
- `region_of[v]`: the one region holding `v` when `v` is in only one region, and
  `0` for a boundary vertex;
- `boundary_vertices`: the vertices that lie in more than one region.
"""
struct RDivision
    regions::Vector{Vector{Int}}
    boundary::Vector{Vector{Int}}
    region_of::Vector{Int}
    boundary_vertices::Vector{Int}
end

# A piece under construction: an induced embedded graph, the original numbers of
# its vertices, and which of them occur in some other piece too.
struct _Piece
    graph::DartGraph
    ids::Vector{Int}
    shared::BitVector
end

"""
    r_division(g, r) -> RDivision

Divide the embedded planar graph `g` into regions such that

1. every vertex is in some region, and every edge has both ends in some region;
2. a vertex is in more than one region exactly when it is a boundary vertex;
3. every region has at most `r` vertices;
4. every region has at most `12 sqrt(r)` boundary vertices.

The number of regions is `O(n/r)` and the number of boundary vertices is
`O(n/sqrt(r))`; both are argued from the recursion and measured, not proved here
(`docs/r_division.md`). `r` must be at least $MIN_R_DIVISION, because the separator
of a piece must be small against the piece for the recursion to shrink it. The graph
may be disconnected, and a region may hold several components. The result is
deterministic and `g` is not modified. The running time is `O(n (1 + log(n/r)))` plus
the cost of cutting the pieces whose boundary is too large.

```julia
g = planar_grid(60, 60)
d = r_division(g, 200)
all(length(R) <= 200 for R in d.regions)      # true
```
"""
function r_division(g::PlanarGraph, r::Integer)
    r >= MIN_R_DIVISION ||
        throw(ArgumentError("r must be at least $MIN_R_DIVISION, got $r"))
    n = Int(nv(g))
    n == 0 && return RDivision(Vector{Int}[], Vector{Int}[], Int[], Int[])
    dg = _planar_dart_graph(g)
    leaves = _divide(_Piece(dg, collect(1:n), falses(n)), Int(r))
    return _finish(_pack(leaves, Int(r), n), n)
end

# --- the two phases -------------------------------------------------------

function _divide(start::_Piece, r::Int)
    leaves = _Piece[]
    stack = _Piece[start]
    while !isempty(stack)
        piece = pop!(stack)
        side = _cut(piece, r)
        side === nothing ? push!(leaves, piece) : append!(stack, _children(piece, side))
    end
    return leaves
end

# The sides of a separator of `piece`, or `nothing` when it is final.
function _cut(piece::_Piece, r::Int)
    q = piece.graph.n
    b = count(piece.shared)
    q > r && return _best_partition(piece.graph, (2 * q) ÷ 3)
    b * b > BOUNDARY_SQUARED * r || return nothing
    return _best_partition(piece.graph, (2 * b) ÷ 3, Int.(piece.shared))
end

# Which separator vertices go into which child: one with a neighbour on a side or
# in the separator goes with that side, or with both for a neighbour in the
# separator; one without neighbours goes into the first child. Every vertex
# and every edge of the piece stays in some child.
function _inclusion(graph::DartGraph, side::Vector{Int8})
    in1 = side .== 1
    in2 = side .== 2
    for x in findall(iszero, side)
        n1 = false
        n2 = false
        for d in incident_darts(graph, x)
            s = side[graph.head[d]]
            n1 |= s != 2
            n2 |= s != 1
        end
        in1[x] = n1 || !n2
        in2[x] = n2
    end
    return in1, in2
end

function _children(piece::_Piece, side::Vector{Int8})
    in1, in2 = _inclusion(piece.graph, side)
    shared = piece.shared .| (in1 .& in2)
    children = _Piece[]
    for mask in (in1, in2)
        any(mask) || continue
        sub, local_ids = induced(piece.graph, mask)
        sub.n < piece.graph.n || error("internal error: a division step made no progress")
        push!(children, _Piece(sub, piece.ids[local_ids], shared[local_ids]))
    end
    return children
end

# --- packing and the result ----------------------------------------------

# Next-fit: a piece joins the current region while the region stays within `r`
# vertices and within the boundary limit (the sum of the shared counts bounds the
# boundary of the union). Two consecutive regions break one of the two limits, so
# there are at most `2 n'/r + 2 B/(12 sqrt(r)) + 1`, with `n'` the total size of the
# pieces and `B` their total shared count.
function _pack(leaves::Vector{_Piece}, r::Int, n::Int)
    seen = zeros(Int, n)
    epoch = 1
    regions = Vector{Int}[]
    current = Int[]
    shared = 0
    for leaf in leaves
        b = count(leaf.shared)
        extra = count(v -> seen[v] != epoch, leaf.ids)
        if !isempty(current) &&
           (length(current) + extra > r || (shared + b)^2 > BOUNDARY_SQUARED * r)
            push!(regions, current)
            current = Int[]
            shared = 0
            epoch += 1
        end
        for v in leaf.ids
            seen[v] == epoch && continue
            seen[v] = epoch
            push!(current, v)
        end
        shared += b
    end
    isempty(current) || push!(regions, current)
    return regions
end

function _finish(regions::Vector{Vector{Int}}, n::Int)
    copies = zeros(Int, n)
    for region in regions, v in region
        copies[v] += 1
    end
    region_of = zeros(Int, n)
    for (i, region) in enumerate(regions), v in region
        copies[v] == 1 && (region_of[v] = i)
    end
    boundary = [Int[v for v in region if copies[v] > 1] for region in regions]
    return RDivision(regions, boundary, region_of, findall(>(1), copies))
end
