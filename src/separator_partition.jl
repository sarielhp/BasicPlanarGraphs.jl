# Recursive separator partition.
#
# Splits an embedded planar graph into one global separator `S` and parts of at
# most `r` vertices with no edge between two different parts. Every recursion
# step takes a group of `q > r` vertices, removes a balanced Lipton-Tarjan
# separator into `S` and continues with the two sides, each at most `2q/3`.
# Groups at one level have disjoint vertex sets, so a level costs `O(n)` and
# the recursion has `O(1 + log(n/r))` levels; the terminal groups are packed
# next-fit into parts of at most `r` vertices.

"""
    SeparatorPartition

Result of [`separator_partition`](@ref). All vertex numbers are those of the
input graph.

- `separator`: the vertices removed by the recursion;
- `parts`: the vertex lists of the parts, each nonempty with at most `r`
  vertices;
- `part_of[v]`: `0` when `v` is in `separator`, else the index of its part;
- `local_index`, `graphs`: `nothing`, or with `subgraphs = true` the position of
  each part vertex in its part (so `parts[i][local_index[v]] == v`) and the
  induced graph of each part, vertex `j` of `graphs[i]` being `parts[i][j]`.
"""
struct SeparatorPartition
    separator::Vector{Int}
    parts::Vector{Vector{Int}}
    part_of::Vector{Int}
    local_index::Union{Nothing,Vector{Int}}
    graphs::Union{Nothing,Vector{PlanarGraph{Int}}}
end

"""
    separator_partition(g, r; subgraphs = false) -> SeparatorPartition

Partition the vertices of the embedded planar graph `g` into a separator `S`
and parts `P₁, …, Pₖ` such that

1. `S, P₁, …, Pₖ` partition `vertices(g)` and every part is nonempty;
2. each part has at most `r` vertices;
3. no edge joins two different parts (edges inside a part, inside `S`, and
   between a part and `S` are allowed);
4. `|S| <= 16 * n / sqrt(r)`;
5. `k <= 2n/r + 1`;
6. `part_of[v] == 0` exactly for `v ∈ S`, else `part_of[v] == i` for `v ∈ Pᵢ`.

The graph may be disconnected and may have isolated vertices. A part may hold
several components of `g - S`. The result is deterministic and `g` is not
modified. For `r >= nv(g)` the separator is empty and there is one part.

The running time is `O(n * (1 + log(n/r)))` and the space is `O(n)`, where the
cost of every conversion, component search and grouping is included. `r` must
be at least `1`. With `subgraphs = true` the induced graph of each part is also
returned, in `O(n)` additional time and space.

```julia
g = planar_grid(30, 30)
p = separator_partition(g, 50)
all(length(P) <= 50 for P in p.parts)     # true
p.part_of[p.separator[1]] == 0            # true
```
"""
function separator_partition(g::PlanarGraph, r::Integer; subgraphs::Bool = false)
    r >= 1 || throw(ArgumentError("r must be at least 1, got $r"))
    n = nv(g)
    n == 0 && return SeparatorPartition(Int[], Vector{Int}[], Int[], nothing, nothing)
    separator, terminals = _recurse(_planar_dart_graph(g), Int(r))
    parts = _next_fit(terminals, Int(r))
    part_of = zeros(Int, n)
    for (i, part) in enumerate(parts), v in part
        part_of[v] = i
    end
    subgraphs || return SeparatorPartition(separator, parts, part_of, nothing, nothing)
    local_index, graphs = _induced_parts(g, parts, part_of)
    return SeparatorPartition(separator, parts, part_of, local_index, graphs)
end

# --- the recursion --------------------------------------------------------

# A group is an induced embedded graph together with the original ids of its
# vertices. A processed group is dropped as soon as its children exist, so the
# live groups are disjoint and their total size stays `O(n)`.
function _recurse(dg::DartGraph, r::Int)
    separator = Int[]
    terminals = Vector{Vector{Int}}()
    stack = Tuple{DartGraph,Vector{Int}}[(dg, collect(1:dg.n))]
    while !isempty(stack)
        graph, ids = pop!(stack)
        if graph.n <= r
            push!(terminals, ids)
            continue
        end
        side = _best_partition(graph, (2 * graph.n) ÷ 3)
        append!(separator, ids[v] for v in findall(iszero, side))
        for s in (Int8(1), Int8(2))
            child = _child(graph, ids, side, s)
            child === nothing || push!(stack, child)
        end
    end
    return separator, terminals
end

function _child(graph::DartGraph, ids::Vector{Int}, side::Vector{Int8}, s::Int8)
    mask = side .== s
    any(mask) || return nothing
    sub, local_ids = induced(graph, mask)
    return (sub, Int[ids[v] for v in local_ids])
end

# Next-fit: the next terminal group joins the current part when it fits. Two
# consecutive parts together hold more than `r` vertices, so there are at most
# `2n/r + 1` of them.
function _next_fit(terminals::Vector{Vector{Int}}, r::Int)
    parts = Vector{Vector{Int}}()
    current = Int[]
    for group in terminals
        if !isempty(current) && length(current) + length(group) > r
            push!(parts, current)
            current = Int[]
        end
        append!(current, group)
    end
    isempty(current) || push!(parts, current)
    return parts
end

# --- induced graphs of the parts ------------------------------------------

# Restricting every rotation to the neighbours in the same part gives the
# inherited embedding; the total work is the number of darts.
function _induced_parts(g::PlanarGraph, parts::Vector{Vector{Int}}, part_of::Vector{Int})
    local_index = zeros(Int, nv(g))
    for part in parts, (j, v) in enumerate(part)
        local_index[v] = j
    end
    graphs = PlanarGraph{Int}[
        PlanarGraph([_local_rotation(g, local_index, part_of, v, i) for v in part]) for
        (i, part) in enumerate(parts)
    ]
    return local_index, graphs
end

function _local_rotation(g::PlanarGraph, local_index, part_of, v::Int, i::Int)
    return Int[local_index[w] for w in neighbor_order(g, v) if part_of[w] == i]
end
