# Incidence-preserving simplification of a multigraph with a rotation system.
#
# Loops are dropped and, for every unordered vertex pair, the edge whose smaller
# incidence ID is least is kept. The supplied rotations are filtered, never reordered.

"""
    IncidenceEdge(u, v, at_u, at_v)

An undirected edge between vertices `u` and `v` (a loop if equal). `at_u` and `at_v` are
its two distinct incidence IDs, the entries naming it in `rotations[u]` and
`rotations[v]`. `IncidenceEdge(v, u, at_v, at_u)` describes the same edge.
"""
struct IncidenceEdge
    u::Int
    v::Int
    at_u::Int
    at_v::Int
end

"""
    SimplifiedEmbedding

Result of [`simplify_embedding`](@ref): for every vertex `v`, `neighbors[v]` and
`incidences[v]` have equal length, and entry `j` of each describes one surviving
incidence (its original ID) and the neighbor it leads to. `neighbors` can be passed to
[`separator_hierarchy`](@ref). The arrays are fresh; treat them as read-only.
"""
struct SimplifiedEmbedding
    neighbors::Vector{Vector{Int}}
    incidences::Vector{Vector{Int}}
end

# Fill owner/opposite (size U) from the edges; reject malformed edge records.
function _index_incidences(n::Int, edges, U::Int)
    owner = zeros(Int, U)
    opposite = zeros(Int, U)
    for e in edges
        (1 <= e.u <= n && 1 <= e.v <= n) ||
            throw(ArgumentError("edge endpoint outside 1:$n"))
        a, b = e.at_u, e.at_v
        (1 <= a <= U && 1 <= b <= U) ||
            throw(ArgumentError("incidence ID outside 1:$U (IDs must be positive)"))
        a != b || throw(ArgumentError("an edge uses the same incidence ID twice"))
        (owner[a] == 0 && owner[b] == 0) ||
            throw(ArgumentError("an incidence ID is used by several edges"))
        owner[a] = e.u
        owner[b] = e.v
        opposite[a] = b
        opposite[b] = a
    end
    return owner, opposite
end

# Every rotation entry must be an unseen incidence of its own vertex, and all must appear.
function _check_rotations(rotations, owner::Vector{Int}, m::Int)
    seen = falses(length(owner))
    count = 0
    for (v, rotation) in enumerate(rotations), id in rotation
        (1 <= id <= length(owner) && owner[id] == v) ||
            throw(ArgumentError("rotation entry $id is foreign or at the wrong vertex"))
        seen[id] && throw(ArgumentError("rotation entry $id is repeated"))
        seen[id] = true
        count += 1
    end
    count == 2 * m ||
        throw(ArgumentError("an edge incidence is missing from the rotations"))
    return nothing
end

# Stable counting sort of `items` by `key[item]` in 1:n (O(n + length(items))).
function _stable_by(items::Vector{Int}, key::Vector{Int}, n::Int)
    start = zeros(Int, n + 1)
    for i in items
        start[key[i]+1] += 1
    end
    start[1] = 1
    for k = 1:n
        start[k+1] += start[k]
    end
    sorted = Vector{Int}(undef, length(items))
    for i in items
        sorted[start[key[i]]] = i
        start[key[i]] += 1
    end
    return sorted
end

# Mark the two incidences of the minimum-key edge of every nonloop endpoint pair.
function _select_edges(edges, n::Int, U::Int)
    lo = Vector{Int}(undef, length(edges))
    hi = similar(lo)
    nonloop = Int[]
    for (i, e) in enumerate(edges)
        lo[i], hi[i] = minmax(e.u, e.v)
        e.u != e.v && push!(nonloop, i)
    end
    order = _stable_by(_stable_by(nonloop, hi, n), lo, n)
    keep = falses(U)
    best = 0
    for (k, i) in enumerate(order)
        if k == 1 || lo[i] != lo[order[k-1]] || hi[i] != hi[order[k-1]]
            best != 0 && _mark!(keep, edges[best])
            best = i
        elseif _key(edges[i]) < _key(edges[best])
            best = i
        end
    end
    best != 0 && _mark!(keep, edges[best])
    return keep
end

_key(e::IncidenceEdge) = min(e.at_u, e.at_v)

function _mark!(keep::BitVector, e::IncidenceEdge)
    keep[e.at_u] = true
    keep[e.at_v] = true
    return keep
end

function _filter_rotations(rotations, owner::Vector{Int}, opposite::Vector{Int}, keep)
    incidences = [Int[id for id in rotation if keep[id]] for rotation in rotations]
    neighbors = [Int[owner[opposite[id]] for id in list] for list in incidences]
    return SimplifiedEmbedding(neighbors, incidences)
end

function _simplify_dense(rotations, edges, U::Int)
    n = length(rotations)
    owner, opposite = _index_incidences(n, edges, U)
    _check_rotations(rotations, owner, length(edges))
    keep = _select_edges(edges, n, U)
    return _filter_rotations(rotations, owner, opposite, keep)
end

"""
    simplify_embedding(rotations, edges; universe_size = largest ID) -> SimplifiedEmbedding

Remove loops and keep one edge per unordered vertex pair of the multigraph on vertices
`1:length(rotations)` with edges `edges::AbstractVector{IncidenceEdge}`;
`rotations[v]` lists the incidence IDs at `v` in cyclic order. The kept edge of a pair
is the one with least `min(at_u, at_v)`, so the result does not depend on edge order or
orientation. Each rotation is filtered, keeping the relative order and starting entry of
the survivors; IDs and vertices are never renumbered. Inputs are not mutated.

Incidence IDs must lie in `1:universe_size`; the default is the largest supplied ID (zero
without edges), and arrays of that size are allocated, so use
[`simplify_embedding_sparse`](@ref) for large, scattered IDs. Deterministic time and
space `O(n + m + U)` (plus `O(R)` for malformed rotation lists of total length `R`), by
counting sorts. Throws `ArgumentError` on malformed input.
"""
function simplify_embedding(
    rotations::AbstractVector{<:AbstractVector{Int}},
    edges::AbstractVector{IncidenceEdge};
    universe_size::Int = maximum((max(e.at_u, e.at_v) for e in edges); init = 0),
)
    universe_size >= 0 || throw(ArgumentError("negative incidence universe"))
    return _simplify_dense(rotations, edges, universe_size)
end

"""
    simplify_embedding_sparse(rotations, edges) -> SimplifiedEmbedding

Same result as [`simplify_embedding`](@ref) for arbitrary positive `Int` incidence IDs.
Storage never depends on the magnitude of an ID: the supplied IDs are sorted (merge
sort), replaced by their ranks, and translated back, with rotation entries located by
binary search. Deterministic `O(n + m log(m + 2))` time and `O(n + m)` space for valid
input. Throws `ArgumentError` on malformed input.
"""
function simplify_embedding_sparse(
    rotations::AbstractVector{<:AbstractVector{Int}},
    edges::AbstractVector{IncidenceEdge},
)
    ids = Int[id for e in edges for id in (e.at_u, e.at_v)]
    all(>(0), ids) || throw(ArgumentError("incidence IDs must be positive"))
    sort!(ids; alg = MergeSort)
    for k = 2:length(ids)
        ids[k] != ids[k-1] ||
            throw(ArgumentError("an incidence ID is used by several edges"))
    end
    rank(id) = _rank(ids, id)
    compact = [IncidenceEdge(e.u, e.v, rank(e.at_u), rank(e.at_v)) for e in edges]
    order = [Int[rank(id) for id in rotation] for rotation in rotations]
    result = _simplify_dense(order, compact, length(ids))
    original = [Int[ids[r] for r in list] for list in result.incidences]
    return SimplifiedEmbedding(result.neighbors, original)
end

function _rank(ids::Vector{Int}, id::Int)
    r = searchsortedfirst(ids, id)
    (r <= length(ids) && ids[r] == id) ||
        throw(ArgumentError("rotation entry $id belongs to no edge"))
    return r
end
