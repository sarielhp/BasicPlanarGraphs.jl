# Five-coloring of planar graphs (Matula, Shiloach and Tarjan, 1980).
#
# Every planar graph has a vertex of degree at most five. Such a vertex `v` is
# removed and remembered; when its degree is five two neighbours `a`, `b` are not
# adjacent (otherwise the five neighbours and `v` would contain a `K5`), and they
# are identified into one vertex, which keeps the graph planar. Undoing the steps
# in reverse order colors the graph: an identified pair gets one color, and `v`
# then sees at most four colors, so a fifth is free.
#
# Adjacency lives in one hash set per vertex, so adjacency tests are `O(1)` and an
# identification moves the smaller set into the larger.

struct _Remove
    v::Int
    neighbors::Vector{Int}
end

struct _Merge
    gone::Int
    kept::Int
end

mutable struct _Elimination
    adj::Vector{Set{Int}}
    alive::BitVector
    remaining::Int
    low::Vector{Int}
    events::Vector{Union{_Remove,_Merge}}
end

function _Elimination(g::SimpleGraph)
    n = nv(g)
    adj = [Set{Int}(neighbors(g, v)) for v = 1:n]
    low = [v for v = 1:n if length(adj[v]) <= 5]
    return _Elimination(adj, trues(n), n, low, Union{_Remove,_Merge}[])
end

function _not_five_colorable()
    return ArgumentError("the graph is not planar: no vertex of degree at most 5 remains")
end

# Queue `v` again when it has become cheap to remove.
function _requeue!(st::_Elimination, v::Int)
    st.alive[v] && length(st.adj[v]) <= 5 && push!(st.low, v)
    return nothing
end

# Two non-adjacent vertices among the five neighbours of a removed vertex.
function _nonadjacent_pair(st::_Elimination, nbrs::Vector{Int})
    for i = 1:4, j = (i+1):5
        nbrs[j] in st.adj[nbrs[i]] || return nbrs[i], nbrs[j]
    end
    throw(ArgumentError("the graph is not planar: it contains K5"))
end

# Identify `a` and `b` (non-adjacent): the vertex with the smaller set disappears.
function _merge!(st::_Elimination, a::Int, b::Int)
    kept, gone = length(st.adj[a]) >= length(st.adj[b]) ? (a, b) : (b, a)
    st.alive[gone] = false
    st.remaining -= 1
    push!(st.events, _Merge(gone, kept))
    for w in st.adj[gone]
        delete!(st.adj[w], gone)
        if w in st.adj[kept]
            _requeue!(st, w)
        else
            push!(st.adj[kept], w)
            push!(st.adj[w], kept)
        end
    end
    empty!(st.adj[gone])
    _requeue!(st, kept)
    return nothing
end

function _remove!(st::_Elimination, v::Int)
    nbrs = collect(st.adj[v])
    st.alive[v] = false
    st.remaining -= 1
    empty!(st.adj[v])
    push!(st.events, _Remove(v, nbrs))
    for w in nbrs
        delete!(st.adj[w], v)
    end
    length(nbrs) == 5 && _merge!(st, _nonadjacent_pair(st, nbrs)...)
    foreach(w -> _requeue!(st, w), nbrs)
    return nothing
end

function _eliminate!(st::_Elimination)
    while !isempty(st.low)
        v = pop!(st.low)
        st.alive[v] && length(st.adj[v]) <= 5 && _remove!(st, v)
    end
    st.remaining == 0 || throw(_not_five_colorable())
    return st
end

# Smallest color in 1:5 absent from `nbrs`, which show at most four colors.
function _free_color(color::Vector{Int}, nbrs::Vector{Int})
    used = 0
    for w in nbrs
        used |= 1 << color[w]
    end
    return findfirst(c -> used & (1 << c) == 0, 1:5)::Int
end

function _unwind(st::_Elimination, n::Int)
    color = zeros(Int, n)
    for event in Iterators.reverse(st.events)
        if event isa _Merge
            color[event.gone] = color[event.kept]
        else
            color[event.v] = _free_color(color, event.neighbors)
        end
    end
    return color
end

"""
    five_coloring(g) -> Vector{Int}

A proper vertex coloring of the planar graph `g` with the colors `1:5`: the
returned vector `c` has `c[u] != c[v]` for every edge `uv`. `g` is a `SimpleGraph`
or a [`PlanarGraph`](@ref) (only its edges are used, not the rotation) and may be
disconnected. Throws an `ArgumentError` when the method gets stuck, which proves that `g`
is not planar; some non-planar graphs are colored properly all the same, and a
returned coloring is always proper.

The method is that of Matula, Shiloach and Tarjan: repeatedly remove a vertex of
degree at most five, identifying two non-adjacent neighbours when there are five,
then color in reverse order. Identifying two vertices moves the smaller adjacency set
into the larger, so the time is near linear in practice (`benchmark/run.jl`), but
linear time is not proved: that needs the linked-list scheme of the paper.

```julia
g = planar_grid(10, 10)
c = five_coloring(g)
is_proper_coloring(g, c)    # true
maximum(c) <= 5             # true
```
"""
function five_coloring(g::SimpleGraph)
    return _unwind(_eliminate!(_Elimination(g)), Int(nv(g)))
end
five_coloring(g::PlanarGraph) = five_coloring(g.graph)
five_coloring(::AbstractGraph) = throw(ArgumentError(_UNDIRECTED_ONLY))

"""
    is_proper_coloring(g, colors) -> Bool

Whether `colors` has one entry per vertex of `g` and no edge of `g` joins two
vertices of the same color. The values themselves are not restricted.
"""
function is_proper_coloring(g::Union{SimpleGraph,PlanarGraph}, colors::AbstractVector)
    length(colors) == nv(g) || return false
    return all(e -> colors[src(e)] != colors[dst(e)], edges(g))
end
