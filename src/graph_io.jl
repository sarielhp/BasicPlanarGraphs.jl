# Reading and writing planar graphs as plain text.
#
# The format is line based and sectioned, so that a file is readable, diffable
# and round-trippable:
#
#     # BasicPlanarGraphs 1
#     vertices 4
#     edges
#     1 2
#     2 3
#     3 4
#     4 1
#     rotation
#     1: 2 4
#     2: 3 1
#     3: 4 2
#     4: 1 3
#     end
#
# `vertices` and `rotation` are optional: without `vertices` the vertex count is
# the largest label seen, and without `rotation` the result is a bare
# `SimpleGraph`. Lines starting with `#` and blank lines are ignored anywhere.

const GRAPH_FILE_MAGIC = "BasicPlanarGraphs 1"

# A vertex count or label above this is refused by the readers unless the caller
# raises `max_vertices`: a file of a few bytes could otherwise ask for terabytes.
const DEFAULT_MAX_VERTICES = 10_000_000

"""
    write_planar_graph(io, g::PlanarGraph)
    write_planar_graph(path::AbstractString, g::PlanarGraph)

Write `g` in the sectioned text format: `vertices`, `edges`, and the rotation
system. The file reads back with [`read_planar_graph`](@ref) into an equal
graph; the same file is also readable with [`read_graph`](@ref) when only the
edges are wanted.
"""
function write_planar_graph(io::IO, g::PlanarGraph)
    n = nv(g)
    println(io, "# ", GRAPH_FILE_MAGIC)
    println(io, "vertices ", n)
    println(io, "edges")
    _write_edges(io, g.graph)
    println(io, "rotation")
    for v = 1:n
        print(io, v, ":")
        for w in g.rotation[v]
            print(io, " ", w)
        end
        println(io)
    end
    println(io, "end")
    return io
end

"""
    write_graph(io, g::SimpleGraph)
    write_graph(path::AbstractString, g::SimpleGraph)

Write the edges of `g` in the same text format, without a rotation section. The
file reads back with [`read_graph`](@ref) as a `SimpleGraph`.
"""
function write_graph(io::IO, g::SimpleGraph)
    println(io, "# ", GRAPH_FILE_MAGIC)
    println(io, "vertices ", nv(g))
    println(io, "edges")
    _write_edges(io, g)
    println(io, "end")
    return io
end

function _write_edges(io::IO, graph::SimpleGraph)
    for e in edges(graph)
        println(io, src(e), " ", dst(e))
    end
    return nothing
end

write_planar_graph(path::AbstractString, g::PlanarGraph) =
    open(io -> write_planar_graph(io, g), path, "w")
write_graph(path::AbstractString, g::SimpleGraph) =
    open(io -> write_graph(io, g), path, "w")

# --- parsing -------------------------------------------------------------

struct _GraphFile
    n::Int
    edges::Vector{Tuple{Int,Int}}
    has_edges::Bool
    rotation::Union{Nothing,Vector{Vector{Int}}}
end

_parse_error(source::AbstractString, lineno::Integer, message) =
    throw(ArgumentError("$source:$lineno: $message"))

function _parse_int(source::AbstractString, lineno::Integer, word::AbstractString)
    value = tryparse(Int, word)
    value === nothing && _parse_error(source, lineno, "`$word` is not an integer")
    return value
end

function _strip_comment(raw::AbstractString)
    cut = findfirst(==('#'), raw)
    cut === nothing && return strip(raw)
    return strip(raw[1:prevind(raw, cut)])
end

# `1: 2 4` or `1:` — the neighbour list may be empty.
function _parse_rotation_line(source::AbstractString, lineno::Integer, line::AbstractString)
    colon = findfirst(':', line)
    colon === nothing && _parse_error(source, lineno, "a rotation line is `v: w1 w2 ...`")
    head = line[1:prevind(line, colon)]
    tail = line[nextind(line, colon):end]
    vertex = _parse_int(source, lineno, strip(head))
    vertex >= 1 || _parse_error(source, lineno, "vertex labels start at 1")
    neighbours = [_parse_int(source, lineno, w) for w in split(tail)]
    all(w -> w != vertex, neighbours) ||
        _parse_error(source, lineno, "vertex $vertex lists itself as a neighbour")
    return vertex, neighbours
end

function _parse_edges_line(source::AbstractString, lineno::Integer, words)
    length(words) == 2 || _parse_error(source, lineno, "an edge is `u v`")
    u = _parse_int(source, lineno, words[1])
    w = _parse_int(source, lineno, words[2])
    (u >= 1 && w >= 1) || _parse_error(source, lineno, "vertex labels start at 1")
    u == w && _parse_error(source, lineno, "self-loops are not supported")
    return min(u, w), max(u, w)
end

# The state of a file being read: the vertex count so far, the edges, whether an
# `edges` section was seen, and the rotation with one flag per known vertex.
mutable struct _GraphParser
    n::Int
    declared::Int
    limit::Int
    edges::Vector{Tuple{Int,Int}}
    edge_set::Set{Tuple{Int,Int}}
    has_edges::Bool
    rotation::Union{Nothing,Vector{Vector{Int}}}
    given::BitVector
    section::Symbol
end

function _GraphParser(limit::Integer)
    return _GraphParser(
        0,
        -1,
        limit,
        Tuple{Int,Int}[],
        Set{Tuple{Int,Int}}(),
        false,
        nothing,
        falses(0),
        :none,
    )
end

function _check_limit(st::_GraphParser, source, lineno, count::Integer)
    count <= st.limit || _parse_error(
        source,
        lineno,
        "vertex count or label $count exceeds max_vertices = $(st.limit)",
    )
    return nothing
end

# Rotation lines may come in any order; every vertex may appear once.
function _parse_rotation!(st::_GraphParser, source, lineno, line)
    rotation = st.rotation::Vector{Vector{Int}}
    given = st.given
    vertex, neighbours = _parse_rotation_line(source, lineno, line)
    _check_limit(st, source, lineno, maximum(neighbours; init = vertex))
    while length(rotation) < vertex
        push!(rotation, Int[])
        push!(given, false)
    end
    given[vertex] && _parse_error(source, lineno, "duplicate rotation for vertex $vertex")
    given[vertex] = true
    rotation[vertex] = neighbours
    return maximum(neighbours; init = vertex)
end

function _parse_graph_file(io::IO, source::AbstractString, limit::Integer)
    st = _GraphParser(limit)
    for (lineno, raw) in enumerate(eachline(io))
        _consume_line!(st, source, lineno, _strip_comment(raw))
    end

    if st.declared >= 0
        st.n <= st.declared || _parse_error(
            source,
            0,
            "vertex $(st.n) exceeds the $(st.declared) vertices that `vertices` declares",
        )
        st.n = st.declared
    end
    rotation = st.rotation
    rotation === nothing && return _GraphFile(st.n, st.edges, st.has_edges, nothing)
    length(rotation) == st.n || _parse_error(
        source,
        0,
        "rotation has $(length(rotation)) lines for $(st.n) vertices; vertex " *
        "$(length(rotation) + 1) has no rotation line",
    )
    return _GraphFile(st.n, st.edges, st.has_edges, rotation)
end

# A keyword sets the state or opens a section; anything else belongs to the
# section that is open.
function _consume_line!(st::_GraphParser, source, lineno, line)
    isempty(line) && return nothing
    words = split(line)
    keyword = String(words[1])
    if keyword in ("vertices", "end", "edges", "rotation")
        return _consume_directive!(st, source, lineno, keyword, words)
    end
    st.section === :none &&
        _parse_error(source, lineno, "unexpected line `$line` (missing a section keyword?)")
    return _consume_data!(st, source, lineno, line, words)
end

function _consume_directive!(st::_GraphParser, source, lineno, keyword, words)
    if keyword == "vertices"
        length(words) == 2 || _parse_error(source, lineno, "expected `vertices <n>`")
        st.declared < 0 || _parse_error(source, lineno, "`vertices` is given twice")
        count = _parse_int(source, lineno, words[2])
        count >= 0 || _parse_error(source, lineno, "vertex count must be non-negative")
        _check_limit(st, source, lineno, count)
        st.declared = count
        return nothing
    end
    if keyword == "end"
        st.section = :none
        return nothing
    end
    length(words) == 1 || _parse_error(source, lineno, "`$keyword` takes no arguments")
    st.section = Symbol(keyword)
    if keyword == "rotation"
        st.rotation === nothing && (st.rotation = Vector{Vector{Int}}())
    else
        st.has_edges = true
    end
    return nothing
end

function _consume_data!(st::_GraphParser, source, lineno, line, words)
    if st.section === :edges
        u, w = _parse_edges_line(source, lineno, words)
        (u, w) in st.edge_set && _parse_error(source, lineno, "duplicate edge $u-$w")
        push!(st.edge_set, (u, w))
        push!(st.edges, (u, w))
        _check_limit(st, source, lineno, w)
        st.n = max(st.n, w)
        return nothing
    end
    st.n = max(st.n, _parse_rotation!(st, source, lineno, line))
    return nothing
end

function _read_file(path::AbstractString, parser; max_vertices::Integer)
    isfile(path) || throw(ArgumentError("no such file: $path"))
    return open(path) do io
        parser(io, path; max_vertices)
    end
end

# Both sections present: they must describe the same graph.
function _check_sections_agree(
    g::PlanarGraph,
    listed::Vector{Tuple{Int,Int}},
    source::AbstractString,
)
    from_rotation = Set((min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(g.graph))
    from_edges = Set(listed)
    from_edges == from_rotation || throw(
        ArgumentError(
            "$source: the `edges` and `rotation` sections describe different graphs",
        ),
    )
    return nothing
end

"""
    read_planar_graph(io; max_vertices = $DEFAULT_MAX_VERTICES) -> PlanarGraph
    read_planar_graph(path; max_vertices = $DEFAULT_MAX_VERTICES) -> PlanarGraph

Read a graph together with its embedding, as written by
[`write_planar_graph`](@ref). The rotation system goes through the `PlanarGraph`
constructor, so a malformed rotation raises an `ArgumentError` instead of
producing a silently broken graph; when the file has both sections they must
describe the same graph. A vertex count or label above `max_vertices` is refused
before anything is allocated, and a label above the count that `vertices`
declares is an error.
"""
function read_planar_graph(
    io::IO,
    source::AbstractString = "<io>";
    max_vertices::Integer = DEFAULT_MAX_VERTICES,
)
    parsed = _parse_graph_file(io, source, max_vertices)
    parsed.rotation === nothing && throw(
        ArgumentError(
            "$source: no `rotation` section; use read_graph for an unembedded graph",
        ),
    )
    g = PlanarGraph(parsed.rotation)
    parsed.has_edges && _check_sections_agree(g, parsed.edges, source)
    return g
end

function read_planar_graph(
    path::AbstractString;
    max_vertices::Integer = DEFAULT_MAX_VERTICES,
)
    return _read_file(path, read_planar_graph; max_vertices)
end

"""
    read_graph(io; max_vertices = $DEFAULT_MAX_VERTICES) -> SimpleGraph
    read_graph(path; max_vertices = $DEFAULT_MAX_VERTICES) -> SimpleGraph

Read the graph of a file written by [`write_graph`](@ref) or
[`write_planar_graph`](@ref), without its embedding. This is the way to obtain a
planar graph's underlying graph, which is the input of [`is_planar`](@ref) and
[`planar_embedding`](@ref).

The edges are those of the `edges` section; a file with only a `rotation`
section describes the graph through its rotations, and then the graph is the one
that rotation induces. `max_vertices` and the `vertices` line are handled as in
[`read_planar_graph`](@ref).
"""
function read_graph(
    io::IO,
    source::AbstractString = "<io>";
    max_vertices::Integer = DEFAULT_MAX_VERTICES,
)
    parsed = _parse_graph_file(io, source, max_vertices)
    if !parsed.has_edges && parsed.rotation !== nothing
        return PlanarGraph(parsed.rotation).graph
    end
    graph = SimpleGraph(parsed.n)
    for (u, w) in parsed.edges
        Graphs.add_edge!(graph, u, w)
    end
    return graph
end

function read_graph(path::AbstractString; max_vertices::Integer = DEFAULT_MAX_VERTICES)
    return _read_file(path, read_graph; max_vertices)
end
