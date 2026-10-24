#!/usr/bin/env julia
# ==============================================================================
#  examples/test_planarity.jl
#
#  Planarity, embedding and separator on a real road network.
#
#  The graph is the Pennsylvania road network from the SNAP collection
#  ("roadNet-PA", built from TIGER/Line): 1,088,092 intersections and
#  1,541,898 road segments. It is a genuine planar graph — roads cross on
#  bridges and overpasses, but a crossing is not an intersection, so no two
#  edges meet at a vertex that is not their own — and it is the kind of input
#  that breaks implementations: the DFS forest has paths hundreds of thousands
#  of vertices long.
#
#  What the program does, in order:
#
#    1. download `<dataset>.txt.gz` into the cache, once, and gunzip it;
#    2. clean it into a pure graph: undirected, no self-loops, no duplicate
#       edges, vertices relabelled 1..n, and *no embedding* — written with
#       `write_graph`, which emits only `vertices` and `edges`;
#    3. read that file back with `read_graph`;
#    4. for each requested size — by default the first 1000, 10,000, 100,000
#       induced vertices and then the whole graph — check planarity, compute an
#       embedding when the graph is planar, write the embedding out, and run the
#       planar separator, verifying the partition and reporting its sizes;
#    5. print one summary table with every running time.
#
#  Sizes are induced subgraphs on vertices 1..k, where labels are assigned in
#  increasing order of the original SNAP ids, so the whole graph is the limit of
#  the sequence as k grows.
#
#  Usage:
#      julia --project=. examples/test_planarity.jl
#      julia --project=. examples/test_planarity.jl --sizes 1000,10000 --dataset roadNet-PA
#      julia --project=. examples/test_planarity.jl --help
#
#  Cache: by default $TMPDIR/<user>/BasicPlanarGraphs, holding the download, the
#  cleaned graph and the written embeddings. Nothing outside the cache is
#  touched.
# ==============================================================================

using Graphs, BasicPlanarGraphs
using Printf: @sprintf

# --- datasets --------------------------------------------------------------

const DATASETS = Dict(
    "roadNet-PA" => (
        url = "https://snap.stanford.edu/data/roadNet-PA.txt.gz",
        what = "Pennsylvania road network (TIGER/Line)",
        nodes = 1_088_092,
        edges = 1_541_898,
    ),
    "roadNet-TX" => (
        url = "https://snap.stanford.edu/data/roadNet-TX.txt.gz",
        what = "Texas road network (TIGER/Line)",
        nodes = 1_379_917,
        edges = 1_921_660,
    ),
    "roadNet-CA" => (
        url = "https://snap.stanford.edu/data/roadNet-CA.txt.gz",
        what = "California road network (TIGER/Line)",
        nodes = 1_965_206,
        edges = 2_766_607,
    ),
)

const VALUE_FLAGS =
    ("--cache", "--dataset", "--sizes", "--embeddings", "--check", "--graph")

default_cache() = joinpath(tempdir(), get(ENV, "USER", "user"), "BasicPlanarGraphs")

# --- options ---------------------------------------------------------------

mutable struct Options
    cache::String
    dataset::String
    sizes::Vector{Union{Int,Symbol}}
    embedding_dir::String
    check_with_graphs::Int
    graph::String   # a local graph file; skips the download when set
    label::String   # names the files written for this input
end

function usage(io::IO = stdout)
    print(
        io,
        """
        Planarity, embedding and separator on a real road network.

          --cache DIR        where downloads, the clean graph and the embeddings live
                             (default $(default_cache()))
          --dataset NAME     $(join(sort(collect(keys(DATASETS))), ", "))
                             (default roadNet-PA)
          --sizes LIST       comma separated vertex counts, `all` for the whole graph
                             (default 1000,10000,100000,all)
          --embeddings DIR   where the rotation systems are written
                             (default <cache>/embeddings)
          --graph FILE       read this graph file (the package format) instead of
                             downloading a dataset; used by the gate's smoke test
          --check N          also run Graphs.is_planar up to N vertices (default 10000,
                             0 disables); that implementation recurses and is
                             quadratic, and it kills the session with a
                             StackOverflowError on larger inputs
          --help             this text
        """,
    )
    return nothing
end

function default_options()
    return Options(
        default_cache(),
        "roadNet-PA",
        Union{Int,Symbol}[1000, 10_000, 100_000, :all],
        "",
        10_000,
        "",
        "roadNet-PA",
    )
end

function set_option!(opts::Options, flag::AbstractString, value::AbstractString)
    if flag == "--cache"
        opts.cache = value
    elseif flag == "--dataset"
        haskey(DATASETS, value) || error("unknown dataset $value (try --help)")
        opts.dataset = value
    elseif flag == "--sizes"
        opts.sizes = Union{Int,Symbol}[_size(s) for s in split(value, ",") if !isempty(s)]
    elseif flag == "--embeddings"
        opts.embedding_dir = value
    elseif flag == "--graph"
        opts.graph = value
        opts.label = replace(basename(value), r"\..*$" => "")
    else
        opts.check_with_graphs = parse(Int, value)
    end
    return opts
end

_size(token::AbstractString) = token == "all" ? :all : parse(Int, token)

function parse_options(args)
    opts = default_options()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--help", "-h")
            usage()
            exit(0)
        end
        arg in VALUE_FLAGS || error("unknown option $arg (try --help)")
        i == length(args) && error("$arg needs a value")
        set_option!(opts, arg, args[i+1])
        i += 2
    end
    isempty(opts.embedding_dir) && (opts.embedding_dir = joinpath(opts.cache, "embeddings"))
    return opts
end

# --- download and clean ----------------------------------------------------

"""Download `url` to `path`, unless it is already there; returns `path`."""
function fetch(url::AbstractString, path::AbstractString)
    isfile(path) && return path
    mkpath(dirname(path))
    println("downloading ", url)
    t0 = time()
    download(url, path)
    println(
        "  -> ",
        path,
        "  ",
        round(filesize(path) / 2^20, digits = 1),
        " MiB in ",
        round(time() - t0, digits = 1),
        "s",
    )
    return path
end

"""Gunzip `gz` into `txt` with the system `gzip`; returns `txt`."""
function gunzip(gz::AbstractString, txt::AbstractString)
    isfile(txt) && return txt
    Sys.which("gzip") === nothing && error("`gzip` is needed to expand $gz")
    println("expanding ", basename(gz))
    t0 = time()
    open(txt, "w") do io
        run(pipeline(`gzip -dc $gz`, stdout = io))
    end
    println(
        "  -> ",
        txt,
        "  ",
        round(filesize(txt) / 2^20, digits = 1),
        " MiB in ",
        round(time() - t0, digits = 1),
        "s",
    )
    return txt
end

# One line of a SNAP edge list: two ids, separated by a tab.
function edge_ids(line::AbstractString)
    parts = split(line)
    length(parts) == 2 || error("malformed edge line: $(repr(line))")
    return parse(Int32, parts[1]), parse(Int32, parts[2])
end

"""
Read a SNAP edge list into the set of its unordered pairs, counting self-loops.
"""
function read_pairs(raw::AbstractString)
    pairs = Set{Tuple{Int32,Int32}}()
    self_loops = 0
    records = 0
    for line in eachline(raw)
        (isempty(line) || startswith(line, '#')) && continue
        u, v = edge_ids(line)
        records += 1
        if u == v
            self_loops += 1
        else
            push!(pairs, u < v ? (u, v) : (v, u))
        end
    end
    return pairs, records, self_loops
end

# Relabel arbitrary ids to 1:n, keeping their order.
function relabelled(pairs)
    ids = sort!(unique!(reduce(vcat, [[u, v] for (u, v) in pairs])))
    index = Dict{Int32,Int}(id => i for (i, id) in enumerate(ids))
    g = SimpleGraph(length(ids))
    for (u, v) in pairs
        add_edge!(g, index[u], index[v])
    end
    return g
end

"""
Clean a SNAP edge list into a pure undirected simple graph.

The file lists both directions of every road and uses arbitrary ids; the clean
version keeps each unordered pair once, drops self-loops, relabels the vertices
`1:n` in increasing order of the original ids, and is written with `write_graph`,
so it carries no embedding at all.
"""
function clean_graph(raw::AbstractString, clean::AbstractString)
    if isfile(clean)
        println("using cached clean graph ", clean)
        return clean
    end
    println("cleaning ", basename(raw), " into a pure graph")
    t0 = time()
    pairs, records, self_loops = read_pairs(raw)
    g = relabelled(pairs)
    write_graph(clean, g)
    println(
        "  ",
        records,
        " edge records -> ",
        ne(g),
        " edges, ",
        nv(g),
        " vertices; ",
        self_loops,
        " self-loops dropped; ",
        round(time() - t0, digits = 1),
        "s",
    )
    println(
        "  -> ",
        clean,
        "  ",
        round(filesize(clean) / 2^20, digits = 1),
        " MiB (pure graph, no rotation)",
    )
    return clean
end

# --- verification ----------------------------------------------------------

# Every vertex of 1:n in exactly one part.
function partition_problems(n::Integer, parts)
    problems = String[]
    seen = falses(n)
    duplicates = 0
    for part in parts, v in part
        (1 <= v <= n) || (push!(problems, "vertex $v out of range"); continue)
        seen[v] && (duplicates += 1)
        seen[v] = true
    end
    duplicates == 0 || push!(problems, "$duplicates vertices listed twice")
    all(seen) || push!(problems, "$(count(!, seen)) vertices listed in no part")
    return problems
end

# Edges with one endpoint in `A` and the other in `B`; there must be none.
function crossing_edges(g, A, B)
    in_b = falses(nv(g))
    for v in B
        1 <= v <= nv(g) && (in_b[v] = true)
    end
    return count(v -> in_b[v], Iterators.flatten(neighbors(g, v) for v in A))
end

"""
Check that `(A, B, C)` really is a separator of `g`, and return what is worth
printing: the three sizes, the number of `A`–`B` edges (which must be zero), the
balance, and any complaint.
"""
function separator_facts(g, A, B, C)
    problems = partition_problems(nv(g), (A, B, C))
    crossing = crossing_edges(g, A, B)
    crossing == 0 || push!(problems, "$crossing edges join A and B")
    balance = maximum((length(A), length(B))) / max(nv(g), 1)
    return (
        a = length(A),
        b = length(B),
        c = length(C),
        crossing = crossing,
        balance = balance,
        problems = problems,
    )
end

fmt_time(t) = @sprintf("%7.2f s", t)
fmt_mem(bytes) = @sprintf("%6.0f MiB", bytes / 2^20)
fmt_ratio(x) = @sprintf("%5.1f%%", 100 * x)
fmt_secs(x) = isnan(x) ? "-" : @sprintf("%.2f", x)

# --- the four stages -------------------------------------------------------

"""Check planarity, and optionally compare with the recursive `Graphs.is_planar`."""
function report_planarity(h, opts)
    check = @timed is_planar_linear(h)
    println(
        "planarity   is_planar_linear      ",
        lpad(string(check.value), 5),
        "   ",
        fmt_time(check.time),
        "  ",
        fmt_mem(check.bytes),
        " allocated",
    )
    flush(stdout)
    if opts.check_with_graphs > 0 && nv(h) <= opts.check_with_graphs
        other = @timed is_planar(h)
        verdict = other.value == check.value ? "" : "   <-- DISAGREES"
        println(
            "            is_planar (Graphs)   ",
            lpad(string(other.value), 5),
            "   ",
            fmt_time(other.time),
            "  ",
            fmt_mem(other.bytes),
            " allocated",
            verdict,
        )
    elseif opts.check_with_graphs > 0
        println(
            "            is_planar (Graphs)   skipped: recursive, overflows the stack ",
            "at this size",
        )
    end
    return check.value, check.time
end

"""Compute the embedding, verify it, and write it next to the cache."""
function report_embedding(h, opts)
    embedded = @timed planar_embedding(h)
    e = embedded.value
    faces = nfaces(e)
    valid = is_planar_embedding(e)
    println(
        "embedding   planar_embedding      ",
        lpad(string(faces), 5),
        " faces  ",
        fmt_time(embedded.time),
        "  ",
        fmt_mem(embedded.bytes),
        " allocated",
    )
    println(
        "            is_planar_embedding   ",
        lpad(string(valid), 5),
        "   (Euler on every component)",
    )
    file = joinpath(opts.embedding_dir, "$(opts.label).n$(nv(h)).embedding.txt")
    written = @timed write_planar_graph(file, e)
    println("output      ", file)
    println(
        "            ",
        round(filesize(file) / 2^20, digits = 1),
        " MiB, ",
        fmt_time(written.time),
        "  ",
        fmt_mem(written.bytes),
        " allocated",
    )
    sample = join(("$(v): $(join(e.rotation[v], " "))" for v = 1:min(3, nv(h))), ", ")
    println("            rotation starts with  ", sample)
    flush(stdout)
    return (
        time = embedded.time,
        faces = faces,
        write = written.time,
        valid = valid,
        embedding = e,
    )
end

"""Run the separator, verify the partition, and print the sizes."""
function report_separator(e)
    separated = @timed separator(e)
    A, B, C = separated.value
    facts = separator_facts(e, A, B, C)
    n = nv(e)
    println(
        "separator   separator(sqrt(n) = ",
        round(Int, sqrt(n)),
        ")   ",
        fmt_time(separated.time),
        "  ",
        fmt_mem(separated.bytes),
        " allocated",
    )
    println(
        @sprintf(
            "            |A| = %d (%s)   |B| = %d (%s)   |C| = %d (%s), A-B edges = %d",
            facts.a,
            fmt_ratio(facts.a / n),
            facts.b,
            fmt_ratio(facts.b / n),
            facts.c,
            fmt_ratio(facts.c / n),
            facts.crossing
        ),
    )
    if isempty(facts.problems)
        println(
            "            verified: A, B, C partition 1:n, no A-B edge, balance ",
            fmt_ratio(facts.balance),
            " <= 66.7%",
        )
    else
        println("            FAILED: ", join(facts.problems, "; "))
    end
    flush(stdout)
    return (time = separated.time, facts = facts, valid = isempty(facts.problems))
end

# --- one size --------------------------------------------------------------

"""
Run the whole pipeline on the induced subgraph of the first `k` vertices:
planarity, embedding, and separator. Returns a row for the summary table.
"""
function run_size(g, k, opts)
    whole = k === :all || k >= nv(g)
    h = whole ? g : first(induced_subgraph(g, 1:k))
    n, m = nv(h), ne(h)

    println()
    println("="^78)
    println(whole ? "the whole graph" : "vertices 1:$n (induced subgraph)")
    println(@sprintf("n = %d, m = %d, average degree %.2f", n, m, 2m / max(n, 1)))
    println("="^78)
    flush(stdout)

    planar, check_time = report_planarity(h, opts)
    row = (
        label = whole ? "all" : string(k),
        n = n,
        m = m,
        planar = planar,
        check = check_time,
        embed = NaN,
        faces = 0,
        write = NaN,
        sep = NaN,
        a = 0,
        b = 0,
        c = 0,
        balance = NaN,
        valid = false,
    )
    if !planar
        println("not planar: no embedding, no separator")
        return row
    end

    embedding = report_embedding(h, opts)
    separation = report_separator(embedding.embedding)
    return merge(
        row,
        (
            embed = embedding.time,
            faces = embedding.faces,
            write = embedding.write,
            sep = separation.time,
            a = separation.facts.a,
            b = separation.facts.b,
            c = separation.facts.c,
            balance = separation.facts.balance,
            valid = separation.valid && embedding.valid,
        ),
    )
end

# --- main ------------------------------------------------------------------

function report_dataset(g, info)
    println(
        "loaded      read_graph            n = ",
        nv(g),
        ", m = ",
        ne(g),
        "   ",
        "expected from SNAP: n = ",
        info.nodes,
        ", m = ",
        info.edges,
        nv(g) == info.nodes && ne(g) == info.edges ? "   (matches)" : "   (differs)",
    )
    return nothing
end

function print_summary(rows, opts)
    println()
    println("="^78)
    println("summary")
    println("="^78)
    println(
        rpad("vertices", 10),
        lpad("n", 10),
        lpad("m", 11),
        lpad("planar", 8),
        lpad("check", 9),
        lpad("embed", 9),
        lpad("write", 9),
        lpad("separator", 11),
        lpad("|C|", 8),
        lpad("balance", 9),
        lpad("verified", 10),
    )
    for r in rows
        println(
            rpad(r.label, 10),
            lpad(r.n, 10),
            lpad(r.m, 11),
            lpad(string(r.planar), 8),
            lpad(fmt_secs(r.check), 9),
            lpad(fmt_secs(r.embed), 9),
            lpad(fmt_secs(r.write), 9),
            lpad(fmt_secs(r.sep), 11),
            lpad(r.c, 8),
            lpad(isnan(r.balance) ? "-" : fmt_ratio(r.balance), 9),
            lpad(string(r.valid), 10),
        )
    end
    println(
        "\ntimes in seconds; peak resident memory ",
        round(Sys.maxrss() / 2^30, digits = 1),
        " GiB",
    )
    println("embeddings written to ", opts.embedding_dir)
    return nothing
end

"""The pure graph to work on: the local file when given, else the cleaned dataset."""
function load_input(opts)
    if !isempty(opts.graph)
        println("graph:   ", opts.graph, " (local file)")
        loaded = @timed read_graph(opts.graph)
        println(
            "loaded      read_graph            ",
            fmt_time(loaded.time),
            "  ",
            fmt_mem(loaded.bytes),
            " allocated",
        )
        return loaded.value
    end

    info = DATASETS[opts.dataset]
    println("dataset: ", opts.dataset, " — ", info.what)
    println("note:    this is a planar graph in the abstract sense; roads cross on")
    println("         bridges, but a crossing is not an intersection.")
    gz = fetch(info.url, joinpath(opts.cache, "$(opts.dataset).txt.gz"))
    raw = gunzip(gz, joinpath(opts.cache, "$(opts.dataset).raw.txt"))
    clean = clean_graph(raw, joinpath(opts.cache, "$(opts.dataset).graph.txt"))

    loaded = @timed read_graph(clean)
    println(
        "loaded      read_graph            ",
        fmt_time(loaded.time),
        "  ",
        fmt_mem(loaded.bytes),
        " allocated",
    )
    report_dataset(loaded.value, info)
    return loaded.value
end

function main(args)
    opts = parse_options(args)
    mkpath(opts.cache)
    mkpath(opts.embedding_dir)
    println("cache:   ", opts.cache)
    g = load_input(opts)

    rows = Any[]
    for k in opts.sizes
        push!(rows, run_size(g, k, opts))
        GC.gc()
    end
    print_summary(rows, opts)
    return rows
end

main(ARGS)
