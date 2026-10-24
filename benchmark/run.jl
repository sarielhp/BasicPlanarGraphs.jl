#!/usr/bin/env julia
# Running-time survey of the public entry points that are documented as linear.
#
#   julia --project=. benchmark/run.jl                  # every suite, default sizes
#   julia --project=. benchmark/run.jl --quick          # small sizes, for a smoke test
#   julia --project=. benchmark/run.jl --suite=operations
#   julia --project=. benchmark/run.jl --suite=separator --reps=5
#   julia --project=. benchmark/run.jl --suite=phases
#
# Every measurement warms the call once and then reports the best of `--reps`
# runs, so compilation and GC pauses do not dominate. Timings are wall-clock
# and machine-dependent; nothing here is asserted by the test suite.

using BasicPlanarGraphs
using Graphs: ne, nv
using Printf
using Random: MersenneTwister

const BPG = BasicPlanarGraphs

const DEFAULT_SIZES = (; random = 20_000, grid = 200, wheel = 50_000)
const QUICK_SIZES = (; random = 5_000, grid = 100, wheel = 12_500)

function print_help()
    println("""
    Usage: julia --project=. benchmark/run.jl [options]

    Options:
      --suite=NAME   all (default), operations, separator, or phases
      --reps=N       repetitions, best one is reported (default 3)
      --seed=N       seed of random_triangulation (default 3)
      --quick        smaller graphs, for a fast smoke test
      -h, --help     show this message

    `operations` compares the public linear-time entry points on the same
    graph; `separator` reports its scaling on three graph families; `phases`
    breaks a single `separator` run into its internal stages.
    """)
    return nothing
end

function parse_args(args)
    suite = "all"
    reps = 3
    seed = 3
    quick = false
    for arg in args
        if arg == "-h" || arg == "--help"
            print_help()
            exit(0)
        end
        if arg == "--quick"
            quick = true
        elseif startswith(arg, "--suite=")
            suite = split(arg, '='; limit = 2)[2]
        elseif startswith(arg, "--reps=")
            reps = parse(Int, split(arg, '='; limit = 2)[2])
        elseif startswith(arg, "--seed=")
            seed = parse(Int, split(arg, '='; limit = 2)[2])
        else
            error("unknown argument: $arg (try --help)")
        end
    end
    suite in ("all", "operations", "separator", "phases") ||
        error("unknown suite: $suite (try --help)")
    reps >= 1 || error("--reps must be at least 1")
    return (; suite, reps, seed, quick)
end

function best_of(f, reps::Int)
    f()
    t = Inf
    for _ = 1:reps
        t = min(t, @elapsed f())
    end
    return t
end

function row(name, g, f, reps)
    t = try
        best_of(f, reps)
    catch err
        err isa ArgumentError || rethrow()
        @printf("  %-26s %10s\n", name, "n/a")
        return
    end
    edges = max(ne(g), 1)
    @printf("  %-26s %9.1f ms  %8.1f ns/edge\n", name, 1e3t, 1e9t / edges)
    return
end

function operations(g, reps, title)
    println("== operations, ", title, " (n=", nv(g), ", m=", ne(g), ")")
    row("faces", g, () -> faces(g), reps)
    row("nfaces", g, () -> nfaces(g), reps)
    row("is_planar_embedding", g, () -> is_planar_embedding(g), reps)
    row("dual", g, () -> dual(g), reps)
    row("separator", g, () -> separator(g), reps)
    row("planar_embedding(PG)", g, () -> planar_embedding(g), reps)
    row("triangulate", g, () -> triangulate(g), reps)
    row("grid_drawing", g, () -> grid_drawing(g), reps)
    row("cycle_separator", g, () -> cycle_separator(g), reps)
    sg = g.graph
    row("planar_embedding(Simple)", sg, () -> planar_embedding(sg), reps)
    row("is_planar_linear(Simple)", sg, () -> is_planar_linear(sg), reps)
    println()
    return
end

function separator_scaling(sizes, seed, reps)
    println("== separator scaling, best of ", reps)
    families = (
        ("random triangulation", n -> random_triangulation(n; rng = MersenneTwister(seed))),
        ("grid", n -> planar_grid(round(Int, sqrt(n)), round(Int, sqrt(n)))),
        ("wheel", n -> planar_wheel(n - 1)),
    )
    for (name, make) in families
        println("-- ", name)
        for n in sizes
            g = make(n)
            t = best_of(() -> separator(g), reps)
            @printf(
                "  n=%6d m=%7d  %9.1f ms  %7.1f ns/vertex\n",
                nv(g),
                ne(g),
                1e3t,
                1e9t / nv(g),
            )
        end
    end
    println()
    return
end

# The band triangulation is the dominant stage on low-diameter graphs. This
# suite reaches into internals on purpose: it explains the numbers above.
function separator_phases(n, seed, reps)
    g = random_triangulation(n; rng = MersenneTwister(seed))
    dg = BPG.DartGraph(g)
    println("== separator phases, random triangulation (n=", dg.n, ")")
    for root in BPG._root_candidates(dg, 1)
        levels = BPG.level_structure(dg, [root])
        l0, l2 = BPG._level_plan(levels.depth, dg.n)
        keep = BitVector(l0 < levels.depth[v] < l2 for v = 1:dg.n)
        band = BPG.restrict(dg, keep)
        ok = Ref(true)
        t_tri = best_of(reps) do
            ok[] = true
            try
                BPG.triangulate(band)
            catch err
                err isa ArgumentError || rethrow()
                ok[] = false
            end
        end
        t_band = best_of(() -> BPG._band_cycles(dg, levels, l0, l2), reps)
        println(
            "  root=",
            root,
            "  |band|=",
            sum(keep),
            "  triangulate=",
            ok[] ? "ok" : "failed",
            "  tri=",
            round(1e3t_tri; digits = 1),
            " ms",
            "  band_cycles=",
            round(1e3t_band; digits = 1),
            " ms",
        )
    end
    println()
    return
end

function main()
    opts = parse_args(ARGS)
    sizes = opts.quick ? QUICK_SIZES : DEFAULT_SIZES
    @printf(
        "BasicPlanarGraphs running-time survey | Julia %s | %s, %d CPU(s)\n\n",
        VERSION,
        Sys.CPU_NAME,
        Sys.CPU_THREADS,
    )
    if opts.suite in ("all", "operations")
        operations(
            random_triangulation(sizes.random; rng = MersenneTwister(opts.seed)),
            opts.reps,
            "random triangulation",
        )
        operations(planar_grid(sizes.grid, sizes.grid), opts.reps, "grid")
        operations(planar_wheel(sizes.wheel), opts.reps, "wheel")
    end
    if opts.suite in ("all", "separator")
        sizes = opts.quick ? (2500, 5000, 10000) : (5000, 10000, 20000, 40000)
        separator_scaling(sizes, opts.seed, opts.reps)
    end
    if opts.suite == "phases"
        separator_phases(opts.quick ? 10000 : 40000, opts.seed, opts.reps)
    end
    return
end

main()
