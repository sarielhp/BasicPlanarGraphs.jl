# The Lipton-Tarjan separator: validity, balance, size, and an independent
# oracle on small graphs.

# Validate a partition using Graphs.jl only, so that this check does not share
# code with the algorithm under test.
function check_separator(g, A, B, C; balance = 2 / 3)
    n = nv(g)
    @test sort(vcat(A, B, C)) == collect(1:n)
    @test isempty(intersect(A, B))
    @test isempty(intersect(A, C))
    @test isempty(intersect(B, C))

    side = zeros(Int8, n)
    side[A] .= 1
    side[B] .= 2
    side[C] .= 3
    cross = 0
    for e in edges(g)
        s, t = side[src(e)], side[dst(e)]
        (s == 1 && t == 2) || (s == 2 && t == 1) ? (cross += 1) : nothing
    end
    @test cross == 0
    @test max(length(A), length(B)) <= balance * n
    @test length(C) <= 2 * sqrt(2 * n)
    return nothing
end

# Smallest separator that leaves every component within `balance * n`, by
# trying every vertex subset. Only for tiny graphs.
function minimum_balanced_separator(g; balance = 2 / 3)
    n = nv(g)
    n <= 16 || throw(ArgumentError("the oracle enumerates subsets, so n <= 16"))
    best = n
    for mask = 1:(2^n-1)
        removed = BitVector([mask >> (v - 1) & 1 == 1 for v = 1:n])
        count(removed) < best || continue
        all_components_fit(g, removed, balance) || continue
        best = count(removed)
    end
    return best
end

function all_components_fit(g, removed::BitVector, balance)
    n = nv(g)
    seen = copy(removed)
    for s = 1:n
        seen[s] && continue
        size = 0
        stack = Int[s]
        seen[s] = true
        while !isempty(stack)
            v = pop!(stack)
            size += 1
            size > balance * n && return false
            for w in neighbors(g, v)
                seen[w] && continue
                seen[w] = true
                push!(stack, w)
            end
        end
    end
    return true
end

@testset "separator" begin
    @testset "degenerate sizes" begin
        for n = 0:3
            g = PlanarGraph([Int[] for _ = 1:n])
            A, B, C = separator(g)
            check_separator(g, A, B, C)
        end
        single = PlanarGraph([[2], [1]])
        check_separator(single, separator(single)...)

        @test_throws ArgumentError separator(PlanarGraph([[2], [1]]); balance = 0.0)
        @test_throws ArgumentError separator(PlanarGraph([[2], [1]]); balance = 1.0)
    end

    @testset "structured graphs" begin
        cases = [
            planar_path(2),
            planar_path(10),
            planar_path(101),
            planar_cycle(3),
            planar_cycle(100),
            planar_grid(1, 9),
            planar_grid(2, 2),
            planar_grid(3, 3),
            planar_grid(12, 12),
            planar_grid(12, 12; diagonals = true),
            planar_grid(3, 40),
            planar_wheel(3),
            planar_wheel(50),
        ]
        for g in cases
            A, B, C = separator(g)
            check_separator(g, A, B, C)
        end
        # A path needs a single vertex, a wheel three, and neither is the whole graph.
        @test length(separator(planar_path(101))[3]) <= 2
        @test length(separator(planar_wheel(50))[3]) == 3
    end

    @testset "random triangulations" begin
        for n in (10, 50, 120, 400), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            A, B, C = separator(g)
            check_separator(g, A, B, C)
        end
    end

    @testset "trees and sparse graphs" begin
        for n in (5, 30, 120), seed = 1:2
            rng = MersenneTwister(seed)
            graph = SimpleGraph(n)
            for v = 2:n
                Graphs.add_edge!(graph, v, rand(rng, 1:(v-1)))
            end
            rot = [collect(neighbors(graph, v)) for v = 1:n]
            g = PlanarGraph(graph, rot)
            A, B, C = separator(g)
            check_separator(g, A, B, C)
        end
    end

    @testset "disconnected graphs" begin
        g = PlanarGraph([Int[] for _ = 1:5])
        check_separator(g, separator(g)...)

        joined = disjoint_union(planar_grid(4, 4), planar_wheel(5))
        check_separator(joined, separator(joined)...)

        # One component dominates, so the separator has to come from inside it.
        dominated = disjoint_union(planar_grid(12, 12), PlanarGraph([Int[] for _ = 1:3]))
        A, B, C = separator(dominated)
        check_separator(dominated, A, B, C)
        @test length(C) <= 12
    end

    @testset "guarantees on every input" begin
        rng = MersenneTwister(2024)
        for trial = 1:200
            n = rand(rng, 3:250)
            g = random_triangulation(n; rng)
            g = trial % 3 == 0 ? g : thin_out(g, rand(rng)^0.5, rng)
            A, B, C = separator(g)
            check_separator(g, A, B, C)
            @test length(C) <= 2 * sqrt(2 * nv(g))
        end
    end

    @testset "minimal separators" begin
        @test length(separator(planar_path(7))[3]) == 1
        @test length(separator(planar_cycle(40))[3]) == 2
        @test length(separator(planar_grid(30, 30))[3]) <= 30
        @test length(separator(planar_path(2))[3]) == 1
    end

    @testset "balance and embedding checks" begin
        @test_throws ArgumentError separator(planar_grid(5, 5); balance = 0.6)
        # Sorted neighbours are not a planar rotation of a grid.
        graph = planar_grid(4, 4).graph
        crossed = PlanarGraph(graph, [collect(neighbors(graph, v)) for v = 1:nv(graph)])
        @test_throws ArgumentError separator(crossed)
    end

    @testset "oracle on small graphs" begin
        for n = 4:8
            g = random_triangulation(n; rng = MersenneTwister(3 * n))
            A, B, C = separator(g)
            check_separator(g, A, B, C)
            # The oracle cannot beat the algorithm: its cut is optimal.
            @test length(C) >= minimum_balanced_separator(g)

            grid = planar_grid(2, n - 2)
            A, B, C = separator(grid)
            check_separator(grid, A, B, C)
            @test length(C) >= minimum_balanced_separator(grid)
        end
    end

    @testset "repeated and typed calls" begin
        g = random_triangulation(60; rng = MersenneTwister(5))
        first_run = separator(g)
        @test separator(g) == first_run

        small = planar_grid(6, 6)
        typed = PlanarGraph([Int32.(r) for r in rotation(small)])
        @test typed isa PlanarGraph{Int32}
        check_separator(small, separator(typed)...)

        looser = separator(planar_grid(9, 9); balance = 0.75)
        check_separator(planar_grid(9, 9), looser...; balance = 0.75)
    end

    @testset "large instance" begin
        g = planar_grid(60, 60; diagonals = true)
        A, B, C = separator(g)
        check_separator(g, A, B, C)

        g = random_triangulation(1000; rng = MersenneTwister(7))
        A, B, C = separator(g)
        check_separator(g, A, B, C)
    end
end

# Balance by weight, checked with Graphs.jl only: the weights of `A` and `B`, and of
# every component of the graph without `C`, are at most `balance` of the total.
function check_weighted_separator(g, (A, B, C), w; balance = 2 / 3)
    n = nv(g)
    @test sort(vcat(A, B, C)) == collect(1:n)
    side = zeros(Int8, n)
    side[A] .= 1
    side[B] .= 2
    side[C] .= 3
    @test all(
        e -> !(side[src(e)] + side[dst(e)] == 3 && side[src(e)] * side[dst(e)] == 2),
        edges(g),
    )
    @test max(sum(w[A]; init = 0), sum(w[B]; init = 0)) <= balance * sum(w)
    @test length(C) <= 2 * sqrt(2 * n)
    keep = trues(n)
    keep[C] .= false
    sub, ids = induced_subgraph(g.graph, findall(keep))
    for comp in connected_components(sub)
        @test sum(w[ids[comp]]) <= balance * sum(w)
    end
    return nothing
end

@testset "weighted separator" begin
    @testset "unit weights change nothing" begin
        for g in (planar_grid(12, 12), random_triangulation(200; rng = MersenneTwister(3)))
            @test separator(g; weights = ones(Int, nv(g))) == separator(g)
            @test separator(g; weights = fill(UInt8(1), nv(g))) == separator(g)
            @test separator(g; weights = ones(Int32, nv(g)), balance = 0.8) ==
                  separator(g; balance = 0.8)
        end
    end

    @testset "random weights" begin
        rng = MersenneTwister(21)
        graphs = PlanarGraph[
            planar_grid(15, 15),
            planar_wheel(40),
            planar_path(60),
            disjoint_union(planar_grid(6, 6), planar_wheel(9), planar_cycle(5)),
        ]
        for n in (5, 20, 100, 500), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            push!(graphs, g)
            push!(graphs, thin_out(g, 0.5, MersenneTwister(seed)))
        end
        for g in graphs, trial = 1:4
            w = rand(rng, trial == 1 ? (0:1) : trial == 2 ? (0:50) : (0:1000), nv(g))
            sum(w) == 0 && (w[1] = 1)
            A, B, C = separator(g; weights = w)
            check_weighted_separator(g, (A, B, C), w)
        end
    end

    @testset "weight concentrated in a corner" begin
        g = planar_grid(20, 20)
        w = [v <= 30 ? 500 : 1 for v = 1:nv(g)]
        A, B, C = separator(g; weights = w)
        check_weighted_separator(g, (A, B, C), w)
        # the unweighted separator is not balanced for these weights
        A0, B0, _ = separator(g)
        @test max(sum(w[A0]), sum(w[B0])) > 2 * sum(w) / 3
        for balance in (0.75, 0.9)
            A, B, C = separator(g; weights = w, balance)
            check_weighted_separator(g, (A, B, C), w; balance)
        end
    end

    @testset "heavy vertices" begin
        g = planar_grid(8, 8)
        w = ones(Int, nv(g))
        w[19] = 1000
        @test separator(g; weights = w)[3] == [19]
        w[40] = 800
        A, B, C = separator(g; weights = w)
        @test C == [19]
        check_weighted_separator(g, (A, B, C), w)
        # a heavy vertex in a small component of a disconnected graph
        d = disjoint_union(planar_grid(10, 10), planar_wheel(6))
        w = ones(Int, nv(d))
        w[nv(d)] = 400
        A, B, C = separator(d; weights = w)
        check_weighted_separator(d, (A, B, C), w)
        # one vertex of weight above two thirds, in a graph of two vertices
        A, B, C = separator(planar_path(2); weights = [5, 1])
        @test sort(vcat(A, B, C)) == [1, 2]
        @test 1 in C
    end

    @testset "zero and tiny cases" begin
        g = planar_grid(5, 5)
        A, B, C = separator(g; weights = zeros(Int, nv(g)))
        check_weighted_separator(g, (A, B, C), zeros(Int, nv(g)))
        @test separator(PlanarGraph([Int[]]); weights = [3]) == (Int[], Int[], [1])
        @test separator(PlanarGraph(Vector{Int}[]); weights = Int[]) ==
              (Int[], Int[], Int[])
    end

    @testset "weights are validated" begin
        g = planar_grid(4, 4)
        @test_throws ArgumentError separator(g; weights = ones(Int, 15))
        @test_throws ArgumentError separator(g; weights = fill(-1, 16))
        @test_throws ArgumentError separator(g; weights = ones(16))
        @test_throws ArgumentError separator(g; weights = fill(typemax(Int) ÷ 4, 16))
        @test_throws ArgumentError separator(g; weights = fill(typemax(UInt64), 16))
    end
end
