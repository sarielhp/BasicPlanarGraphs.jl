# The Har-Peled / Nayyeri cycle separator: the result must be a closed walk of
# the input graph (a simple cycle when the graph is a triangulation) whose two
# sides are balanced. The checks use Graphs.jl primitives, so they do not share
# code with the implementation.

# Check that `C` is a closed walk of `g` (simple when `simple` is set), that the
# two sides form a partition with no edge crossing the walk, and that every
# component of the graph without the walk is at most `balance * n`.
function check_cycle(g, C; balance = 2 / 3, simple = false)
    n = nv(g)
    @test C isa Vector{Int}
    @test length(C) >= 3
    simple && @test allunique(C)
    @test all(v -> 1 <= v <= n, C)
    for i in eachindex(C)
        @test has_edge(g, C[i], C[mod1(i + 1, length(C))])
    end

    inside, outside = cycle_sides(g, C)
    @test sort(vcat(inside, outside, unique(C))) == collect(1:n)
    @test isempty(intersect(inside, outside))
    @test isempty(intersect(inside, C))
    @test isempty(intersect(outside, C))

    side = zeros(Int8, n)
    side[inside] .= 1
    side[outside] .= 2
    side[C] .= 3
    for e in edges(g)
        s, t = side[src(e)], side[dst(e)]
        @test !((s == 1 && t == 2) || (s == 2 && t == 1))
    end

    @test max(length(inside), length(outside)) <= balance * n
    keep = trues(n)
    keep[C] .= false
    sub, _ = induced_subgraph(g.graph, findall(keep))
    components = connected_components(sub)
    @test isempty(components) || maximum(length, components) <= balance * n
    # cycle_separator reports the smaller side as the inside
    @test length(inside) <= length(outside)
    return inside, outside
end

@testset "cycle separator" begin
    @testset "structured graphs" begin
        graphs = PlanarGraph[
            planar_cycle(3),
            planar_cycle(4),
            planar_cycle(10),
            planar_cycle(100),
            planar_wheel(3),
            planar_wheel(10),
            planar_wheel(100),
            planar_grid(2, 2),
            planar_grid(2, 9),
            planar_grid(3, 3),
            planar_grid(4, 7),
            planar_grid(10, 10),
            planar_grid(20, 20),
            planar_grid(3, 40),
            planar_grid(12, 12; diagonals = true),
        ]
        for g in graphs
            check_cycle(g, cycle_separator(g))
        end
    end

    @testset "random triangulations" begin
        for n in (5, 6, 10, 30, 100, 400), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            C = cycle_separator(g)
            check_cycle(g, C; simple = true)
            # On a triangulation the walk is a cycle and the paper's bound
            # applies to it.
            @test length(C) <= 2 * sqrt(2 * n) + 1
        end
    end

    @testset "phases 2 and 3 with a small step" begin
        # With the default step the fundamental cycle of a small triangulation is
        # already short, so the ladder is never reached; a small step forces it.
        for n in (30, 100, 400), delta in (2, 3, 5), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            C = BasicPlanarGraphs._cycle_separator(BasicPlanarGraphs.DartGraph(g); delta)
            check_cycle(g, C; simple = true)
        end
    end

    @testset "a heavy ring outside the last rung" begin
        # Every rung of the ladder is light, because more than two thirds of the
        # faces lie beyond the last one; the paper cuts the outermost ring. Without
        # that case the answer was the fundamental cycle, of length 2L + 4.
        for (L, q) in ((50, 125), (200, 500), (800, 2000))
            g = tube_blob(L, q)
            @test ne(g) == 3nv(g) - 6
            C = cycle_separator(g)
            check_cycle(g, C; simple = true)
            @test length(C) <= 2 * sqrt(2 * nv(g)) + 1
        end
    end

    @testset "cycle_sides convention" begin
        g = planar_wheel(10)
        C = cycle_separator(g)
        inside, outside = cycle_sides(g, C)
        # A counter-clockwise cycle has the inside on its left; reversing the
        # traversal swaps the two sides.
        @test cycle_sides(g, reverse(C)) == (outside, inside)
        @test cycle_sides(g, collect(C)) == (inside, outside)
        @test cycle_sides(g, Int32.(C)) == (inside, outside)
        @test_throws ArgumentError cycle_sides(g, Float64.(C))
        @test_throws ArgumentError cycle_sides(g, Any[])
    end

    @testset "cycle_sides validation" begin
        g = planar_grid(4, 4)
        n = nv(g)
        @test_throws ArgumentError cycle_sides(g, [1, 2])
        @test_throws ArgumentError cycle_sides(g, [1, 1, 2])
        @test_throws ArgumentError cycle_sides(g, [1, 2, 3])     # 3-1 is not an edge
        @test_throws ArgumentError cycle_sides(g, [1, 2, n + 1])
        @test_throws ArgumentError cycle_sides(g, Int[])
        @test_throws ArgumentError cycle_sides(PlanarGraph([[2], [1], Int[]]), [1, 2, 3])
    end

    @testset "preconditions" begin
        @test_throws ArgumentError cycle_separator(planar_cycle(2))
        # not connected
        @test_throws ArgumentError cycle_separator(PlanarGraph([[2], [1], Int[]]))
        # not planar: sorted neighbours are not a planar rotation of a grid
        graph = planar_grid(4, 4).graph
        crossed = PlanarGraph(graph, [collect(neighbors(graph, v)) for v = 1:16])
        @test_throws ArgumentError cycle_separator(crossed)
    end

    @testset "graphs with bridges, cut vertices and large faces" begin
        rng = MersenneTwister(31)
        tested = 0
        for trial = 1:300
            g = random_triangulation(rand(rng, 3:120); rng)
            g = thin_out(g, rand(rng) * 0.85, rng)
            is_connected(g.graph) || continue
            check_cycle(g, cycle_separator(g))
            tested += 1
        end
        @test tested > 100
        # a path, a star-like tree and a triangle with a pendant edge
        for g in (
            planar_path(3),
            planar_path(40),
            PlanarGraph([[4], [3, 4], [4, 2], [1, 2, 3]]),
            planar_grid(2, 60),
        )
            check_cycle(g, cycle_separator(g))
        end
    end

    @testset "closed walks" begin
        # A triangle with a pendant edge: the walk goes out along the spike and
        # back, so every vertex is on it and both sides are empty.
        g = PlanarGraph([[4], [3, 4], [4, 2], [1, 2, 3]])
        @test cycle_sides(g, [1, 4, 3, 2, 4]) == (Int[], Int[])
        # A spike into a region puts its end on the walk and leaves the rest alone.
        grid = planar_grid(3, 3)
        ring = [1, 2, 3, 6, 9, 8, 7, 4]
        @test cycle_sides(grid, ring) == (Int[], [5])
        @test cycle_sides(grid, [ring; 5; 4]) == (Int[], Int[])
        # A walk that goes around twice crosses itself and is rejected.
        @test_throws ArgumentError cycle_sides(grid, [1, 2, 5, 4, 1, 2, 5, 4])
    end

    @testset "repeated calls" begin
        g = planar_grid(9, 9)
        C = cycle_separator(g)
        @test cycle_separator(g) == C
        @test cycle_separator(g) == C
    end

    @testset "forced ladder step" begin
        # Exercise the ladder branch of phase 3, which the default step skips
        # on small graphs; the ladder only runs on triangulations.
        for g in (
            triangulate(planar_grid(12, 12)),
            planar_wheel(60),
            random_triangulation(150; rng = MersenneTwister(2)),
        )
            dg = BasicPlanarGraphs.DartGraph(g)
            for delta in (3, 5, 8)
                C = BasicPlanarGraphs._cycle_separator(dg; delta = delta)
                check_cycle(g, C; simple = true)
            end
        end
    end
end

@testset "walks that return to a hub many times" begin
    # The darts of a walk used to be found by rescanning the hub at every visit:
    # a star, and triangles around one vertex, took quadratic time.
    star(k) = PlanarGraph(vcat([collect(2:k)], [[1] for _ = 2:k]))
    windmill(k) = begin
        rot = [Int[] for _ = 1:(2k+1)]
        for i = 1:k
            a, b = 2i, 2i + 1
            append!(rot[1], (a, b))
            rot[a] = [1, b]
            rot[b] = [a, 1]
        end
        PlanarGraph(rot)
    end
    for make in (star, windmill)
        for k in (6, 7, 50)
            g = make(k)
            c = cycle_separator(g)
            inside, outside = cycle_sides(g, c)
            @test max(length(inside), length(outside)) <= 2nv(g) / 3
        end
        small, large = make(2_000), make(32_000)
        @test fastest(cycle_separator, large) <= 80 * fastest(cycle_separator, small) + 0.05
    end
end
