# The Har-Peled / Nayyeri cycle separator: the result must be a genuine simple
# cycle of the input graph whose two sides are balanced. All checks use
# Graphs.jl primitives, so they do not share code with the implementation.

# Check that `C` is a simple cycle of `g` and that the two sides of the cycle
# form a partition with no edge crossing the cycle.
function check_cycle(g, C; balance = 2 / 3)
    n = nv(g)
    @test C isa Vector{Int}
    @test length(C) >= 3
    @test allunique(C)
    @test all(v -> 1 <= v <= n, C)
    for i in eachindex(C)
        @test has_edge(g, C[i], C[mod1(i + 1, length(C))])
    end

    inside, outside = cycle_sides(g, C)
    @test sort(vcat(inside, outside, C)) == collect(1:n)
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
            check_cycle(g, C)
            # On a triangulation no shortcut is needed, so the paper's bound
            # applies to the returned cycle.
            @test length(C) <= 2 * sqrt(2 * n) + 1
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
        # not 2-connected
        @test_throws ArgumentError cycle_separator(PlanarGraph([[2], [1, 3], [2]]))
        @test_throws ArgumentError cycle_separator(
            PlanarGraph([[2, 3], [3], [1, 2], Int[]]),
        )
    end

    @testset "repeated calls" begin
        g = planar_grid(9, 9)
        C = cycle_separator(g)
        @test cycle_separator(g) == C
        @test cycle_separator(g) == C
    end

    @testset "forced ladder step" begin
        # Exercise the ladder branch of phase 3, which the default step skips
        # on small graphs.
        for g in (planar_grid(12, 12), planar_wheel(60))
            dg = BasicPlanarGraphs.DartGraph(g)
            for delta in (3, 5, 8)
                C = BasicPlanarGraphs._cycle_separator(dg; delta = delta)
                @test length(C) >= 3
                @test allunique(C)
                for i in eachindex(C)
                    @test has_edge(g, C[i], C[mod1(i + 1, length(C))])
                end
            end
        end
    end
end
