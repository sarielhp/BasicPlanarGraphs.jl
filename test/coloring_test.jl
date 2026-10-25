# Five-coloring: every result is checked edge by edge with Graphs.jl primitives.

function check_coloring(g, c)
    @test c isa Vector{Int}
    @test length(c) == nv(g)
    @test all(x -> 1 <= x <= 5, c)
    @test all(e -> c[src(e)] != c[dst(e)], edges(g))
    @test is_proper_coloring(g, c)
    return c
end

# The icosahedron: planar and 5-regular, so every first step identifies two vertices.
function icosahedron()
    g = SimpleGraph(12)
    for i = 2:6
        add_edge!(g, 1, i)
        add_edge!(g, 12, i + 5)
        add_edge!(g, i, i == 6 ? 2 : i + 1)
        add_edge!(g, i + 5, i == 6 ? 7 : i + 6)
        add_edge!(g, i, i + 5)
        add_edge!(g, i, i == 2 ? 11 : i + 4)
    end
    return g
end

@testset "five coloring" begin
    @testset "structured graphs" begin
        ico = icosahedron()
        @test ne(ico) == 30 && all(v -> degree(ico, v) == 5, vertices(ico))
        @test is_planar_linear(ico)
        check_coloring(ico, five_coloring(ico))
        for g in (
            planar_cycle(3),
            planar_cycle(9),
            planar_wheel(5),
            planar_wheel(8),
            planar_grid(7, 9),
            planar_grid(12, 12; diagonals = true),
            planar_path(30),
            triangulate(planar_grid(10, 10)),
        )
            check_coloring(g, five_coloring(g))
            check_coloring(g.graph, five_coloring(g.graph))
        end
        @test length(unique(five_coloring(planar_grid(5, 5)))) <= 5
    end

    @testset "random triangulations" begin
        for n in (4, 5, 12, 50, 300, 2000), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            check_coloring(g, five_coloring(g))
        end
    end

    @testset "sparse and disconnected graphs" begin
        rng = MersenneTwister(11)
        for trial = 1:200
            g = thin_out(random_triangulation(rand(rng, 3:150); rng), rand(rng), rng)
            check_coloring(g, five_coloring(g))
        end
        g = disjoint_union(planar_wheel(7), PlanarGraph([Int[]]), planar_cycle(5))
        check_coloring(g, five_coloring(g))
    end

    @testset "tiny graphs" begin
        @test five_coloring(SimpleGraph(0)) == Int[]
        @test five_coloring(SimpleGraph(1)) == [1]
        check_coloring(SimpleGraph(2), five_coloring(SimpleGraph(2)))
        k2 = path_graph(2)
        check_coloring(k2, five_coloring(k2))
        k4 = complete_graph(4)
        @test sort(five_coloring(k4)) == [1, 2, 3, 4]
    end

    @testset "other vertex types" begin
        g = SimpleGraph{Int32}(planar_grid(6, 6).graph)
        check_coloring(g, five_coloring(g))
    end

    @testset "non-planar and unsupported input" begin
        for g in (complete_graph(6), complete_graph(8), complete_graph(7))
            @test_throws ArgumentError five_coloring(g)
        end
        # not planar, but the method never gets stuck: the coloring is still proper
        for g in (complete_graph(5), complete_bipartite_graph(3, 3), smallgraph(:petersen))
            check_coloring(g, five_coloring(g))
        end
        @test_throws ArgumentError five_coloring(path_digraph(3))
    end

    @testset "is_proper_coloring" begin
        g = planar_cycle(4)
        @test is_proper_coloring(g, [1, 2, 1, 2])
        @test !is_proper_coloring(g, [1, 2, 2, 1])
        @test !is_proper_coloring(g, [1, 2, 1])
        @test is_proper_coloring(g, [10, 20, 10, 20])
    end

    @testset "running time" begin
        small = random_triangulation(5_000; rng = MersenneTwister(1)).graph
        large = random_triangulation(40_000; rng = MersenneTwister(1)).graph
        @test fastest(five_coloring, large) <= 40 * fastest(five_coloring, small) + 0.05
        many = disjoint_union([planar_wheel(6) for _ = 1:20_000]...).graph
        check_coloring(many, five_coloring(many))
        path = planar_path(500_000).graph
        check_coloring(path, five_coloring(path))
    end
end
