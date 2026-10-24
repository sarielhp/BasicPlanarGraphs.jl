# Planarity testing and embedding computation: the Left-Right test.

function graph_edge_keys(g)
    return Set((min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(g))
end

function cube_graph()
    g = SimpleGraph(8)
    for i = 1:4
        add_edge!(g, i, mod1(i + 1, 4))
        add_edge!(g, i + 4, mod1(i + 1, 4) + 4)
        add_edge!(g, i, i + 4)
    end
    return g
end

# Replace edge `index` of `g` by a path of length two.
function subdivide_edge(g::SimpleGraph, index::Integer)
    h = SimpleGraph(nv(g) + 1)
    middle = nv(g) + 1
    for (i, e) in enumerate(collect(edges(g)))
        u, v = src(e), dst(e)
        if i == index
            add_edge!(h, u, middle)
            add_edge!(h, v, middle)
        else
            add_edge!(h, u, v)
        end
    end
    return h
end

# Random graphs, mostly sparse so that both outcomes occur.
function random_graph(n::Integer, p::Float64, rng)
    g = SimpleGraph(n)
    for u = 1:n, v = (u+1):n
        rand(rng) < p && add_edge!(g, u, v)
    end
    return g
end

@testset "planarity" begin
    @testset "planar graphs are accepted" begin
        planar = [
            SimpleGraph(0),
            SimpleGraph(1),
            SimpleGraph(4),
            path_graph(7),
            cycle_graph(9),
            complete_graph(4),
            complete_bipartite_graph(2, 3),
            cube_graph(),
            smallgraph(:dodecahedral),
            smallgraph(:octahedral),
            planar_path(8).graph,
            planar_cycle(12),
            planar_grid(4, 5),
            planar_grid(4, 4; diagonals = true),
            planar_wheel(10),
        ]
        for n in (3, 4, 5, 12, 40)
            push!(planar, random_triangulation(n; rng = MersenneTwister(n)).graph)
        end
        for g in planar
            @test is_planar(g)
            embedded = planar_embedding(g)
            @test embedded isa PlanarGraph
            @test nv(embedded) == nv(g)
            @test graph_edge_keys(embedded) == graph_edge_keys(g)
            @test is_planar_embedding(embedded)
        end
    end

    @testset "non-planar graphs are rejected" begin
        non_planar = [
            complete_graph(5),
            complete_graph(6),
            complete_bipartite_graph(3, 3),
            complete_bipartite_graph(3, 4),
            smallgraph(:petersen),
            subdivide_edge(complete_graph(5), 1),
            subdivide_edge(complete_bipartite_graph(3, 3), 2),
            subdivide_edge(subdivide_edge(complete_graph(5), 3), 5),
        ]
        for g in non_planar
            @test !is_planar(g)
            @test planar_embedding(g) === nothing
        end

        # A dense graph is rejected by the edge-count bound alone.
        dense = complete_graph(30)
        @test !is_planar(dense)
        @test planar_embedding(dense) === nothing
    end

    @testset "the linear decision procedure" begin
        # Same verdict as Graphs' implementation wherever that one survives.
        agree = [
            SimpleGraph(0),
            SimpleGraph(1),
            path_graph(6),
            cycle_graph(7),
            complete_graph(4),
            complete_graph(5),
            complete_bipartite_graph(3, 3),
            smallgraph(:dodecahedral),
            smallgraph(:petersen),
            planar_grid(3, 3),
            planar_wheel(7),
            random_triangulation(30; rng = MersenneTwister(2)).graph,
        ]
        for n in (6, 8, 10), p in (0.2, 0.5, 0.8)
            push!(agree, random_graph(n, p, MersenneTwister(n)))
        end
        for g in agree
            @test is_planar_linear(g) == is_planar(g)
        end

        # A PlanarGraph with a valid rotation is planar without re-running the
        # test; a broken rotation is not evidence of anything, so it is recomputed.
        embedded = planar_wheel(9)
        @test is_planar_linear(embedded)
        broken = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])
        @test !is_planar_embedding(broken)
        @test is_planar_linear(broken)

        # The regression that matters: Graphs' recursive test overflows its stack
        # on a path of this size, because the DFS depth is the path length.
        long_path = path_graph(100_000)
        @test is_planar_linear(long_path)

        # ... and the depth grows with n, so the running time must not.
        short = @elapsed is_planar_linear(path_graph(25_000))
        long = @elapsed is_planar_linear(path_graph(100_000))
        @test long <= 12 * short + 0.5
    end

    @testset "the embedding is a rotation system" begin
        g = smallgraph(:dodecahedral)
        embedded = planar_embedding(g)
        @test all(
            v -> sort(collect(neighbors(embedded, v))) == sort(embedded.rotation[v]),
            1:nv(g),
        )
        @test nfaces(embedded) == ne(embedded) - nv(embedded) + 2
        # The embedding is usable downstream.
        A, B, C = separator(embedded)
        @test length(A) + length(B) + length(C) == nv(g)
        @test isempty(intersect(A, B))

        # Embeddings of disconnected graphs keep the components apart.
        g = SimpleGraph(9)
        for (u, v) in ((1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4), (7, 8))
            add_edge!(g, u, v)
        end
        embedded = planar_embedding(g)
        @test nv(embedded) == 9
        @test graph_edge_keys(embedded) == graph_edge_keys(g)
        @test is_planar_embedding(embedded)
        @test isempty(neighbors(embedded, 9))
    end

    @testset "agreement with Graphs.jl" begin
        rng = MersenneTwister(21)
        for trial = 1:60
            n = rand(rng, 5:11)
            p = rand(rng) * 0.7
            g = random_graph(n, p, rng)
            expected = Graphs.is_planar(g)
            @test is_planar(g) == expected
            embedded = planar_embedding(g)
            @test (embedded === nothing) == !expected
            if embedded !== nothing
                @test is_planar_embedding(embedded)
                @test graph_edge_keys(embedded) == graph_edge_keys(g)
            end
        end
    end

    @testset "planar graphs stay planar when edges are removed" begin
        rng = MersenneTwister(5)
        for trial = 1:40
            n = rand(rng, 4:16)
            g = copy(random_triangulation(n; rng = rng).graph)
            for e in collect(edges(g))
                rand(rng) < 0.5 && rem_edge!(g, src(e), dst(e))
            end
            @test is_planar(g)
            embedded = planar_embedding(g)
            @test embedded !== nothing
            @test is_planar_embedding(embedded)
            @test graph_edge_keys(embedded) == graph_edge_keys(g)
        end
    end

    @testset "the PlanarGraph method" begin
        g = planar_wheel(7)
        @test is_planar(g)
        @test planar_embedding(g) === g

        # A stored rotation that cannot be an embedding: K5's rotation lists
        # every neighbour, so the constructor accepts it.
        k5 = PlanarGraph([
            [2, 3, 4, 5],
            [1, 3, 4, 5],
            [1, 2, 4, 5],
            [1, 2, 3, 5],
            [1, 2, 3, 4],
        ])
        @test !is_planar_embedding(k5)
        @test !is_planar(k5)
        @test planar_embedding(k5) === nothing

        # A broken rotation of a planar graph: two K4s pinned at vertex 1 can be
        # embedded, but not with this rotation. The embedding is recomputed.
        broken = PlanarGraph([
            [2, 3, 4, 5],
            [1, 3, 4],
            [1, 2, 4],
            [1, 2, 3],
            [6, 7, 8, 1],
            [5, 7, 8],
            [5, 6, 8],
            [5, 6, 7],
        ])
        @test is_planar(broken)
        @test !is_planar_embedding(broken)
        recomputed = planar_embedding(broken)
        @test recomputed isa PlanarGraph
        @test recomputed !== broken
        @test is_planar_embedding(recomputed)
        @test graph_edge_keys(recomputed) == graph_edge_keys(broken)
    end

    @testset "linear scaling" begin
        small = random_triangulation(1000; rng = MersenneTwister(1)).graph
        large = random_triangulation(4000; rng = MersenneTwister(1)).graph
        planar_embedding(small)
        t_small = minimum(@elapsed(planar_embedding(small)) for _ = 1:2)
        t_large = minimum(@elapsed(planar_embedding(large)) for _ = 1:2)
        # n grows by 4: allow a factor of 30 before calling it superlinear.
        @test t_large <= 30 * t_small + 0.05
        @test is_planar(large)
    end
end
