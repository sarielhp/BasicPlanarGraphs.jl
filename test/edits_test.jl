# Embedding-preserving edits and the Graphs.jl interface of PlanarGraph.

# The graph and the rotation describe the same edges (the constructor checks it) and
# the embedding is planar.
function consistent(g::PlanarGraph)
    PlanarGraph(g.graph, rotation(g))
    return is_planar_embedding(g)
end

# Edges of the contraction of `u v` into `u`, computed on plain edge lists.
function contracted_edges(g, u, v)
    out = Set{Tuple{Int,Int}}()
    for e in edges(g)
        a, b = src(e) == v ? u : src(e), dst(e) == v ? u : dst(e)
        a == b || push!(out, (min(a, b), max(a, b)))
    end
    return out
end

edge_set(g) = Set((min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(g))

@testset "embedding edits" begin
    @testset "delete_edge!" begin
        rng = MersenneTwister(4)
        for g in (random_triangulation(60; rng), planar_grid(6, 6), planar_wheel(12))
            es = shuffle(rng, collect(edges(g)))
            m = ne(g)
            for e in es[1:(end÷2)]
                @test delete_edge!(g, src(e), dst(e)) === g
                m -= 1
                @test ne(g) == m
                @test !has_edge(g, src(e), dst(e))
                @test consistent(g)
            end
        end
        g = planar_cycle(5)
        @test_throws ArgumentError delete_edge!(g, 1, 3)
        @test_throws ArgumentError delete_edge!(g, 1, 99)
        @test rem_edge!(g, 1, 3) == false
        @test rem_edge!(g, 1, 2) == true
        @test !has_edge(g, 1, 2) && consistent(g)
        @test ne(g) == 4
    end

    @testset "contract_edge!" begin
        rng = MersenneTwister(5)
        for g in (
            random_triangulation(80; rng),
            planar_grid(7, 7),
            planar_wheel(15),
            thin_out(random_triangulation(80; rng), 0.4, rng),
        )
            for _ = 1:30
                es = collect(edges(g))
                isempty(es) && break
                e = rand(rng, es)
                u, v = src(e), dst(e)
                expected = contracted_edges(g, u, v)
                @test contract_edge!(g, u, v) === g
                @test edge_set(g) == expected
                @test degree(g, v) == 0
                @test consistent(g)
                active = count(w -> degree(g, w) > 0, vertices(g))
                active >= 3 && @test ne(g) <= 3 * active - 6
            end
        end
        g = planar_wheel(6)
        @test_throws ArgumentError contract_edge!(g, 2, 4)
        # a triangle of the triangulation disappears with its two parallel edges
        t = random_triangulation(30; rng = MersenneTwister(9))
        e = first(edges(t))
        common = length(intersect(neighbors(t, src(e)), neighbors(t, dst(e))))
        before = ne(t)
        contract_edge!(t, src(e), dst(e))
        @test ne(t) == before - 1 - common
    end

    @testset "add_edge_in_face!" begin
        g = planar_cycle(6)
        add_edge_in_face!(g, 1, 4; after_u = 2, after_v = 5)
        @test consistent(g) && nfaces(g) == 3 && has_edge(g, 1, 4)
        # a path: both ends have one corner
        p = planar_path(4)
        add_edge_in_face!(p, 1, 4)
        @test consistent(p) && nfaces(p) == 2
        # isolated endpoints go into any face
        s = PlanarGraph([[2], [1], Int[]])
        add_edge_in_face!(s, 3, 2)
        @test consistent(s) && ne(s) == 2
        z = PlanarGraph([Int[], Int[]])
        add_edge_in_face!(z, 1, 2)
        @test consistent(z) && ne(z) == 1
        # adding then deleting restores the embedding
        h = planar_grid(4, 4)
        copy_h = PlanarGraph(h.graph, rotation(h))
        add_edge_in_face!(h, 1, 6; after_u = 5, after_v = 2)
        delete_edge!(h, 1, 6)
        @test h == copy_h
    end

    @testset "add_edge_in_face! accepts exactly the corners on one face" begin
        for g0 in (planar_wheel(5), planar_grid(3, 3), planar_cycle(5))
            corners = [(v, a) for v = 1:nv(g0) for a in neighbor_order(g0, v)]
            for (u, a) in corners, (v, b) in corners
                u == v && continue
                has_edge(g0, u, v) && continue
                rot = rotation(g0)
                insert!(rot[u], findfirst(==(a), rot[u]) + 1, v)
                insert!(rot[v], findfirst(==(b), rot[v]) + 1, u)
                planar = is_planar_embedding(PlanarGraph(rot))
                g = PlanarGraph(g0.graph, rotation(g0))
                done = try
                    add_edge_in_face!(g, u, v; after_u = a, after_v = b)
                    true
                catch err
                    err isa ArgumentError || rethrow()
                    false
                end
                @test done == planar
                done && @test g == PlanarGraph(rot)
            end
        end
    end

    @testset "add_edge_in_face! validation" begin
        g = planar_grid(3, 3)
        @test_throws ArgumentError add_edge_in_face!(g, 1, 1)
        @test_throws ArgumentError add_edge_in_face!(g, 1, 2; after_u = 4, after_v = 1)
        @test_throws ArgumentError add_edge_in_face!(g, 1, 99)
        @test_throws ArgumentError add_edge_in_face!(g, 1, 5)            # corners needed
        @test_throws ArgumentError add_edge_in_face!(g, 1, 5; after_u = 7, after_v = 2)
        @test_throws ArgumentError add_edge_in_face!(
            PlanarGraph([Int[], [3], [2]]),
            1,
            2;
            after_u = 3,
        )
    end
end

@testset "PlanarGraph is a Graphs.jl graph" begin
    g = thin_out(
        random_triangulation(120; rng = MersenneTwister(2)),
        0.3,
        MersenneTwister(2),
    )
    s = g.graph
    @test g isa AbstractGraph{Int}
    @test !is_directed(g) && !is_directed(typeof(g))
    @test edgetype(g) == edgetype(s)
    @test has_vertex(g, 5) && !has_vertex(g, nv(g) + 1)
    @test outneighbors(g, 3) == neighbors(s, 3) == inneighbors(g, 3)
    @test connected_components(g) == connected_components(s)
    @test bfs_parents(g, 1) == bfs_parents(s, 1)
    @test dijkstra_shortest_paths(g, 1).dists == dijkstra_shortest_paths(s, 1).dists
    @test adjacency_matrix(g) == adjacency_matrix(s)
    @test is_bipartite(g) == is_bipartite(s)
    @test core_number(g) == core_number(s)
    @test degree_histogram(g) == degree_histogram(s)
    @test is_connected(planar_grid(4, 4)) && diameter(planar_grid(4, 4)) == 6
    @test global_clustering_coefficient(g) == global_clustering_coefficient(s)
    g32 = PlanarGraph([Int32[2, 3], Int32[3, 1], Int32[1, 2]])
    @test g32 isa AbstractGraph{Int32} && connected_components(g32) == [Int32[1, 2, 3]]
    # the package's own functions still take the graph
    @test is_planar(g) && five_coloring(g) isa Vector{Int}
    @test_throws ArgumentError five_coloring(path_digraph(3))
end
