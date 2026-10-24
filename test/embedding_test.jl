# The public type, the rotation system, faces and the dual.

@testset "PlanarGraph" begin
    triangle = PlanarGraph([[2, 3], [3, 1], [1, 2]])
    @test nv(triangle) == 3
    @test ne(triangle) == 3
    @test vertices(triangle) == 1:3
    @test length(collect(edges(triangle))) == 3
    @test degree(triangle, 1) == 2
    @test neighbors(triangle, 1) == [2, 3]
    @test has_edge(triangle, 1, 2)
    @test !has_edge(triangle, 1, 1)

    @test neighbor_order(triangle, 1) == [2, 3]
    @test rotation(triangle) == [[2, 3], [3, 1], [1, 2]]
    @test rotation(triangle) !== triangle.rotation
    @test planar_graph([[2, 3], [3, 1], [1, 2]]) == triangle

    @testset "faces" begin
        walks = faces(triangle)
        @test length(walks) == 2
        @test all(length(w) == 3 for w in walks)
        @test nfaces(triangle) == 2
        @test is_planar_embedding(triangle)

        square = planar_cycle(4)
        @test nfaces(square) == 2
        @test nfaces(planar_path(7)) == 1
        @test nfaces(planar_wheel(4)) == 5
        @test nfaces(planar_grid(3, 4)) == 7
    end

    @testset "validation" begin
        @test_throws ArgumentError PlanarGraph([[2, 3], [3, 1]])
        @test_throws ArgumentError PlanarGraph([[2, 2], [1, 1]])
        @test_throws ArgumentError PlanarGraph([[2, 9], [1, 1]])
        @test_throws ArgumentError PlanarGraph([[2, 3], [1, 1], [1, 2]])
        g = SimpleGraph(3)
        Graphs.add_edge!(g, 1, 2)
        @test PlanarGraph(g, [[2], [1], Int[]]) isa PlanarGraph
        @test_throws ArgumentError PlanarGraph(g, [[2, 3], [1], [1]])
        @test_throws ArgumentError PlanarGraph(g, [[2], [1]])
        @test_throws ArgumentError PlanarGraph([[2], [1.5]])
    end

    @testset "element type" begin
        g = PlanarGraph([Int32[2, 3], Int32[3, 1], Int32[1, 2]])
        @test g isa PlanarGraph{Int32}
        @test neighbor_order(g, 2) == Int32[3, 1]
        @test nfaces(g) == 2
    end
end

@testset "is_planar_embedding" begin
    @test is_planar_embedding(planar_path(1))
    @test is_planar_embedding(planar_path(9))
    @test is_planar_embedding(planar_cycle(7))
    @test is_planar_embedding(planar_grid(5, 5))
    @test is_planar_embedding(planar_grid(5, 5; diagonals = true))
    @test is_planar_embedding(planar_wheel(8))

    # K5 cannot be embedded: whatever the rotation, Euler's formula fails.
    k5 = PlanarGraph([[2, 3, 4, 5], [1, 3, 4, 5], [1, 2, 4, 5], [1, 2, 3, 5], [1, 2, 3, 4]])
    @test !is_planar_embedding(k5)
end

@testset "dual" begin
    for g in (
        triangulate(planar_grid(3, 4)),
        triangulate(planar_wheel(6)),
        triangulate(planar_cycle(7)),
        random_triangulation(30; rng = MersenneTwister(3)),
    )
        d = dual(g)
        @test nv(d) == nfaces(g)
        @test ne(d) == ne(g)
        @test nfaces(d) == nv(g)
        @test is_planar_embedding(d)
        @test all(degree(d, v) == 3 for v = 1:nv(d))

        # The double dual is the primal graph again, up to relabelling by face.
        back = dual(d)
        @test nv(back) == nv(g)
        @test ne(back) == ne(g)
        @test nfaces(back) == nfaces(g)
        @test sort(collect(degree(back, v) for v = 1:nv(back))) ==
              sort(collect(degree(g, v) for v = 1:nv(g)))
        @test is_planar_embedding(back)
    end

    # A bridge would become a self-loop in the dual, which has no PlanarGraph.
    @test_throws ArgumentError dual(planar_path(5))
end

@testset "internal dart layer" begin
    dg = BasicPlanarGraphs.DartGraph(planar_grid(2, 3))
    @test dg.n == 6
    @test BasicPlanarGraphs.ndarts(dg) == 2 * ne(planar_grid(2, 3))
    @test BasicPlanarGraphs.nedges(dg) == ne(planar_grid(2, 3))
    @test all(p -> dg.twin[dg.twin[p]] == p, 1:BasicPlanarGraphs.ndarts(dg))
    @test BasicPlanarGraphs.rotation_lists(dg) == rotation(planar_grid(2, 3))
    @test BasicPlanarGraphs.PlanarGraph(dg) == planar_grid(2, 3)
    @test copy(dg).tail == dg.tail
    @test copy(dg).head == dg.head

    keep = BitVector([true, true, true, false, true, true])
    sub = BasicPlanarGraphs.restrict(dg, keep)
    @test is_planar_embedding(BasicPlanarGraphs.PlanarGraph(sub))
    @test BasicPlanarGraphs.darts_of(sub, 4) == Int[]

    small, ids = BasicPlanarGraphs.induced(dg, keep)
    @test ids == [1, 2, 3, 5, 6]
    @test small.n == 5
    @test BasicPlanarGraphs.nedges(small) == BasicPlanarGraphs.nedges(sub)
    @test is_planar_embedding(BasicPlanarGraphs.PlanarGraph(small))
    # `induced` renumbers, so vertex 4 of `small` is vertex 5 of `dg`.
    @test sort(collect(ids[small.head[d]] for d in BasicPlanarGraphs.darts_of(small, 4))) ==
          sort(
        collect(dg.head[d] for d in BasicPlanarGraphs.darts_of(dg, 5) if keep[dg.head[d]]),
    )

    label, ncomps = BasicPlanarGraphs.component_labels(dg)
    @test ncomps == 1
    @test all(==(1), label)

    parent, _, depth, order = BasicPlanarGraphs.bfs_tree(dg, [1])
    @test length(order) == 6
    @test depth[1] == 0
    @test maximum(depth) == 3
    @test parent[1] == 0
    @test all(depth[parent[v]] == depth[v] - 1 for v = 2:6 if parent[v] != 0)
end
