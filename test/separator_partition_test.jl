# The recursive separator partition: the seven guarantees of the contract, the
# threshold extremes, immutability, and the optional induced graphs.

function check_partition(g, p, r)
    n = nv(g)
    @test sort(vcat(p.separator, p.parts...)) == collect(1:n)
    @test all(!isempty, p.parts)
    @test all(length(part) <= r for part in p.parts)
    @test length(p.part_of) == n
    for v in p.separator
        @test p.part_of[v] == 0
    end
    for (i, part) in enumerate(p.parts), v in part
        @test p.part_of[v] == i
    end
    for e in edges(g)
        a, b = p.part_of[src(e)], p.part_of[dst(e)]
        @test a == 0 || b == 0 || a == b
    end
    @test length(p.separator) <= 16 * n / sqrt(r)
    @test length(p.parts) <= 2 * n / r + 1
    return nothing
end

@testset "separator partition" begin
    @testset "families and thresholds" begin
        rng = MersenneTwister(77)
        graphs = [
            planar_grid(20, 20),
            planar_grid(7, 60; diagonals = true),
            planar_wheel(200),
            planar_path(150),
            planar_cycle(99),
            random_triangulation(500; rng),
            thin_out(random_triangulation(600; rng), 0.5, rng),
            disjoint_union(planar_grid(9, 9), planar_wheel(30), planar_path(17)),
            PlanarGraph([Int[] for _ = 1:25]),
        ]
        for g in graphs, r in (1, 2, 3, 10, 64, nv(g) - 1, nv(g), nv(g) + 5)
            r >= 1 || continue
            check_partition(g, separator_partition(g, r), r)
        end
    end

    @testset "extremes and errors" begin
        empty_graph = PlanarGraph([Int[] for _ = 1:0])
        p = separator_partition(empty_graph, 3)
        @test isempty(p.separator) && isempty(p.parts) && isempty(p.part_of)

        g = planar_grid(6, 6)
        p = separator_partition(g, nv(g))
        @test isempty(p.separator) && p.parts == [collect(1:nv(g))]
        @test all(==(1), p.part_of)

        single = separator_partition(PlanarGraph([Int[]]), 1)
        @test isempty(single.separator) && single.parts == [[1]]

        @test_throws ArgumentError separator_partition(g, 0)
        crossed = PlanarGraph(g.graph, [collect(neighbors(g.graph, v)) for v = 1:nv(g)])
        @test_throws ArgumentError separator_partition(crossed, 5)
        # Also when no separator would be needed (r >= n), and when the bad
        # component is not the one that is cut.
        @test_throws ArgumentError separator_partition(crossed, nv(crossed))
        k7 = PlanarGraph([[mod1(i + d, 7) for d in (1, 3, 2, 6, 4, 5)] for i = 1:7])
        @test !is_planar_embedding(k7)
        @test_throws ArgumentError separator_partition(k7, 7)
        @test_throws ArgumentError separator_partition(k7, 100)
        @test_throws ArgumentError separator(k7)
        @test_throws ArgumentError separator(disjoint_union(k7, planar_grid(6, 6)))
        @test_throws ArgumentError separator(
            disjoint_union(k7, PlanarGraph([Int[] for _ = 1:7])),
        )
    end

    @testset "determinism and immutability" begin
        g = random_triangulation(300; rng = MersenneTwister(4))
        before = rotation(g)
        a = separator_partition(g, 20)
        b = separator_partition(g, 20)
        @test a.separator == b.separator && a.parts == b.parts && a.part_of == b.part_of
        @test rotation(g) == before
    end

    @testset "induced graphs" begin
        g = planar_grid(15, 15)
        p = separator_partition(g, 30; subgraphs = true)
        check_partition(g, p, 30)
        @test separator_partition(g, 30).graphs === nothing
        for (i, part) in enumerate(p.parts)
            h = p.graphs[i]
            @test nv(h) == length(part)
            for (j, v) in enumerate(part)
                @test p.local_index[v] == j
                kept = [w for w in neighbor_order(g, v) if p.part_of[w] == i]
                @test [part[w] for w in neighbor_order(h, j)] == kept
            end
        end
    end

    @testset "separator size shrinks like 1/sqrt(r)" begin
        g = planar_grid(60, 60)
        sizes = [length(separator_partition(g, r).separator) for r in (4, 16, 64, 256)]
        @test issorted(sizes; rev = true)
        @test sizes[1] <= 16 * nv(g) / 2
    end
end
