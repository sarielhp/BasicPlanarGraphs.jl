# r-division: every property is checked with Graphs.jl primitives and plain counting.

# The guarantees of `r_division` that hold for every input.
function check_division(g, d, r)
    n = nv(g)
    copies = zeros(Int, n)
    for R in d.regions
        @test !isempty(R) && allunique(R)
        @test length(R) <= r
        for v in R
            copies[v] += 1
        end
    end
    @test all(>=(1), copies)
    @test d.boundary_vertices == findall(>(1), copies)
    for (i, R) in enumerate(d.regions)
        @test sort(d.boundary[i]) == sort([v for v in R if copies[v] > 1])
        @test length(d.boundary[i])^2 <= 144 * r
    end
    @test all(
        v -> d.region_of[v] == 0 ? copies[v] > 1 : v in d.regions[d.region_of[v]],
        1:n,
    )
    # every edge has both ends in one region
    index = [Set(R) for R in d.regions]
    holder = [Int[] for _ = 1:n]
    for (i, R) in enumerate(d.regions), v in R
        push!(holder[v], i)
    end
    for e in edges(g)
        @test any(i -> dst(e) in index[i], holder[src(e)])
    end
    return nothing
end

@testset "r-division" begin
    @testset "structured graphs" begin
        for (g, r) in (
            (planar_grid(40, 40), 100),
            (planar_grid(60, 60), 200),
            (planar_grid(30, 30; diagonals = true), 150),
            (planar_wheel(3000), 100),
            (planar_path(2000), 100),
            (planar_cycle(2000), 100),
            (triangulate(planar_grid(35, 35)), 120),
        )
            check_division(g, r_division(g, r), r)
        end
    end

    @testset "random and sparse graphs" begin
        rng = MersenneTwister(8)
        for n in (500, 3000), r in (100, 250), seed = 1:2
            g = random_triangulation(n; rng = MersenneTwister(seed))
            check_division(g, r_division(g, r), r)
            h = thin_out(g, 0.5, rng)
            check_division(h, r_division(h, r), r)
        end
    end

    @testset "hubs, components and trees" begin
        star(k) = PlanarGraph(vcat([collect(2:k)], [[1] for _ = 2:k]))
        k2n = PlanarGraph(
            vcat([collect(3:402)], [reverse(collect(3:402))], [[1, 2] for _ = 3:402]),
        )
        many = disjoint_union([planar_wheel(30) for _ = 1:200]...)
        mixed = disjoint_union(
            planar_grid(30, 30),
            PlanarGraph([Int[]]),
            planar_wheel(500),
            planar_path(300),
        )
        for g in (star(2000), k2n, many, mixed)
            check_division(g, r_division(g, 100), 100)
        end
        @test isempty(r_division(many, 100).boundary_vertices)
    end

    @testset "sizes of the result" begin
        # measured, with room; these are test constants, not proved bounds
        for (g, r) in (
            (planar_grid(150, 150), 200),
            (random_triangulation(8000; rng = MersenneTwister(2)), 200),
        )
            d = r_division(g, r)
            n = nv(g)
            @test length(d.regions) <= 4n / r + 4
            @test length(d.boundary_vertices) <= 6n / sqrt(r) + 10
        end
    end

    @testset "small graphs and options" begin
        g = planar_grid(8, 8)
        d = r_division(g, 100)
        @test d.regions == [collect(1:64)] && isempty(d.boundary_vertices)
        @test d.region_of == fill(1, 64)
        @test r_division(g, 100).regions == d.regions
        e = r_division(PlanarGraph(Vector{Int}[]), 100)
        @test isempty(e.regions) && isempty(e.region_of)
        @test_throws ArgumentError r_division(g, 99)
        @test_throws ArgumentError r_division(g, 0)
        graph = planar_grid(4, 4).graph
        crossed = PlanarGraph(graph, [collect(neighbors(graph, v)) for v = 1:16])
        @test_throws ArgumentError r_division(crossed, 100)
        g32 = PlanarGraph(
            SimpleGraph{Int32}(planar_grid(30, 30).graph),
            rotation(planar_grid(30, 30)),
        )
        check_division(g32, r_division(g32, 100), 100)
    end

    @testset "pieces with too many shared vertices are cut" begin
        # Natural inputs seldom reach the second phase, so start from a piece whose
        # vertices are all shared.
        g = planar_grid(15, 15)
        dg = BasicPlanarGraphs.DartGraph(g)
        start = BasicPlanarGraphs._Piece(dg, collect(1:225), trues(225))
        leaves = BasicPlanarGraphs._divide(start, 225)
        @test length(leaves) > 1
        covered = falses(225)
        for leaf in leaves
            @test leaf.graph.n <= 225
            @test count(leaf.shared)^2 <= 144 * 225
            covered[leaf.ids] .= true
        end
        @test all(covered)
    end

    @testset "running time" begin
        small = random_triangulation(4_000; rng = MersenneTwister(1))
        large = random_triangulation(64_000; rng = MersenneTwister(1))
        @test fastest(g -> r_division(g, 200), large; repeats = 3) <=
              80 * fastest(g -> r_division(g, 200), small; repeats = 3) + 0.2
    end
end
