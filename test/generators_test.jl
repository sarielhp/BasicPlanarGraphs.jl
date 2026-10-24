# The constructors used by the examples and by the separator tests.

@testset "generators" begin
    @testset "planar_path" begin
        for n = 1:6
            g = planar_path(n)
            @test nv(g) == n
            @test ne(g) == n - 1
            @test nfaces(g) == (n == 1 ? 0 : 1)
            @test is_planar_embedding(g)
            @test degree(g, 1) == (n == 1 ? 0 : 1)
        end
        @test_throws ArgumentError planar_path(0)
    end

    @testset "planar_cycle" begin
        for n = 3:9
            g = planar_cycle(n)
            @test nv(g) == n
            @test ne(g) == n
            @test all(degree(g, v) == 2 for v = 1:n)
            @test nfaces(g) == 2
            @test is_planar_embedding(g)
        end
        @test_throws ArgumentError planar_cycle(2)
    end

    @testset "planar_grid" begin
        for (r, c) in ((1, 1), (1, 5), (2, 2), (3, 4), (5, 5))
            g = planar_grid(r, c)
            n = r * c
            @test nv(g) == n
            @test ne(g) == 2r * c - r - c
            @test nfaces(g) == (ne(g) == 0 ? 0 : ne(g) - n + 2)
            @test is_planar_embedding(g)
        end
        for (r, c) in ((2, 2), (3, 4), (5, 5))
            g = planar_grid(r, c; diagonals = true)
            n = r * c
            @test ne(g) == 2r * c - r - c + (r - 1) * (c - 1)
            @test is_planar_embedding(g)
        end
        @test_throws ArgumentError planar_grid(0, 3)
        @test_throws ArgumentError planar_grid(3, 0)
    end

    @testset "planar_wheel" begin
        for rim = 3:7
            g = planar_wheel(rim)
            @test nv(g) == rim + 1
            @test ne(g) == 2rim
            @test degree(g, 1) == rim
            @test nfaces(g) == rim + 1
            @test is_planar_embedding(g)
        end
        @test_throws ArgumentError planar_wheel(2)
    end

    @testset "random_triangulation" begin
        for n in (3, 4, 8, 25), seed = 1:3
            g = random_triangulation(n; rng = MersenneTwister(seed))
            @test nv(g) == n
            @test ne(g) == 3n - 6
            @test is_planar_embedding(g)
            same = random_triangulation(n; rng = MersenneTwister(seed))
            @test rotation(same) == rotation(g)
        end
        @test_throws ArgumentError random_triangulation(2)
    end
end
