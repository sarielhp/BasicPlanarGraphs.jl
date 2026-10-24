# The algorithm is linear, so quadrupling the number of vertices must not blow
# the running time up. Generous factors keep this stable on loaded machines.

function measure(g; repeats = 2)
    separator(g)                     # warm up compilation and caches
    return minimum(@elapsed(separator(g)) for _ = 1:repeats)
end

@testset "linear scaling" begin
    small = planar_grid(60, 60; diagonals = true)
    large = planar_grid(120, 120; diagonals = true)

    t_small = measure(small)
    t_large = measure(large)

    # n grows by 4: allow a factor of 30 before calling it superlinear.
    @test t_large <= 30 * t_small + 0.05

    _, _, c_small = separator(small)
    _, _, c_large = separator(large)
    @test length(c_large) <= 4 * length(c_small)
    @test length(c_large) <= 2 * sqrt(2 * nv(large)) + 1

    tri_small = random_triangulation(2000; rng = MersenneTwister(11))
    tri_large = random_triangulation(8000; rng = MersenneTwister(11))
    t_small = measure(tri_small)
    t_large = measure(tri_large)
    @test t_large <= 30 * t_small + 0.05

    # One huge face used to make the cycle separator quadratic (a ladder has a
    # single outer face of length n): an 8-fold size is allowed 25 times the time.
    function measure_cycle(g)
        cycle_separator(g)
        return minimum(@elapsed(cycle_separator(g)) for _ = 1:2)
    end
    ladder_small = planar_grid(2, 5_000)
    ladder_large = planar_grid(2, 40_000)
    @test measure_cycle(ladder_large) <= 25 * measure_cycle(ladder_small) + 0.05

    # The same single huge face made triangulate quadratic as well.
    function measure_triangulate(g)
        triangulate(g)
        return minimum(@elapsed(triangulate(g)) for _ = 1:2)
    end
    @test measure_triangulate(planar_grid(2, 40_000)) <=
          25 * measure_triangulate(planar_grid(2, 5_000)) + 0.05
end
