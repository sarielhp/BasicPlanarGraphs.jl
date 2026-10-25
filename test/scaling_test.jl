# The algorithm is linear, so quadrupling the number of vertices must not blow
# the running time up. Generous factors keep this stable on loaded machines.

# The fastest of several runs after a warm-up, with a garbage collection before
# each, so that a pause or a busy core does not decide the comparison.
function fastest(f, g; repeats = 5)
    f(g)
    return minimum(begin
        GC.gc()
        @elapsed f(g)
    end for _ = 1:repeats)
end

measure(g) = fastest(separator, g)

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
    # single outer face of length n) and triangulate as well. A 16-fold size would
    # cost 256 times as much if quadratic and about 16 times if linear; allow 80.
    ladder_small = planar_grid(2, 5_000)
    ladder_large = planar_grid(2, 80_000)
    for f in (cycle_separator, triangulate)
        @test fastest(f, ladder_large) <= 80 * fastest(f, ladder_small) + 0.05
    end

    # Many components used to make the drawing and the planarity test quadratic:
    # each component cost a pass over all the vertices.
    isolated_small = SimpleGraph(20_000)
    isolated_large = SimpleGraph(320_000)
    for f in (grid_drawing, planar_embedding)
        @test fastest(f, isolated_large) <= 80 * fastest(f, isolated_small) + 0.05
    end
end
