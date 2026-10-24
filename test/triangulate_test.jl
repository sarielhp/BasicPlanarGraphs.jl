# Triangulation: chords inside faces, maximal planar counts, Euler.

# Every face walk is a triangle.
function all_faces_triangular(g::PlanarGraph)
    return all(length(w) == 3 for w in faces(g))
end

@testset "triangulate" begin
    @testset "planar graph inputs" begin
        for g in (
            planar_cycle(5),
            planar_cycle(9),
            planar_wheel(5),
            planar_wheel(12),
            planar_grid(3, 4),
            planar_grid(6, 6),
            planar_grid(4, 7; diagonals = true),
        )
            t = triangulate(g)
            n = nv(g)
            @test t isa PlanarGraph
            @test nv(t) == n
            @test ne(t) == 3n - 6
            @test nfaces(t) == 2n - 4
            @test all_faces_triangular(t)
            @test is_planar_embedding(t)
            @test triangulate(t) == t
        end
    end

    @testset "graphs with bridges" begin
        # A path can still be triangulated: its faces are simple cycles.
        path = triangulate(planar_path(6))
        @test ne(path) == 3 * 6 - 6
        @test all_faces_triangular(path)
        @test is_planar_embedding(path)

        # Bridges are fine: the faces are still simple cycles.
        two_k4 = [
            [2, 3, 4, 5],
            [1, 4, 3],
            [1, 2, 4],
            [1, 3, 2],
            [6, 7, 8, 1],
            [5, 8, 7],
            [5, 6, 8],
            [5, 7, 6],
        ]
        bridged = PlanarGraph(two_k4)
        @test is_planar_embedding(bridged)
        @test ne(triangulate(bridged)) == 3 * 8 - 6
    end

    @testset "rotation systems that are not planar" begin
        # Same graph, different rotations: this one is K4 plus K4 pinned at a
        # vertex, which no embedding of that edge set can realise. Chord
        # insertion must refuse instead of returning a broken triangulation.
        rot = [
            [2, 3, 4, 5],
            [1, 3, 4],
            [1, 2, 4],
            [1, 2, 3],
            [6, 7, 8, 1],
            [5, 7, 8],
            [5, 6, 8],
            [5, 6, 7],
        ]
        broken = PlanarGraph(rot)
        @test !is_planar_embedding(broken)
        @test_throws ArgumentError triangulate(broken)
    end

    @testset "random triangulations" begin
        for n in (3, 4, 5, 10, 40), seed = 1:2
            g = random_triangulation(n; rng = MersenneTwister(seed))
            @test ne(g) == 3n - 6
            @test nfaces(g) == 2n - 4
            @test all_faces_triangular(g)
            @test is_planar_embedding(g)
            @test triangulate(g) == g
        end
    end
end
