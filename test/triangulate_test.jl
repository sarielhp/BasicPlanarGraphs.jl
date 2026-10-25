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
        # A path's single face walk repeats every internal vertex; each repeat
        # is cut at its second appearance and the rest becomes a simple cycle.
        path = triangulate(planar_path(6))
        @test ne(path) == 3 * 6 - 6
        @test all_faces_triangular(path)
        @test is_planar_embedding(path)

        for n in (4, 5, 7, 12)
            p = triangulate(planar_path(n))
            @test ne(p) == 3n - 6
            @test all_faces_triangular(p)
            @test is_planar_embedding(p)
        end

        # Two triangles sharing a cut vertex: the shared vertex repeats on the
        # only non-triangular face, so it is peeled before that face is fanned.
        bowtie = PlanarGraph([[2, 3, 4, 5], [1, 3], [1, 2], [1, 5], [1, 4]])
        t = triangulate(bowtie)
        @test ne(t) == 3 * 5 - 6
        @test all_faces_triangular(t)
        @test is_planar_embedding(t)

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

    @testset "sparse graphs with large faces" begin
        # Properties of the result on connected graphs with bridges, cut
        # vertices and big faces: a simple maximal planar supergraph that keeps
        # the original edges and their rotation order.
        rng = MersenneTwister(5)
        tested = 0
        for trial = 1:300
            g = thin_out(random_triangulation(rand(rng, 3:120); rng), rand(rng) * 0.95, rng)
            is_connected(g.graph) || continue
            t = triangulate(g)
            n = nv(g)
            @test ne(t) == 3n - 6
            @test all_faces_triangular(t)
            @test is_planar_embedding(t)
            rg, rt = rotation(g), rotation(t)
            for v = 1:n
                kept = [w for w in rt[v] if has_edge(g, v, w)]
                @test BasicPlanarGraphs.same_cyclic_order(kept, rg[v])
            end
            tested += 1
        end
        @test tested > 50
    end

    @testset "faces whose vertices are adjacent elsewhere" begin
        # K_{2,n} has hubs that share all of its faces; a ladder has a rung
        # between every pair of opposite vertices of its outer face. A chord
        # between such vertices would be a parallel edge.
        for g in (planar_grid(2, 40), planar_grid(3, 25), planar_wheel(40))
            t = triangulate(g)
            @test ne(t) == 3nv(g) - 6
            @test all_faces_triangular(t)
        end
        rot = vcat([collect(3:12)], [reverse(collect(3:12))], [[1, 2] for _ = 3:12])
        k2n = PlanarGraph(rot)
        @test is_planar_embedding(k2n)
        t = triangulate(k2n)
        @test ne(t) == 3nv(k2n) - 6
        @test all_faces_triangular(t)
    end

    @testset "chords are spread over the vertices of a face" begin
        # A ladder cuts the face into quadrangles, so no vertex collects the
        # chords of the whole face (a fan from one vertex would give a degree of
        # n - 1): a vertex gains a constant number of edges per face, also when
        # the rungs are blocked by edges outside (the ladder graph, a path).
        for g in (planar_cycle(4000), planar_path(4000), planar_grid(2, 2000))
            t = triangulate(g)
            @test maximum(degree(t, v) - degree(g, v) for v = 1:nv(g)) <= 6
        end
    end
end
