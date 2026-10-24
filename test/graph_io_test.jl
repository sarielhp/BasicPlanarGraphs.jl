# Plain-text graph files: the format, the round trip, and what is rejected.

function edge_keys(g)
    return Set((min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(g))
end

@testset "graph files" begin
    @testset "round trip" begin
        graphs = (
            PlanarGraph([[2, 3], [3, 1], [1, 2]]),
            planar_cycle(6),
            planar_wheel(6),
            planar_grid(3, 4),
            planar_grid(3, 3; diagonals = true),
            random_triangulation(25; rng = MersenneTwister(4)),
            planar_path(5),
        )
        for g in graphs
            io = IOBuffer()
            write_planar_graph(io, g)
            text = String(take!(io))

            # A header, both sections, and the terminator.
            @test startswith(text, "# BasicPlanarGraphs 1\n")
            @test occursin("\nvertices $(nv(g))\n", text)
            @test occursin("\nedges\n", text)
            @test occursin("\nrotation\n", text)
            @test endswith(text, "end\n")

            back = read_planar_graph(IOBuffer(text))
            @test back isa PlanarGraph
            @test back == g
            @test rotation(back) == rotation(g)
            @test nv(back) == nv(g)
            @test edge_keys(back) == edge_keys(g)

            bare = read_graph(IOBuffer(text))
            @test bare isa SimpleGraph
            @test nv(bare) == nv(g)
            @test edge_keys(bare) == edge_keys(g)
        end
    end

    @testset "the format" begin
        text = """
        # written by hand
        vertices 4

        edges
        1 2
        2 3   # a trailing comment
        3 4
        4 1
        rotation
        1: 2 4
        2: 3 1
        3: 4 2
        4: 1 3
        end
        """
        g = read_planar_graph(IOBuffer(text))
        @test nv(g) == 4
        @test ne(g) == 4
        @test rotation(g) == [[2, 4], [3, 1], [4, 2], [1, 3]]
        @test is_planar_embedding(g)
        @test nfaces(g) == 2

        bare = read_graph(IOBuffer(text))
        @test nv(bare) == 4
        @test ne(bare) == 4
        @test has_edge(bare, 1, 2) && has_edge(bare, 3, 4)

        # Without a `vertices` line the count is the largest label.
        g = read_graph(IOBuffer("edges\n1 3\n"))
        @test nv(g) == 3
        @test ne(g) == 1

        # A rotation-only file describes the graph through its rotations.
        triangle = "vertices 3\nrotation\n1: 2 3\n2: 3 1\n3: 1 2\nend\n"
        derived = read_graph(IOBuffer(triangle))
        @test nv(derived) == 3
        @test ne(derived) == 3
        @test edge_keys(derived) == edge_keys(read_planar_graph(IOBuffer(triangle)))

        # Empty graphs and isolated vertices survive the round trip.
        for g in (SimpleGraph(0), SimpleGraph(1), SimpleGraph(4))
            io = IOBuffer()
            write_graph(io, g)
            back = read_graph(IOBuffer(String(take!(io))))
            @test nv(back) == nv(g)
            @test ne(back) == 0
        end

        # A graph with an isolated vertex: the rotation of that vertex is empty.
        g = PlanarGraph([[2], [1], Int[]])
        io = IOBuffer()
        write_planar_graph(io, g)
        back = read_planar_graph(IOBuffer(String(take!(io))))
        @test back == g
        @test rotation(back)[3] == Int[]
    end

    @testset "written files" begin
        mktempdir() do dir
            path = joinpath(dir, "graph.txt")
            g = random_triangulation(12; rng = MersenneTwister(11))
            write_planar_graph(path, g)
            @test isfile(path)
            @test read_planar_graph(path) == g
            @test edge_keys(read_graph(path)) == edge_keys(g)

            simple_path = joinpath(dir, "simple.txt")
            s = SimpleGraph(5)
            for (u, v) in ((1, 2), (2, 3), (4, 5))
                add_edge!(s, u, v)
            end
            write_graph(simple_path, s)
            @test nv(read_graph(simple_path)) == 5
            @test edge_keys(read_graph(simple_path)) == edge_keys(s)
            @test_throws ArgumentError read_planar_graph(simple_path)
        end
    end

    @testset "rejected input" begin
        # No rotation section at all.
        @test_throws ArgumentError read_planar_graph(IOBuffer("vertices 3\nedges\n1 2\n"))
        # Data outside a section.
        @test_throws ArgumentError read_graph(IOBuffer("1 2\n"))
        # Not a pair of integers.
        @test_throws ArgumentError read_graph(IOBuffer("edges\n1\n"))
        @test_throws ArgumentError read_graph(IOBuffer("edges\n1 x\n"))
        # Self-loops and label 0.
        @test_throws ArgumentError read_graph(IOBuffer("edges\n1 1\n"))
        @test_throws ArgumentError read_graph(IOBuffer("edges\n0 1\n"))
        # The same edge twice.
        @test_throws ArgumentError read_graph(IOBuffer("edges\n1 2\n2 1\n"))
        # Vertex count is not an integer, or negative.
        @test_throws ArgumentError read_graph(IOBuffer("vertices many\n"))
        @test_throws ArgumentError read_graph(IOBuffer("vertices -1\n"))
        # Rotation problems.
        @test_throws ArgumentError read_planar_graph(
            IOBuffer("vertices 3\nrotation\n1: 2\n2: 1\n3:\n3: 2\nend\n"),
        )
        @test_throws ArgumentError read_planar_graph(
            IOBuffer("vertices 3\nrotation\n1: 1 2\n2: 3\n3: 2\nend\n"),
        )
        @test_throws ArgumentError read_planar_graph(
            IOBuffer("vertices 3\nrotation\n1: 2\n2: 3\n3: 2\nend\n"),
        )
        @test_throws ArgumentError read_planar_graph(
            IOBuffer("vertices 2\nrotation\n1: 2 3\n2: 1\nend\n"),
        )
        # The two sections disagree.
        @test_throws ArgumentError read_planar_graph(
            IOBuffer("vertices 3\nedges\n1 2\nrotation\n1: 2 3\n2: 1\n3: 1\nend\n"),
        )
        # A rotation that is not a planar embedding is still a valid file: the
        # format stores the embedding as it is, it does not check planarity.
        k5 = PlanarGraph([
            [2, 3, 4, 5],
            [1, 3, 4, 5],
            [1, 2, 4, 5],
            [1, 2, 3, 5],
            [1, 2, 3, 4],
        ])
        io = IOBuffer()
        write_planar_graph(io, k5)
        @test read_planar_graph(IOBuffer(String(take!(io)))) == k5
        # Missing files.
        @test_throws ArgumentError read_graph(joinpath(tempdir(), "no-such-graph.txt"))
    end
end
