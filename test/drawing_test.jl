# Straight-line grid drawings: no crossings, integer coordinates, linear area.

_orientation(a, b, c) = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

# `p` lies on the closed segment `a`-`b`.
function on_segment(p, a, b)
    _orientation(a, b, p) == 0 || return false
    return min(a[1], b[1]) <= p[1] <= max(a[1], b[1]) &&
           min(a[2], b[2]) <= p[2] <= max(a[2], b[2])
end

# The two segments cross properly, in their interiors.
function segments_cross(p1, p2, p3, p4)
    d1 = _orientation(p3, p4, p1)
    d2 = _orientation(p3, p4, p2)
    d3 = _orientation(p1, p2, p3)
    d4 = _orientation(p1, p2, p4)
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
           ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
end

# Everything a straight-line drawing has to satisfy, checked edge pair by edge
# pair: distinct integer positions, no vertex inside an edge, no two edges
# crossing.
function vertex_on_edge_problem(positions, es, i, n)
    a, b = src(es[i]), dst(es[i])
    for v = 1:n
        v in (a, b) && continue
        on_segment(positions[v], positions[a], positions[b]) &&
            return "vertex $v lies on edge $a-$b"
    end
    return nothing
end

function edge_crossing_problem(positions, es, i)
    a, b = src(es[i]), dst(es[i])
    for j = (i+1):length(es)
        c, e = src(es[j]), dst(es[j])
        (a == c || a == e || b == c || b == e) && continue
        segments_cross(positions[a], positions[b], positions[c], positions[e]) &&
            return "edges $a-$b and $c-$e cross"
    end
    return nothing
end

function drawing_problem(g::SimpleGraph, d::GridDrawing)
    n = nv(g)
    length(d) == n || return "wrong number of positions"
    positions = coordinates(d)
    allunique(positions) || return "two vertices share a position"
    all(p -> p[1] isa Int && p[2] isa Int, positions) ||
        return "coordinates are not integers"
    all(p -> 0 <= p[1] <= width(d) && 0 <= p[2] <= height(d), positions) ||
        return "a vertex is outside the bounding box"

    es = collect(edges(g))
    for i in eachindex(es)
        problem = vertex_on_edge_problem(positions, es, i, n)
        problem === nothing || return problem
        problem = edge_crossing_problem(positions, es, i)
        problem === nothing || return problem
    end
    return nothing
end

drawing_problem(g::PlanarGraph, d::GridDrawing) = drawing_problem(g.graph, d)

@testset "grid drawings" begin
    @testset "known graphs" begin
        graphs = [
            PlanarGraph([[2, 3], [3, 1], [1, 2]]),
            planar_path(6),
            planar_cycle(7),
            planar_wheel(8),
            planar_grid(4, 4),
            planar_grid(3, 5; diagonals = true),
        ]
        for n in (4, 5, 8, 30)
            push!(graphs, random_triangulation(n; rng = MersenneTwister(n)))
        end
        for g in graphs
            d = grid_drawing(g)
            @test d isa GridDrawing
            @test length(d) == nv(g)
            @test drawing_problem(g, d) === nothing
        end
    end

    @testset "small graphs" begin
        empty = SimpleGraph(0)
        d = grid_drawing(empty)
        @test length(d) == 0
        @test width(d) == 0 && height(d) == 0

        single = SimpleGraph(1)
        d = grid_drawing(single)
        @test length(d) == 1
        @test d[1] == (0, 0)

        pair = SimpleGraph(2)
        add_edge!(pair, 1, 2)
        d = grid_drawing(pair)
        @test d[1] != d[2]
        @test drawing_problem(pair, d) === nothing

        triangles = SimpleGraph(6)
        for (u, v) in ((1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4))
            add_edge!(triangles, u, v)
        end
        d = grid_drawing(triangles)
        @test drawing_problem(triangles, d) === nothing
    end

    @testset "trees and disconnected graphs" begin
        rng = MersenneTwister(9)
        for n in (2, 3, 10, 40)
            g = SimpleGraph(n)
            for v = 2:n
                add_edge!(g, v, rand(rng, 1:(v-1)))
            end
            d = grid_drawing(g)
            @test drawing_problem(g, d) === nothing
        end

        # Several components plus isolated vertices.
        g = SimpleGraph(15)
        for (u, v) in ((1, 2), (2, 3), (3, 1), (4, 5), (5, 6), (6, 4), (7, 8), (12, 13))
            add_edge!(g, u, v)
        end
        d = grid_drawing(g)
        @test drawing_problem(g, d) === nothing
        @test all(p -> 0 <= p[1] <= width(d), coordinates(d))
    end

    @testset "planar graphs with removed edges" begin
        rng = MersenneTwister(13)
        for trial = 1:15
            n = rand(rng, 4:12)
            g = copy(random_triangulation(n; rng = rng).graph)
            for e in collect(edges(g))
                rand(rng) < 0.4 && rem_edge!(g, src(e), dst(e))
            end
            d = grid_drawing(g)
            @test drawing_problem(g, d) === nothing
        end
    end

    @testset "area is linear in both directions" begin
        for n in (5, 20, 60, 200)
            g = random_triangulation(n; rng = MersenneTwister(n))
            d = grid_drawing(g)
            @test width(d) <= 2 * (n - 2)
            @test height(d) <= n - 2
            @test width(d) <= 2 * n
            @test height(d) <= n
            # A triangulation fills the box exactly.
            @test width(d) == 2 * (n - 2)
            @test height(d) == n - 2
        end

        # Many small components pack side by side, not on top of each other.
        g = SimpleGraph(60)
        for v = 1:3:58
            add_edge!(g, v, v + 1)
            add_edge!(g, v + 1, v + 2)
            add_edge!(g, v + 2, v)
        end
        d = grid_drawing(g)
        @test drawing_problem(g, d) === nothing
        @test width(d) <= 2 * nv(g)
        @test height(d) <= nv(g)
    end

    @testset "non-planar input" begin
        @test_throws ArgumentError grid_drawing(complete_graph(5))
        @test_throws ArgumentError grid_drawing(complete_bipartite_graph(3, 3))
        @test_throws ArgumentError grid_drawing(smallgraph(:petersen))
        k5 = PlanarGraph([
            [2, 3, 4, 5],
            [1, 3, 4, 5],
            [1, 2, 4, 5],
            [1, 2, 3, 5],
            [1, 2, 3, 4],
        ])
        @test_throws ArgumentError grid_drawing(k5)

        # A rotation that is not an embedding of a planar graph is not an error:
        # the embedding is recomputed and the graph is drawn.
        broken = PlanarGraph([[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]])
        @test !is_planar_embedding(broken)
        d = grid_drawing(broken)
        @test drawing_problem(complete_graph(4), d) === nothing
    end

    @testset "the drawing type" begin
        g = planar_grid(3, 3)
        d = grid_drawing(g)
        @test length(d) == 9
        @test coordinates(d) == [d[v] for v = 1:9]
        @test all(p -> 0 <= p[1] <= width(d) && 0 <= p[2] <= height(d), coordinates(d))
        @test occursin("GridDrawing(9 vertices", sprint(show, d))
        # The same input gives the same drawing.
        @test coordinates(grid_drawing(g)) == coordinates(d)
        # Vertex labels are untouched: every original vertex has a position.
        @test length(d.x) == 9 && length(d.y) == 9
    end

    @testset "drawing from an unembedded graph" begin
        g = smallgraph(:dodecahedral)
        @test drawing_problem(g, grid_drawing(g)) === nothing
        # A planar graph with a broken rotation: the embedding is recomputed.
        broken = PlanarGraph([
            [2, 3, 4, 5],
            [1, 3, 4],
            [1, 2, 4],
            [1, 2, 3],
            [6, 7, 8, 1],
            [5, 7, 8],
            [5, 6, 8],
            [5, 6, 7],
        ])
        d = grid_drawing(broken)
        @test drawing_problem(broken, d) === nothing
    end
end
