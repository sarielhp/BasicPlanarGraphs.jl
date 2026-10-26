# The separator hierarchy: ownership, sizes, balance, separation, intervals, routing.

function subtree_sets(h)
    sets = [Int[] for _ in h.nodes]
    for t = length(h.nodes):-1:1
        sets[t] = vcat(h.nodes[t].owned, (sets[c] for c in h.nodes[t].children)...)
    end
    return sets
end

function check_node(h, sets, t)
    nd = h.nodes[t]
    @test nd.size == length(sets[t])
    @test length(nd.owned) <= 2 * sqrt(2 * nd.size) + 1e-9
    @test length(nd.children) <= 2
    isempty(nd.children) && @test length(nd.owned) == nd.size
    isempty(nd.owned) && @test length(nd.children) == 2
    next = t + 1
    for c in nd.children
        @test c == next
        @test 3 * h.nodes[c].size <= 2 * nd.size
        @test h.nodes[c].size >= 1
        next = h.nodes[c].last + 1
    end
    @test next == nd.last + 1
    @test sort(sets[t]) == sort(vcat(nd.owned, (sets[c] for c in nd.children)...))
    return nothing
end

function check_separation(g, h, sets)
    for nd in h.nodes
        length(nd.children) == 2 || continue
        a, b = sets[nd.children[1]], sets[nd.children[2]]
        for e in edges(g)
            @test !(src(e) in a && dst(e) in b) && !(src(e) in b && dst(e) in a)
        end
    end
    for e in edges(g)
        lo, hi = minmax(h.owner[src(e)], h.owner[dst(e)])
        @test hi <= h.nodes[lo].last
    end
    return nothing
end

function check_routing(h)
    for t in eachindex(h.nodes), label = t:h.nodes[t].last
        c = child_containing(h, t, label)
        label == t && @test c == 0
        label != t && @test c in h.nodes[t].children && c <= label <= h.nodes[c].last
    end
    return nothing
end

function check_hierarchy(g, h)
    n = nv(g)
    if n == 0
        @test h.root == 0 && isempty(h.nodes) && isempty(h.owner)
        return nothing
    end
    @test h.root == 1
    @test length(h.nodes) <= 2n - 1
    @test sort(vcat((collect(nd.owned) for nd in h.nodes)...)) == collect(1:n)
    for (t, nd) in enumerate(h.nodes), v in nd.owned
        @test h.owner[v] == t
    end
    sets = subtree_sets(h)
    @test sort(sets[1]) == collect(1:n)
    foreach(t -> check_node(h, sets, t), eachindex(h.nodes))
    check_separation(g, h, sets)
    check_routing(h)
    return nothing
end

@testset "separator hierarchy" begin
    @testset "families" begin
        rng = MersenneTwister(5)
        graphs = [
            planar_path(1),
            planar_path(2),
            planar_path(40),
            planar_cycle(30),
            planar_grid(12, 12),
            planar_grid(6, 20; diagonals = true),
            planar_wheel(50),
            random_triangulation(150; rng = rng),
            PlanarGraph([[2, 3, 4, 5, 6, 7], [1], [1], [1], [1], [1], [1]]),
            PlanarGraph([Int[] for _ = 1:30]),
        ]
        for g in graphs
            h = separator_hierarchy(g)
            check_hierarchy(g, h)
            @test separator_hierarchy(g).owner == h.owner
        end
    end

    @testset "empty and singleton" begin
        h = separator_hierarchy(PlanarGraph(Vector{Int}[]))
        check_hierarchy(PlanarGraph(Vector{Int}[]), h)
        h = separator_hierarchy([Int[]])
        @test h.root == 1 && collect(h.nodes[1].owned) == [1]
        @test h.nodes[1].size == 1 && h.nodes[1].last == 1 && isempty(h.nodes[1].children)
    end

    @testset "disconnected, cut vertices, rotations" begin
        # two triangles joined at nothing, plus a path through a cut vertex
        rot = [[2, 3], [3, 1], [1, 2], [5, 6], [6, 4], [4, 5], [8], [7, 9], [8]]
        g = PlanarGraph(rot)
        check_hierarchy(g, separator_hierarchy(g))
        # unsorted cyclic starts and a graph that has an empty-owned node
        shifted = [circshift(r, 1) for r in rot]
        copy_before = deepcopy(shifted)
        h = separator_hierarchy(shifted)
        @test shifted == copy_before
        check_hierarchy(PlanarGraph(shifted), h)
        @test any(nd -> isempty(nd.owned), h.nodes)
        h2 = separator_hierarchy([Int[] for _ = 1:7])
        check_hierarchy(PlanarGraph([Int[] for _ = 1:7]), h2)
    end

    @testset "both input forms and immutability" begin
        g = planar_grid(9, 9)
        rot = [copy(r) for r in g.rotation]
        edges_before = collect(edges(g.graph))
        h1 = separator_hierarchy(g)
        h2 = separator_hierarchy(rot)
        check_hierarchy(g, h1)
        check_hierarchy(g, h2)
        @test h1.owner == h2.owner
        @test rot == g.rotation && collect(edges(g.graph)) == edges_before
    end

    @testset "routing errors" begin
        h = separator_hierarchy(planar_grid(5, 5))
        c = h.nodes[1].children[1]
        @test child_containing(h, 1, 1) == 0
        @test_throws ArgumentError child_containing(h, c, 1)
        @test_throws ArgumentError child_containing(h, c, h.nodes[c].last + 1)
        @test_throws ArgumentError child_containing(h, 1, length(h.nodes) + 1)
        @test_throws BoundsError child_containing(h, length(h.nodes) + 1, 1)
    end

    @testset "invalid input" begin
        @test_throws ArgumentError separator_hierarchy([[2], [3]])
        @test_throws ArgumentError separator_hierarchy([[1]])
        @test_throws ArgumentError separator_hierarchy([[2], [1, 1]])
        # K5 rotation of higher genus is rejected only with validate
        k5 = [[j for j = 1:5 if j != i] for i = 1:5]
        @test_throws ArgumentError separator_hierarchy(k5; validate = true)
    end
end
