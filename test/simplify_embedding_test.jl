# Incidence-preserving simplification: oracle, invariances, validation, sparse IDs.

const IE = IncidenceEdge

# Independent pairwise selector: an edge survives unless a same-pair edge has a smaller key.
function pairwise_survivors(edges)
    pair(e) = minmax(e.u, e.v)
    key(e) = min(e.at_u, e.at_v)
    kept = Set{Int}()
    for e in edges
        e.u == e.v && continue
        any(f -> pair(f) == pair(e) && key(f) < key(e), edges) && continue
        push!(kept, e.at_u, e.at_v)
    end
    return kept
end

function oracle(rotations, edges)
    kept = pairwise_survivors(edges)
    other = Dict{Int,Int}()
    owner = Dict{Int,Int}()
    for e in edges
        other[e.at_u] = e.at_v
        other[e.at_v] = e.at_u
        owner[e.at_u] = e.u
        owner[e.at_v] = e.v
    end
    inc = [[id for id in r if id in kept] for r in rotations]
    return [[owner[other[id]] for id in l] for l in inc], inc
end

# Random valid multigraph with scattered or compact IDs and shuffled rotations.
function random_multigraph(rng, n, m; spread = 1)
    ids = randperm(rng, spread * 2m)[1:(2m)]
    edges = [IE(rand(rng, 1:n), rand(rng, 1:n), ids[2i-1], ids[2i]) for i = 1:m]
    rotations = [Int[] for _ = 1:n]
    for e in edges
        push!(rotations[e.u], e.at_u)
        push!(rotations[e.v], e.at_v)
    end
    return shuffle!.(rng, rotations), edges
end

flip(e) = IE(e.v, e.u, e.at_v, e.at_u)

@testset "simplify_embedding" begin
    @testset "specification example" begin
        edges = [IE(1, 2, 7, 4), IE(1, 2, 2, 9), IE(1, 1, 5, 8)]
        rot = [[7, 5, 2, 8], [4, 9]]
        for g in (simplify_embedding(rot, edges), simplify_embedding_sparse(rot, edges))
            @test g.neighbors == [[2], [1]]
            @test g.incidences == [[2], [9]]
        end
    end

    @testset "empty, isolates, loops only" begin
        g = simplify_embedding(Vector{Int}[], IE[])
        @test isempty(g.neighbors) && isempty(g.incidences)
        g = simplify_embedding([Int[], Int[], Int[]], IE[]; universe_size = 0)
        @test g.neighbors == [Int[], Int[], Int[]]
        g = simplify_embedding_sparse([[3, 1], Int[]], [IE(1, 1, 1, 3)])
        @test g.neighbors == [Int[], Int[]] && g.incidences == [Int[], Int[]]
    end

    @testset "random vs pairwise oracle, invariances" begin
        rng = Random.MersenneTwister(7)
        for trial = 1:300
            n = rand(rng, 1:7)
            m = rand(rng, 0:14)
            rot, edges = random_multigraph(rng, n, m; spread = rand(rng, 1:3))
            nb, inc = oracle(rot, edges)
            g = simplify_embedding(rot, edges)
            @test (g.neighbors, g.incidences) == (nb, inc)
            s = simplify_embedding_sparse(rot, edges)
            @test (s.neighbors, s.incidences) == (nb, inc)
            perm = [flip(e) for e in shuffle(rng, edges)]
            p = simplify_embedding(rot, perm)
            @test (p.neighbors, p.incidences) == (nb, inc)
            # rotating a list rotates its filtered result
            v = rand(rng, 1:n)
            r2 = [copy(r) for r in rot]
            circshift!(r2[v], 1)
            nb2, inc2 = oracle(r2, edges)
            @test simplify_embedding(r2, edges).incidences == inc2
        end
    end

    @testset "many parallel edges and stars" begin
        edges = [IE(1, 2, 10i + 3, 10i + 4) for i = 1:50]
        push!(edges, IE(2, 1, 7, 5))
        rot = [[e.at_u for e in edges[1:50]], [e.at_v for e in edges[1:50]]]
        push!(rot[2], 7)
        push!(rot[1], 5)
        g = simplify_embedding(rot, edges)
        @test g.incidences == [[5], [7]]
        star = [IE(1, k, 2k, 2k + 1) for k = 2:6]
        rot = [[2k for k = 2:6], [[2k + 1] for k = 2:6]...]
        g = simplify_embedding(rot, star)
        @test g.neighbors[1] == 2:6 && g.neighbors[3] == [1]
    end

    @testset "dense equals sparse, IDs near typemax" begin
        top = typemax(Int)
        edges = [IE(1, 2, top, top - 5), IE(2, 1, top - 7, top - 1), IE(2, 2, 3, top - 2)]
        rot = [[top, top - 1], [top - 7, top - 5, 3, top - 2]]
        s = simplify_embedding_sparse(rot, edges)
        @test s.incidences == [[top - 1], [top - 7]]
        @test s.neighbors == [[2], [1]]
        rng = Random.MersenneTwister(3)
        rot, edges = random_multigraph(rng, 6, 12; spread = 2)
        d = simplify_embedding(rot, edges; universe_size = 48)
        s = simplify_embedding_sparse(rot, edges)
        @test d.neighbors == s.neighbors && d.incidences == s.incidences
    end

    @testset "malformed input" begin
        ok = [[1], [2]]
        e = [IE(1, 2, 1, 2)]
        for f in (simplify_embedding, simplify_embedding_sparse)
            @test_throws ArgumentError f(ok, [IE(1, 3, 1, 2)])
            @test_throws ArgumentError f(ok, [IE(0, 2, 1, 2)])
            @test_throws ArgumentError f([[0], [2]], [IE(1, 2, 0, 2)])
            @test_throws ArgumentError f([[-1], [2]], [IE(1, 2, -1, 2)])
            @test_throws ArgumentError f([[1], [1]], [IE(1, 2, 1, 1)])
            @test_throws ArgumentError f([[1, 3], [2, 4]], [IE(1, 2, 1, 2), IE(1, 2, 1, 4)])
            @test_throws ArgumentError f([[1], Int[]], e)
            @test_throws ArgumentError f([[1, 1], [2]], e)
            @test_throws ArgumentError f([[2], [1]], e)
            @test_throws ArgumentError f([[1, 9], [2]], e)
            @test_throws ArgumentError f([[1], [2, 5]], e)
        end
        @test_throws ArgumentError simplify_embedding(ok, e; universe_size = -1)
        @test_throws ArgumentError simplify_embedding(ok, e; universe_size = 1)
    end

    @testset "immutability and output independence" begin
        rot = [[7, 5, 2, 8], [4, 9]]
        edges = [IE(1, 2, 7, 4), IE(1, 2, 2, 9), IE(1, 1, 5, 8)]
        rot0, edges0 = deepcopy(rot), copy(edges)
        g = simplify_embedding(rot, edges)
        @test rot == rot0 && edges == edges0
        rot[1][1] = 99
        @test g.incidences == [[2], [9]]
        @test_throws ArgumentError simplify_embedding([[1], [3]], [IE(1, 2, 1, 2)])
    end

    @testset "feeds separator_hierarchy" begin
        w = planar_wheel(9)
        rot = [neighbor_order(w, v) for v = 1:nv(w)]
        ids = Dict{Tuple{Int,Int},Int}()
        edges = IE[]
        for v = 1:nv(w), u in rot[v]
            v < u || continue
            k = length(ids) + 1
            ids[(v, u)] = 2k - 1
            ids[(u, v)] = 2k
            push!(edges, IE(v, u, 2k - 1, 2k))
            push!(edges, IE(u, v, 4 * nv(w)^2 + 2k, 4 * nv(w)^2 + 2k + 1))
        end
        rotations = [Int[ids[(v, u)] for u in rot[v]] for v = 1:nv(w)]
        for e in edges[2:2:end]  # parallel copies, appended after the originals
            push!(rotations[e.u], e.at_u)
            push!(rotations[e.v], e.at_v)
        end
        g = simplify_embedding(rotations, edges)
        @test g.neighbors == rot
        h = separator_hierarchy(g.neighbors)
        @test h.root == 1
    end
end
