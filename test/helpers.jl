# Helpers shared by the test files.

# Rotation systems of several graphs side by side, vertex numbers shifted.
function disjoint_union(parts::PlanarGraph...)
    rot = Vector{Vector{Int}}()
    offset = 0
    for part in parts
        append!(rot, [Int[w + offset for w in r] for r in rotation(part)])
        offset += nv(part)
    end
    return PlanarGraph(rot)
end

# Delete each edge with probability `p`, keeping the rotation order.
function thin_out(g::PlanarGraph, p::Real, rng)
    gone = Set{Tuple{Int,Int}}()
    for e in edges(g)
        rand(rng) < p && push!(gone, (src(e), dst(e)))
    end
    kept(v, w) = !((v, w) in gone || (w, v) in gone)
    rot = rotation(g)
    return PlanarGraph([Int[w for w in rot[v] if kept(v, w)] for v = 1:nv(g)])
end

# The fastest of several runs after a warm-up, with a garbage collection before
# each, so that a pause or a busy core does not decide the comparison.
function fastest(f, g; repeats = 5)
    f(g)
    return minimum(begin
        GC.gc()
        @elapsed f(g)
    end for _ = 1:repeats)
end

# A triangulated tube of circumference 4 and length L whose last ring is fanned to
# a rim of 4q vertices capped by one vertex; vertex 1 is the apex at the other end.
# Once 4q is above about 8L most faces are beyond the last BFS levels, so every rung
# of the ladder is light, and the fundamental cycle runs down the whole tube. The
# rotations are written out: a triangulation has one embedding, but the BFS and the
# dual tree, and with them the bug, depend on which one the rotation starts from.
tube_ring(i, c) = i == 0 ? 1 : 1 + 4 * (i - 1) + mod1(c, 4)
tube_rim(L, q, j) = 1 + 4L + mod1(j, 4q)

# The rim vertices that ring vertex `c` of the last ring is joined to.
tube_fan(L, q, c) = [tube_rim(L, q, j) for j = ((c-1)*q):(c*q)]

# Counter-clockwise neighbours of ring vertex `c` of ring `i`: odd columns turn
# one way, even columns the other; the entries beyond the last ring are the fan.
function tube_ring_rotation(L, q, i, c)
    up(col) = i < L ? [tube_ring(i + 1, col)] : Int[]
    fan = i < L ? Int[] : tube_fan(L, q, c)
    same(col) = tube_ring(i, col)
    below(col) = i > 1 ? [tube_ring(i - 1, col)] : Int[]
    prev = tube_ring(i - 1, c)
    if isodd(c)
        return vcat(up(c), fan, up(c + 1), same(c + 1), prev, same(c - 1), up(c - 1))
    end
    return vcat(up(c), fan, same(c + 1), below(c + 1), prev, below(c - 1), same(c - 1))
end

# The rim vertex `j`: the cap, the next rim vertex, the ring vertex or two it is
# joined to, and the previous rim vertex.
function tube_rim_rotation(L, q, j)
    cap = 2 + 4L + 4q
    c = (j - 1) ÷ q + 1
    joined = j % q == 0 ? [tube_ring(L, c + 1), tube_ring(L, c)] : [tube_ring(L, c)]
    return vcat(cap, tube_rim(L, q, j + 1), joined, tube_rim(L, q, j - 1))
end

function tube_blob(L, q)
    cap = 2 + 4L + 4q
    rot = [Int[] for _ = 1:cap]
    rot[1] = [tube_ring(1, c) for c = 1:4]
    for i = 1:L, c = 1:4
        rot[tube_ring(i, c)] = tube_ring_rotation(L, q, i, c)
    end
    for j = 1:(4q)
        rot[tube_rim(L, q, j)] = tube_rim_rotation(L, q, j)
    end
    rot[cap] = [tube_rim(L, q, j) for j = (4q):-1:1]
    return PlanarGraph(rot)
end
