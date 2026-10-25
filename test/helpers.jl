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
