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

# Tube of circumference 4 and length L from vertex 1 (r) to a blob: rim of
# m = 4q vertices fanned from the last ring, capped by one vertex z. All faces
# are triangles. Rotations are counter-clockwise for rings drawn concentrically.
function tube_blob(L, q; zrev = true)
    m = 4q
    r = 1
    a(i) = 1 + i
    b(i) = 1 + L + i
    c(i) = 1 + 2L + i
    d(i) = 1 + 3L + i
    w(j) = 1 + 4L + mod1(j, m)
    z = 2 + 4L + m
    n = z
    rot = [Int[] for _ = 1:n]
    rot[r] = [a(1), b(1), c(1), d(1)]
    for i = 1:L
        ap = i == 1 ? r : a(i - 1)
        bp = i == 1 ? r : b(i - 1)
        cp = i == 1 ? r : c(i - 1)
        dp = i == 1 ? r : d(i - 1)
        if i < L
            rot[a(i)] =
                i == 1 ? [a(2), b(2), b(1), r, d(1), d(2)] :
                [a(i + 1), b(i + 1), b(i), ap, d(i), d(i + 1)]
            rot[b(i)] =
                i == 1 ? [b(2), c(1), r, a(1)] :
                [b(i + 1), c(i), c(i - 1), bp, a(i - 1), a(i)]
            rot[c(i)] =
                i == 1 ? [c(2), d(2), d(1), r, b(1), b(2)] :
                [c(i + 1), d(i + 1), d(i), cp, b(i), b(i + 1)]
            rot[d(i)] =
                i == 1 ? [d(2), a(1), r, c(1)] :
                [d(i + 1), a(i), a(i - 1), dp, c(i - 1), c(i)]
        else
            rot[a(L)] = vcat([w(4q)], [w(j) for j = 1:q], [b(L), ap, d(L)])
            rot[b(L)] = vcat([w(j) for j = q:2q], [c(L), c(L - 1), bp, a(L - 1), a(L)])
            rot[c(L)] = vcat([w(j) for j = 2q:3q], [d(L), cp, b(L)])
            rot[d(L)] = vcat([w(j) for j = 3q:4q], [a(L), a(L - 1), dp, c(L - 1), c(L)])
        end
    end
    corner =
        Dict(q => (a(L), b(L)), 2q => (b(L), c(L)), 3q => (c(L), d(L)), 4q => (d(L), a(L)))
    for j = 1:m
        if haskey(corner, j)
            X, Y = corner[j]
            rot[w(j)] = [z, w(j + 1), Y, X, w(j - 1)]
        else
            X = j < q ? a(L) : j < 2q ? b(L) : j < 3q ? c(L) : d(L)
            rot[w(j)] = [z, w(j + 1), X, w(j - 1)]
        end
    end
    rot[z] = zrev ? [w(j) for j = m:-1:1] : [w(j) for j = 1:m]
    return PlanarGraph(rot)
end
