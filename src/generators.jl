# Small embedded graph constructors, mainly for examples and tests.

"""
    planar_path(n) -> PlanarGraph

Path on `n` vertices, embedded in the plane.
"""
function planar_path(n::Integer)
    n >= 1 || throw(ArgumentError("n must be positive, got $n"))
    rot = [Int[] for _ = 1:n]
    for v = 1:n
        v > 1 && push!(rot[v], v - 1)
        v < n && push!(rot[v], v + 1)
    end
    return PlanarGraph(rot)
end

"""
    planar_cycle(n) -> PlanarGraph

Cycle on `n >= 3` vertices, counter-clockwise.
"""
function planar_cycle(n::Integer)
    n >= 3 || throw(ArgumentError("a simple cycle needs at least 3 vertices, got $n"))
    rot = [[mod1(v - 1, n), mod1(v + 1, n)] for v = 1:n]
    return PlanarGraph(rot)
end

"""
    planar_grid(rows, cols; diagonals=false) -> PlanarGraph

`rows x cols` grid embedded in the plane. With `diagonals = true` every cell
gets its top-left to bottom-right diagonal, so every face is a triangle.
"""
function planar_grid(rows::Integer, cols::Integer; diagonals::Bool = false)
    rows >= 1 || throw(ArgumentError("rows must be positive, got $rows"))
    cols >= 1 || throw(ArgumentError("cols must be positive, got $cols"))
    n = rows * cols
    rot = [Int[] for _ = 1:n]
    for r = 1:rows, c = 1:cols
        v = _grid_id(r, c, cols)
        for (dr, dc) in _grid_directions(diagonals)
            rr = r + dr
            cc = c + dc
            (1 <= rr <= rows && 1 <= cc <= cols) || continue
            push!(rot[v], _grid_id(rr, cc, cols))
        end
    end
    return PlanarGraph(rot)
end

_grid_id(r::Integer, c::Integer, cols::Integer) = (r - 1) * cols + c

# Counter-clockwise neighbour directions.
function _grid_directions(diagonals::Bool)
    diagonals || return ((0, 1), (-1, 0), (0, -1), (1, 0))
    return ((0, 1), (-1, 0), (-1, -1), (0, -1), (1, 0), (1, 1))
end

"""
    planar_wheel(rim) -> PlanarGraph

Wheel with a `rim`-cycle and a hub vertex `1` in the centre.
"""
function planar_wheel(rim::Integer)
    rim >= 3 || throw(ArgumentError("a wheel rim needs at least 3 vertices, got $rim"))
    n = rim + 1
    rot = [Int[] for _ = 1:n]
    rot[1] = collect(2:n)
    for i = 1:rim
        v = i + 1
        nxt = mod1(i + 1, rim) + 1
        prv = mod1(i - 1, rim) + 1
        rot[v] = [1, prv, nxt]
    end
    return PlanarGraph(rot)
end

"""
    random_triangulation(n; rng=Random.default_rng()) -> PlanarGraph

Random maximal planar graph on `n >= 3` vertices, built by repeatedly inserting
a vertex into a random face.
"""
function random_triangulation(n::Integer; rng = Random.default_rng())
    n >= 3 || throw(ArgumentError("need at least 3 vertices, got $n"))
    dg = _embedded_triangle(n)
    faces, _ = face_orbits(dg)
    for v = 4:n
        index = rand(rng, eachindex(faces))
        face = faces[index]
        faces[index] = faces[end]
        pop!(faces)
        append!(faces, _split_triangle!(dg, face, v))
    end
    return PlanarGraph(dg)
end

# The triangle 1-2-3 on a graph with `n` vertices, counter-clockwise.
function _embedded_triangle(n::Int)
    g = DartGraph(n)
    a = add_dart!(g, 1, 2, 0)
    c2 = add_dart!(g, 1, 3, a)
    b = add_dart!(g, 2, 3, 0)
    a2 = add_dart!(g, 2, 1, b)
    c = add_dart!(g, 3, 1, 0)
    b2 = add_dart!(g, 3, 2, c)
    _link!(g, a, a2)
    _link!(g, b, b2)
    _link!(g, c, c2)
    return g
end

# Insert vertex `v` inside the triangular face `face` (three darts) and return
# the three new face orbits.
function _split_triangle!(g::DartGraph, face::Vector{Int}, v::Int)
    a = tail(g, face[1])
    b = tail(g, face[2])
    c = tail(g, face[3])
    va = add_dart!(g, v, a, 0)
    vb = add_dart!(g, v, b, va)
    vc = add_dart!(g, v, c, vb)
    av = add_dart!(g, a, v, face[1])
    bv = add_dart!(g, b, v, face[2])
    cv = add_dart!(g, c, v, face[3])
    _link!(g, va, av)
    _link!(g, vb, bv)
    _link!(g, vc, cv)
    return Vector{Int}[[face[1], bv, va], [face[2], cv, vb], [face[3], av, vc]]
end

function _link!(g::DartGraph, d::Int, twin::Int)
    g.twin[d] = twin
    g.twin[twin] = d
    return g
end
