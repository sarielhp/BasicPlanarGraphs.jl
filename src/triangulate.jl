# Triangulation of an embedded graph by inserting chords inside faces.

"""
    triangulate(g) -> same type as g

Add edges to `g` until every face is a triangle, preserving the embedding and
the vertex set. The result is a maximal planar graph when the input is
2-connected.

Chords are inserted inside faces of the given rotation system, so the result is
verified before it is returned: when a face is not a simple cycle and the
insertion cannot go on, or when the rotation system does not describe a planar
embedding in the first place, an `ArgumentError` is raised instead of a graph
that is not a triangulation. For a [`PlanarGraph`](@ref) the result is returned
as a `PlanarGraph`.
"""
function triangulate(g::DartGraph)
    h = _add_chords!(copy(g))
    _verify_triangulation(h)
    return h
end

triangulate(g::PlanarGraph) = PlanarGraph(triangulate(DartGraph(g)))

function _add_chords!(h::DartGraph)
    orbits, _ = face_orbits(h)
    stack = orbits
    adjacency = [Set{Int}() for _ = 1:h.n]
    for d = 1:ndarts(h)
        push!(adjacency[tail(h, d)], h.head[d])
    end
    while !isempty(stack)
        orbit = pop!(stack)
        length(orbit) <= 3 && continue
        chord = _find_chord(orbit, h, adjacency)
        chord === nothing && continue
        faces_added = _insert_chord!(h, orbit, chord[1], chord[2], adjacency)
        push!(stack, faces_added[1])
        push!(stack, faces_added[2])
    end
    return h
end

# Chord insertion is only sound for faces that are simple cycles; make sure the
# result really is a triangulation of a valid embedding.
function _verify_triangulation(h::DartGraph)
    orbits, _ = face_orbits(h)
    all(orbit -> length(orbit) == 3, orbits) ||
        throw(ArgumentError("cannot triangulate: some face is not a simple cycle"))
    euler_holds(h) ||
        throw(ArgumentError("cannot triangulate: the rotation system is not planar"))
    return nothing
end

# First chord that splits the face `orbit` into two proper faces.
function _find_chord(orbit::Vector{Int}, h::DartGraph, adjacency)
    k = length(orbit)
    for i = 1:k
        for step = 2:(k-2)
            j = mod1(i + step, k)
            vi = tail(h, orbit[i])
            vj = tail(h, orbit[j])
            vi == vj && continue
            vj in adjacency[vi] && continue
            return (i, j)
        end
    end
    return nothing
end

# Insert the chord (v_i, v_j) into `orbit`; return the two new face orbits.
function _insert_chord!(h::DartGraph, orbit::Vector{Int}, i::Int, j::Int, adjacency)
    k = length(orbit)
    di = orbit[i]
    dj = orbit[j]
    vi = tail(h, di)
    vj = tail(h, dj)
    forward = add_dart!(h, vi, vj, di)
    backward = add_dart!(h, vj, vi, dj)
    h.twin[forward] = backward
    h.twin[backward] = forward
    push!(adjacency[vi], vj)
    push!(adjacency[vj], vi)
    return _split_face(orbit, i, j, forward, backward)
end

function _split_face(orbit::Vector{Int}, i::Int, j::Int, forward::Int, backward::Int)
    k = length(orbit)
    first_length = mod(i - j, k)
    second_length = mod(j - i, k)
    first_face = Vector{Int}(undef, first_length + 1)
    first_face[1] = forward
    for t = 1:first_length
        first_face[t+1] = orbit[mod1(j + t - 1, k)]
    end
    second_face = Vector{Int}(undef, second_length + 1)
    for t = 1:second_length
        second_face[t] = orbit[mod1(i + t - 1, k)]
    end
    second_face[end] = backward
    return (first_face, second_face)
end
