module BasicPlanarGraphs

using Graphs: Graphs, SimpleGraph, degree, edges, has_edge, ne, neighbors, nv, vertices
using Random: Random

export PlanarGraph, planar_graph
export faces, nfaces, dual, is_planar_embedding, neighbor_order, rotation
export triangulate, separator
export planar_path, planar_cycle, planar_grid, planar_wheel, random_triangulation

include("embedded.jl")
include("planar_graph.jl")
include("triangulate.jl")
include("separator.jl")
include("generators.jl")

end
