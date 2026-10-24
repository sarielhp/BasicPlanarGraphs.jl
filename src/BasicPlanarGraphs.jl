module BasicPlanarGraphs

using Graphs:
    Graphs, SimpleGraph, degree, dst, edges, has_edge, ne, neighbors, nv, src, vertices
import Graphs: is_planar
using Random: Random

export PlanarGraph, planar_graph
export faces, nfaces, dual, is_planar_embedding, neighbor_order, rotation
export triangulate, separator, cycle_separator, cycle_sides
export separator_partition, SeparatorPartition
export planar_path, planar_cycle, planar_grid, planar_wheel, random_triangulation
export write_planar_graph, read_planar_graph, write_graph, read_graph
export is_planar, is_planar_linear, planar_embedding
export grid_drawing, GridDrawing, coordinates, width, height

include("embedded.jl")
include("planar_graph.jl")
include("triangulate.jl")
include("separator.jl")
include("separator_partition.jl")
include("cycle_separator.jl")
include("cycle_walk.jl")
include("generators.jl")
include("graph_io.jl")
include("planarity.jl")
include("grid_drawing.jl")

end
