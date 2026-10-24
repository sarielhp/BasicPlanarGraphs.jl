using Aqua
using Graphs
using BasicPlanarGraphs
using Random
using Test

@testset "BasicPlanarGraphs" begin
    include("embedding_test.jl")
    include("triangulate_test.jl")
    include("generators_test.jl")
    include("separator_test.jl")
    include("cycle_separator_test.jl")
    include("scaling_test.jl")
    include("graph_io_test.jl")
    include("planarity_test.jl")
    include("drawing_test.jl")

    @testset "Aqua" begin
        Aqua.test_all(BasicPlanarGraphs)
    end
end
