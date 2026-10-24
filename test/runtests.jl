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
    include("scaling_test.jl")

    @testset "Aqua" begin
        Aqua.test_all(BasicPlanarGraphs)
    end
end
