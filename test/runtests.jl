using SimpleGraphs
using ExactPredicates
using Test
using Random
using DelaunayTriangulation
const DT = DelaunayTriangulation

#=
@testset "Primitives" begin
    include("primitives.jl")
end
@testset "Collections" begin
    include("collections.jl")
end
@testset "Struct definitions" begin
    include("struct_definitions.jl")
end
@testset "Struct updates" begin
    include("struct_updates.jl")
end
@testset "Predicates" begin
    include("predicates.jl")
end
@testset "Triangulation" begin
    include("triangulate.jl")
end
@testset "Utility functions" begin
    include("utils.jl")
end
@testset "Points with labels" begin
    include("cells.jl")
end
=#

@testset "All" begin include("algorithm.jl") end

