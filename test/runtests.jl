using SimpleGraphs
using ExactPredicates
using Test
using Random
using DelaunayTriangulation
const DT = DelaunayTriangulation
using StaticArrays
import DataStructures: DefaultDict

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

# @testset "All" begin include("algorithm.jl") end
include("functions.jl")
include("interface.jl")
include("data_structures.jl")
include("predicates.jl")
include("add_triangle.jl")
include("delete_triangle.jl")
include("flip_edge.jl")
include("split_triangle.jl")
include("utils.jl")
include("deberg.jl")
include("split_edge.jl")
include("legalise_edge.jl")
include("point_location.jl")
include("add_triangle_ghosts.jl")
include("delete_triangle_ghosts.jl")
include("incircle_ghosts.jl")
include("centroid.jl")
include("intriangle_ghosts.jl")
include("point_location_ghosts.jl")
include("bowyerwatson.jl")
include("addremove_ghost_triangles_to_triangulation.jl")
include("point_exclusion.jl")