using SimpleGraphs
using ExactPredicates
using Test
using Random
using DelaunayTriangulation
const DT = DelaunayTriangulation
using StaticArrays
import DataStructures: DefaultDict

include("functions.jl")
@testset "Interface" begin
    include("interface.jl")
end
@testset "Interface" begin
    include("data_structures.jl")
end
@testset "Geometrical Predicates" begin
    include("predicates.jl")
end
@testset "Operations" begin
    include("operations.jl")
end
@testset "Utility Functions" begin
    include("utils.jl")
end
@testset "Triangulation Algorithms" begin
    include("triangulate.jl")
end
@testset "Point Location" begin
    include("point_location.jl")
end
@testset "Gmsh" begin
    include("gmsh.jl")
end
@testset "Convex Hull" begin
    include("convex_hull.jl")
end