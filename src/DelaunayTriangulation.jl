module DelaunayTriangulation

using SimpleGraphs 
using ExactPredicates
using ExactPredicates.Codegen 
using Random

#=
const TriangleType = NTuple{3, Int64}
const EdgeType = NTuple{2, Int64}
const LargeRightIdx = 0 # p₋₁
const LargeLeftIdx = -1 # p₋₂
const BoundaryIdx = -2  # ∅
const Adjacent2VertexVector = Set{EdgeType}

include("data_structures.jl")
include("predicates.jl")
include("update_structure.jl")
include("triangulate.jl")
=#

include("algorithm.jl")

export Edge, Triangle, Point
export initial, terminal
export indices, geti, getj, getk 
export coords, getx, gety
export Triangles, triangles 
export Points, points, add_point!
export Triangulation 
export adjacent, adjacent2vertex, graph, history, triangles, points, neighbours
export triangulate
export num_points, num_triangles, edges, num_edges
end
