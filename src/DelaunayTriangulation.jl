module DelaunayTriangulation

using SimpleGraphs 
using ExactPredicates
using Random

include("primitives.jl")
include("constants.jl")
include("collections.jl")
include("struct_definitions.jl")
include("struct_updates.jl")
include("predicates.jl")
include("triangulate.jl")
include("utils.jl")

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
