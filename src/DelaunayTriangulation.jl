module DelaunayTriangulation

using SimpleGraphs
using ExactPredicates
using Random
import DataStructures: DefaultDict

#include("primitives.jl")
include("constants.jl")
include("interface.jl")
include("data_structures.jl")
include("predicates.jl")
include("add_triangle.jl")
include("utils.jl")
include("delete_triangle.jl")
include("split_triangle.jl")
include("flip_edge.jl")
include("legalise_edge.jl")
include("split_edge.jl")
include("locate_triangle.jl")
include("deberg.jl")
#include("algorithm.jl")
#include("operations.jl")
#include("collections.jl")
#include("struct_definitions.jl")
#include("struct_updates.jl")
#include("predicates.jl")
#include("triangulate.jl")
#include("utils.jl")

#export Edge, Triangle, Point
#export initial, terminal
#export indices, geti, getj, getk 
#export coords, getx, gety
#export Triangles, triangles 
#export Points, points, add_point!
#export Triangulation 
#export adjacent, adjacent2vertex, graph, history, triangles, points, neighbours
#export triangulate, triangulate!, initialise_triangulation
#export num_points, num_triangles, edges, num_edges
#export get_point

export adjacent, adjacent2vertex, graph, pointlocation, triangles, points, neighbours
export triangulate, triangulate!
export UnconstrainedTriangulation
export num_points, num_triangles, edges, num_edges
export add_point!
export geti, getj, getk, indices
export getx, gety
export get_edge


end
