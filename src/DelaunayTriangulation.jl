module DelaunayTriangulation

using SimpleGraphs
using ExactPredicates
using Random
using MutableNamedTuples
import DataStructures: DefaultDict

using Scanf
using Printf
using ElasticArrays

include("constants.jl")
include("interface.jl")
include("data_structures.jl")
include("predicates.jl")
include("utils.jl")
include("operations.jl")
include("bowyerwatson.jl")
include("locate_triangle.jl")
include("deberg.jl")
include("gmsh.jl")

export adjacent, adjacent2vertex, graph, pointlocation, triangles, points, neighbours
export triangulate, triangulate!
export UnconstrainedTriangulation
export edges
export add_point!
export geti, getj, getk, indices
export getx, gety
export get_edge
export locate_triangle, jump_and_march
export generate_mesh, triangulate 
export num_triangles, num_points
export number_type

end
