module DelaunayTriangulation

using SimpleGraphs
using ExactPredicates
using Random
import DataStructures: DefaultDict

include("constants.jl")
include("interface.jl")
include("data_structures.jl")
include("predicates.jl")
include("utils.jl")
include("operations.jl")
include("bowyerwatson.jl")
include("locate_triangle.jl")
include("deberg.jl")

export adjacent, adjacent2vertex, graph, pointlocation, triangles, points, neighbours
export triangulate, triangulate!
export UnconstrainedTriangulation
export num_points, num_triangles, edges, num_edges
export add_point!
export geti, getj, getk, indices
export getx, gety
export get_edge
export locate_triangle, jump_and_march

end
