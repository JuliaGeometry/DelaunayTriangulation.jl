module DelaunayTriangulation

####################################################
##
## CONSTANTS 
##
####################################################
const DefaultAdjacentValue = 0
const BoundaryIndex = -1
const FirstPointIndex = DefaultAdjacentValue + 1

####################################################
##
## DEPENDENCIES
##
####################################################
using DataStructures
using SimpleGraphs
using ExactPredicates
using EnumX
using Random

####################################################
##
## FILES AND EXPORTS
##
####################################################
include("interfaces/points.jl")
include("interfaces/triangles.jl")
include("interfaces/edges.jl")
include("interfaces/boundary_nodes.jl")

export indices
export num_triangles
export each_triangle
export geti
export getj
export getk
export initial
export terminal
export num_edges
export each_edge
export getx
export gety
export getxy
export each_point
export each_point_index
export get_point
export num_points
export get_boundary_nodes
export num_boundary_edges

include("data_structures/adjacent.jl")
include("data_structures/adjacent2vertex.jl")
include("data_structures/graph.jl")
include("data_structures/convex_hull.jl")
include("data_structures/triangulation/definition.jl")
include("data_structures/triangulation/constructors.jl")
include("data_structures/triangulation/adjacent.jl")
include("data_structures/triangulation/adjacent2vertex.jl")
include("data_structures/triangulation/graph.jl")
include("data_structures/triangulation/convex_hull.jl")
include("data_structures/triangulation/points.jl")
include("data_structures/triangulation/triangles.jl")
include("data_structures/triangulation/edges.jl")
include("data_structures/triangulation/boundary_nodes.jl")
include("data_structures/triangulation/predicates.jl")
include("data_structures/triangulation/representative_points.jl")
include("data_structures/representative.jl")
include("data_structures/statistics.jl")
include("data_structures/refinement/refinement_targets.jl")
include("data_structures/refinement/refinement_queue.jl")
include("data_structures/refinement/event_history.jl")
include("data_structures/point_location_history.jl")
include("data_structures/polylabel/cell.jl")
include("data_structures/polylabel/cell_queue.jl")
include("data_structures/voronoi/voronoi.jl")

export get_adjacent
export get_adjacent2vertex
export get_graph
export get_edges
export get_neighbours
export get_points
export get_triangles
export get_boundary_map
export get_constrained_edges
export get_boundary_nodes
export get_all_constrained_edges
export get_convex_hull
export get_boundary_edge_map
export get_boundary_index_ranges
export Triangulation
export ConvexHull
export convex_hull
export convex_hull!
export each_solid_triangle
export each_ghost_triangle
export num_solid_triangles 
export num_ghost_triangles
export get_vertices
export clear_empty_features!
export get_indices
export get_convex_hull_indices
export each_vertex
export num_vertices
export num_solid_vertices 
export num_ghost_vertices
export each_solid_edge
export each_ghost_edge
export num_solid_edges
export num_ghost_edges
export each_solid_vertex
export each_ghost_vertex
export each_constrained_edge
export statistics
export get_total_area
export get_all_stat
export VoronoiTessellation
export num_polygons
export get_polygon
export each_polygon
export get_polygon_point
export get_polygon_vertex
export get_area
export each_generator
export get_generator
export each_polygon_index
export each_polygon_vertex
export num_polygon_vertices

include("predicates/certificate.jl")
include("predicates/boundaries_and_ghosts.jl")
include("predicates/general.jl")
include("predicates/index_and_ghost_handling.jl")

export Certificate

include("operations/add_triangle.jl")
include("operations/add_boundary_information.jl")
include("operations/add_ghost_triangles.jl")
include("operations/delete_triangle.jl")
include("operations/delete_ghost_triangles.jl")
include("operations/add_point.jl")
include("operations/flip_edge.jl")
include("operations/split_edge.jl")
include("operations/split_triangle.jl")
include("operations/legalise_edge.jl")
include("operations/delete_point.jl")
include("operations/add_edge.jl")
include("operations/delete_holes.jl")
include("operations/clear_empty_features.jl")
include("operations/lock_convex_hull.jl")
include("operations/unlock_convex_hull.jl")

export add_ghost_triangles!
export delete_ghost_triangles!
export add_point!
export add_triangle!
export delete_triangle!
export flip_edge!
export add_boundary_information!
export split_edge!
export split_triangle!
export legalise_edge!
export delete_point!
export add_edge!
export lock_convex_hull!
export unlock_convex_hull!

include("triangulation/rectangle.jl")
include("triangulation/bowyer_watson.jl")
include("triangulation/triangulate.jl")
include("triangulation/convex_triangulation.jl")
include("triangulation/triangulate_constrained.jl")

export triangulate_rectangle
export triangulate
export triangulate_convex

include("point_location/brute_force.jl")
include("point_location/initialisers.jl")
include("point_location/jump_and_march.jl")

export brute_force_search
export jump_and_march

include("constrained_triangulation/segment_location.jl")
include("constrained_triangulation/segment_insertion.jl")

include("utils.jl")

export convert_boundary_points_to_indices

include("geometry_utils/polygons.jl")
include("geometry_utils/polylabel.jl")
include("geometry_utils/intersections.jl")
include("geometry_utils/sutherland_hodgman.jl")

const polylabel = pole_of_inaccessibility

include("refinement/encroachment.jl")
include("refinement/quality_assessment.jl")
include("refinement/refinement_operations.jl")
include("refinement/refine.jl")

export refine!

include("voronoi/main.jl")
include("voronoi/unbounded_construction.jl")
include("voronoi/clipped_construction.jl")
include("voronoi/lloyd.jl")
include("voronoi/coordinates.jl")

export voronoi
export centroidal_smooth
export get_polygon_coordinates

export edge_indices

end