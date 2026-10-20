module DelaunayTriangulation

include("setup.jl")

import ExactPredicates as EP
import AdaptivePredicates as AP
import EnumX
import Random

abstract type AbstractPredicateKernel end # needs to be defined early for use in the data structures below
const PredicateCacheType = Union{Nothing, <:Tuple}

# Data structures. These are loaded first since the predicates and algorithms use their types.
include("utils/max_priority_queue.jl")
include("utils/queue.jl")
include("trees/rtree.jl")
include("trees/polygon_hierarchy.jl")
include("triangulation/adjacent.jl")
include("triangulation/adjacent2vertex.jl")
include("triangulation/graph.jl")
include("curves/abstract.jl")
include("curves/linesegment.jl")
include("curves/piecewiselinear.jl")
include("curves/circulararc.jl")
include("curves/ellipticalarc.jl")
include("curves/beziercurve.jl")
include("curves/bspline.jl")
include("curves/catmullromspline.jl")
include("triangulation/representative_coordinates.jl")
include("geometry/pole_of_inaccessibility/cell.jl")
include("geometry/pole_of_inaccessibility/cell_queue.jl")
include("geometry/convex_hull/convex_hull.jl")
include("triangulation/triangulation.jl")
include("triangulation/triangulation_cache.jl")
include("curves/boundary_enricher.jl")
include("triangulation/methods/adjacent.jl")
include("triangulation/methods/adjacent2vertex.jl")
include("triangulation/methods/boundary_curves.jl")
include("triangulation/methods/boundary_edge_map.jl")
include("triangulation/methods/boundary_nodes.jl")
include("triangulation/methods/convex_hull.jl")
include("triangulation/methods/exterior_curve_indices.jl")
include("triangulation/methods/ghost_vertex_map.jl")
include("triangulation/methods/ghost_vertex_ranges.jl")
include("triangulation/methods/graph.jl")
include("triangulation/methods/iterators.jl")
include("triangulation/methods/points.jl")
include("triangulation/methods/representative_point_list.jl")
include("triangulation/methods/segments.jl")
include("triangulation/methods/triangles.jl")
include("triangulation/methods/weights.jl")
include("triangulation/methods/checks.jl")
include("point_location/point_location_history.jl")
include("statistics/individual_triangle_statistics.jl")
include("statistics/triangulation_statistics.jl")
include("refinement/insertion_event_history.jl")
include("refinement/refinement_constraints.jl")
include("refinement/refinement_queue.jl")
include("refinement/refinement_arguments.jl")
include("voronoi/tessellation.jl")
include("geometry/polygon_clipping/polygon.jl")
include("triangulation/construction/shuffled_polygon_linked_list.jl")

# Predicates
include("predicates/kernels.jl")
include("predicates/certificate.jl")
include("predicates/predicate_definitions.jl")
include("predicates/predicates.jl")
include("predicates/boundaries_and_ghosts.jl")

# Geometric primitives: the interface for points, edges, triangles and boundary nodes
include("primitives/boundary_nodes.jl")
include("primitives/edges.jl")
include("primitives/points.jl")
include("primitives/triangles.jl")

# Utilities
include("geometry/geometry_utils.jl")
include("utils/utils.jl")

# Algorithms
include("trees/rtree_intersections.jl")
include("geometry/convex_hull/monotone_chain.jl")
include("point_location/brute_force.jl")
include("point_location/jump_and_march.jl")
include("point_location/nearest_neighbour.jl")
include("point_location/find_polygon.jl")
include("geometry/polygon_clipping/liang_barsky.jl")
include("geometry/polygon_clipping/sutherland_hodgman.jl")
include("geometry/pole_of_inaccessibility/pole_of_inaccessibility.jl")
include("triangulation/construction/constrained.jl")
include("triangulation/construction/check_args.jl")
include("triangulation/construction/triangulate.jl")
include("triangulation/construction/curve_bounded.jl")
include("refinement/mesh_refinement.jl")
include("triangulation/construction/convex.jl")
include("triangulation/construction/rectangle.jl")
include("triangulation/construction/unconstrained.jl")
include("triangulation/operations/add_boundary_information.jl")
include("triangulation/operations/add_ghost_triangles.jl")
include("triangulation/operations/add_point.jl")
include("triangulation/operations/add_segment.jl")
include("triangulation/operations/add_triangle.jl")
include("triangulation/operations/clear_empty_features.jl")
include("triangulation/operations/delete_ghost_triangles.jl")
include("triangulation/operations/delete_holes.jl")
include("triangulation/operations/delete_point.jl")
include("triangulation/operations/delete_triangle.jl")
include("triangulation/operations/flip_edge.jl")
include("triangulation/operations/legalise_edge.jl")
include("triangulation/operations/lock_convex_hull.jl")
include("triangulation/operations/split_edge.jl")
include("triangulation/operations/split_triangle.jl")
include("triangulation/operations/unlock_convex_hull.jl")
include("voronoi/centroidal.jl")
include("voronoi/clipped_coordinates.jl")
include("voronoi/clipped.jl")
include("voronoi/construction.jl")
include("voronoi/unbounded.jl")

# Validation
include("triangulation/validation.jl")

include("exports.jl")
include("public.jl")

end
