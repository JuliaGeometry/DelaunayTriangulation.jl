module DelaunayTriangulation

"""
    DefaultAdjacentValue = 0 

Default value used for representing an empty result 
from an adjacency query.
"""
const DefaultAdjacentValue = 0

"""
    ∅ = DefaultAdjacentValue

Alias for [`DefaultAdjacentValue`](@ref).
"""
const ∅ = DefaultAdjacentValue

"""
    GhostVertex = -1

Number used for representing initial ghost vertices. 
All other ghost vertices are derived from subtracting from 
this number. See https://juliageometry.github.io/DelaunayTriangulation.jl/stable/manual/ghost_triangles/.
"""
const GhostVertex = -1

"""
    𝒢 = GhostVertex

Alias for [`𝒢`](@ref).
"""
const 𝒢 = GhostVertex

"""
    ε = sqrt(eps(Float64))

Number used as a tolerance in certain functions, e.g. 
for mesh refinement when using [`check_precision`](@ref) to 
avoid degenerate circumcenters.
"""
const ε = sqrt(eps(Float64))

const INF_WARN = Ref(true)
"""
    toggle_inf_warn!()

Toggle the warning for infinite circumcenters in the Voronoi tessellation.
By default, this warning is enabled.
"""
toggle_inf_warn!() = (INF_WARN[] = !INF_WARN[])

@static if VERSION ≥ v"1.6"
    using Preferences
end

@static if VERSION ≥ v"1.6"
    const USE_EXACTPREDICATES = @load_preference("USE_EXACTPREDICATES", true)::Bool
else 
    const USE_EXACTPREDICATES = true 
end
@doc """
    USE_EXACTPREDICATES

Whether to use ExactPredicates.jl for computing predicates. By default, 
this is true, but a user can change this by defining a preference with Preferences.jl, i.e. 
you could do the following 

```julia-repl
julia> using Preferences: set_preferences!

julia> set_preferences!("DelaunayTriangulation", "USE_EXACTPREDICATES" => false)

julia> using DelaunayTriangulation # load only after setting the preference

julia> DelaunayTriangulation.USE_EXACTPREDICATES
false
```

You have set `USE_EXACTPREDICATES = $USE_EXACTPREDICATES`. 

!!! note "Precision"

    Even if you have disabled ExactPredicates.jl, the predicates 
    are still computed in Float64 precision.
"""
USE_EXACTPREDICATES

using ExactPredicates
using EnumX
using Random

include("data_structures/queue/max_priority_queue.jl")
include("data_structures/queue/queue.jl")
include("data_structures/trees/bst.jl")
include("data_structures/trees/rtree.jl")
include("data_structures/trees/polygon_hierarchy.jl")
include("data_structures/triangulation/adjacent.jl")
include("data_structures/triangulation/adjacent2vertex.jl")
include("data_structures/triangulation/graph.jl")
include("data_structures/mesh_refinement/curves.jl")
include("data_structures/representative_coordinates/representative_coordinates.jl")
include("data_structures/representative_coordinates/cell.jl")
include("data_structures/representative_coordinates/cell_queue.jl")
include("data_structures/convex_hull.jl")
include("data_structures/triangulation/triangulation.jl")
include("data_structures/triangulation/triangulation_cache.jl")
include("data_structures/mesh_refinement/boundary_enricher.jl")
include("data_structures/triangulation/methods/adjacent.jl")
include("data_structures/triangulation/methods/adjacent2vertex.jl")
include("data_structures/triangulation/methods/boundary_curves.jl")
include("data_structures/triangulation/methods/boundary_edge_map.jl")
include("data_structures/triangulation/methods/boundary_nodes.jl")
include("data_structures/triangulation/methods/convex_hull.jl")
include("data_structures/triangulation/methods/exterior_curve_indices.jl")
include("data_structures/triangulation/methods/ghost_vertex_map.jl")
include("data_structures/triangulation/methods/ghost_vertex_ranges.jl")
include("data_structures/triangulation/methods/graph.jl")
include("data_structures/triangulation/methods/iterators.jl")
include("data_structures/triangulation/methods/points.jl")
include("data_structures/triangulation/methods/representative_point_list.jl")
include("data_structures/triangulation/methods/segments.jl")
include("data_structures/triangulation/methods/triangles.jl")
include("data_structures/triangulation/methods/weights.jl")
include("data_structures/triangulation/methods/checks.jl")
include("data_structures/point_location_history.jl")
include("data_structures/statistics/individual_triangle_statistics.jl")
include("data_structures/statistics/triangulation_statistics.jl")
include("data_structures/mesh_refinement/insertion_event_history.jl")
include("data_structures/mesh_refinement/refinement_constraints.jl")
include("data_structures/mesh_refinement/refinement_queue.jl")
include("data_structures/mesh_refinement/refinement_arguments.jl")
include("data_structures/voronoi.jl")
include("data_structures/polygon.jl")
include("data_structures/shuffled_polygon_linked_list.jl")

include("geometric_primitives/boundary_nodes.jl")
include("geometric_primitives/edges.jl")
include("geometric_primitives/points.jl")
include("geometric_primitives/triangles.jl")

include("predicates/certificate.jl")
include("predicates/exactpredicates_definitions.jl")
include("predicates/predicates.jl")
include("predicates/boundaries_and_ghosts.jl")

include("utils/geometry_utils.jl")
include("utils/utils.jl")

include("algorithms/intersections/rtree.jl")
include("algorithms/convex_hull.jl")
include("algorithms/point_location/brute_force.jl")
include("algorithms/point_location/jump_and_march.jl")
include("algorithms/point_location/nearest_neighbour.jl")
include("algorithms/point_location/find_polygon.jl")
include("algorithms/polygon_clipping/liang_barsky.jl")
include("algorithms/polygon_clipping/sutherland_hodgman.jl")
include("algorithms/pole_of_inaccessibility.jl")
include("algorithms/triangulation/constrained_triangulation.jl")
include("algorithms/triangulation/check_args.jl")
include("algorithms/triangulation/main.jl")
include("algorithms/triangulation/triangulate_curve_bounded.jl")
include("algorithms/triangulation/mesh_refinement.jl")
include("algorithms/triangulation/triangulate_convex.jl")
include("algorithms/triangulation/triangulate_rectangle.jl")
include("algorithms/triangulation/unconstrained_triangulation.jl")
include("algorithms/triangulation/basic_operations/add_boundary_information.jl")
include("algorithms/triangulation/basic_operations/add_ghost_triangles.jl")
include("algorithms/triangulation/basic_operations/add_point.jl")
include("algorithms/triangulation/basic_operations/add_segment.jl")
include("algorithms/triangulation/basic_operations/add_triangle.jl")
include("algorithms/triangulation/basic_operations/clear_empty_features.jl")
include("algorithms/triangulation/basic_operations/delete_ghost_triangles.jl")
include("algorithms/triangulation/basic_operations/delete_holes.jl")
include("algorithms/triangulation/basic_operations/delete_point.jl")
include("algorithms/triangulation/basic_operations/delete_triangle.jl")
include("algorithms/triangulation/basic_operations/flip_edge.jl")
include("algorithms/triangulation/basic_operations/legalise_edge.jl")
include("algorithms/triangulation/basic_operations/lock_convex_hull.jl")
include("algorithms/triangulation/basic_operations/split_edge.jl")
include("algorithms/triangulation/basic_operations/split_triangle.jl")
include("algorithms/triangulation/basic_operations/unlock_convex_hull.jl")
include("algorithms/voronoi/centroidal.jl")
include("algorithms/voronoi/clipped_coordinates.jl")
include("algorithms/voronoi/clipped.jl")
include("algorithms/voronoi/main.jl")
include("algorithms/voronoi/unbounded.jl")

include("validation.jl")

export
    each_triangle,
    each_solid_triangle,
    each_ghost_triangle,
    each_edge,
    each_solid_edge,
    each_ghost_edge,
    each_vertex,
    each_solid_vertex,
    each_ghost_vertex,
    num_triangles,
    num_solid_triangles,
    num_ghost_triangles,
    num_edges,
    num_solid_edges,
    num_ghost_edges,
    num_vertices,
    num_solid_vertices,
    num_ghost_vertices,
    get_points,
    get_boundary_nodes,
    get_interior_segments,
    get_all_segments,
    get_adjacent,
    get_adjacent2vertex,
    get_ghost_vertex_map,
    get_ghost_vertex_ranges,
    get_ghost_vertex_range,
    get_triangles,
    get_boundary_edge_map,
    statistics,
    each_segment,
    get_area,
    get_all_stat,
    num_polygons,
    get_polygon,
    each_polygon,
    get_polygon_point,
    each_generator,
    get_generator,
    each_polygon_index,
    each_polygon_vertex,
    num_polygon_vertices,
    get_neighbours,
    get_convex_hull,
    get_convex_hull_vertices,
    triangle_vertices,
    edge_vertices,
    getx,
    gety,
    getxy,
    get_point,
    num_boundary_edges,
    add_ghost_triangles!,
    delete_ghost_triangles!,
    add_point!,
    add_triangle!,
    delete_triangle!,
    flip_edge!,
    split_edge!,
    split_triangle!,
    legalise_edge!,
    delete_point!,
    add_segment!,
    lock_convex_hull!,
    unlock_convex_hull!,
    triangulate,
    triangulate_rectangle,
    triangulate_convex,
    brute_force_search,
    jump_and_march,
    convert_boundary_points_to_indices,
    refine!,
    centroidal_smooth,
    get_polygon_coordinates,
    get_nearest_neighbour,
    convex_hull,
    get_convex_hull,
    get_graph,
    voronoi,
    retriangulate,
    Triangulation,
    VoronoiTessellation,
    ConvexHull,
    TriangulationStatistics,
    convex_hull!,
    add_weight!,
    get_weight,
    get_weights,
    CircularArc,
    LineSegment,
    EllipticalArc,
    BezierCurve,
    BSpline,
    CatmullRomSpline,
    find_polygon,
    get_boundary_curves,
    map_ghost_vertex,
    polygon_bounds,
    num_neighbours

end