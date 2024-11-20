struct Triangulation
    # Geometry
    points::P
    boundary_nodes::BN
    interior_segments::S
    all_segments::S
    weights::W
    triangle_counts::TriangleCounts
    # Topology 
    adjacent::Adjacent{I,E}
    adjacent2vertex_candidates::Adjacent2VertexCandidates{I}
    # Boundary Handling 
    boundary_curves::BC
    ghost_vertex_map::GhostVertexMap{GVM}
    boundary_edge_map::BoundaryEdgeMap{E,BEM}
    ghost_vertex_ranges::GhostVertexRanges{I}
    convex_hull::ConvexHull{P,I}
    representative_point_list::RepresentativeCoordinatesList{T}
    polygon_hierarchy::PolygonHierarchy{I}
    boundary_enricher::BE
    # Diagnostics
    predicate_diagnostics::PD # could be Nothing, meaning no diagnostics
    point_location_diagnostics::PLD # could be Nothing, meaning no diagnostics
    # Other 
    cache::C
end

get_points(tri::Triangulation) = tri.points
get_boundary_nodes(tri::Triangulation) = tri.boundary_nodes
get_interior_segments(tri::Triangulation) = tri.interior_segments
get_all_segments(tri::Triangulation) = tri.all_segments
get_weights(tri::Triangulation) = tri.weights
get_triangle_counts(tri::Triangulation) = tri.triangle_counts
get_adjacent(tri::Triangulation) = tri.adjacent
get_adjacent2vertex_candidates(tri::Triangulation) = tri.adjacent2vertex_candidates
get_boundary_curves(tri::Triangulation) = tri.boundary_curves
get_ghost_vertex_map(tri::Triangulation) = tri.ghost_vertex_map
get_boundary_edge_map(tri::Triangulation) = tri.boundary_edge_map
get_ghost_vertex_ranges(tri::Triangulation) = tri.ghost_vertex_ranges
get_convex_hull(tri::Triangulation) = tri.convex_hull
get_representative_point_list(tri::Triangulation) = tri.representative_point_list
get_polygon_hierarchy(tri::Triangulation) = tri.polygon_hierarchy
get_boundary_enricher(tri::Triangulation) = tri.boundary_enricher
get_predicate_diagnostics(tri::Triangulation) = tri.predicate_diagnostics
get_point_location_diagnostics(tri::Triangulation) = tri.point_location_diagnostics
get_cache(tri::Triangulation) = tri.cache

@inline function validate_incircle_cache(tri::Triangulation, cache)
    isnothing(cache) && return true 
    incircle_cache = get_incircle_cache(tri) 
    return typeof(incircle_cache) === typeof(cache)  
end

@inline function validate_orient3_cache(tri::Triangulation, cache)
    isnothing(cache) && return true 
    orient3_cache = get_orient3_cache(tri) 
    return typeof(orient3_cache) === typeof(cache)
end

@inline fix_incircle_cache(tri::Triangulation, cache) = validate_incircle_cache(tri, cache) ? cache : nothing

@inline fix_orient3_cache(tri::Triangulation, cache) = validate_orient3_cache(tri, cache) ? cache : nothing