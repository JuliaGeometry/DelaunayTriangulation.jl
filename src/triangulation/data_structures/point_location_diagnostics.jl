mutable struct PointLocationDiagnostics
    interior_vertex_starts::Int
    boundary_vertex_starts::Int
    restarts::Int
    degenerate_arrangements::Int
    edges_sampled::Int
    exterior_entries::Int
    exterior_rotations::Int
    adjacent_boundary_edge_searches::Int
    jumps::Int
end
PointLocationDiagnostics() = PointLocationDiagnostics(0, 0, 0, 0, 0, 0, 0, 0, 0)

@inline add_interior_vertex_start!(pld::PointLocationDiagnostics) = (pld.interior_vertex_starts += 1)
@inline add_boundary_vertex_start!(pld::PointLocationDiagnostics) = (pld.boundary_vertex_starts += 1)
@inline add_restart!(pld::PointLocationDiagnostics) = (pld.restarts += 1)
@inline add_degenerate_arrangement!(pld::PointLocationDiagnostics) = (pld.degenerate_arrangements += 1)
@inline add_edges_sampled!(pld::PointLocationDiagnostics) = (pld.edges_sampled += 1)
@inline add_exterior_entry!(pld::PointLocationDiagnostics) = (pld.exterior_entries += 1)
@inline add_exterior_rotation!(pld::PointLocationDiagnostics) = (pld.exterior_rotations += 1)
@inline add_adjacent_boundary_edge_search!(pld::PointLocationDiagnostics) = (pld.adjacent_boundary_edge_searches += 1)
@inline add_jump!(pld::PointLocationDiagnostics) = (pld.jumps += 1)

@inline num_interior_vertex_starts(pld::PointLocationDiagnostics) = pld.interior_vertex_starts
@inline num_boundary_vertex_starts(pld::PointLocationDiagnostics) = pld.boundary_vertex_starts
@inline num_restarts(pld::PointLocationDiagnostics) = pld.restarts
@inline num_degenerate_arrangements(pld::PointLocationDiagnostics) = pld.degenerate_arrangements
@inline num_edges_sampled(pld::PointLocationDiagnostics) = pld.edges_sampled
@inline num_exterior_entries(pld::PointLocationDiagnostics) = pld.exterior_entries
@inline num_exterior_rotations(pld::PointLocationDiagnostics) = pld.exterior_rotations
@inline num_adjacent_boundary_edge_searches(pld::PointLocationDiagnostics) = pld.adjacent_boundary_edge_searches
@inline num_jumps(pld::PointLocationDiagnostics) = pld.jumps