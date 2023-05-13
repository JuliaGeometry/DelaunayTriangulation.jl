"""
    add_edge!(tri::Triangulation, segment; rng::AbstractRNG=Random.default_rng())
    add_edge!(tri::Triangulation, i, j; rng::AbstractRNG=Random.default_rng())

Adds the edge `segment = (i, j)` into the triangulation `tri`.
"""
function add_edge!(tri::Triangulation, segment; rng::AbstractRNG=Random.default_rng())
    e = contains_boundary_edge(tri, segment) ? segment : sort_edge_by_degree(tri, segment) # We want to preserve the orientation of the boundary edges as it makes it easier to split them later if we run into collinearities
    constrained_edges = get_constrained_edges(tri)
    if initial(e) ≠ initial(segment) && contains_edge(segment, constrained_edges) # If we switch the edge, let's also change it in constrained_edges
        delete_edge!(constrained_edges, segment)
        add_edge!(constrained_edges, e)
    end
    all_constrained_edges = get_all_constrained_edges(tri)
    if !(contains_edge(e, constrained_edges) || contains_edge(reverse_edge(e), constrained_edges)) && !contains_boundary_edge(tri, segment) # In case it isn't already in the list, let's make sure we add it
        add_edge!(constrained_edges, e)
    end
    if !contains_edge(reverse_edge(e), all_constrained_edges)
        add_edge!(all_constrained_edges, e)
    end
    if edge_exists(tri, e) || edge_exists(tri, reverse_edge(e)) # Don't need to add edges that already appear
        return nothing
    end
    intersecting_triangles, collinear_segments, left_cavity, right_cavity = locate_intersecting_triangles(tri, e, !contains_boundary_edge(tri, segment), rng)
    flag = process_collinear_segments!(all_constrained_edges, constrained_edges, e, collinear_segments, tri; rng)
    flag && return nothing
    delete_intersected_triangles!(tri, intersecting_triangles)
    triangulated_left_cavity = triangulate_cavity_cdt(tri, left_cavity; rng)
    triangulated_right_cavity = triangulate_cavity_cdt(tri, right_cavity; rng)
    add_new_triangles!(tri, triangulated_left_cavity, triangulated_right_cavity)
    return nothing
end
add_edge!(tri::Triangulation, i, j; rng=Random.default_rng()) =
    let E = edge_type(tri)
        add_edge!(tri, construct_edge(E, i, j); rng)
    end

function process_collinear_segments!(all_constrained_edges, constrained_edges, e, collinear_segments, tri::Triangulation; rng=Random.default_rng())
    if !is_empty(collinear_segments)
        delete_edge!(all_constrained_edges, e)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
        split_constrained_edge!(constrained_edges, e, collinear_segments)
        if contains_boundary_edge(tri, e)
            split_boundary_edge_at_collinear_segments!(tri, collinear_segments)
        end
        for η in each_edge(collinear_segments)
            add_edge!(tri, η; rng)
        end
        return true
    end
    return false
end

