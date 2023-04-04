"""
    add_edge!(tri::Triangulation, e; rng::AbstractRNG=Random.default_rng())

Adds the edge `e` into the triangulation `tri`.
"""
function add_edge!(tri::Triangulation, e; rng::AbstractRNG=Random.default_rng())
    e = sort_edge_by_degree(tri, e)
    constrained_edges = get_constrained_edges(tri)
    all_constrained_edges = get_all_constrained_edges(tri)
    add_edge!(all_constrained_edges, e)
    if edge_exists(tri, e) || edge_exists(tri, reverse_edge(e))
        return nothing 
    end
    intersecting_triangles, collinear_segments, left_cavity, right_cavity = locate_intersecting_triangles(tri, e; rng)
    if !is_empty(collinear_segments)
        delete_edge!(all_constrained_edges, e)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
        split_constrained_edge!(constrained_edges, e, collinear_segments)
        for η in each_edge(collinear_segments)
            add_edge!(tri, η; rng)
        end
        return nothing 
    end
    delete_intersected_triangles!(tri, intersecting_triangles)
    triangulated_left_cavity = triangulate_cavity_cdt(tri, left_cavity; rng)
    triangulated_right_cavity = triangulate_cavity_cdt(tri, right_cavity; rng)
    add_new_triangles!(tri, triangulated_left_cavity, triangulated_right_cavity)
    return nothing 
end

