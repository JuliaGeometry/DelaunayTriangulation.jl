"""
    unlock_convex_hull!(tri::Triangulation)

Unlocks the convex hull of the triangulation `tri` by removing it from the constrained edges
of `tri`. If `has_boundary_nodes(tri)` is already false or if `get_boundary_nodes(tri) ≠ get_convex_hull_indices(tri)`, 
an error is thrown.
"""
function unlock_convex_hull!(tri::Triangulation)
    if !has_boundary_nodes(tri)
        throw(ArgumentError("Cannot unlock the convex hull of a triangulation without boundary nodes."))
    end
    if !circular_equality(get_boundary_nodes(tri), get_convex_hull_indices(tri))
        throw(ArgumentError("Cannot unlock the convex hull of a triangulation with boundary nodes that are not the convex hull."))
    end
    I = integer_type(tri)
    bn = get_boundary_nodes(tri)
    bn_map = get_boundary_map(tri)
    bnn_map = get_boundary_edge_map(tri)    
    empty!(bn)
    bn_map[I(BoundaryIndex)] = bn 
    edges = get_constrained_edges(tri)
    all_edges  = get_all_constrained_edges(tri)
    for e in keys(bnn_map)
        e′ = reverse_edge(e)
        if !contains_edge(e, edges) && !contains_edge(e′, edges)
            delete_edge!(edges, e)
            delete_edge!(edges, e′)
            delete_edge!(all_edges, e)
            delete_edge!(all_edges, e′)
        end
    end
    empty!(bnn_map)
    return nothing
end