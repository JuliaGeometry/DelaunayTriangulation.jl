"""
    lock_convex_hull!(tri::Triangulation)

Locks the convex hull of the triangulation `tri` by adding it to the constrained edges 
of `tri` in place of `boundary_nodes`. If `has_boundary_nodes(tri)` is already true, 
an error is thrown.
"""
function lock_convex_hull!(tri::Triangulation)
    if has_boundary_nodes(tri)
        throw(ArgumentError("Cannot lock the convex hull of a triangulation with boundary nodes."))
    end
    I = integer_type(tri)
    idx = get_convex_hull_indices(tri)
    bn = get_boundary_nodes(tri)
    resize!(bn, num_boundary_edges(idx) + 1)
    copyto!(bn, idx)
    bn_map = get_boundary_map(tri)
    bnn_map = get_boundary_edge_map(tri)
    bn_map[I(BoundaryIndex)] = bn
    ne = num_boundary_edges(bn)
    E = edge_type(tri)
    for i in 1:ne
        u = get_boundary_nodes(bn, i)
        v = get_boundary_nodes(bn, i + 1)
        e = construct_edge(E, u, v)
        bnn_map[e] = (bn, i)
    end
    for e in keys(bnn_map)
        add_edge!(tri, e)
    end
    return nothing
end