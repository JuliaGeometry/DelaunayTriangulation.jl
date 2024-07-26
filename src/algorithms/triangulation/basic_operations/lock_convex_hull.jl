"""
    lock_convex_hull!(tri::Triangulation; rng=Random.default_rng(), predicates::AbstractPredicateType=def_alg222())

Locks the convex hull of the unconstrained triangulation `tri` so that it is now treated as a constrained triangulation 
with boundary given by its convex hull.

The random number generator (used inside [`add_segment!`](@ref) can be provided with the `rng` keyword argument,
and similarly for `predicates`.

!!! warning 

    If an edge is encountered along the convex hull that contains a segment from `tri.interior_segments`,
    then this edge will be deleted from `tri.interior_segments`; this will be undone from `unlock_convex_hull!`, 
    possibly splitting the segments in case they were split before unlocking.
"""
function lock_convex_hull!(tri::Triangulation; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateType=def_alg222())
    if has_boundary_nodes(tri)
        throw(ArgumentError("Cannot lock the convex hull of a triangulation with boundary nodes."))
    end
    I = integer_type(tri)
    idx = get_convex_hull_vertices(tri)
    bn = get_boundary_nodes(tri)
    resize!(bn, num_boundary_edges(idx) + 1)
    copyto!(bn, idx)
    bn_map = get_ghost_vertex_map(tri)
    bnn_map = get_boundary_edge_map(tri)
    bn_map[I(𝒢)] = bn
    ne = num_boundary_edges(bn)
    E = edge_type(tri)
    interior_segments = get_interior_segments(tri)
    interior_segments_on_hull = get_interior_segments_on_hull(get_cache(tri))
    empty!(interior_segments_on_hull)
    for i in 1:ne
        u = get_boundary_nodes(bn, i)
        v = get_boundary_nodes(bn, i + 1)
        e = construct_edge(E, u, v)
        bnn_map[e] = (bn, i)
        if contains_unoriented_edge(e, interior_segments)
            delete_unoriented_edge!(interior_segments, e)
            add_edge!(interior_segments_on_hull, e)
        end
    end
    for e in keys(bnn_map)
        add_segment!(tri, e; rng, predicates)
    end
    return tri
end