"""
    unlock_convex_hull!(tri::Triangulation; reconstruct=false)

Unlocks the convex hull of the constrained triangulation `tri` so that it is now treated as an unconstrained triangulation,
assuming that it was locked using [`lock_convex_hull!`](@ref). If `reconstruct = true`, then the 
convex hull of `tri` will be reconstructed from the boundary nodes of `tri`. This is useful if, for example, 
you have split some of the boundary edges during mesh refinement.
"""
function unlock_convex_hull!(tri::Triangulation; reconstruct = false)
    if !has_boundary_nodes(tri)
        throw(ArgumentError("Cannot unlock the convex hull of a triangulation without boundary nodes."))
    end
    if !reconstruct && !circular_equality(get_boundary_nodes(tri), get_convex_hull_vertices(tri))
        throw(ArgumentError("Cannot unlock the convex hull of a triangulation with boundary nodes that are not the convex hull. If the boundary nodes have been split, consider setting reconstruct=true."))
    end
    bn = get_boundary_nodes(tri)
    if reconstruct
        idx = get_convex_hull_vertices(tri)
        resize!(idx, num_boundary_edges(bn) + 1)
        copyto!(idx, bn)
    end
    I = integer_type(tri)
    E = edge_type(tri)
    bn_map = get_ghost_vertex_map(tri)
    bnn_map = get_boundary_edge_map(tri)
    empty!(bn)
    bn_map[I(𝒢)] = bn
    segments = get_interior_segments(tri)
    all_segments = get_all_segments(tri)
    for e in keys(bnn_map)
        if !contains_unoriented_edge(e, segments)
            delete_unoriented_edge!(segments, e)
            delete_unoriented_edge!(all_segments, e)
        end
    end
    empty!(bnn_map)
    # repair any interior segments that appeared on the convex hull 
    cache = get_cache(tri)
    interior_segments = get_interior_segments(tri)
    interior_segments_on_hull = get_interior_segments_on_hull(cache)
    for e in each_edge(interior_segments_on_hull)
        u, v = edge_vertices(e)
        chain = get_boundary_chain(tri, u, v, I(𝒢))
        ne = length(chain) - 1
        for ℓ in 1:ne
            i, j = chain[ℓ], chain[ℓ + 1]
            e_int = construct_edge(E, i, j)
            add_edge!(all_segments, e_int)
            add_edge!(interior_segments, e_int)
        end
    end
    empty!(interior_segments_on_hull)
    return tri
end

"""
    get_boundary_chain(tri::Triangulation, i, j) -> Edges 

Given two boundary vertices `i` and `j` on a boundary with ghost vertex `ghost_vertex`, 
walks counter-clockwise from `i` to `j` along the boundary and returns the collection of all vertices encountered in 
counter-clockwise order.
"""
function get_boundary_chain(tri::Triangulation, i, j, ghost_vertex)
    I = integer_type(tri)
    chain = Vector{I}()
    push!(chain, i)
    w = I(∅)
    while w ≠ j
        w = get_adjacent(tri, i, ghost_vertex)
        add_edge!(chain, w)
        i = w
    end
    return chain
end