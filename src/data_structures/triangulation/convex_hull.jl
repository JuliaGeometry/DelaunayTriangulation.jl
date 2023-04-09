"""
    get_convex_hull_indices(tri::Triangulation)

Calls `get_indices(get_convex_hull(tri))`.
"""
@inline get_convex_hull_indices(tri::Triangulation) = get_indices(get_convex_hull(tri))

"""
    convex_hull!(tri::Triangulation; reconstruct=is_constrained(tri))

Computes the convex hull of the points included in the triangulation `tri`. If `reconstruct`, 
then the method in [`convex_hull`](@ref) will be used, whereas if `!reconstruct` then the triangulation's 
boundary is extracted to get the convex hull directly. Note that this latter method fails 
if there are any constrained edges on the boundary of `tri`.
"""
function convex_hull!(tri::Triangulation; reconstruct=is_constrained(tri))
    I = integer_type(tri)
    if reconstruct
        convex_hull!(get_convex_hull(tri))
        return nothing
    elseif has_ghost_triangles(tri)
        outer_boundary_edges = get_adjacent2vertex(tri, I(BoundaryIndex)) # Not exactly all of them, e.g. the outer boundary could be made up of multiple segments, but this will be close
        indices = get_convex_hull_indices(tri)
        empty!(indices)
        sizehint!(indices, num_edges(outer_boundary_edges))
        e = first(each_edge(outer_boundary_edges))
        u = initial(e) # Need to start somewhere 
        push!(indices, u)
        start_index = u
        v = get_right_boundary_node(tri, u, I(BoundaryIndex))
        while v ≠ start_index
            push!(indices, v)
            v = get_right_boundary_node(tri, v, I(BoundaryIndex))
        end
        push!(indices, start_index)
        return nothing
    else
        add_ghost_triangles!(tri)
        convex_hull!(tri; reconstruct=false)
        delete_ghost_triangles!(tri)
    end
    return nothing
end