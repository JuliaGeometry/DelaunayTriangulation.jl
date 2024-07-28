"""
    get_convex_hull_vertices(tri::Triangulation) -> Vector{Vertex}

Returns the vertices on the convex hull of `tri`, in counter-clockwise order. 

See also [`ConvexHull`](@ref).
"""
get_convex_hull_vertices(tri::Triangulation) = (get_vertices ∘ get_convex_hull)(tri)

"""
    convex_hull!(tri::Triangulation; reconstruct=has_boundary_nodes(tri), predicates::AbstractPredicateKernel=AdaptiveKernel())

Updates the `convex_hull` field of `tri` to match the current triangulation. 

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).

# Keyword Arguments
- `reconstruct=has_boundary_nodes(tri)`: If `true`, then the convex hull is reconstructed from scratch, using [`convex_hull`](@ref) on the points. Otherwise,
   computes the convex hull using the ghost triangles of `tri`. If there are no ghost triangles but `reconstruct=true`, then the convex hull is reconstructed from scratch.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
"""
function convex_hull!(tri::Triangulation; reconstruct=has_boundary_nodes(tri), predicates::AbstractPredicateKernel=AdaptiveKernel())
    I = integer_type(tri)
    if reconstruct
        convex_hull!(get_convex_hull(tri); predicates)
        return tri
    elseif has_ghost_triangles(tri)
        outer_boundary_edges = get_adjacent2vertex(tri, I(𝒢)) # Not exactly all of them, e.g. the outer boundary could be made up of multiple sections, but this will be close
        indices = get_convex_hull_vertices(tri)
        empty!(indices)
        sizehint!(indices, num_edges(outer_boundary_edges))
        e = first(each_edge(outer_boundary_edges))
        u = initial(e) # Need to start somewhere 
        push!(indices, u)
        start_index = u
        v = get_right_boundary_node(tri, u, I(𝒢))
        while v ≠ start_index
            push!(indices, v)
            v = get_right_boundary_node(tri, v, I(𝒢))
        end
        push!(indices, start_index)
        return tri
    else
        add_ghost_triangles!(tri)
        convex_hull!(tri; predicates, reconstruct=false)
        delete_ghost_triangles!(tri)
    end
    return tri
end