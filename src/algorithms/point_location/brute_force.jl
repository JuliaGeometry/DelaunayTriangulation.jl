struct PointNotFoundError{T, P} <: Exception
    tri::T
    q::P
end
function Base.showerror(io::IO, e::PointNotFoundError)
    err = "The point, $(e.q), could not be located."
    if !has_ghost_triangles(e.tri)
        err *= " This may be due to the point being outside of the triangulation and your triangulation not containing ghost triangles. Consider using `add_ghost_triangles!` on your triangulation and try again."
    end
    return print(io, err)
end

"""
    brute_force_search(tri::Triangulation, q; itr = each_triangle(tri), predicates::AbstractPredicateKernel=AdaptiveKernel())

Searches for the triangle containing the point `q` by brute force. An exception will be 
raised if no triangle contains the point.

See also [`find_triangle`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The point to be located.

# Keyword Arguments 
- `itr = each_triangle(tri)`: The iterator over the triangles of the triangulation.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output 
- `V`: The triangle containing the point `q`.
"""
function brute_force_search(tri::Triangulation, q; itr = each_triangle(tri), predicates::AbstractPredicateKernel = AdaptiveKernel())
    for V in itr
        cert = point_position_relative_to_triangle(predicates, tri, V, q)
        !is_outside(cert) && return V
    end
    return throw(PointNotFoundError(tri, q))
end

"""
    brute_force_search_enclosing_circumcircle(tri::Triangulation, i, predicates::AbstractPredicateKernel=AdaptiveKernel(); cache = nothing) -> Triangle 

Searches for a triangle in `tri` containing the vertex `i` in its circumcircle using brute force. If 
`tri` is a weighted Delaunay triangulation, the triangle returned instead has the lifted vertex `i` 
below its witness plane. If no such triangle exists, `($∅, $∅, $∅)` is returned. You can control 
the method used for computing predicates via the `predicates` argument.

The `cache` argument is passed to [`point_position_relative_to_circumcircle`] and should be one of 
- `nothing`: No cache is used.
- `get_incircle_cache(tri)`: The cache stored inside `tri`.
- `AdaptivePredicates.incircleadapt_cache(number_type(tri))`: Compute a new cache.
The cache is only needed if an `AdaptiveKernel()` is used.
"""
function brute_force_search_enclosing_circumcircle(tri::Triangulation, i, predicates::AbstractPredicateKernel = AdaptiveKernel(); cache::PredicateCacheType = nothing)
    for V in each_triangle(tri)
        cert = point_position_relative_to_circumcircle(predicates, tri, V, i; cache)
        !is_outside(cert) && return V
    end
    tri_type = triangle_type(tri)
    V = construct_triangle(tri_type, ∅, ∅, ∅)
    return V
end
