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
    brute_force_search(tri::Triangulation, q; itr = each_triangle(tri))

Searches for the triangle containing the point `q` by brute force. An exception will be 
raised if no triangle contains the point.

See also [`find_triangle`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The point to be located.

# Keyword Arguments 
- `itr = each_triangle(tri)`: The iterator over the triangles of the triangulation.

# Output 
- `V`: The triangle containing the point `q`.
"""
function brute_force_search(tri::Triangulation, q; itr = each_triangle(tri))
    for V in itr
        cert = point_position_relative_to_triangle(tri, V, q)
        !is_outside(cert) && return V
    end
    return throw(PointNotFoundError(tri, q))
end


"""
    brute_force_search_enclosing_circumcircle(tri::Triangulation, i) -> Triangle 

Searches for a triangle in `tri` containing the vertex `i` in its circumcircle using brute force. If 
`tri` is a weighted Delaunay triangulation, the triangle returned instead has the lifted vertex `i` 
below its witness plane. If no such triangle exists, `($∅, $∅, $∅)` is returned.
"""
function brute_force_search_enclosing_circumcircle(tri::Triangulation, i)
    for V in each_triangle(tri)
        cert = point_position_relative_to_circumcircle(tri, V, i)
        !is_outside(cert) && return V
    end
    tri_type = triangle_type(tri)
    V = construct_triangle(tri_type, ∅, ∅, ∅)
    return V
end