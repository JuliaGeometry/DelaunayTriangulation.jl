"""
    is_exterior_ghost_vertex(tri::Triangulation, i) -> Bool 
    
Returns `true` if the ghost vertex `i` in `tri` is an exterior ghost vertex, and `false` otherwise. 

See also [`is_ghost_vertex`](@ref) and [`is_interior_ghost_vertex`](@ref).

# Extended help 
An exterior ghost vertex is a ghost vertex corresponding to a curve or section appearing on the exterior boundary.
"""
function is_exterior_ghost_vertex(tri::Triangulation, i)
    !is_ghost_vertex(i) && return false
    curve_index = get_curve_index(tri, i)
    return is_exterior_curve(tri, curve_index)
end


"""
    is_interior_ghost_vertex(tri::Triangulation, i) -> Bool

Returns `true` if the ghost vertex `i` in `tri` is an interior ghost vertex, and `false` otherwise.

See also [`is_ghost_vertex`](@ref) and [`is_exterior_ghost_vertex`](@ref).

# Extended help
An interior ghost vertex is a ghost vertex corresponding to a curve or section appearing on the interior boundary.
"""
is_interior_ghost_vertex(tri::Triangulation, i) = !is_exterior_ghost_vertex(tri, i)


"""
    get_curve_index(tri::Triangulation, ℓ) -> Integer

Returns the curve index corresponding to the ghost vertex `ℓ` in `tri`.
"""
get_curve_index(tri::Triangulation, ℓ) = get_curve_index(get_ghost_vertex_map(tri), ℓ)


"""
    get_section_index(tri::Triangulation, ℓ) -> Integer

Returns the section index corresponding to the ghost vertex `ℓ` in `tri`.
"""
get_section_index(tri::Triangulation, ℓ) = get_section_index(ghost_vertex_map(tri), ℓ)


"""
    map_ghost_vertex(tri::Triangulation, ℓ) -> Vertex

Given a ghost vertex `ℓ` in `tri`, returns the corresponding section in the 
`boundary_nodes` of `tri`. See also [`get_ghost_vertex_map`](@ref).
"""
map_ghost_vertex(tri::Triangulation, ℓ) = get_ghost_vertex_map(tri)[ℓ]