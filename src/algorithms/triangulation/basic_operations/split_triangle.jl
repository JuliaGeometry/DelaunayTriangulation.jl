"""
    split_triangle!(tri::Triangulation, i, j, k, r)

Splits the triangle `(i, j, k)` at the vertex `r`, assumed to be inside the triangle. 

See also [`legalise_split_triangle!`](@ref) and [`complete_split_triangle_and_legalise!`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the triangle.
- `j`: The second vertex of the triangle.
- `k`: The third vertex of the triangle.
- `r`: The vertex to split the triangle at.

# Outputs
There is no output, but `tri` will be updated so that it now contains the triangles `(i, j, r)`, `(j, k, r)`, and `(k, i, r)`.
"""
function split_triangle!(tri::Triangulation, i, j, k, r)
    delete_triangle!(tri, i, j, k; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, i, j, r; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, j, k, r; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, k, i, r; protect_boundary=true, update_ghost_edges=false)
    return tri
end

"""
    legalise_split_triangle!(tri::Triangulation, i, j, k, r; predicates::AbstractPredicateKernel=AdaptiveKernel())

Legalises the newly added edges in `tri` after the triangle `(i, j, k)` was split using [`split_triangle!`](@ref).

See also [`complete_split_triangle_and_legalise!`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the triangle.
- `j`: The second vertex of the triangle.
- `k`: The third vertex of the triangle.
- `r`: The vertex to split the triangle at.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, as `tri` is updated in-place.
"""
function legalise_split_triangle!(tri::Triangulation, i, j, k, r; predicates::AbstractPredicateKernel=AdaptiveKernel())
    legalise_edge!(tri, i, j, r; predicates)
    legalise_edge!(tri, j, k, r; predicates)
    legalise_edge!(tri, k, i, r; predicates)
    return tri
end

"""
    complete_split_triangle_and_legalise!(tri::Triangulation, i, j, k, r)

Splits the triangle `(i, j, k)` at the vertex `r`, assumed to be inside the triangle, and legalises the newly added edges in `tri`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the triangle.
- `j`: The second vertex of the triangle.
- `k`: The third vertex of the triangle.
- `r`: The vertex to split the triangle at.

# Outputs
There is no output, as `tri` is updated in-place.
"""
function complete_split_triangle_and_legalise!(tri::Triangulation, i, j, k, r; predicates::AbstractPredicateKernel=AdaptiveKernel())
    split_triangle!(tri, i, j, k, r)
    legalise_split_triangle!(tri, i, j, k, r; predicates)
    return tri
end
