"""
    split_triangle!(tri::Triangulation, i, j, k, r) 

Given a triangulation `tri`, a triangle `(i, j, k)`, and a 
point `r` inside the triangle, splits the triangle at `r` 
so that `(i, j, k)` is replaced by the three triangles 
`(i, j, r)`, `(j, k, r)`, and `(k, i, r)`, respectively.
"""
function split_triangle!(tri::Triangulation, i, j, k, r)
    delete_triangle!(tri, i, j, k; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, i, j, r; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, j, k, r; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, k, i, r; protect_boundary=true, update_ghost_edges=false)
    return nothing
end

"""
    legalise_split_triangle!(tri::Triangulation, i, j, k, r)

Given a triangulation `tri`, a triangle `(i, j, k)` that has 
already been split by [`split_triangle!`](@ref) at the point `r`,
legalises the new edges using [`legalise_edge!`](@ref).
"""
function legalise_split_triangle!(tri::Triangulation, i, j, k, r)
    legalise_edge!(tri, i, j, r)
    legalise_edge!(tri, j, k, r)
    legalise_edge!(tri, k, i, r)
    return nothing
end

"""
    complete_split_triangle_and_legalise!(tri::Triangulation, i, j, k, r)

Given a triangulation `tri`, a triangle `(i, j, k)`, and a point `r`,
splits `(i, j, k)` at `r` using [`split_triangle!`](@ref) and then legalises
the new edges using [`legalise_split_triangle!`](@ref).
"""
function complete_split_triangle_and_legalise!(tri::Triangulation, i, j, k, r)
    split_triangle!(tri, i, j, k, r)
    legalise_split_triangle!(tri, i, j, k, r)
    return nothing
end
