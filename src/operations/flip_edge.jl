"""
    flip_edge!(tri::Triangulation, i, j)

Given a triangulation `tri` and an edge `(i, j)` appearing in the triangulation, 
flips the edge `(i, j)` so that it becomes `(ℓ, k)`, where `ℓ = get_adjacent(tri, i, j)`
and `k = get_adjacent(tri, j, i)`. Note that if `(i, j, ℓ, k)` is not a convex 
quadrilateral, than this edge flip makes the triangulation non-planar.
"""
function flip_edge!(tri::Triangulation, i, j)
    ℓ = get_adjacent(tri, i, j)
    k = get_adjacent(tri, j, i)
    flip_edge!(tri, i, j, k, ℓ)
    return nothing
end
function flip_edge!(tri::Triangulation, i, j, k, ℓ)
    delete_triangle!(tri, i, k, j; protect_boundary=true)
    delete_triangle!(tri, i, j, ℓ; protect_boundary=true)
    add_triangle!(tri, ℓ, k, j)
    add_triangle!(tri, ℓ, i, k)
    return nothing
end