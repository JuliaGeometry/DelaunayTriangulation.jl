"""
    split_edge!(tri::Triangulation, i, j, r)

Given a triangulation `tri` and an edge `(i, j)`, splits the edge at 
the point `r` so that the edges `(i, r)` and `(r, j)` now appear in `tri` 
(with the triangles updated accordingly). Note that this assumes that 
the points corresponding to `i`, `j`, and `r` are collinear.
"""
function split_edge!(tri::Triangulation, i, j, r)
    k = get_adjacent(tri, i, j)
    delete_triangle!(tri, i, j, k; protect_boundary = !is_boundary_edge(tri, j, i))
    add_triangle!(tri, i, r, k)
    add_triangle!(tri, r, j, k)
    return nothing
end 