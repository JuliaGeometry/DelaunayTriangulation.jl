"""
    split_triangle!(tri::Triangulation, i, j, k, r) 

Given a triangulation `tri`, a triangle `(i, j, k)`, and a 
point `r` inside the triangle, splits the triangle at `r` 
so that `(i, j, k)` is replaced by the three triangles 
`(i, j, r)`, `(j, k, r)`, and `(k, i, r)`, respectively.
"""
function split_triangle!(tri::Triangulation, i, j, k, r)
    delete_triangle!(tri, i, j, k; protect_boundary=true)
    add_triangle!(tri, i, j, r)
    add_triangle!(tri, j, k, r)
    add_triangle!(tri, k, i, r)
    return nothing
end