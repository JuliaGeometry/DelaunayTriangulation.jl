"""
    delete_triangle!(tri::Triangulation, u, v, w)
    delete_triangle!(tri::Triangulation, T)

Deletes the triangle `T = (u, v, w)` from `tri`. This won't delete the points from `tri`, it will just update the fields 
so that its non-existence in the triangulation is known.

# Arguments 
- `tri::Triangulation`: The triangulation to delete the triangle from.
- `u, v, w`: The vertices of the triangle to delete.

# Outputs 
There are no outputs as `tri` is updated in-place.
"""
function delete_triangle!(tri::Ts, u::Integer, v::Integer, w::Integer) 
    adj = get_adjacent(tri)
    delete_triangle!(adj, u, v, w)
    return tri
end
function delete_triangle!(tri::Triangulation, T)
    u, v, w = triangle_vertices(T)
    delete_triangle!(tri, u, v, w)
    return tri
end