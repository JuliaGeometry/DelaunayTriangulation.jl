"""
    add_triangle!(tri::Triangulation, u, v, w)
    add_triangle!(tri::Triangulation, T)

Adds the triangle `T = (u, v, w)` into `tri`. This won't add the points into `tri`, it will just update the fields 
so that its existence in the triangulation is known.

# Arguments 
- `tri::Triangulation`: The triangulation to add the triangle to.
- `u, v, w`: The vertices of the triangle to add.

# Outputs 
There are no outputs as `tri` is updated in-place.
"""
function add_triangle!(tri::Ts, u::Integer, v::Integer, w::Integer) where {Ts <: Triangulation}
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    add_triangle!(adj, u, v, w)
    add_candidates!(adj2v, u, v, w)
    return tri
end
function add_triangle!(tri::Triangulation, T; protect_boundary = false, update_ghost_edges = false)
    u, v, w = triangle_vertices(T)
    add_triangle!(tri, u, v, w; protect_boundary, update_ghost_edges)
    return tri
end
