"""
    add_ghost_triangles!(tri::Triangulation)

Given a [`Triangulation`](@ref) `tri`, adds the ghost triangles into `tri`. In particular, 
the `triangles, [`Adjacent`](@ref), and [`Adjacent2Vertex`](@ref) fields are updated so that 
ghost triangles are stored in them.

A ghost triangle is a triangle of the form `(i, j, k)` where only one of the indices, say `i`, 
satisfies `is_boundary_index(i)`.

No values are returned.
"""
function add_ghost_triangles!(tri::Triangulation)
    boundary_indices = all_boundary_indices(tri)
    T = get_triangles(tri)
    for g in boundary_indices
        for e in each_edge(get_adjacent2vertex(tri, g))
            u, v = edge_indices(e)
            add_adjacent!(tri, v, g, u)
            add_adjacent!(tri, g, u, v)
            add_adjacent2vertex!(tri, u, v, g)
            add_adjacent2vertex!(tri, v, g, u)
            add_triangle!(T, u, v, g)
        end
    end
    return nothing
end
