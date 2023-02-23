"""
    delete_ghost_triangles!(tri::Triangulation)

Given a [`Triangulation`](@ref) `tri`, deletes the ghost triangles from `tri`. In particular, 
the `triangles`, [`Adjacent`](@ref), and [`Adjacent2Vertex`](@ref) fields are updated so that 
ghost triangles are no longer stored in them.

A ghost triangle is a triangle of the form `(i, j, k)` where one of the indices, say `i`, 
satisfies `is_boundary_index(i)`.

No values are returned.
"""
function delete_ghost_triangles!(tri::Triangulation)
    boundary_indices = all_boundary_indices(tri)
    T = get_triangles(tri)
    for g in boundary_indices
        for uv in each_edge(get_adjacent2vertex(tri, g))
            u = initial(uv)
            v = terminal(uv)
            delete_adjacent!(tri, v, g)
            delete_adjacent!(tri, g, u)
            delete_adjacent2vertex!(tri, u, v, g)
            delete_adjacent2vertex!(tri, v, g, u)
            delete_triangle!(T, u, v, g)
        end
    end
    return nothing
end
