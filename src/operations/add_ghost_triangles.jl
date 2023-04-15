"""
    add_ghost_triangles!(tri::Triangulation; boundary_indices=all_boundary_indices(tri), add_neighbours=false

Given a [`Triangulation`](@ref) `tri`, adds the ghost triangles into `tri`. In particular, 
the `triangles`, [`Adjacent`](@ref), and [`Adjacent2Vertex`](@ref) fields are updated so that 
ghost triangles are stored in them.

!!! note 

    A ghost triangle is a triangle of the form `(i, j, k)` where only one of the indices, say `i`, 
    satisfies `is_boundary_index(i)`.
"""
function add_ghost_triangles!(tri::Triangulation; boundary_indices=all_boundary_indices(tri), add_neighbours=false)
    T = get_triangles(tri)
    for g in boundary_indices
        for e in each_edge(get_adjacent2vertex(tri, g))
            u, v = edge_indices(e)
            add_adjacent!(tri, v, g, u)
            add_adjacent!(tri, g, u, v)
            add_adjacent2vertex!(tri, u, v, g)
            add_adjacent2vertex!(tri, v, g, u)
            add_triangle!(T, u, v, g)
            if add_neighbours
                add_neighbour!(tri, u, v, g)
                add_neighbour!(tri, v, g)
            end
        end
    end
    return nothing
end
