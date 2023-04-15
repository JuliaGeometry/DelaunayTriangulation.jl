"""
    delete_ghost_triangles!(tri::Triangulation; boundary_indices = all_boundary_indices(tri), delete_neighbours=false)

Deletes the ghost triangles from the triangulation `tri`.

!!! note 

    A ghost triangle is a triangle of the form `(i, j, k)` where only one of the indices, say `i`, 
    satisfies `is_boundary_index(i)`.
"""
function delete_ghost_triangles!(tri::Triangulation; boundary_indices=all_boundary_indices(tri), delete_neighbours=false)
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
            if delete_neighbours
                delete_neighbour!(tri, u, v, g)
                delete_neighbour!(tri, v, g)
            end
        end
    end
    return nothing
end
