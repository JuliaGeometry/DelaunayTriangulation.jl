"""
    add_ghost_triangles!(tri::Triangulation)

Adds all the ghost triangles to `tri`.
"""
function add_ghost_triangles!(tri::Triangulation)
    T = get_triangles(tri)
    # If we have boundary nodes but the boundary_vertex_to_ghost map is empty,
    # we need to populate it. This can happen when ghost triangles were deleted
    # and are now being re-added.
    if has_boundary_nodes(tri) && isempty(get_boundary_vertex_to_ghost(tri))
        add_boundary_information!(tri)
    end
    for g in each_ghost_vertex(tri)
        for e in each_edge(get_adjacent2vertex(tri, g))
            u, v = edge_vertices(e)
            add_adjacent!(tri, v, g, u)
            add_adjacent!(tri, g, u, v)
            add_adjacent2vertex!(tri, u, v, g)
            add_adjacent2vertex!(tri, v, g, u)
            add_triangle!(T, u, v, g)
        end
    end
    set_has_ghosts!(tri, true)
    return tri
end
