"""
    add_ghost_triangles!(tri::Triangulation)

Adds all the ghost triangles to `tri`.
"""
function add_ghost_triangles!(tri::Triangulation)
    T = get_triangles(tri)
    for g in each_ghost_vertex(tri)
        for e in each_edge(get_adjacent2vertex(tri, g))
            u, v = edge_vertices(e)
            add_adjacent!(tri, v, g, u)
            add_adjacent!(tri, g, u, v)
            add_adjacent2vertex!(tri, u, v, g)
            add_adjacent2vertex!(tri, v, g, u)
            # add_triangle!(T, u, v, g)
        end
    end
    return tri
end
