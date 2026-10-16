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
    # Use all_ghost_vertices instead of each_ghost_vertex because each_ghost_vertex
    # depends on has_ghost_vertices(graph) which may be false if ghost vertices
    # were removed from the graph by clear_empty_vertices!.
    adj2v_dict = get_adjacent2vertex(get_adjacent2vertex(tri))
    for g in all_ghost_vertices(tri)
        # Skip ghost vertices that don't have any edges in adjacent2vertex.
        # This can happen when ghost_vertex_ranges has a default entry but no
        # actual boundary data was set up.
        haskey(adj2v_dict, g) || continue
        for e in each_edge(get_adjacent2vertex(tri, g))
            u, v = edge_vertices(e)
            add_adjacent!(tri, v, g, u)
            add_adjacent!(tri, g, u, v)
            add_adjacent2vertex!(tri, u, v, g)
            add_adjacent2vertex!(tri, v, g, u)
            add_neighbour!(tri, g, u, v)
            add_triangle!(T, u, v, g)
        end
    end
    set_has_ghosts!(tri, true)
    return tri
end
