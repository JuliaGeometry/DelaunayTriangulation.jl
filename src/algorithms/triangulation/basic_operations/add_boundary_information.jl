"""
    add_boundary_information!(tri::Triangulation)

Updates `tri` so that the ghost triangle information defined by the boundary nodes in `tri` is added to the triangulation. 
"""
function add_boundary_information!(tri::Triangulation)
    I = integer_type(tri)
    ghost_vertex= I(𝒢)
    bn = get_boundary_nodes(tri)
    if has_multiple_curves(tri)
        add_boundary_curve_information!(tri, bn, ghost_vertex)
    elseif has_multiple_sections(tri)
        add_boundary_segment_information!(tri, bn, ghost_vertex)
    else
        add_boundary_node_information!(tri, bn, ghost_vertex)
    end
    return tri
end
function add_boundary_node_information!(tri::Triangulation, bn, ghost_vertex)
    n_edge = num_boundary_edges(bn)
    u = get_boundary_nodes(bn, n_edge + 1)
    for j in n_edge:-1:1
        v = get_boundary_nodes(bn, j)
        add_adjacent!(tri, u, v, ghost_vertex)
        add_adjacent2vertex!(tri, ghost_vertex, u, v)
        add_neighbour!(tri, ghost_vertex, u, v)
        u = v
    end
    return tri
end
function add_boundary_segment_information!(tri::Triangulation, bn, ghost_vertex)
    for n in 1:num_sections(bn)
        bn_n = get_boundary_nodes(bn, n)
        add_boundary_node_information!(tri, bn_n, ghost_vertex)
        ghost_vertex -= 1
    end
    return ghost_vertex
end
function add_boundary_curve_information!(tri::Triangulation, bn, ghost_vertex)
    for m in 1:num_curves(bn)
        bn_m = get_boundary_nodes(bn, m)
        ghost_vertex = add_boundary_segment_information!(tri, bn_m, ghost_vertex)
    end
    return tri
end
