"""
    add_boundary_information!(tri::Triangulation)

Given a [`Triangulation`](@ref) `tri`, adds boundary information into `tri`. In particular, 
the [`Adjacent`](@ref), [`Adjacent2Vertex`](@ref), and [`Graph`](@ref) fields are updated so that e.g. 
boundary edges map to a corresponding boundary index from [`Adjacent`](@ref), boundary indices map 
to sets of boundary edges from [`Adjacent2Vertex`](@ref), and boundary indices map to boundary 
nodes from [`Graph`](@ref). 

No values are returned.
"""
function add_boundary_information!(tri::Triangulation)
    I = integer_type(tri)
    boundary_index = I(BoundaryIndex)
    bn = get_boundary_nodes(tri)
    if has_multiple_curves(tri)
        add_boundary_curve_information!(tri, bn, boundary_index)
    elseif has_multiple_segments(tri)
        add_boundary_segment_information!(tri, bn, boundary_index)
    else
        add_boundary_node_information!(tri, bn, boundary_index)
    end
    return nothing
end
function add_boundary_node_information!(tri::Triangulation, bn, boundary_index)
    n_edge = num_boundary_edges(bn)
    u = get_boundary_nodes(bn, n_edge + 1)
    for j in n_edge:-1:1
        v = get_boundary_nodes(bn, j)
        add_adjacent!(tri, u, v, boundary_index)
        add_adjacent2vertex!(tri, boundary_index, u, v)
        add_neighbour!(tri, boundary_index, u, v)
        u = v
    end
    return nothing
end
function add_boundary_segment_information!(tri::Triangulation, bn, boundary_index)
    for n in 1:num_segments(bn)
        bn_n = get_boundary_nodes(bn, n)
        add_boundary_node_information!(tri, bn_n, boundary_index)
        boundary_index -= 1
    end
    return boundary_index
end
function add_boundary_curve_information!(tri::Triangulation, bn, boundary_index)
    for m in 1:num_curves(bn)
        bn_m = get_boundary_nodes(bn, m)
        boundary_index = add_boundary_segment_information!(tri, bn_m, boundary_index)
    end
    return nothing
end
