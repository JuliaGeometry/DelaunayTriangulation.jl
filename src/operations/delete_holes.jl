function get_initial_seed(tri::Triangulation, curve_index, segment_index=1, init=1)
    if has_multiple_curves(tri)
        nodes = get_boundary_nodes(tri, (curve_index, segment_index))
        u = get_boundary_nodes(nodes, init)
        v = get_boundary_nodes(nodes, init + 1)
    elseif has_multiple_segments(tri)
        nodes = get_boundary_nodes(tri, segment_index)
        u = get_boundary_nodes(nodes, init)
        v = get_boundary_nodes(nodes, init + 1)
    else
        u = get_boundary_nodes(tri, init)
        v = get_boundary_nodes(tri, init + 1)
    end
    return get_adjacent(tri, v, u)
end

function find_all_points_to_delete(tri::Triangulation)
    points_to_delete = Set{integer_type(tri)}()
    all_bn = get_all_boundary_nodes(tri)
    if has_multiple_curves(tri)
        find_all_points_to_delete_curve!(points_to_delete, tri, all_bn)
    elseif has_multiple_segments(tri)
        find_all_points_to_delete_segment!(points_to_delete, tri, all_bn, get_boundary_nodes(tri), 1)
    else
        find_all_points_to_delete_node!(points_to_delete, tri, all_bn)
    end
    filter!(∉(all_bn), points_to_delete)
    return points_to_delete
end
function find_all_points_to_delete_curve!(points_to_delete, tri::Triangulation, all_bn)
    nc = num_curves(tri)
    for curve_index in 1:nc
        curve_boundary_nodes = get_boundary_nodes(tri, curve_index)
        find_all_points_to_delete_segment!(points_to_delete, tri, all_bn, curve_boundary_nodes, curve_index)
    end
    return nothing
end
function find_all_points_to_delete_segment!(points_to_delete, tri::Triangulation, all_bn, curve_boundary_nodes, curve_index)
    ns = num_segments(curve_boundary_nodes)
    for segment_index in 1:ns
        seed = get_initial_seed(tri, curve_index, segment_index)
        seed ∉ all_bn && push!(points_to_delete, seed) # Need seed ∉ all_bn in case the seed is not inside, but on the boundary
        find_all_points_to_delete!(points_to_delete, tri, seed, all_bn, curve_boundary_nodes, curve_index == 1)
    end
    return nothing
end
function find_all_points_to_delete_node!(points_to_delete, tri::Triangulation, all_bn)
    seed = get_initial_seed(tri, nothing)
    seed ∉ all_bn && push!(points_to_delete, seed) # Need seed ∉ all_bn in case the seed is not inside, but on the boundary
    find_all_points_to_delete!(points_to_delete, tri, seed, all_bn, get_boundary_nodes(tri), true)
    return nothing
end

function find_all_points_to_delete!(points_to_delete, tri::Triangulation, seed, all_bn, curve_boundary_nodes, check_outside=false) # if check_outside=false, then check_inside
    for new_seed in get_neighbours(tri, seed)
        q = get_point(tri, new_seed)
        δ = distance_to_polygon(q, get_points(tri), curve_boundary_nodes) # Need to make sure that we only look inside the boundary
        δ = !check_outside ? δ : -δ
        if new_seed ∉ points_to_delete && new_seed ∉ all_bn && δ > zero(δ)
            push!(points_to_delete, new_seed)
            find_all_points_to_delete!(points_to_delete, tri, new_seed, all_bn, curve_boundary_nodes, check_outside)
        end
    end
    return nothing
end

function find_all_triangles_to_delete(tri::Triangulation, points_to_delete)
    T = triangle_type(tri)
    triangles = Set{T}()
    for point_to_delete in points_to_delete
        for e in each_edge(get_adjacent2vertex(tri, point_to_delete))
            u, v = edge_indices(e)
            V = construct_triangle(T, u, v, point_to_delete)
            if !contains_triangle(V, triangles)[2]
                push!(triangles, V)
            end
        end
    end
    return triangles
end

function delete_all_exterior_triangles(tri::Triangulation, triangles)
    for T in each_triangle(triangles)
        delete_triangle!(tri, T; protect_boundary=true)
    end
    return nothing
end

function clear_deleted_points!(tri::Triangulation, points_to_delete)
    for i in points_to_delete
        delete_adjacent2vertex!(tri, i)
        delete_vertex!(tri, i)
    end
    return nothing
end

function delete_remaining_triangles_connecting_boundary_edges!(tri::Triangulation)
    # There is almost certainly a better way to do this. We are looking to delete 
    # any triangles that have edges that join two boundary nodes, since those won't be
    # detected in the work above. Our approach for this is to just look over all 
    # adjacent2vertex results. We need to do this per curve rather than checking with 
    # all_bn, since there could be an edge connecting two separate boundaries. 
    T = triangle_type(tri)
    triangles_to_delete = Set{T}()
    if has_multiple_curves(tri)
        get_triangles_connecting_boundary_edges_curve!(triangles_to_delete, tri)
    elseif has_multiple_segments(tri)
        curve_boundary_nodes = get_boundary_nodes(tri)
        get_triangles_connecting_boundary_edges_segment!(triangles_to_delete, tri, true, curve_boundary_nodes)
    else
        get_triangles_connecting_boundary_edges_node!(triangles_to_delete, tri)
    end
    for V in each_triangle(triangles_to_delete)
        delete_triangle!(tri, V; protect_boundary=true)
    end
    return nothing
end
function get_triangles_connecting_boundary_edges_curve!(triangles, tri::Triangulation)
    nc = num_curves(tri)
    for curve_index in 1:nc
        curve_boundary_nodes = get_boundary_nodes(tri, curve_index)
        get_triangles_connecting_boundary_edges_segment!(triangles, tri, curve_index == 1, curve_boundary_nodes)
    end
    return nothing
end
function get_triangles_connecting_boundary_edges_segment!(triangles, tri::Triangulation, check_inside, curve_boundary_nodes)
    I = integer_type(tri)
    flattened_nodes = Set{I}()
    ns = num_segments(curve_boundary_nodes)
    for segment_index in 1:ns
        segment_boundary_nodes = get_boundary_nodes(curve_boundary_nodes, segment_index)
        ne = num_boundary_edges(segment_boundary_nodes)
        for node_index in 1:ne
            node = get_boundary_nodes(segment_boundary_nodes, node_index)
            push!(flattened_nodes, node)
        end
    end
    search_flattened_nodes!(triangles, tri, flattened_nodes, check_inside, curve_boundary_nodes)
    return nothing
end
function get_triangles_connecting_boundary_edges_node!(triangles, tri::Triangulation)
    I = integer_type(tri)
    flattened_nodes = Set{I}()
    boundary_nodes = get_boundary_nodes(tri)
    for node in each_boundary_node(boundary_nodes)
        push!(flattened_nodes, node)
    end
    search_flattened_nodes!(triangles, tri, flattened_nodes, true, boundary_nodes)
end
function search_flattened_nodes!(triangles, tri::Triangulation, flattened_nodes, check_inside=false, outer_boundary_nodes=nothing)
    T = triangle_type(tri)
    for node in flattened_nodes
        S = get_adjacent2vertex(tri, node)
        for e in each_edge(S)
            u, v = edge_indices(e)
            if u ∈ flattened_nodes && v ∈ flattened_nodes && !contains_constrained_edge(tri, u, v)
                V = construct_triangle(T, u, v, node)
                if !check_inside
                    push!(triangles, V)
                else
                    p, q = get_point(tri, u, v)
                    m = ((getx(p) + getx(q)) / 2, (gety(p) + gety(q)) / 2)
                    δ = distance_to_polygon(m, get_points(tri), outer_boundary_nodes)
                    δ < 0 && push!(triangles, V)
                end
            end
        end
    end
end

"""
    delete_holes!(tri::Triangulation)

Deletes all the exterior faces to the boundary nodes specified in the 
triangulation `tri`.
"""
function delete_holes!(tri::Triangulation)
    points_to_delete = find_all_points_to_delete(tri)
    triangles = find_all_triangles_to_delete(tri, points_to_delete)
    delete_all_exterior_triangles(tri, triangles)
    delete_remaining_triangles_connecting_boundary_edges!(tri)
    clear_deleted_points!(tri, points_to_delete)
    return nothing
end