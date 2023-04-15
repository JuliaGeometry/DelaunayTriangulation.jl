## I'm really not convinced that I need to be relying so heavily on 
## point-in-polygon querying in these functions. There should be something 
## much smarter and faster that I could do - not pleased with any of this at all.
## There's some repeated code used here also, not going to make it more 
## confusing until I figure out the smarter way.

"""
    delete_holes!(tri::Triangulation)

Deletes all the exterior faces to the boundary nodes specified in the 
triangulation `tri`.
"""
function delete_holes!(tri::Triangulation)
    points_to_process = find_all_points_to_delete(tri)
    triangles_to_delete = find_all_triangles_to_delete(tri, points_to_process)
    delete_all_exterior_triangles!(tri, triangles_to_delete)
    return nothing
end

"""
    has_interiors_within_interiors(tri::Triangulation)

Returns `true` if the triangulation has multiple curves and the first curve has a positive area and all other curves have negative areas,
meaning there are some interior curves that are inside other interior curves. Returns `false` otherwise.
"""
function has_interiors_within_interiors(tri::Triangulation)
    points = get_points(tri)
    boundary_nodes = get_boundary_nodes(tri)
    return has_interiors_within_interiors(points, boundary_nodes)
end
function has_interiors_within_interiors(points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        nc = num_curves(boundary_nodes)
        areas = zeros(nc)
        for curve_index in 1:nc
            curve_boundary_nodes = get_boundary_nodes(boundary_nodes, curve_index)
            areas[curve_index] = polygon_features(points, curve_boundary_nodes)[1]
        end
        if areas[1] ≥ 0 && all(<(0), @views areas[2:end])
            return Val(false)
        else
            return Val(true)
        end
    else
        return Val(false)
    end
end

"""
    delete_all_exterior_triangles!(tri::Triangulation, triangles)

Given `triangles` in the exterior face of `tri`, deletes them from `tri`.
"""
function delete_all_exterior_triangles!(tri::Triangulation, triangles)
    for T in each_triangle(triangles)
        delete_triangle!(tri, T; protect_boundary=true)
    end
    return nothing
end

function find_all_points_to_delete(tri::Triangulation)
    ## Find all points that might qualify for deletion
    points_to_process = Set{integer_type(tri)}()
    all_bn = get_all_boundary_nodes(tri)
    # Now spread the seed
    if has_multiple_curves(tri)
        nc = num_curves(tri)
        for curve_index in 1:nc
            curve_boundary_nodes = get_boundary_nodes(tri, curve_index)
            ns = num_segments(curve_boundary_nodes)
            for segment_index in 1:ns
                segment_boundary_nodes = get_boundary_nodes(tri, curve_index, segment_index)
                ne = num_boundary_edges(segment_boundary_nodes)
                for edge_index in 1:ne
                    u = get_boundary_nodes(segment_boundary_nodes, edge_index)
                    v = get_boundary_nodes(segment_boundary_nodes, edge_index + 1)
                    seed = get_adjacent(tri, v, u)
                    !edge_exists(seed) && break
                    seed ∉ all_bn && push!(points_to_process, seed) # boundary nodes are processed separately 
                    find_all_points_to_delete!(points_to_process, tri, seed, all_bn)
                end
            end
        end
    elseif has_multiple_segments(tri)
        ns = num_segments(tri)
        for segment_index in 1:ns
            segment_boundary_nodes = get_boundary_nodes(tri, segment_index)
            ne = num_boundary_edges(segment_boundary_nodes)
            for edge_index in 1:ne
                u = get_boundary_nodes(segment_boundary_nodes, edge_index)
                v = get_boundary_nodes(segment_boundary_nodes, edge_index + 1)
                seed = get_adjacent(tri, v, u)
                !edge_exists(seed) && break
                seed ∉ all_bn && push!(points_to_process, seed)
                find_all_points_to_delete!(points_to_process, tri, seed, all_bn)
            end
        end
    else
        boundary_nodes = get_boundary_nodes(tri)
        ne = num_boundary_edges(boundary_nodes)
        for edge_index in 1:ne
            u = get_boundary_nodes(boundary_nodes, edge_index)
            v = get_boundary_nodes(boundary_nodes, edge_index + 1)
            seed = get_adjacent(tri, v, u)
            !edge_exists(seed) && break
            seed ∉ all_bn && push!(points_to_process, seed)
            find_all_points_to_delete!(points_to_process, tri, seed, all_bn)
        end
    end
    return points_to_process
end

function find_all_points_to_delete!(points_to_process, tri::Triangulation, seed, all_bn)
    points = get_points(tri)
    complete_boundary_nodes = get_boundary_nodes(tri)
    for new_seed in get_neighbours(tri, seed)
        q = get_point(tri, new_seed)
        δ = distance_to_polygon(q, points, complete_boundary_nodes)
        if new_seed ∈ points_to_process && δ > zero(δ) && !is_boundary_index(new_seed) # want to make sure we look at the ghost triangles
            delete!(points_to_process, new_seed)
        elseif new_seed ∉ points_to_process && new_seed ∉ all_bn && δ < zero(δ)
            push!(points_to_process, new_seed)
            find_all_points_to_delete!(points_to_process, tri, new_seed, all_bn)
        end
    end
    return nothing
end

function find_all_triangles_to_delete(tri::Triangulation, points_to_process)
    ## Process the non-boundary nodes for deletion
    T = triangle_type(tri)
    triangles_to_delete = Set{T}()
    for point in points_to_process
        S = get_adjacent2vertex(tri, point)
        for e in each_edge(S)
            u, v = edge_indices(e)
            V = construct_triangle(T, u, v, point)
            !contains_triangle(V, triangles_to_delete)[2] && add_triangle!(triangles_to_delete, V)
        end
    end
    ## Now process the boundary nodes
    all_bn = get_all_boundary_nodes(tri)
    points = get_points(tri)
    complete_boundary_nodes = get_boundary_nodes(tri)
    for node in all_bn
        S = get_adjacent2vertex(tri, node)
        for e in each_edge(S)
            u, v = edge_indices(e)
            V = construct_triangle(T, u, v, node)
            if !contains_triangle(V, triangles_to_delete)[2]
                p, q, r = get_point(tri, u, v, node)
                c = triangle_centroid(p, q, r)
                δ = distance_to_polygon(c, points, complete_boundary_nodes)
                if δ < zero(δ)
                    add_triangle!(triangles_to_delete, V)
                end
            end
        end
    end
    return triangles_to_delete
end