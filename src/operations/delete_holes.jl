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

"""
    has_interiors_within_interiors(tri::Triangulation)

Returns `true` if the triangulation has multiple curves and the first curve has a positive area and all other curves have negative areas,
meaning there are some interior curves that are inside other interior curves. Returns `false` otherwise.
"""
function has_interiors_within_interiors(tri::Triangulation)
    if has_multiple_curves(tri)
        nc = num_curves(tri)
        areas = zeros(nc)
        for curve_index in 1:nc
            curve_boundary_nodes = get_boundary_nodes(tri, curve_index)
            areas[curve_index] = polygon_features(get_points(tri), curve_boundary_nodes)[1]
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

"""
    find_all_points_to_delete(tri::Triangulation)

Finds all the points that are outside the boundary specified in the 
triangulation `tri`.
"""
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
    counterclockwise_orientation = polygon_features(get_points(tri), curve_boundary_nodes)[1] > 0 # rather than dealing with the case where we have interiors within interiors separately, let's just compute the orientation at each stage
    for segment_index in 1:ns
        ## This originally didn't iterate over every seed, but it becomes problematic when the boundary is blocked 
        ## into separate pieces by solid triangles. In that case, we need to iterate over every seed to make sure
        ## that we delete all the points that are outside the boundary.
        segment_nodes = get_boundary_nodes(curve_boundary_nodes, segment_index)
        ne = num_boundary_edges(segment_nodes)
        for i in 1:ne
            seed = get_initial_seed(tri, curve_index, segment_index, i)
            !edge_exists(seed) && break
            seed ∉ all_bn && push!(points_to_delete, seed) # Need seed ∉ all_bn in case the seed is not inside, but on the boundary
            find_all_points_to_delete!(points_to_delete, tri, seed, all_bn, curve_boundary_nodes, counterclockwise_orientation)
        end
    end
    return nothing
end
function find_all_points_to_delete_node!(points_to_delete, tri::Triangulation, all_bn)
    ne = length(all_bn)  # not using num_boundary_edges(all_bn) because all_bn is a Set and not a Vector
    for i in 1:ne
        seed = get_initial_seed(tri, nothing, 1, i)
        !edge_exists(seed) && break
        seed ∉ all_bn && push!(points_to_delete, seed) # Need seed ∉ all_bn in case the seed is not inside, but on the boundary
        find_all_points_to_delete!(points_to_delete, tri, seed, all_bn, get_boundary_nodes(tri), true)
    end
    return nothing
end

function find_all_points_to_delete!(points_to_delete, tri::Triangulation, seed, all_bn, curve_boundary_nodes, check_outside=false) # if check_outside=false, then check_inside
    for new_seed in get_neighbours(tri, seed)
        q = get_point(tri, new_seed)
        δ = distance_to_polygon(q, get_points(tri), curve_boundary_nodes) # Need to make sure that we only look inside the boundary
        δ = !check_outside ? δ : -δ
        if new_seed ∈ points_to_delete && δ < zero(δ)
            delete!(points_to_delete, new_seed)
        elseif new_seed ∉ points_to_delete && new_seed ∉ all_bn && δ > zero(δ)
            push!(points_to_delete, new_seed)
            find_all_points_to_delete!(points_to_delete, tri, new_seed, all_bn, curve_boundary_nodes, check_outside)
        end
    end
    return nothing
end

"""
    find_all_triangles_to_delete(tri::Triangulation, points_to_delete)

Given a set of `points_to_delete`, returns all triangles that have it as a vertex.
If `tri` has interiors within interiors, then point-in-polygon location has to be 
performed to verify that each triangle should actually be deleted.
"""
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
    if is_true(has_interiors_within_interiors(tri))
        complete_boundary = get_boundary_nodes(tri)
        # Things can go wrong if we have interiors within interiors. Let's just 
        # double check each triangle. 
        for V in each_triangle(triangles)
            u, v, w = indices(V)
            p, q, r = get_point(tri, u, v, w)
            m = ((getx(p) + getx(q) + getx(r)) / 3, (gety(p) + gety(q) + gety(r)) / 3)
            δ = distance_to_polygon(m, get_points(tri), complete_boundary)
            if δ > zero(δ)
                delete!(triangles, V)
            end
        end
    end
    return triangles
end

"""
    delete_all_exterior_triangles(tri::Triangulation, triangles)

Given `triangles` in the exterior face of `tri`, deletes them from `tri`.
"""
function delete_all_exterior_triangles(tri::Triangulation, triangles)
    for T in each_triangle(triangles)
        delete_triangle!(tri, T; protect_boundary=true)
    end
    return nothing
end

"""
    clear_deleted_points!(tri::Triangulation, points_to_delete)

Clears the empty neighbours and sets left behind when deleting the points 
in `points_to_delete`.
"""
function clear_deleted_points!(tri::Triangulation, points_to_delete)
    for i in points_to_delete
        delete_adjacent2vertex!(tri, i)
        delete_vertex!(tri, i)
    end
    return nothing
end

"""
    delete_remaining_triangles_connecting_boundary_edges!(tri::Triangulation)

Deletes triangles from `tri` that have a boundary edge for an edge and lie in the 
exterior face of `tri`. This covers all triangles that were not found by 
`find_all_points_to_delete!` as these triangles might not have had any points in the 
exterior as faces, and these edges could cross over the exteror. Point-in-polygon 
location may be used for the exterior faces but not the interior faces. If there 
are interiors within interiors, then point-in-polygon location is used for all
exterior faces.
"""
function delete_remaining_triangles_connecting_boundary_edges!(tri::Triangulation)
    # There is almost certainly a better way to do this. We are looking to delete 
    # any triangles that have edges that join two boundary nodes, since those won't be
    # detected in the work above. Our approach for this is to just look over all 
    # adjacent2vertex results. We need to do this per curve rather than checking with 
    # all_bn, since there could be an edge connecting two separate boundaries. 
    # In another case where there are interiors within interiors, there can be 
    # boundary nodes connecting the interior its within. So, once we are done, we scan 
    # for this case, using distance_to_polygon to perform point-in-polygon location.
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
    if is_true(has_interiors_within_interiors(tri))
        complete_boundary = get_boundary_nodes(tri)
        # Start by validating each triangle in triangles_to_delete. 
        for V in each_triangle(triangles_to_delete)
            u, v, w = indices(V)
            p, q, r = get_point(tri, u, v, w)
            m = ((getx(p) + getx(q) + getx(r)) / 3, (gety(p) + gety(q) + gety(r)) / 3)
            δ = distance_to_polygon(m, get_points(tri), complete_boundary)
            if δ > zero(δ)
                delete!(triangles_to_delete, V)
            end
        end
        # Now sure for any missed triangles
        all_bn = get_all_boundary_nodes(tri)
        for i in all_bn
            S = get_adjacent2vertex(tri, i)
            for e in each_edge(S)
                u, v = edge_indices(e)
                V = construct_triangle(T, u, v, i)
                if !contains_triangle(V, triangles_to_delete)[2]
                    if all(∈(all_bn), (u, v))
                        p, q, r = get_point(tri, u, v, i)
                        m = ((getx(p) + getx(q) + getx(r)) / 3, (gety(p) + gety(q) + gety(r)) / 3)
                        δ = distance_to_polygon(m, get_points(tri), complete_boundary)
                        if δ < zero(δ)
                            push!(triangles_to_delete, construct_triangle(T, u, v, i))
                        end
                    end
                end
            end
        end
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
