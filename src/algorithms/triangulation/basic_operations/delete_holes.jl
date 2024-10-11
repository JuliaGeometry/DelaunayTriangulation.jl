## I'm really not convinced that I need to be relying so heavily on 
## point-in-polygon querying in these functions. There should be something 
## much smarter and faster that I could do - not pleased with any of this at all.
## There's some repeated code used here also, not going to make it more 
## confusing until I figure out the smarter way.

"""
    delete_holes!(tri::Triangulation)

Deletes all the exterior faces to the boundary nodes specified in the triangulation `tri`.

# Extended help
This function works in several stages:

1. First, [`find_all_points_to_delete`](@ref) is used to identify all points in the exterior faces.
2. Once all the points to delete have been found, all the associated triangles are found using [`find_all_triangles_to_delete`](@ref), taking care for any incorrectly identified triangles and points. 
3. Once the correct set of triangles to delete has been found, they are deleted using [`delete_all_exterior_triangles!`](@ref).
"""
function delete_holes!(tri::Triangulation)
    points_to_process = find_all_points_to_delete(tri)
    triangles_to_delete = find_all_triangles_to_delete(tri, points_to_process)
    delete_all_exterior_triangles!(tri, triangles_to_delete)
    return tri
end

"""
    delete_all_exterior_triangles!(tri::Triangulation, triangles)

Deletes all the triangles in the set `triangles` from the triangulation `tri`.
"""
function delete_all_exterior_triangles!(tri::Triangulation, triangles)
    for T in each_triangle(triangles)
        delete_triangle!(tri, T; protect_boundary = true, update_ghost_edges = false)
    end
    return tri
end

"""
    find_all_points_to_delete(tri::Triangulation) -> Set{Int}

Returns a set of all the points that are in the exterior faces of the triangulation `tri`.

# Extended help
This function works by 'spreading' from some initial vertex. In particular, starting at each boundary node, we spread outwards towards adjacent vertices, recursively spreading so that all exterior points are identified
with the help of [`find_all_points_to_delete!`](@ref).
"""
function find_all_points_to_delete(tri::Triangulation)
    ## Find all points that might qualify for deletion
    I = integer_type(tri)
    points_to_process = Set{I}()
    all_bn = get_all_boundary_nodes(tri)
    # Now spread the seed
    if has_multiple_curves(tri)
        nc = num_curves(tri)
        for curve_index in 1:nc
            curve_boundary_nodes = get_boundary_nodes(tri, curve_index)
            ns = num_sections(curve_boundary_nodes)
            for section_index in 1:ns
                segment_boundary_nodes = get_boundary_nodes(tri, curve_index, section_index)
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
    elseif has_multiple_sections(tri)
        ns = num_sections(tri)
        for section_index in 1:ns
            segment_boundary_nodes = get_boundary_nodes(tri, section_index)
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
    I(𝒢) ∉ points_to_process && push!(points_to_process, I(𝒢)) # need to make sure we include this to get rid of any lingering ghost triangles
    return points_to_process
end

"""
    find_all_points_to_delete!(points_to_process, tri::Triangulation, seed, all_bn)

Starting at `seed`, finds more points to spread to and mark for deletion. 

# Arguments 
- `points_to_process`: The current list of points marked for deletion. 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `seed`: The seed vertex to start spreading from.
- `all_bn`: All the boundary nodes in the triangulation, obtained from [`get_all_boundary_nodes`](@ref).

# Outputs 
There are no outputs, as `points_to_process` is updated in-place.

# Extended help 
This function works by considering the neighbours around the vertex `seed`. For each neighbouring vertex, we designate that as a new seed, 
and consider if it needs to be added into `points_to_process` according to its distance from the triangulation computed from [`distance_to_polygon`](@ref).
We then call [`find_all_points_to_delete!`](@ref) recursively again on the new seed.
"""
function find_all_points_to_delete!(points_to_process, tri::Triangulation, seed, all_bn)
    for new_seed in get_neighbours(tri, seed)
        q = get_point(tri, new_seed)
        δ = dist(tri, q)
        # δ = distance_to_polygon(q, points, complete_boundary_nodes)
        if new_seed ∈ points_to_process && δ > zero(δ) && !is_ghost_vertex(new_seed) # want to make sure we look at the ghost triangles
            delete!(points_to_process, new_seed)
        elseif new_seed ∉ points_to_process && new_seed ∉ all_bn && δ < zero(δ)
            push!(points_to_process, new_seed)
            find_all_points_to_delete!(points_to_process, tri, new_seed, all_bn)
        end
    end
    return points_to_process
end

"""
    find_all_triangles_to_delete(tri::Triangulation, points_to_process) -> Set{Triangle}

Returns a set of all the triangles that are in the exterior faces of the triangulation `tri`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `points_to_process`: The set of points that are in the exterior faces of the triangulation `tri`, obtained from [`find_all_points_to_delete`](@ref).

# Outputs 
- `triangles_to_delete`: The set of triangles that are in the exterior faces of the triangulation `tri`.

# Extended help
This function works in two stages. 
1. Firstly, all the non-boundary vertices, i.e. those from `points_to_process`, are processed. For each vertex `v`, the triangles adjoining it, given by `get_adjacent2vertex(tri, v)`, aremarked for deletion. 
2. Next, all the boundary vertices need to be processed and carefully analysed to determine if any other triangles need to be deleted since, for example, a triangle may be adjoining three vertices that are all 
   boundary vertices, and it might not be obvious if it is inside or outside of the triangulation. By applying [`dist`](@ref) to compute the distance between the triangle's centroid and the triangulation, 
   the triangle can be accurately marked for deletion if it is outside of the triangulation.
"""
function find_all_triangles_to_delete(tri::Triangulation, points_to_process)
    ## Process the non-boundary nodes for deletion
    I = integer_type(tri)
    triangles_to_delete = Set{NTuple{3,I}}()
    for point in points_to_process
        S = get_adjacent2vertex(tri, point)
        for e in each_edge(S)
            u, v = edge_vertices(e)
            V = (u, v, point)
            !contains_triangle(V, triangles_to_delete)[2] && add_triangle!(triangles_to_delete, V)
        end
    end
    ## Now process the boundary nodes
    all_bn = get_all_boundary_nodes(tri)
    for node in all_bn
        S = get_adjacent2vertex(tri, node)
        for e in each_edge(S)
            u, v = edge_vertices(e)
            V = (u, v, node)
            if !contains_triangle(V, triangles_to_delete)[2]
                p, q, r = get_point(tri, u, v, node)
                c = triangle_centroid(p, q, r)
                δ = dist(tri, c)
                if δ < zero(δ)
                    add_triangle!(triangles_to_delete, V)
                end
            end
        end
    end
    return triangles_to_delete
end
