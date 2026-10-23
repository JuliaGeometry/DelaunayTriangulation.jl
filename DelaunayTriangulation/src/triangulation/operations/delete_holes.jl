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
    distances = Dict{I, number_type(tri)}() # shared between the seeds, so that each distance is only computed once
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
                    seed ∈ points_to_process && continue # already spread from, e.g. a ghost vertex that is the seed for every edge of its curve
                    seed ∉ all_bn && push!(points_to_process, seed) # boundary nodes are processed separately 
                    find_all_points_to_delete!(points_to_process, tri, seed, all_bn, distances)
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
                seed ∈ points_to_process && continue
                seed ∉ all_bn && push!(points_to_process, seed)
                find_all_points_to_delete!(points_to_process, tri, seed, all_bn, distances)
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
            seed ∈ points_to_process && continue
            seed ∉ all_bn && push!(points_to_process, seed)
            find_all_points_to_delete!(points_to_process, tri, seed, all_bn, distances)
        end
    end
    I(𝒢) ∉ points_to_process && push!(points_to_process, I(𝒢)) # need to make sure we include this to get rid of any lingering ghost triangles
    return points_to_process
end

"""
    find_all_points_to_delete!(points_to_process, tri::Triangulation, seed, all_bn, distances = Dict{integer_type(tri), number_type(tri)}())

Starting at `seed`, finds more points to spread to and mark for deletion.

# Arguments
- `points_to_process`: The current list of points marked for deletion.
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `seed`: The seed vertex to start spreading from.
- `all_bn`: All the boundary nodes in the triangulation, obtained from [`get_all_boundary_nodes`](@ref).
- `distances`: The distances from the boundary of the vertices considered so far, computed from [`dist`](@ref) and updated in-place.

# Outputs
There are no outputs, as `points_to_process` is updated in-place.

# Extended help
This function works by considering the neighbours around the vertex `seed`. For each neighbouring vertex, we designate that as a new seed,
and consider if it needs to be added into `points_to_process` according to its distance from the triangulation computed from [`distance_to_polygon`](@ref).
We then do the same for each new seed. This uses a stack of the seeds left to consider rather than recursion, since the search can go through
every point outside of the domain, which would overflow the stack if there are many of them.

Computing a distance takes time proportional to the number of boundary nodes, so it is only computed for the vertices that could be marked
or unmarked, and at most once for each vertex. Otherwise, for a seed that is a ghost vertex, whose neighbours are all the boundary nodes on its
curve, the distance would be computed for every boundary node each time the seed is considered.
"""
function find_all_points_to_delete!(points_to_process, tri::Triangulation, seed, all_bn, distances = Dict{integer_type(tri), number_type(tri)}())
    seeds = [seed]
    while !isempty(seeds)
        seed = pop!(seeds)
        for new_seed in get_neighbours(tri, seed)
            if new_seed ∈ points_to_process
                is_ghost_vertex(new_seed) && continue # want to make sure we look at the ghost triangles
                δ = get!(() -> dist(tri, get_point(tri, new_seed)), distances, new_seed)
                δ > zero(δ) && delete!(points_to_process, new_seed)
            elseif new_seed ∉ all_bn
                δ = get!(() -> dist(tri, get_point(tri, new_seed)), distances, new_seed)
                if δ < zero(δ)
                    push!(points_to_process, new_seed)
                    push!(seeds, new_seed)
                end
            end
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
   boundary vertices, and it might not be obvious if it is inside or outside of the triangulation. This is decided by [`is_exterior_boundary_triangle`](@ref).
"""
function find_all_triangles_to_delete(tri::Triangulation, points_to_process)
    ## Process the non-boundary nodes for deletion
    T = triangle_type(tri)
    triangles_to_delete = Set{T}()
    for point in points_to_process
        S = get_adjacent2vertex(tri, point)
        for e in each_edge(S)
            u, v = edge_vertices(e)
            V = construct_triangle(T, u, v, point)
            !contains_triangle(V, triangles_to_delete)[2] && add_triangle!(triangles_to_delete, V)
        end
    end
    ## Now process the boundary nodes
    all_bn = get_all_boundary_nodes(tri)
    for node in all_bn
        S = get_adjacent2vertex(tri, node)
        for e in each_edge(S)
            u, v = edge_vertices(e)
            V = construct_triangle(T, u, v, node)
            if !contains_triangle(V, triangles_to_delete)[2] && is_exterior_boundary_triangle(tri, V)
                add_triangle!(triangles_to_delete, V)
            end
        end
    end
    return triangles_to_delete
end

"""
    is_exterior_boundary_triangle(tri::Triangulation, T) -> Bool

Tests if the positively oriented triangle `T`, adjoining a boundary node of `tri`, is outside of the domain of `tri`.

# Extended help
The domain is always to the left of its boundary edges. So, if `T` has a boundary edge `(u, v)` as one of its edges, `T` is inside the domain, and if
it has the reverse `(v, u)` of a boundary edge, `T` is outside (see [`boundary_edge_side`](@ref)). If `T` has no boundary edges, then the triangles that can be
reached from `T` without crossing a boundary edge are either all inside or all outside of the domain, so we search through them for one that has a boundary edge,
or for an edge of the convex hull, beyond which everything is outside of the domain.

This only uses the boundary edges rather than the coordinates, so it also works for thin triangles along the boundary, e.g. between boundary nodes that
are only just inside the convex hull and the edge of the convex hull over them. The centroid of such a triangle can be so close to the boundary that its
computed distance to the boundary is zero or has the wrong sign. The distance from the centroid of `T` to the boundary, computed with [`dist`](@ref), is only used
for ghost triangles, and if the search finds neither a boundary edge nor an edge of the convex hull.
"""
function is_exterior_boundary_triangle(tri::Triangulation, T)
    is_ghost_triangle(T) && return centroid_is_exterior(tri, T)
    side = boundary_edge_side(tri, T)
    !iszero(side) && return side < 0
    V = triangle_type(tri)
    visited = Set{V}((sort_triangle(T),))
    queue = [T]
    while !isempty(queue)
        S = popfirst!(queue) # S has no boundary edges, so we can cross any of its edges
        i, j, k = triangle_vertices(S)
        for (u, v) in ((i, j), (j, k), (k, i))
            w = get_adjacent(tri, v, u)
            edge_exists(w) || continue
            is_ghost_vertex(w) && return true # (u, v) is an edge of the convex hull
            R = construct_triangle(V, v, u, w)
            sort_triangle(R) ∈ visited && continue
            side = boundary_edge_side(tri, R)
            !iszero(side) && return side < 0
            push!(visited, sort_triangle(R))
            push!(queue, R)
        end
    end
    return centroid_is_exterior(tri, T)
end

"""
    boundary_edge_side(tri::Triangulation, T) -> Int

Returns `1` if the positively oriented triangle `T` has a boundary edge of `tri`, so that `T` is to the left of that edge and hence inside the domain,
`-1` if `T` has the reverse of a boundary edge, so that `T` is outside of the domain, and `0` if `T` has neither.
"""
function boundary_edge_side(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    for (u, v) in ((i, j), (j, k), (k, i))
        contains_boundary_edge(tri, v, u) && return -1
        contains_boundary_edge(tri, u, v) && return 1
    end
    return 0
end

"""
    centroid_is_exterior(tri::Triangulation, T) -> Bool

Tests if the centroid of the triangle `T` is outside of the domain of `tri`, using its distance to the boundary computed with [`dist`](@ref).
"""
function centroid_is_exterior(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    c = triangle_centroid(p, q, r)
    δ = dist(tri, c)
    return δ < zero(δ)
end
