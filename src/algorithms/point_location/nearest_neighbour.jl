@doc """
    get_nearest_neighbour(tri_or_vor, q; kwargs...)

Get the index of the nearest neighbour of `q` in `tri_or_vor`. 

For power diagrams, distance is measured using [`get_power_distance`](@ref) (with `q` being assigned zero weight).

# Arguments 
- `tri_or_vor`: A [`Triangulation`](@ref) or [`VoronoiTessellation`](@ref).
- `q`: The point to be located.

# Keyword Arguments
- `kwargs...`: Keyword arguments passed to [`find_triangle`](@ref).

# Output 
- `i`: The index of the nearest neighbour. This is a point of the triangulation if `tri_or_vor` is a [`Triangulation`](@ref) or of a generator if `tri_or_vor` is a [`VoronoiTessellation`](@ref).
"""
get_nearest_neighbour
get_nearest_neighbour(vor::VoronoiTessellation, q; kwargs...) = find_triangle(vor, q; kwargs...)
get_nearest_neighbour(tri::Triangulation, q; kwargs...) = jump_to_voronoi_polygon(tri, q; kwargs...)

function find_triangle(vor::VoronoiTessellation, q; kwargs...)
    return jump_to_voronoi_polygon(get_triangulation(vor), q; kwargs...)
end
function jump_to_voronoi_polygon(tri::Triangulation, q; kwargs...)
    V = find_triangle(tri, q; kwargs...)
    qx, qy = getxy(q)
    V = sort_triangle(V)
    i, j, k = triangle_vertices(V)
    a, b = get_point(tri, i, j)
    daq² = is_weighted(tri) ? get_power_distance(tri, i, q) : dist_sqr(a, q)
    dbq² = is_weighted(tri) ? get_power_distance(tri, j, q) : dist_sqr(b, q)
    if !is_ghost_vertex(k)
        c = get_point(tri, k)
        dcq² = is_weighted(tri) ? get_power_distance(tri, k, q) : dist_sqr(c, q)
    else
        dcq² = typemax(number_type(tri))
    end
    min_dist² = min(daq², dbq², dcq²)
    u = min_dist² == daq² ? i :
        min_dist² == dbq² ? j : k
    current_idx = u
    current_dist = min_dist²
    # The code below checks the polygon surrounding u. It's essentially 
    # just the get_surrounding_polygon code, but we don't store S.
    neighbouring_vertices = get_neighbours(tri, u)
    points = get_points(tri)
    v = first(neighbouring_vertices)
    if !is_ghost_vertex(v)
        # current_dist, current_idx = compare_distance(current_dist, current_idx, points, v, qx, qy)
        # Can't use compare_distance it only uses the Euclidean distance
        r = get_point(tri, v)
        pv_dist = is_weighted(tri) ? get_power_distance(tri, v, q) : dist_sqr(r, q)
        if pv_dist < current_dist
            current_dist = pv_dist
            current_idx = v
        end
    end
    k = num_neighbours(tri, u)
    nbnd = count(is_ghost_vertex, neighbouring_vertices)
    if nbnd > 0
        k = k - nbnd + 1
    end
    for i in 2:k
        v = get_adjacent(tri, u, v)
        if !is_ghost_vertex(v)
            # current_dist, current_idx = compare_distance(current_dist, current_idx, points, v, qx, qy)
            # Can't use compare_distance it only uses the Euclidean distance
            r = get_point(tri, v)
            pv_dist = is_weighted(tri) ? get_power_distance(tri, v, q) : dist_sqr(r, q)
            if pv_dist < current_dist
                current_dist = pv_dist
                current_idx = v
            end
        end
    end
    return current_idx
end
