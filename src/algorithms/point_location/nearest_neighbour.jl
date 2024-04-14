@doc """
    get_nearest_neighbour(tri_or_vor, q; kwargs...)

Get the index of the nearest neighbour of `q` in `tri_or_vor`. 

# Arguments 
- `tri_or_vor`: A [`Triangulation`](@ref) or [`VoronoiTessellation`](@ref).
- `q`: The point to be located.

# Keyword Arguments
- `kwargs...`: Keyword arguments passed to [`jump_and_march`](@ref).

# Output 
- `i`: The index of the nearest neighbour. This is a point of the triangulation if `tri_or_vor` is a [`Triangulation`](@ref) or of a generator if `tri_or_vor` is a [`VoronoiTessellation`](@ref).
"""
get_nearest_neighbour
get_nearest_neighbour(vor::VoronoiTessellation, q; kwargs...) = jump_and_march(vor, q; kwargs...)
get_nearest_neighbour(tri::Triangulation, q; kwargs...) = jump_to_voronoi_polygon(tri, q; kwargs...)

function jump_and_march(vor::VoronoiTessellation, q; kwargs...)
    return jump_to_voronoi_polygon(get_triangulation(vor), q; kwargs...)
end
function jump_to_voronoi_polygon(tri::Triangulation, q; kwargs...)
    V = jump_and_march(tri, q; kwargs...)
    qx, qy = _getxy(q)
    V = sort_triangle(V)
    i, j, k = triangle_vertices(V)
    a, b = get_point(tri, i, j)
    #ax, ay = _getxy(a)
    #bx, by = _getxy(b)
    #daq = (qx - ax)^2 + (qy - ay)^2
    #dbq = (qx - bx)^2 + (qy - by)^2
    daq² = dist_sqr(a, q)
    dbq² = dist_sqr(b, q)
    if !is_ghost_vertex(k)
        c = get_point(tri, k)
        #cx, cy = _getxy(c)
        #dcq = (qx - cx)^2 + (qy - cy)^2
        dcq² = dist_sqr(c, q)
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
        current_dist, current_idx = compare_distance(current_dist, current_idx, points, v, qx, qy)
    end
    k = num_neighbours(tri, u)
    nbnd = count(is_ghost_vertex, neighbouring_vertices)
    if nbnd > 0
        k = k - nbnd + 1
    end
    for i in 2:k
        v = get_adjacent(tri, u, v)
        if !is_ghost_vertex(v)
            current_dist, current_idx = compare_distance(current_dist, current_idx, points, v, qx, qy)
        end
    end
    return current_idx
end