"""
    is_first_ghost_vertex(cell, i) -> Bool 

Assuming that the circular vector `cell` is such that ghost vertices only appear next to eachother in `cell` and there are only two, tests if `i` is the first ghost vertex in `cell`.

See also [`is_last_ghost_vertex`](@ref).

# Arguments 
- `cell`: The circular vector.
- `i`: The index of the vertex in `cell`.

# Outputs
- `flag`: `true` if `i` is the first ghost vertex in `cell`, and `false` otherwise.
"""
function is_first_ghost_vertex(cell, i)
    prev = previndex_circular(cell, i)
    return !is_ghost_vertex(cell[prev])
end

"""
    is_last_ghost_vertex(cell, i) -> Bool

Assuming that the circular vector `cell` is such that ghost vertices only appear next to eachother in `cell`, tests if `i` is the last ghost vertex in `cell`.

See also [`is_first_ghost_vertex`](@ref).

# Arguments
- `cell`: The circular vector.
- `i`: The index of the vertex in `cell`.

# Outputs
- `flag`: `true` if `i` is the last ghost vertex in `cell`, and `false` otherwise.
"""
function is_last_ghost_vertex(cell, i)
    prev = previndex_circular(cell, i)
    return is_ghost_vertex(cell[prev])
end

"""
    get_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box=nothing; predicates::AbstractPredicateType=Adaptive()) -> Vector{NTuple{2,Number}}

Returns the coordinates of the polygon with index `i` in `vorn`. If `bounding_box` is provided, then the polygon is clipped to the bounding box. If the polygon is unbounded, then `bounding_box` must be provided.

See also [`get_unbounded_polygon_coordinates`](@ref) and [`get_bounded_polygon_coordinates`](@ref).

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The index of the polygon.
- `bounding_box=nothing`: The bounding box to clip the polygon to. If `nothing`, then the polygon is not clipped. If the polygon is unbounded, then `bounding_box` must be provided.

# Keyword Arguments 
- `predicates::AbstractPredicateType=Adaptive()`: Method to use for computing predicates. Can be one of [`Fast`](@ref), [`Exact`](@ref), and [`Adaptive`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `coords`: The coordinates of the polygon. This is a circular vector.
"""
function get_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box=nothing; predicates::AbstractPredicateType=Adaptive())
    if !isnothing(bounding_box)
        a, b, c, d = bounding_box
        @assert a < b && c < d "The bounding box must be of the form (xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax."
    end
    if i ∈ get_unbounded_polygons(vorn)
        isnothing(bounding_box) && throw(ArgumentError("The polygon is unbounded, so a bounding box must be provided. See DelaunayTriangulation.polygon_bounds for a reasonable default."))
        return get_unbounded_polygon_coordinates(vorn, i, bounding_box; predicates)
    else
        return get_bounded_polygon_coordinates(vorn, i, bounding_box; predicates)
    end
end

"""
    get_clipping_poly_structs(vorn::VoronoiTessellation, i, bounding_box) -> (Polygon, Polygon)

Returns the polygons used for clipping the `i`th polygon of `vorn` to `bounding_box`.

See also [`clip_polygon`](@ref).

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The index of the polygon.
- `bounding_box`: The bounding box to clip the polygon to.

# Outputs
- `poly`: The polygon to clip.
- `clip_poly`: The polygon to clip to.
"""
function get_clipping_poly_structs(vorn::VoronoiTessellation, i, bounding_box)
    vertices = get_polygon(vorn, i)
    points = get_polygon_points(vorn)
    clip_vertices = (1, 2, 3, 4)
    a, b, c, d = bounding_box
    clip_points = ((a, c), (b, c), (b, d), (a, d))
    return Polygon(vertices, points), Polygon(clip_vertices, clip_points)
end

"""
    clip_bounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive()) -> Vector{NTuple{2,Number}}

Clips the `i`th polygon of `vorn` to `bounding_box`.

See also [`clip_polygon`](@ref).

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The index of the polygon.
- `bounding_box`: The bounding box to clip the polygon to.

# Keyword Arguments 
- `predicates::AbstractPredicateType=Adaptive()`: Method to use for computing predicates. Can be one of [`Fast`](@ref), [`Exact`](@ref), and [`Adaptive`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `coords`: The coordinates of the clipped polygon. This is a circular vector.
"""
function clip_bounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive())
    poly, clip_poly = get_clipping_poly_structs(vorn, i, bounding_box)
    return clip_polygon(poly, clip_poly; predicates)
end

"""
    _get_ray(vorn, i, ghost_vertex, predicates::AbstractPredicateType=Adaptive()) -> (Point, Point)

Extracts the ray from the `i`th polygon of `vorn` corresponding to the `ghost_vertex`, where `ghost_vertex`
here means that `get_polygon(vorn, i)[ghost_vertex]` is a ghost vertex.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The index of the polygon.
- `ghost_vertex`: The index of the ghost vertex in the polygon.
- `predicates::AbstractPredicateType=Adaptive()`: Method to use for computing predicates. Can be one of [`Fast`](@ref), [`Exact`](@ref), and [`Adaptive`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `p`: The first point of the ray.
- `q`: A second point of the ray, so that `pq` gives the direction of the ray (which extends to infinity).
"""
function _get_ray(vorn::VoronoiTessellation, i, ghost_vertex, predicates::AbstractPredicateType=Adaptive())
    C = get_polygon(vorn, i)
    ghost_tri = get_circumcenter_to_triangle(vorn, C[ghost_vertex])
    u, v, _ = triangle_vertices(ghost_tri) # w is the ghost vertex
    p, q = get_generator(vorn, u, v)
    qx, qy = getxy(q)
    mx, my = midpoint(p, q)
    m = (mx, my)
    is_first = is_first_ghost_vertex(C, ghost_vertex)
    if is_first
        prev_index = previndex_circular(C, ghost_vertex)
        r = get_polygon_point(vorn, C[prev_index])
    else
        next_index = nextindex_circular(C, ghost_vertex)
        r = get_polygon_point(vorn, C[next_index])
    end
    if r == m # It's possible for the circumcenter to lie on the edge and exactly at the midpoint (e.g. [(0.0,1.0),(-1.0,2.0),(-2.0,-1.0)]). In this case, just rotate 
        dx, dy = qx - mx, qy - my
        if is_right(point_position_relative_to_line(predicates, p, q, r))
            rotated_dx, rotated_dy = dy, -dx
            r = mx + rotated_dx, my + rotated_dy
        else
            rotated_dx, rotated_dy = -dy, dx
            r = mx + rotated_dx, my + rotated_dy
        end
    end
    if is_right(point_position_relative_to_line(predicates, p, q, r)) # in this case, the circumcenter is inside and so mr points inwards rather than outwards 
        m, r = r, m
    end
    return m, r
end

"""
    maximum_distance_to_box(a, b, c, d, p) -> Number

Computes the maximum squared distance from the point `p` to the box with corners `(a, c)`, `(b, c)`, `(b, d)`, `(a, d)`.

# Arguments
- `a`: The minimum `x`-coordinate of the box. 
- `b`: The maximum `x`-coordinate of the box.
- `c`: The minimum `y`-coordinate of the box.
- `d`: The maximum `y`-coordinate of the box.
- `p`: The point.

# Outputs
- `dist`: The maximum squared distance from `p` to the box.
"""
function maximum_distance_to_box(a, b, c, d, p)
    F = number_type(p)
    p1 = (F(a), F(c))
    p2 = (F(b), F(c))
    p3 = (F(b), F(d))
    p4 = (F(a), F(d))
    d1 = dist_sqr(p1, p)
    d2 = dist_sqr(p2, p)
    d3 = dist_sqr(p3, p)
    d4 = dist_sqr(p4, p)
    return max(d1, d2, d3, d4)
end

"""
    grow_polygon_outside_of_box(vorn::VoronoiTessellation, i, bounding_box)

Truncates unbounded edges of the `i`th polygon of `vorn`, assumed to be unbounded,
so that the line connecting the truncated unbounded edges is entirely outside 
of the polygon. The method of growth is iterative, utilising the Liang-Barsky algorithm 
at each stage while we translate the line. The returned polygon does not satisfy 
`P[begin] == P[end]`.
"""

"""
    grow_polygon_outside_of_box(vorn::VoronoiTessellation, i, bounding_box, predicates::AbstractPredicateType=Adaptive()) -> (Vector{Int}, Vector{NTuple{2,Number}})

Truncates the unbounded edges of the `i`th polygon of `vorn` so that the line connecting the truncated unbounded edges is entirely outside of `bounding_box`.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `i`: The index of the polygon. The polygon must be unbounded.
- `bounding_box`: The bounding box to clip the polygon to. See also [`polygon_bounds`](@ref).
- `predicates::AbstractPredicateType=Adaptive()`: Method to use for computing predicates. Can be one of [`Fast`](@ref), [`Exact`](@ref), and [`Adaptive`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
- `new_vertices`: The new vertices of the polygon. This is not a circular vector.
- `new_points`: The new points of the polygon. This is not a circular vector.
"""
function grow_polygon_outside_of_box(vorn::VoronoiTessellation, i, bounding_box, predicates::AbstractPredicateType=Adaptive())
    a, b, c, d = bounding_box
    vertices = get_polygon(vorn, i)
    new_vertices, new_points, ghost_vertices = get_new_polygon_indices(vorn, vertices)
    inside = true
    t = 1.0 # don't do 0.5 so we get t = 1 later, else we get duplicated vertices for polygons completely outside of the box
    u, v = ghost_vertices
    u_m, u_r = _get_ray(vorn, i, u, predicates)
    v_m, v_r = _get_ray(vorn, i, v, predicates)
    u_mx, u_my = getxy(u_m)
    u_rx, u_ry = getxy(u_r)
    v_mx, v_my = getxy(v_m)
    v_rx, v_ry = getxy(v_r)
    p = (0.0, 0.0)
    q = (0.0, 0.0)
    dist_to_box = maximum_distance_to_box(a, b, c, d, u_m) # this is a squared distance
    dist_to_box = max(dist_to_box, maximum_distance_to_box(a, b, c, d, v_m))
    while inside
        t *= 2.0
        p = (u_mx + t * (u_rx - u_mx), u_my + t * (u_ry - u_my))
        q = (v_mx + t * (v_rx - v_mx), v_my + t * (v_ry - v_my))
        int1, int2 = liang_barsky(a, b, c, d, p, q)
        outside = all(isnan, int1) && all(isnan, int2)
        # We need to be careful of the case where the generator is outside of the bounding box. In this case, 
        # the unbounded edge might start initially outside of the box but then find it's way inside. 
        # So, to avoid this, we also apply a conservative check that the length of each ray is greater than 
        # the maximum distance from the generators to the bounding box.
        # See the example with [(-3,7),(1,6),(-1,3),(-2,4),(3,-2),(5,5),(-4,-3),(3,8)] and bb = (0,5,-15,15) with the 7th polygon.
        p_length = dist_sqr(p, (u_mx, u_my))
        q_length = dist_sqr(q, (v_mx, v_my))
        might_be_inside = min(p_length, q_length) < dist_to_box
        outside = outside && !might_be_inside
        inside = !outside
    end
    new_points[u] = p
    new_points[v] = q
    new_vertices[u] = u
    new_vertices[v] = v
    return new_vertices, new_points
end

"""
    get_new_polygon_indices(vorn, vertices) -> (Vector{Int}, Vector{NTuple{2,Float64}}, Tuple{Int, Int})

Returns the new vertices and points of the polygon, as well as the indices of the ghost vertices in the polygon.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `vertices`: The vertices of the polygon.

# Outputs 
- `new_vertices`: The new vertices of the polygon. This is not a circular vector. The vertices corresponding to a ghost vertex will be given by the ghost vertex itself.
- `new_points`: The new points of the polygon. This is not a circular vector. The points corresponding to a ghost vertex will be given by by `(NaN, NaN)`.
- `ghost_vertices`: The indices of the ghost vertices in `new_vertices`.
"""
function get_new_polygon_indices(vorn, vertices)
    new_points = NTuple{2,Float64}[]
    sizehint!(new_points, length(vertices))
    new_vertices = similar(vertices, length(vertices) - 1)
    ghost_vertices = (0, 0)
    for i in firstindex(vertices):(lastindex(vertices)-1)
        v = vertices[i]
        if is_ghost_vertex(v)
            is_first = is_first_ghost_vertex(vertices, i)
            if is_first
                ghost_vertices = (i, ghost_vertices[2])
            else
                ghost_vertices = (ghost_vertices[1], i)
            end
            push!(new_points, (NaN, NaN))
            new_vertices[i] = v
        else
            push!(new_points, getxy(get_polygon_point(vorn, v)))
            new_vertices[i] = length(new_points)
        end
    end
    return new_vertices, new_points, ghost_vertices
end

"""
    get_bounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive()) -> Vector{NTuple{2,Number}}

Returns the coordinates of the `i`th polygon of `vorn`, clipped to `bounding_box`.

Use the keyword arguments `predicates` to determine how predicates are computed. Should be one of [`Exact`](@ref),
[`Adaptive`](@ref), and [`Fast`](@ref). See the documentation for more information about these choices.
"""
function get_bounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive())
    if isnothing(bounding_box)
        C = get_polygon(vorn, i)
        F = number_type(vorn)
        coords = Vector{NTuple{2,F}}(undef, length(C))
        for j in eachindex(C)
            coords[j] = get_polygon_point(vorn, C[j])
        end
        return coords
    else
        return clip_bounded_polygon_to_bounding_box(vorn, i, bounding_box; predicates)
    end
end

"""
    get_unbounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive()) -> Vector{NTuple{2,Number}}

Returns the coordinates of the `i`th polygon of `vorn`, clipped to `bounding_box`. The polygon is assumed to be unbounded.

Use the keyword arguments `predicates` to determine how predicates are computed. Should be one of [`Exact`](@ref),
[`Adaptive`](@ref), and [`Fast`](@ref). See the documentation for more information about these choices.
"""
function get_unbounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive())
    return clip_unbounded_polygon_to_bounding_box(vorn, i, bounding_box; predicates)
end

"""
    clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box; predicates::AbsractPredicateType=Adaptive()) -> Vector{NTuple{2,Number}}

Clips the `i`th polygon of `vorn` to `bounding_box`. The polygon is assumed to be unbounded. See also [`clip_polygon`](@ref).

Use the keyword arguments `predicates` to determine how predicates are computed. Should be one of [`Exact`](@ref),
[`Adaptive`](@ref), and [`Fast`](@ref). See the documentation for more information about these choices.
"""
function clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box; predicates::AbstractPredicateType=Adaptive())
    new_vertices, new_points = grow_polygon_outside_of_box(vorn, i, bounding_box, predicates)
    clip_vertices = (1, 2, 3, 4)
    a, b, c, d = bounding_box
    clip_points = ((a, c), (b, c), (b, d), (a, d))
    return clip_polygon(new_vertices, new_points, clip_vertices, clip_points; predicates)
end