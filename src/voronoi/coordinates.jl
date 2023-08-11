"""
    get_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box = nothing)

Returns a vector for the coordinates of the `i`th polygon in `vorn`. If `bounding_box` 
is provided, the polygon will be clipped to the bounding box, assuming that it takes 
the form `(xmin, xmax, ymin, ymax)`. Some specific cases:

- If the polygon is unbounded but `bounding_box` is `nothing`, then an error will be thrown. 
- If the polygon is bounded and `bounding_box` is `nothing`, then the polygon coordinates will be returned as if without any clipping.
- If the polygon is outside of the bounding box entirely, then an empty vector will be returned.

If you do need to consider clipping your polygon to an arbitrary polygon, see the 
[`polygon_clip`](@ref) function; this function (`get_polygon_coordinates`) uses 
`polygon_clip` for rectangular clipping when `bounding_box` is considered.

See also [`polygon_bounds`](@ref) for a good default for `bounding_box`.
"""
function get_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box=nothing)
    if !isnothing(bounding_box)
        a, b, c, d = bounding_box
        @assert a < b && c < d "The bounding box must be of the form (xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax."
    end
    if i ∈ get_unbounded_polygons(vorn)
        isnothing(bounding_box) && throw(ArgumentError("The polygon is unbounded, so a bounding box must be provided. See DelaunayTriangulation.polygon_bounds for a reasonable default."))
        return get_unbounded_polygon_coordinates(vorn, i, bounding_box)
    else
        return get_bounded_polygon_coordinates(vorn, i, bounding_box)
    end
end

#=
"""
    polygon_position_relative_to_box(vorn::VoronoiTessellation, bounding_box, i)

Tests the position of the `i`th polygon of `vorn` relative to the bounding box`, ignoring 
boundary indices. Returns:

- `Cert.Inside`: The polygon is contained entirely within the box.
- `Cert.Outside`: The polygon is entirely outside the box.
- `Cert.Touching`: The polygon is contained entirely within the box, but touches the boundary of the box. Note that if a polygon touches the boundary of the box but is entirely outside otherwise, `Cert.Outside` will be returned instead. 
- `Cert.Multiple`: The polygon intersects the box in multiple places.
"""
function polygon_position_relative_to_box(vorn::VoronoiTessellation, bounding_box, i)
    a, b, c, d = bounding_box
    vertices = get_polygon(vorn, i)
    points = get_polygon_points(vorn)
    flag = polygon_position_relative_to_box(a, b, c, d, vertices, points)
    return flag
end
=#

function get_clipping_poly_structs(vorn::VoronoiTessellation, i, bounding_box)
    vertices = get_polygon(vorn, i)
    points = get_polygon_points(vorn)
    clip_vertices = (1, 2, 3, 4)
    a, b, c, d = bounding_box
    clip_points = ((a, c), (b, c), (b, d), (a, d))
    return Polygon(vertices, points), Polygon(clip_vertices, clip_points)
end

"""
    clip_bounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)

Clips the `i`th polygon of `vorn` to the bounding box, assuming it is bounded.
The Sutherland-Hodgman algorithm is used. 
"""
function clip_bounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)
    poly, clip_poly = get_clipping_poly_structs(vorn, i, bounding_box)
    return clip_polygon(poly, clip_poly)
end

"""
    _get_ray(vorn, i, boundary_index)

Extracts the ray from the `i`th polygon of `vorn` corresponding to the `boundary_index`, where `boundary_index` 
here means that `get_polygon(vorn, i)[boundary_index]` is a boundary index.
The returned points are given in the form `(p, q)`, defining the oriented line `pq` such that the line 
is in the direction of infinity.
"""
function _get_ray(vorn, i, boundary_index)
    C = get_polygon(vorn, i)
    ghost_tri = get_circumcenter_to_triangle(vorn, C[boundary_index])
    u, v, _ = indices(ghost_tri) # w is the ghost vertex
    p, q = get_generator(vorn, u, v)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    mx, my = (px + qx) / 2, (py + qy) / 2
    m = (mx, my)
    is_first = is_first_boundary_index(C, boundary_index)
    if is_first
        prev_index = previndex_circular(C, boundary_index)
        r = get_polygon_point(vorn, C[prev_index])
    else
        next_index = nextindex_circular(C, boundary_index)
        r = get_polygon_point(vorn, C[next_index])
    end
    if r == m # It's possible for the circumcenter to lie on the edge and exactly at the midpoint (e.g. [(0.0,1.0),(-1.0,2.0),(-2.0,-1.0)]). In this case, just rotate 
        dx, dy = qx - mx, qy - my
        if (is_right ∘ point_position_relative_to_line)(p, q, r)
            rotated_dx, rotated_dy = dy, -dx
            r = mx + rotated_dx, my + rotated_dy
        else
            rotated_dx, rotated_dy = -dy, dx
            r = mx + rotated_dx, my + rotated_dy
        end
    end
    if (is_right ∘ point_position_relative_to_line)(p, q, r) # in this case, the circumcenter is inside and so mr points inwards rather than outwards 
        m, r = r, m
    end
    return m, r
end

"""
    maximum_distance_to_box(a, b, c, d, p)

Given a bounding box `[a, b] × [c, d]`, returns the maximum distance 
between `p` and the boundary of the box. The returned distance is squared.
"""
function maximum_distance_to_box(a, b, c, d, p)
    p1x, p1y = a, c 
    p2x, p2y = b, c
    p3x, p3y = b, d
    p4x, p4y = a, d
    px, py = _getxy(p)
    d1 = (p1x - px)^2 + (p1y - py)^2
    d2 = (p2x - px)^2 + (p2y - py)^2
    d3 = (p3x - px)^2 + (p3y - py)^2
    d4 = (p4x - px)^2 + (p4y - py)^2
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
function grow_polygon_outside_of_box(vorn::VoronoiTessellation, i, bounding_box)
    a, b, c, d = bounding_box
    vertices = get_polygon(vorn, i)
    new_vertices, new_points, boundary_indices = get_new_polygon_indices(vorn, vertices)
    inside = true
    t = 1.0 # don't do 0.5 so we get t = 1 later, else we get duplicated vertices for polygons completely outside of the box
    u, v = boundary_indices
    u_m, u_r = _get_ray(vorn, i, u)
    v_m, v_r = _get_ray(vorn, i, v)
    u_mx, u_my = _getxy(u_m)
    u_rx, u_ry = _getxy(u_r)
    v_mx, v_my = _getxy(v_m)
    v_rx, v_ry = _getxy(v_r)
    p = (0.0, 0.0)
    q = (0.0, 0.0)
    dist_to_box = maximum_distance_to_box(a, b, c, d, u_m)
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
        px, py = _getxy(p)
        qx, qy = _getxy(q)
        p_length = (px - u_mx)^2 + (py - u_my)^2
        q_length = (qx - v_mx)^2 + (qy - v_my)^2
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
    get_new_polygon_indices(vorn, vertices)

Given an unbounded Voronoi polygon from `vorn` with `vertices`, returns 
`(new_vertices, new_points, boundary_indices)`, where `new_vertices` are vertices 
mapping to the points in `new_points`, which is just a vector of all the polygon points, 
and `boundary_indices` is a `Tuple` of the indices of the points in `new_points` that correspond 
to points out at infinity.
"""
function get_new_polygon_indices(vorn, vertices)
    new_points = NTuple{2,Float64}[]
    sizehint!(new_points, length(vertices))
    new_vertices = similar(vertices, length(vertices) - 1)
    boundary_indices = (0, 0)
    for i in firstindex(vertices):(lastindex(vertices)-1)
        v = vertices[i]
        if is_boundary_index(v)
            is_first = is_first_boundary_index(vertices, i)
            if is_first
                boundary_indices = (i, boundary_indices[2])
            else
                boundary_indices = (boundary_indices[1], i)
            end
            push!(new_points, (NaN, NaN))
            new_vertices[i] = v
        else
            push!(new_points, _getxy(get_polygon_point(vorn, v)))
            new_vertices[i] = length(new_points)
        end
    end
    return new_vertices, new_points, boundary_indices
end

"""
    get_bounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box)

Returns the coordinates of the `i`th polygon of `vorn`, assuming it is bounded. If 
`bounding_box` is `nothing`, then the polygon coordinates will be returned as if without
any clipping. If `bounding_box` is provided, the polygon will be clipped to the bounding box 
using [`clip_to_bounding_box`](@ref).
"""
function get_bounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box)
    if isnothing(bounding_box)
        C = get_polygon(vorn, i)
        F = number_type(vorn)
        coords = Vector{NTuple{2,F}}(undef, length(C) - 1)
        for j in firstindex(C):(lastindex(C)-1)
            coords[j] = get_polygon_point(vorn, C[j])
        end
        return coords
    else
        return clip_bounded_polygon_to_bounding_box(vorn, i, bounding_box)
    end
end

"""
    get_unbounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box)

Returns the coordinates of the `i`th polygon of `vorn`, assuming it is unbounded, clipping 
the polygon to the bounding box.
"""
function get_unbounded_polygon_coordinates(vorn::VoronoiTessellation, i, bounding_box)
    return clip_unbounded_polygon_to_bounding_box(vorn, i, bounding_box)
end

"""
    clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)

Clips the `i`th polygon of `vorn` to the bounding box, assuming it is unbounded.
The unbounded polygon is truncated so that the line connecting the unbounded edges 
is outside of the bounded box, and then the Sutherland-Hodgman algorithm 
is used to clip the resulting polygon to the bounding box.
"""
function clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)
    new_vertices, new_points = grow_polygon_outside_of_box(vorn, i, bounding_box)
    clip_vertices = (1, 2, 3, 4)
    a, b, c, d = bounding_box
    clip_points = ((a, c), (b, c), (b, d), (a, d))
    return clip_polygon(new_vertices, new_points, clip_vertices, clip_points)
end