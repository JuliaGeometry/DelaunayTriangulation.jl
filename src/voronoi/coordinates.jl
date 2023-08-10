"""
    get_polygon_coordinates(vorn::VoronoiTessellation, i; bounding_box = nothing)

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
function get_polygon_coordinates(vorn::VoronoiTessellation, i; bounding_box=nothing)
    if !isnothing(bounding_box)
        a, b, c, d = bounding_box
        @assert a < b && c < d "The bounding box must be of the form (xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax."
    end
    if i ∈ get_unbounded_polygons(vorn)
        isnothing(bounding_box) && throw(ArgumentError("The polygon is unbounded, so a bounding box must be provided."))
        return get_unbounded_polygon_coordinates(vorn, i, bounding_box)
    else
        return get_bounded_polygon_coordinates(vorn, i, bounding_box)
    end
end

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

Extracts the ray from the `i`th polygon of `vorn` corresponding to the `boundary_index`th boundary index. 
The returned points are given in the form `(p, q)`, defining the oriented line `pq` such that the line 
is in the direction of infinity.
"""
function _get_ray(vorn, i, boundary_index)
    ghost_tri = get_circumcenter_to_triangle(vorn, boundary_index)
    u, v, _ = indices(ghost_tri) # w is the ghost vertex
    p, q = get_generator(vorn, u, v)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    mx, my = (px + qx) / 2, (py + qy) / 2
    m = (mx, my)
    C = get_polygon(vorn, i)
    is_first = is_first_boundary_index(C, i)
    if is_first
        prev_index = previndex_circular(C, i)
        r = get_polygon_point(vorn, C[prev_index])
    else
        next_index = nextindex_circular(C, i)
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

# wrap into a simpler vector so that the vertices are [1, 2, ...], 
# making it easier to push intersection points into the points vector
function get_new_polygon_indices(vorn, vertices)
    nv = count(!is_boundary_index, vertices)
    F = Float64
    T = NTuple{2,F}
    new_points = Vector{T}(undef, nv)
    new_vertices = zero(vertices)
    local ctr = 0
    for (i, v) in enumerate(vertices)
        if is_boundary_index(v)
            new_vertices[i] = v
        else
            ctr += 1
            new_vertices[i] = ctr
            new_points[ctr] = _getxy(get_polygon_point(vorn, v))
        end
    end
    return new_vertices, new_points
end

"""
    clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)

Clips the `i`th polygon of `vorn` to the bounding box, assuming it is unbounded.
A modified form of the Sutherland-Hodgman algorithm is used. 
"""
function clip_unbounded_polygon_to_bounding_box(vorn::VoronoiTessellation, i, bounding_box)
    poly, clip_poly = get_clipping_poly_structs(vorn, i, bounding_box)
    vertices, clip_vertices, clip_points = poly.vertices, clip_poly.vertices, clip_poly.points
    new_vertices, new_point_list = get_new_polygon_indices(vorn, vertices)
    output_vertices = new_vertices
    output_points = deepcopy(new_point_list)
    q = clip_points[end]
    for p in clip_points
        input_vertices = output_vertices
        output_vertices = _clip_unbounded_polygon_to_edge(vorn, i, input_vertices, output_points, q, p)
        q = p
    end
    return output_vertices
end
function _clip_unbounded_polygon_to_edge(vorn::VoronoiTessellation, i, input_vertices, new_point_list, q, p)
    T = typeof(q)
    I = eltype(input_vertices)
    output_vertices = I[]
    output_points = T[]
    s_vertex = input_vertices[end]
    for vertex in input_vertices
        if is_segment_between_two_ghosts(s_vertex, vertex)
            s_vertex = vertex
            continue
        end
        if is_ray_going_in(s_vertex, vertex) # is_boundary_index(s_vertex)
            _clip_unbounded_polygon_edge_to_ray_going_in!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
        elseif is_ray_going_out(s_vertex, vertex) # is_boundary_index(vertex)
            _clip_unbounded_polygon_edge_to_ray_going_out!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
        else # is_finite_segment(s, vert). This is just the standard Sutherland-Hodgman algorithm
            _clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
        end
        s_vertex = vertex
    end
    return output_vertices
end
function _clip_unbounded_polygon_edge_to_ray_going_in!(vorn::VoronoiTessellation, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
    m, r = _get_ray(vorn, i, s_vertex)
    point = new_point_list[vertex]
    # Here, we have a ray that's coming from infinity into vertex. 
    # The intersection in this case requires that vertex is 
    # is to the right of the edge. So, first we check if vertex
    # is to the left.
    if (is_left ∘ point_position_relative_to_line)(q, p, point)
        # Now, vertex is to the left of the line, so there is nothing to be done. 
        push!(output_vertices, vertex)
        push!(output_points, point)
    else
        # We know that vert is to the right of the edge, so we have an intersection. 
        r = intersection_of_ray_with_edge(q, p, m, r)
        push!(new_point_list, r)
        push!(output_vertices, length(new_point_list))
        push!(output_points, r)
    end
end
function _clip_unbounded_polygon_edge_to_ray_going_out!(vorn::VoronoiTessellation, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
    m, r = _get_ray(vorn, i, s_vertex)
    point = new_point_list[s_vertex]
    # Here, we have a ray that's going from s_vertex to infinity. An intersection 
    # in this case requires that vertex is to the right of the edge.
    if (is_left ∘ point_position_relative_to_line)(q, p, point)
        push!(output_vertices, vertex)
        push!(output_points, point)
    else
        r = intersection_of_ray_with_edge(q, p, m, r)
        push!(new_point_list, r)
        push!(output_vertices, length(new_point_list))
        push!(output_points, r)
    end
end
function _clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
    s_point = new_point_list[s_vertex]
    point = new_point_list[vertex]
    if (is_left ∘ point_position_relative_to_line)(q, p, point)
        flag = point_position_relative_to_line(q, p, s_point)
        if !is_left(flag)
            r = segment_intersection_coordinates(q, p, s_point, point)
            push!(new_point_list, r)
            push!(output_vertices, length(new_point_list))
            push!(output_points, r)
        end
        push!(output_vertices, vertex)
        push!(output_points, point)
    elseif (is_left ∘ point_position_relative_to_line)(q, p, s_point)
        r = segment_intersection_coordinates(q, p, s_point, point)
        push!(new_point_list, r)
        push!(output_vertices, length(new_point_list))
        push!(output_points, r)
    end
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
    F = number_type(vorn)
    if (is_outside ∘ polygon_position_relative_to_box)(vorn, bounding_box, i)
        return NTuple{2,F}[]
    else
        return clip_unbounded_polygon_to_bounding_box(vorn, i, bounding_box)
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