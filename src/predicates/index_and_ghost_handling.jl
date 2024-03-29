
"""
    triangle_orientation(i, j, k, pts, representative_point_list, boundary_map)
    triangle_orientation(T, pts, representative_point_list, boundary_map)

Given a triangle `T = (i, j, k)`, with indices corresponding to points in `pts`, computes the orientation of the triangle,
using `boundary_map` from [`construct_boundary_map`](@ref) to map boundary indices to their corresponding points in
`representative_point_list`. We return:

- `Certificate.PositivelyOriented`: The triangle is positively oriented.
- `Certificate.Degenerate`: The triangle is degenerate, meaning the coordinates are collinear. 
- `Certificate.NegativelyOriented`: The triangle is negatively oriented.

!!! note

    A test is also made for the case that `is_outer_ghost_triangle(T)`: If  `T` 
    is a ghost triangle, then the index corresponding to a boundary index 
    points to a centroid, in which case one of the edges has its orientation 
    flipped. This case will also be handled correctly. In case the boundary 
    index corresponds to an interior curve, this flip is not necessary.
"""
@inline function triangle_orientation(i, j, k, pts, representative_point_list, boundary_map)
    p, q, r = get_point(pts, representative_point_list, boundary_map, i, j, k)
    if is_outer_ghost_triangle(i, j, k, boundary_map)
        return triangle_orientation(r, q, p) # Exterior ghost triangles have the boundary index represented at a centroid which is inwards, so flip the orientation
    end
    return triangle_orientation(p, q, r)
end
@inline function triangle_orientation(T, pts, representative_point_list, boundary_map::AbstractDict)
    i, j, k = indices(T)
    return triangle_orientation(i, j, k, pts, representative_point_list, boundary_map)
end

"""
    point_position_relative_to_circumcircle(i, j, k, ℓ, pts, representative_point_list, boundary_map)
    point_position_relative_to_circumcircle(T, ℓ, pts, representative_point_list, boundary_map)

Tests if the `ℓ`th point of `pts` is inside the circumcircle of the triangle `T = (i, j, k)`, using the `boundary_map` to map 
boundary indices to their corresponding points in `representative_point_list`, returning:

- `Certificate.Outside`: `pₗ` is outside of the circumcircle.
- `Certificate.On`: `pₗ` is on the circumcircle.
- `Certificate.Inside`: `pₗ` is inside the circumcircle.

!!! note

    A test is also made for the case that `is_ghost_triangle(T)`: When `T` is a ghost triangle, one of its indices is a boundary index, say `i`. Since this vertex 
    is treated as being out at infinity, the circumcircle degenerates into the line through the other two vertices and out to infinity in that direction. 
    Thus, we test that the `ℓ`th point is inside this circumcircle by seeing if it is in the oriented outer halfplane defined by the other two vertices, 
    accomplished via [`point_position_relative_to_oriented_outer_halfplane`](@ref).
"""
function point_position_relative_to_circumcircle(i, j, k, ℓ, pts, representative_point_list, boundary_map)
    a, b, c, p = get_point(pts, representative_point_list, boundary_map, i, j, k, ℓ)
    if is_boundary_index(i)
        return point_position_relative_to_oriented_outer_halfplane(b, c, p)
    elseif is_boundary_index(j)
        return point_position_relative_to_oriented_outer_halfplane(c, a, p)
    elseif is_boundary_index(k)
        return point_position_relative_to_oriented_outer_halfplane(a, b, p)
    end
    return point_position_relative_to_circle(a, b, c, p)
end
function point_position_relative_to_circumcircle(T, ℓ, pts, representative_point_list, boundary_map)
    i, j, k = indices(T)
    return point_position_relative_to_circumcircle(i, j, k, ℓ, pts, representative_point_list, boundary_map)
end

"""
    point_position_relative_to_line(i, j, u, pts, representative_point_list, boundary_map)

Computes the position of the `u`th point of `pts` relative to the line through the `i`th and `j`th points, 
respectively, of `pts`. Boundary indices are mapped to their corresponding points in `representative_point_list` via 
the `boundary_map` argument from [`construct_boundary_map`](@ref). The returned values are:

- `Certificate.Left`: `p` is to the left of the line. 
- `Certificate.Collinear`: `p` is on the line.
- `Certificate.Right`: `p` is to the right of the line,

where `p` is the `u`th point of `pts`.

!!! note

    If `is_outer_ghost_edge(i, j, boundary_map)`, the orientation of the line is flipped as the point corresponding to the boundary index will be a
    centroid which swaps the orientation.
"""
function point_position_relative_to_line(i, j, u, pts, representative_point_list, boundary_map)
    a, b, p = get_point(pts, representative_point_list, boundary_map, i, j, u)
    if !is_outer_ghost_edge(i, j, boundary_map)
        return point_position_relative_to_line(a, b, p)
    else
        return point_position_relative_to_line(b, a, p) # Exterior ghost edges have the boundary index represented at a centroid which is inwards, so flip the orientation
    end
end

"""
    point_closest_to_line(i, j, u, v, pts)

Let `a, b, p, q` be the points corresponding to the indices `i, j, u, v`, respectively, in `pts`, and let `ℓ`
be the oriented line through `a` and `b`. This function tests if `p` is closer to `ℓ` than `q` is, returning:

- `Certificate.Closer`: `p` is closer to `ℓ`.
- `Certificate:Further`: `q` is closer to `ℓ`.
- `Certificate.Equidistant`: `p` and `q` are the same distance from `ℓ`.

!!! note 

    It is assumed that `p` and `q` are to the left of `ℓ`.
"""
function point_closest_to_line(i, j, u, v, pts)
    a, b, p, q = get_point(pts, i, j, u, v)
    return point_closest_to_line(a, b, p, q)
end

"""
    point_position_on_line_segment(i, j, u, pts)

Given indices `i`, `j`, and `u` corresponding to points `a`, `b`, and `p` in `pts`, respectively, computes the position 
of `p` relative to the oriented line segment `(a, b)`, assuming that the three points are collinear. The returned values are:

- `Certificate.On`: `p` is on the line segment, meaning between `a` and `b`.
- `Certificate.Degenerate`: Either `p == a` or `p == b`, i.e. `p` is one of the endpoints. 
- `Certificate.Left`: `p` is off and to the left of the line segment.
- `Certificate.Right`: `p` is off and to the right of the line segment.
"""
function point_position_on_line_segment(i, j, u, pts)
    a, b, p = get_point(pts, i, j, u)
    return point_position_on_line_segment(a, b, p)
end

"""
    line_segment_intersection_type(u, v, i, j, pts)

Let `u`, `v`, `i`, and `j` be indices corresponding to points `p`, `q`, `a`, and `b`, respectively, in `pts`. This 
function tests the number of intersections between the two line segments `(p, q)` and `(a, b)`, returning:

- `Certificate.None`: The line segments do not meet at any points. 
- `Certificate.Multiple`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `Certificate.Single`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `Certificate.On`: One of the endpoints is on `[a, b]`, but there are no other intersections.
"""
function line_segment_intersection_type(u, v, i, j, pts)
    p, q, a, b = get_point(pts, u, v, i, j)
    return line_segment_intersection_type(p, q, a, b)
end

"""
    point_position_relative_to_triangle(i, j, k, u, pts, representative_point_list, boundary_map)
    point_position_relative_to_triangle(T, u, pts, representative_point_list, boundary_map)

Given a triangle `T = (i, j, k)`, with indices referring to points in `pts`,
computes the position of `u`, corresponding to a point `p`, relative to `T`, with 
any boundary indices mapped to their corresponding representative points in `representative_point_list` 
via the `boundary_map` argument from [`construct_boundary_map`](@ref). The returned values are:

- `Certificate.Outside`: `p` is outside of the triangle. 
- `Certificate.On`: `p` is on one of the edges. 
- `Certificate.Inside`: `p` is inside the triangle.
"""
function point_position_relative_to_triangle(i, j, k, u, pts, representative_point_list, boundary_map::AbstractDict)
    if !is_outer_ghost_triangle(i, j, k, boundary_map)
        a, b, c, p = get_point(pts, representative_point_list, boundary_map, i, j, k, u)
        return point_position_relative_to_triangle(a, b, c, p)
    else
        i, j, k = rotate_ghost_triangle_to_standard_form(i, j, k)
        a, b, c, p = get_point(pts, representative_point_list, boundary_map, i, j, k, u) # a and b are solid vertices, but c is a ghost (mapped to a centroid in the opposite direction)
        edge_ab = point_position_relative_to_line(a, b, p)
        is_right(edge_ab) && return Cert.Outside
        if is_collinear(edge_ab)
            cert = point_position_on_line_segment(a, b, p)
            return (is_left(cert) || is_right(cert)) ? Cert.Outside : Cert.On
        end
        edge_bc = point_position_relative_to_line(c, b, p) # Flipped to match centroid location
        is_right(edge_bc) && return Cert.Outside
        #is_collinear(edge_bc) && return Cert.On # Don't need to check that it's not on the (c, b) part, since we already know we're to the left of (a, b) at this point
        edge_ca = point_position_relative_to_line(c, a, p)
        is_left(edge_ca) && return Cert.Outside
        #is_collinear(edge_ca) && return Cert.On
        return Cert.Inside

        # The collinear tests were deleted for the ghost edges. It doesn't really make much sense to see if a point is 
        # on the ghost edges. It's not like we can do anything with that information, and if we are using it then 
        # there's no point distinguishing between the two adjacent ghost triangles in that case.
    end
end
function point_position_relative_to_triangle(T, u, pts, representative_point_list, boundary_map::AbstractDict)
    i, j, k = indices(T)
    return point_position_relative_to_triangle(i, j, k, u, pts, representative_point_list, boundary_map)
end

"""
    triangle_line_segment_intersection(i, j, k, u, v, pts)

Given a triangle `(i, j, k)` and a line segment `(u, v)`,
with indices corresponding to points in `pts`, tests if `(u, v)` 
intersects the triangle's interior. Letting `(p, q, r)` be the coordinates 
corresponding to the triangle's vertices, and `(a, b)` those for the edge's 
vertices, returns:

- `Cert.Inside`: `(a, b)` is entirely inside `(p, q, r)`.
- `Cert.Single`: `(a, b)` has one endpoint inside `(p, q, r)`, and the other is outside.
- `Cert.Outside`: `(a, b)` is entirely outside `(p, q, r)`.
- `Cert.Touching`: `(a, b)` is on `(p, q, r)`'s boundary, but not in its interior.
- `Cert.Multiple`: `(a, b)` passes entirely through `(p, q, r)`. This includes the case where a point is on the boundary of `(p, q, r)`.
"""
@inline function triangle_line_segment_intersection(i, j, k, u, v, pts)
    p, q, r, a, b = get_point(pts, i, j, k, u, v)
    return triangle_line_segment_intersection(p, q, r, a, b)
end