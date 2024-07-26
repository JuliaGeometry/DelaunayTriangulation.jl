"""
    intersection_of_ray_with_bounding_box(p, q, a, b, c, d) -> NTuple{2, Number}

Compute the intersection of the ray emanating from `p` and passing through `q` with the box `[a, b] × [c, d]`. It is assumed that 
`p` is inside of the box.
"""
function intersection_of_ray_with_bounding_box(p, q, a, b, c, d)
    px, py = getxy(p)
    qx, qy = getxy(q)
    pℓbx, pℓby = a, c
    pℓtx, pℓty = a, d
    prtx, prty = b, d
    prbx, prby = b, c
    θ = mod(atan(qy - py, qx - px), 2π)
    θlb = mod(atan(pℓby - py, pℓbx - px), 2π)
    θlt = mod(atan(pℓty - py, pℓtx - px), 2π)
    θrt = mod(atan(prty - py, prtx - px), 2π)
    θrb = mod(atan(prby - py, prbx - px), 2π)
    if θ == 0.0
        return b, py
    elseif θ == π / 2
        return px, d
    elseif θ == π
        return a, py
    elseif θ == 3π / 2
        return px, c
    elseif θlb ≤ θ ≤ θrb # y = c, a ≤ x ≤ b
        # y = py + Rsinθ = c ⟹ R = (c - py) / sinθ 
        # x = px + Rcosθ = px + (c - py)cotθ
        return px + (c - py) * cot(θ), c
    elseif θrt ≤ θ ≤ θlt # y = d, a ≤ x ≤ b 
        # y = py + Rsinθ = d ⟹ R = (d - py) / sinθ
        # x = px + Rcosθ = px + (d - py)cotθ
        return px + (d - py) * cot(θ), d
    elseif θlt ≤ θ ≤ θlb # x = a, c ≤ y ≤ d
        # x = px + Rcosθ = a ⟹ R = (a - px) / cosθ
        # y = py + Rsinθ = py + (a - px)tanθ
        return a, py + (a - px) * tan(θ)
    else # x = b, c ≤ y ≤ d
        # x = px + Rcosθ = b ⟹ R = (b - px) / cosθ
        # y = py + Rsinθ = py + (b - px)tanθ
        return b, py + (b - px) * tan(θ)
    end
end

"""
    segment_intersection_coordinates(a, b, c, d) -> NTuple{2, Number}

Given two segments `(a, b)` and `(c, d)` that are assumed to intersect, computes the coordinates of the intersection.
"""
function segment_intersection_coordinates(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    num = (cx - ax) * (dy - ay) - (cy - ay) * (dx - ax)
    den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)
    α = num / den
    return ax + α * (bx - ax), ay + α * (by - ay)
end

"""
    intersection_of_edge_and_bisector_ray([method::AbstractPredicateType=Adaptive(),] a, b, c) -> (Certificate, NTuple{2, Number})

Given an edge `(a, b)` and a ray emanating from `c` perpendicular
with the edge and collinear with its midpoint, tests if `c` intersects the edge. The returned value is `(cert, p)`, where:

- `cert`: A [`Certificate`](@ref) indicating the position of `c` relative to the line through `(a, b)`.
- `p`: The intersection point (which is the midpoint) if `c` intersects the edge, `(NaN, NaN)` otherwise.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function intersection_of_edge_and_bisector_ray(method::AbstractPredicateType, a, b, c)
    cert = point_position_relative_to_line(method, a, b, c)
    if !is_left(cert)
        ax, ay = getxy(a)
        bx, by = getxy(b)
        m = midpoint((ax, ay), (bx, by))
        return cert, m
    else
        F = number_type(a)
        return cert, (F(NaN), F(NaN))
    end
end
intersection_of_edge_and_bisector_ray(a, b, c) = intersection_of_edge_and_bisector_ray(Adaptive(), a, b, c)

"""
    classify_and_compute_segment_intersection([method::AbstractPredicateType,] a, b, c, d) -> (Certificate, Certificate, Certificate, NTuple{2, Number})

Given two line segments `(a, b)` and `(c, d)`, classifies the intersection of the two segments. The returned value is `(cert, cert_c, cert_d, p)`, where:

- `cert`: A [`Certificate`](@ref) indicating the intersection type.
- `cert_c`: A [`Certificate`](@ref) indicating the position of `c` relative to the line through `(a, b)`.
- `cert_d`: A [`Certificate`](@ref) indicating the position of `d` relative to the line through `(a, b)`.
- `p`: The intersection point if `cert` is `Cert.Single` or `Cert.Touching`, and `(NaN, NaN)` otherwise.
"""
function classify_and_compute_segment_intersection(method::AbstractPredicateType, a, b, c, d)
    F = number_type(a)
    if any(!isfinite, getxy(a)) || any(!isfinite, getxy(b)) || any(!isfinite, getxy(c)) || any(!isfinite, getxy(d))
        return Cert.None, Cert.None, Cert.None, (F(NaN), F(NaN))
    end
    cert = line_segment_intersection_type(method, a, b, c, d)
    cert_c = point_position_relative_to_line(method, a, b, c)
    cert_d = point_position_relative_to_line(method, a, b, d)
    if !is_none(cert)
        p = segment_intersection_coordinates(a, b, c, d)
        return cert, cert_c, cert_d, p
    else
        F = number_type(a)
        return cert, cert_c, cert_d, (F(NaN), F(NaN))
    end
end
classify_and_compute_segment_intersections(a, b, c, p) = classify_and_compute_segment_intersections(Adaptive(), a, b, c, p)

"""
    polygon_features(points, boundary_nodes) -> (Number, NTuple{2, Number})

Computes the signed area and centroid of the polygon defined by `(points, boundary_nodes)`. The `boundary_nodes` must match the specification in the documentation 
and in [`check_args`](@ref).
"""
function polygon_features(points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return polygon_features_multiple_curves(points, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return polygon_features_multiple_segments(points, boundary_nodes)
    else
        return polygon_features_single_segment(points, boundary_nodes)
    end
end
function polygon_features_single_segment(points, boundary_nodes; scale=Val(true))
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(points, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge+1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(points, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        area_contribution = xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ
        cx += (xᵢ + xᵢ₊₁) * area_contribution
        cy += (yᵢ + yᵢ₊₁) * area_contribution
        a += area_contribution
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    if is_true(scale)
        return a / 2, (cx / (3a), cy / (3a))
    else
        return a, (cx, cy)
    end
end
function polygon_features_multiple_segments(points, boundary_nodes)
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_single_segment(points, bn;
            scale=Val(false))
        a += sub_a
        cx += sub_cx
        cy += sub_cy
    end
    return a / 2, (cx / (3a), cy / (3a))
end
function polygon_features_multiple_curves(points, boundary_nodes)
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_multiple_segments(points, bn)
        cx += sub_a * sub_cx
        cy += sub_a * sub_cy
        a += sub_a
    end
    return a, (cx / a, cy / a)
end

"""
    squared_distance_to_segment(x₁, y₁, x₂, y₂, x, y) -> Number 

Given a line segment `(x₁, y₁) → (x₂, y₂)` and a query point `(x, y)`, returns the
squared distance from `(x, y)` to the line segment.
"""
function squared_distance_to_segment(x₁, y₁, x₂, y₂, x, y)
    qp₁_x = x - x₁
    qp₁_y = y - y₁
    p₁p₂_x = x₂ - x₁
    p₁p₂_y = y₂ - y₁
    denom = dist_sqr((x₁, y₁), (x₂, y₂))
    t = (qp₁_x * p₁p₂_x + qp₁_y * p₁p₂_y) / denom
    ts = min(max(t, zero(t)), one(t)) # https://math.stackexchange.com/a/330329/861404
    intersect_x = x₁ + ts * p₁p₂_x
    intersect_y = y₁ + ts * p₁p₂_y
    δ² = dist_sqr((x, y), (intersect_x, intersect_y))
    return δ²
end

"""
    distance_to_polygon(q, points, boundary_nodes) -> Number

Given a query point `q` and a polygon defined by `(points, boundary_nodes)`, returns the signed distance from `q` to the polygon. The `boundary_nodes` must match the specification in the documentation
and in [`check_args`](@ref).

See also [`dist`](@ref).
"""
function distance_to_polygon(q, points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return distance_to_polygon_multiple_curves(q, points, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return distance_to_polygon_multiple_segments(q, points, boundary_nodes)
    else
        return distance_to_polygon_single_segment(q, points, boundary_nodes)
    end
end
function distance_to_polygon_single_segment(q, points, boundary_nodes, is_in_outer=false; return_sqrt=Val(true))
    x, y = getxy(q)
    F = number_type(points)
    dist = typemax(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(points, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge+1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(points, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        if (yᵢ₊₁ > y) ≠ (yᵢ > y)
            intersect_x = (xᵢ - xᵢ₊₁) * (y - yᵢ₊₁) / (yᵢ - yᵢ₊₁) + xᵢ₊₁
            if x < intersect_x
                is_in_outer = !is_in_outer
            end
        end
        new_dist = squared_distance_to_segment(xᵢ, yᵢ, xᵢ₊₁, yᵢ₊₁, x, y)
        dist = new_dist < dist ? new_dist : dist
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer ? dist : -dist
end
function distance_to_polygon_multiple_segments(q, points, boundary_nodes, is_in_outer=-one(number_type(points)); return_sqrt=Val(true))
    F = number_type(points)
    dist = typemax(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_single_segment(q, points, bn, is_in_outer == one(F); return_sqrt=Val(false))
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer * dist
end
function distance_to_polygon_multiple_curves(q, points, boundary_nodes)
    F = number_type(points)
    is_in_outer = -one(F)
    dist = typemax(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_multiple_segments(q, points, bn, is_in_outer == one(F);
            return_sqrt=Val(false))
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    return is_in_outer * sqrt(dist)
end

"""
    polygon_bounds(points, boundary_nodes, check_all_curves=Val(false)) -> (Number, Number, Number, Number)

Computes the bounding box of the polygon defined by `(points, boundary_nodes)`. The `boundary_nodes` must match the specification in the documentation
and in [`check_args`](@ref). If `check_all_curves` is `true`, then the bounding box of the union of all curves of the `polygon` is computed instead of just the first curve.
"""
function polygon_bounds(points, boundary_nodes, check_all_curves=Val(false))
    if has_multiple_curves(boundary_nodes)
        if !is_true(check_all_curves)
            return polygon_bounds_multiple_segments(points, get_boundary_nodes(boundary_nodes, 1)) # 1 is the outermost boundary, unless you have a multiple polygon. PolygonHierarchy would be better for this but too bad 
        else
            F = number_type(points)
            xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
            nc = num_curves(boundary_nodes)
            for i in 1:nc
                bn = get_boundary_nodes(boundary_nodes, i)
                xminᵢ, xmaxᵢ, yminᵢ, ymaxᵢ = polygon_bounds_multiple_segments(points, bn)
                xmin = min(xminᵢ, xmin)
                xmax = max(xmaxᵢ, xmax)
                ymin = min(yminᵢ, ymin)
                ymax = max(ymaxᵢ, ymax)
            end
            return xmin, xmax, ymin, ymax
        end
    elseif has_multiple_sections(boundary_nodes)
        return polygon_bounds_multiple_segments(points, boundary_nodes)
    else
        return polygon_bounds_single_segment(points, boundary_nodes)
    end
end
function polygon_bounds_single_segment(points, boundary_nodes)
    F = number_type(points)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    n_edge = num_boundary_edges(boundary_nodes)
    for i in 1:n_edge
        vᵢ = get_boundary_nodes(boundary_nodes, i)
        pᵢ = get_point(points, vᵢ)
        xᵢ, yᵢ = getxy(pᵢ)
        xmin = min(xᵢ, xmin)
        xmax = max(xᵢ, xmax)
        ymin = min(yᵢ, ymin)
        ymax = max(yᵢ, ymax)
    end
    return xmin, xmax, ymin, ymax
end
function polygon_bounds_multiple_segments(points, boundary_nodes)
    F = number_type(points)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        _xmin, _xmax, _ymin, _ymax = polygon_bounds_single_segment(points, bn)
        xmin = xmin > _xmin ? _xmin : xmin
        xmax = xmax < _xmax ? _xmax : xmax
        ymin = ymin > _ymin ? _ymin : ymin
        ymax = ymax < _ymax ? _ymax : ymax
    end
    return xmin, xmax, ymin, ymax
end

"""
    sort_convex_polygon!(vertices, points) 

Sorts the vertices of a convex polygon in counter-clockwise order. The polygon is defined by
`(points, vertices)`, and the vertices are sorted in-place. It is 
assumed that the vertices are not circular, i.e. `vertices[begin] ≠ vertices[end]`.
"""
function sort_convex_polygon!(vertices, points)
    cx, cy = mean_points(points, vertices)
    to_angle = p -> atan(gety(p) - cy, getx(p) - cx)
    vert_to_angle = v -> to_angle(get_point(points, v))
    sort!(vertices, by=vert_to_angle)
    return vertices
end

"""
    get_plane_through_three_points(a, b, c) -> NTuple{4, Number}

Given three points `(a, b, c)` in `ℝ³` represented as `Tuple`s, computes the equation of the plane 
through the points. The result is given in the form `(α, β, γ, δ)`, so that the plane 
is given by 

    αx + βy + γz + δ = 0.

# Extended help
The equation of the plane is computed by expanding the equation 

```math 
\\det \\begin{bmatrix} x & y & z & 1 \\\\ a_x & a_y & a_z & 1 \\\\ b_x & b_y & b_z & 1 \\\\ c_x & c_y & c_z & 1 \\end{bmatrix} = 0.
```

From this, we find:

```math 
\\begin{align*}
\\alpha &= a_y b_z - a_z b_y - a_y c_z + a_z c_y + b_y c_z - b_z c_y, \\\\
\\beta &= a_z b_x - a_x b_z + a_x c_z - a_z c_x - b_x c_z + b_z c_x, \\\\
\\gamma &= a_x b_y - a_y b_x - a_x c_y + a_y c_x + b_x c_y - b_y c_x, \\\\
\\delta &= a_x b_z c_y - a_x b_y c_z + a_y b_x c_z - a_y b_z c_x - a_z b_x c_y + a_z b_y c_x.
\\end{align*}
```
"""
@inline function get_plane_through_three_points(a::Tuple, b::Tuple, c::Tuple)
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    α = ay * bz - az * by - ay * cz + az * cy + by * cz - bz * cy
    β = az * bx - ax * bz + ax * cz - az * cx - bx * cz + bz * cx
    γ = ax * by - ay * bx - ax * cy + ay * cx + bx * cy - by * cx
    δ = ax * bz * cy - ax * by * cz + ay * bx * cz - ay * bz * cx - az * bx * cy + az * by * cx
    return α, β, γ, δ
end

"""
    get_steepest_descent_direction(a, b, c) -> NTuple{2, Number}

Given three points in `ℝ³` defining a plane, returns the direction `(x, y)`
of the steepest descent along the plane. In particular, if 

    αx + βy + γz + δ = 0 

is the plane, then the steepest descent direction is `(α, β)/γ`. The returned 
value is given by `(x, y) = sign(γ)(α, β)`.

See also [`get_plane_through_three_points`](@ref).
"""
function get_steepest_descent_direction(a::Tuple, b::Tuple, c::Tuple)
    α, β, γ, _ = get_plane_through_three_points(a, b, c)
    sγ = sign(γ)
    return sγ * α, sγ * β
end

"""
    get_vertical_distance_to_plane(a, b, c, p) -> Number 

Returns the vertical distance from the point `p` to the plane defined by the points 
`(a, b, c)`. The distance is positive if `p` is above the plane.
"""
function get_vertical_distance_to_plane(a::Tuple, b::Tuple, c::Tuple, p::Tuple)
    α, β, γ, δ = get_plane_through_three_points(a, b, c)
    x, y, z = p
    plane_z = -(α * x + β * y + δ) / γ
    return z - plane_z
end

"""
    get_distance_to_plane(a, b, c, p) -> Number 

Returns the distance from the point `p` to the plane defined by the points 
`(a, b, c)`. The distance is positive if `p` is above the plane.
""" # https://doi.org/10.1016/B978-0-08-050755-2.50050-6
function get_distance_to_plane(a::Tuple, b::Tuple, c::Tuple, p::Tuple)
    α, β, γ, δ = get_plane_through_three_points(a, b, c)
    if γ < 0 # so that nz > 0
        α, β, γ, δ = -α, -β, -γ, -δ
    end
    nmag = sqrt(α^2 + β^2 + γ^2)
    nx, ny, nz, Jd = α / nmag, β / nmag, γ / nmag, δ / nmag
    x, y, z = p
    return x * nx + y * ny + z * nz + Jd
end

"""
    angle_between(p, q) -> Float64

Returns the angle between the vectors `p` and `q` in radians, treating `q`
as the base. See [this article](https://straypixels.net/angle-between-vectors/).
The returned angle is in `[0, 2π)`.
"""
function angle_between(p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    a = px * qx + py * qy
    b = px * -qy + py * qx
    return mod(atan(b, a), 2π)
end