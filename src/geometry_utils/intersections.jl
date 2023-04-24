"""
    intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)

Finds the intersection of the ray through `p` and `q`, oriented from `p` to `q`,
with the boundary of the polygon defined by `points` and `boundary_nodes`. 

It is assumed that `p` is inside the polygon, but `q` could be outside or inside.

Currently, this function has only been tested on rectangular boundaries.
"""
function intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)
    p == q && throw(ArgumentError("p and q must be distinct."))
    px, py = getxy(p)
    qx, qy = getxy(q)

    ## Start by finding a bracketing for the intersection point 
    t1 = zero(px)
    t2 = one(px)
    δ1 = distance_to_polygon(p, points, boundary_nodes)
    δ2 = distance_to_polygon(q, points, boundary_nodes)
    sign(δ1) == -1 && throw(ArgumentError("p must be inside the polygon."))
    sign(δ1) == 0 && return p
    sign(δ2) == 0 && return q
    iters = 0
    while sign(δ2) == 1
        t2 *= 2
        r = px + t2 * (qx - px), py + t2 * (qy - py)
        δ2 = distance_to_polygon(r, points, boundary_nodes)
        iters += 1 
        iters ≥ 1000 && throw(ArgumentError("Something went wrong. Could not find a bracketing interval."))
    end

    ## Perform bisection 
    t = (t1 + t2) / 2
    r = px + t * (qx - px), py + t * (qy - py)
    δ = distance_to_polygon(r, points, boundary_nodes)
    iters = 0
    while (t2 - t1) / 2 > tol && abs(δ) > tol
        if sign(δ) == sign(δ1)
            t1 = t
            δ1 = δ
        else
            t2 = t
            δ2 = δ
        end
        t = (t1 + t2) / 2
        r = px + t * (qx - px), py + t * (qy - py)
        δ = distance_to_polygon(r, points, boundary_nodes)
        iters += 1
        iters ≥ 1000 && throw(ArgumentError("Something went wrong. Could not find the intersection point."))
    end
    return r
end

"""
    segment_intersection_coordinates(a, b, c, d)

Finds the coordinates of the intersection of the line segment from `a` to `b` 
with the line segment from `c` to `d`.
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
    intersection_of_edge_and_bisector_ray(a, b, c)

Given an edge `(a, b)` and a ray emanating from `c` perpendicular  
with the edge and collinear with its midpoint, tests if `c` intersects 
the edge, and if so, returns the intersection point which is the midpoint.
If there is no intersection, `(NaN, NaN)` is returned. The ray should be directed
to the left of the edge.
"""
function intersection_of_edge_and_bisector_ray(a, b, c)
    cert = point_position_relative_to_line(a, b, c)
    if !is_left(cert)
        ax, ay = getxy(a)
        bx, by = getxy(b)
        m = (ax + bx) / 2, (ay + by) / 2
        return m 
    else
        F = number_type(a) 
        return F(NaN), F(NaN)
    end
end

"""
    classify_and_compute_segment_intersection(a, b, c, d)

Given two line segments `(a, b)` and `(c, d)`, classifies the intersection
of the two segments, returning `(cert, cert_c, cert_d, p)`, where `p` is the intersection point
or `(NaN, NaN)` if there is no intersection. The certificate `cert` determines the intersection 
type, giving 

- `Cert.None`: No intersections.
- `Cert.Single`: There is an intersection point.
- `Cert.Touching`: There is an intersection point, and one of `c` and `d` is the intersection point.
- `Cert.Multiple`: The closed segments meet in one or several points.

The certificates `cert_c` and `cert_d` similarly return the positions of `c` and `d` relative to `(a, b)`,
respectively.
"""
function classify_and_compute_segment_intersection(a, b, c, d)
    cert = line_segment_intersection_type(a, b, c, d)
    cert_c = point_position_relative_to_line(a, b, c)
    cert_d = point_position_relative_to_line(a, b, d)
    if !is_none(cert)
        p = segment_intersection_coordinates(a, b, c, d)
        return cert, cert_c, cert_d, p
    else
        F = number_type(a)
        return cert, cert_c, cert_d, (F(NaN), F(NaN))
    end
end

