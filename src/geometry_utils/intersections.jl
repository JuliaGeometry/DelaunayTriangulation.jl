"""
    intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)

Finds the intersection of the ray through `p` and `q`, oriented from `p` to `q`,
with the boundary of the polygon defined by `points` and `boundary_nodes`. 

It is assumed that `p` is inside the polygon, but `q` could be outside or inside.

Currently, this function has only been tested on rectangular boundaries.
"""
function intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)
    # TODO: Write this in terms of an angle θ rather than q, computing θ from pq. 
    # Probably don't need bisection if we do that.
    p == q && throw(ArgumentError("p and q must be distinct."))
    px, py = _getxy(p)
    qx, qy = _getxy(q)

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
    identify_side(r, a, b, c, d)

Given a rectangle `[a, b] × [c, d]` and a point `r` 
on the rectangle, returns a symbol denoting the side 
of the rectangle that the point is on:

- `:left`: `r` is on the left side of the rectangle.
- `:right`: `r` is on the right side of the rectangle.
- `:bottom`: `r` is on the bottom side of the rectangle.
- `:top`: `r` is on the top side of the rectangle.
"""
function identify_side(r, a, b, c, d)
    rx, ry = _getxy(r)
    if rx == a
        return :left
    elseif rx == b
        return :right
    elseif ry == c
        return :bottom
    elseif ry == d
        return :top
    else
        throw(ArgumentError("Something went wrong."))
    end
end

"""
    intersection_of_ray_with_bounding_box(p, q, a, b, c, d)

Given a ray starting at `p` and in the direction of `q`, finds the intersection 
of the ray with the bounding box `[a, b] × [c, d]`. It is assumed that `p` is inside 
the bounding box, but `q` can be inside or outside.
"""
function intersection_of_ray_with_bounding_box(p, q, a, b, c, d)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
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
    segment_intersection_coordinates(a, b, c, d)

Finds the coordinates of the intersection of the line segment from `a` to `b` 
with the line segment from `c` to `d`.
"""
function segment_intersection_coordinates(a, b, c, d)
    ax, ay = _getxy(a)
    bx, by = _getxy(b)
    cx, cy = _getxy(c)
    dx, dy = _getxy(d)
    num = (cx - ax) * (dy - ay) - (cy - ay) * (dx - ax)
    den = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)
    α = num / den
    return ax + α * (bx - ax), ay + α * (by - ay)
end

"""
    intersection_of_edge_and_bisector_ray(a, b, c)

Given an edge `(a, b)` and a ray emanating from `c` perpendicular  
with the edge and collinear with its midpoint, tests if `c` intersects 
the edge, and if so, returns the `(cert, p)`, where `p` is the intersection point (which is the midpoint)
and `c` is the position of `c` relative to `(a, b)`.
If there is no intersection, `p = (NaN, NaN)` is returned (together with `cert`). The ray should be directed
to the left of the edge.
"""
function intersection_of_edge_and_bisector_ray(a, b, c)
    cert = point_position_relative_to_line(a, b, c)
    if !is_left(cert)
        ax, ay = _getxy(a)
        bx, by = _getxy(b)
        m = (ax + bx) / 2, (ay + by) / 2
        return cert, m
    else
        F = number_type(a)
        return cert, (F(NaN), F(NaN))
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
    if any(isinf, a) || any(isinf, b) || any(isinf, c) || any(isinf, d)
        return Cert.None, Cert.None, Cert.None, (NaN, NaN)
    end
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

"""
    intersection_of_ray_with_edge(p, q, a, b)

Given a ray starting at `p` and in the direction of `q` out to 
infinity, finds the intersection of the ray with the edge from `a` to `b`.
If no such intersection exists, then the coordinates of the intersection are returned 
as `(NaN, NaN)`.
"""
function intersection_of_ray_with_edge(p, q, a, b)
    p1, p2 = _getxy(p)
    q1, q2 = _getxy(q)
    a1, a2 = _getxy(a)
    b1, b2 = _getxy(b)
    den = a2 * q1 - a2 * p1 - b1 * p2 + b2 * p1 + b1 * q2 - b2 * q1 + a1 * (p2 - q2)
    t = -(a1 * b2 - a2 * b1 - a1 * p2 + a2 * p1 + b1 * p2 - b2 * p1) / den
    u = (a2 * q1 - a2 * p1 + p1 * q2 - p2 * q1 + a1 * (p2 - q2)) / den
    if t < 0 || u < 0 || u > 1
        return (NaN, NaN)
    else
        return (p1 + t * (q1 - p1), p2 + t * (q2 - p2))
    end
end

"""
    liang_barsky(a, b, c, d, p, q)

Applies the Liang-Barsky algorithm to find the intersection of the line segment from `p` to `q`
with the rectangle from `[a, b] × [c, d]`. Returns `(u, v)`, where:

- If there is an intersection, then `u` and `v` are the coordinates of the two intersections. 
- If there is no intersection, then `u = v = (NaN, NaN)`.
"""
function liang_barsky(a, b, c, d, p, q)
    t1 = 0.0
    t2 = 1.0
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    Δx = qx - px
    t1, t2, inside = liang_barsky_clipper(-Δx, px - a, t1, t2)
    if inside
        t1, t2, inside = liang_barsky_clipper(Δx, b - px, t1, t2)
        if inside
            Δy = qy - py
            t1, t2, inside = liang_barsky_clipper(-Δy, py - c, t1, t2)
            if inside
                t1, t2, inside = liang_barsky_clipper(Δy, d - py, t1, t2)
                if inside
                    if t2 < 1
                        qx = px + t2 * Δx
                        qy = py + t2 * Δy
                    end
                    if t1 > 0
                        px = px + t1 * Δx
                        py = py + t1 * Δy
                    end
                end
            end
        end
    end
    if inside
        return (px, py), (qx, qy)
    else
        return (NaN, NaN), (NaN, NaN)
    end
end
function liang_barsky_clipper(p, q, t1, t2)
    inside = true
    if p < 0
        r = q / p
        if r > t2
            inside = false
        elseif r > t1
            t1 = r
        end
    elseif p > 0
        r = q / p
        if r < t1
            inside = false
        elseif r < t2
            t2 = r
        end
    elseif q < 0
        inside = false
    end
    return t1, t2, inside
end
