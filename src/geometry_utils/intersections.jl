"""
    intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)

Finds the intersection of the ray through `p` and `q`, oriented from `p` to `q`,
with the boundary of the polygon defined by `points` and `boundary_nodes`. 

It is assumed that `p` is inside the polygon, but `q` could be outside or inside.

Currently, this function has only been tested on rectangular boundaries.
"""
function intersection_of_ray_with_boundary(points, boundary_nodes, p, q, tol=1e-9)
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
    while sign(δ2) == 1
        t2 *= 2
        r = px + t2 * (qx - px), py + t2 * (qy - py)
        δ2 = distance_to_polygon(r, points, boundary_nodes)
    end

    ## Perform bisection 
    t = (t1 + t2) / 2
    r = px + t * (qx - px), py + t * (qy - py)
    δ = distance_to_polygon(r, points, boundary_nodes)
    while (t2 - t1) / 2 > tol || abs(δ) > tol
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
    end
    return r
end

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