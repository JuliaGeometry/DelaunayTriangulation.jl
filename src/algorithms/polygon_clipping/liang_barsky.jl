"""
    liang_barsky(a, b, c, d, p, q) -> (Point, Point)

Applies the Liang-Barsky algorithm to find the intersection of the line segment `pq` with the rectangle 
`[a, b] × [c, d]`.

# Arguments 
- `p`: The first point of the line segment.
- `q`: The second point of the line segment.
- `a`: The minimum x-coordinate of the rectangle.
- `b`: The maximum x-coordinate of the rectangle.
- `c`: The minimum y-coordinate of the rectangle.
- `d`: The maximum y-coordinate of the rectangle.

# Output
- `u`: The first coordinate of the intersection, or `(NaN, NaN)` if there is no intersection. 
- `v`: The second coordinate of the intersection, or `(NaN, NaN)` if there is no intersection.
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
