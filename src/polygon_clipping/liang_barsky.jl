struct LiangBarskyClipper{T}
    p::T 
    q::T
    t1::T 
    t2::T
end 

function liang_barsky(a, b, c, d, p, q)
    F = number_type(p)
    t1 = zero(F)
    t2 = one(F)
    px, py = getxy(p)
    qx, qy = getxy(q)
    Δx = qx - px
    clipper = LiangBarskyClipper(-Δx, px - a, t1, t2)
    t1, t2, inside = liang_barsky_clipper(clipper)
    if inside
        clipper = LiangBarskyClipper(Δx, b - px, t1, t2)
        t1, t2, inside = liang_barsky_clipper(clipper)
        if inside
            Δy = qy - py
            clipper = LiangBarskyClipper(-Δy, py - c, t1, t2)
            t1, t2, inside = liang_barsky_clipper(clipper)
            if inside
                clipper = LiangBarskyClipper(Δy, d - py, t1, t2)
                t1, t2, inside = liang_barsky_clipper(clipper)
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
        return getxy(p), getxy(q)
    else
        return (F(NaN), F(NaN)), (F(NaN), F(NaN))
    end
end
function liang_barsky_clipper(clipper::LiangBarskyClipper)
    p, q = get_p(clipper), get_q(clipper)
    t1, t2 = get_t1(clipper), get_t2(clipper)
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
