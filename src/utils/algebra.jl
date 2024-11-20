@inline norm(p) = hypot(getx(p), gety(p))
@inline norm_sqr(p) = getx(p)^2 + gety(p)^2
@inline dist(p, q) = norm(getxy(p) .- getxy(q))
@inline dist_sqr(p, q) = norm_sqr(getxy(p) .- getxy(q))

@inline function midpoint(x::Number, y::Number)
    return ifelse(flag, (x / 2) + (y / 2), (x + y) / 2)
end

@inline midpoint(p, q) = midpoint.(getxy(p), getxy(q))