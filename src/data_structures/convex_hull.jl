"""
    ConvexHull{P,I}

Struct storing the results for a convex hull.

# Fields 
- `points::P`

The complete set of points.
- `indices::I`

Indices of points in `points` corresponding to the convex hull, in counter-clockwise order, 
and `indices[begin] == indices[end]`.
"""
struct ConvexHull{P,I}
    points::P
    indices::I
end
get_points(ch::ConvexHull) = ch.points
get_indices(ch::ConvexHull) = ch.indices
num_points(ch::ConvexHull) = num_points(get_points(ch))

function Base.:(==)(ch::ConvexHull, ch2::ConvexHull)
    get_points(ch) ≠ get_points(ch2) && return false
    i1 = get_indices(ch)
    i2 = get_indices(ch2)
    return circular_equality(i1, i2)
end
function Base.show(io::IO, m::MIME"text/plain", ch::ConvexHull)
    println(io, "Convex hull.")
    println(io, "    Indices:")
    show(io,m,get_indices(ch))
end

"""
    convex_hull(points; IntegerType::Type{I}=Int64) where {I}

Computes the convex hull of `points` using Graham's scan. Returns a [`ConvexHull`](@ref) object.

Note that if there are a trio of points on the convex hull that are collinear, they will 
all be included, instead of only taking the endpoints of the collinear points.
"""
function convex_hull(points; IntegerType::Type{I}=Int64) where {I}
    num_points(points) ≤ 2 && throw("Need at least 3 points to compute a convex hull.")
    ch = ConvexHull(points, I[])
    sizehint!(ch, num_points(points))
    convex_hull!(ch)
    return ch
end
function convex_hull!(ch::ConvexHull)
    indices = get_indices(ch)
    points = get_points(ch)
    num_points(points) ≤ 2 && throw("Need at least 3 points to compute a convex hull.")
    empty!(indices)
    sizehint!(indices, num_points(ch))
    point_order = lexicographic_order(points)
    u = point_order[begin]
    v = point_order[begin + 1]
    p, q = get_point(points, u, v)
    upper = [u, v]
    sizehint!(upper, num_points(points) ÷ 2)
    for i in (firstindex(point_order) + 2):lastindex(point_order)
        p, q, u, v = _add_to_hull!(upper, i, points, point_order, p, q, u, v)
    end
    u = point_order[end]
    v = point_order[end - 1]
    p, q = get_point(points, u, v)
    lower = [u, v]
    sizehint!(lower, num_points(points) ÷ 2)
    for i in (lastindex(point_order) - 2):-1:firstindex(point_order)
        p, q, u, v = _add_to_hull!(lower, i, points, point_order, p, q, u, v)
    end
    popfirst!(lower)
    append!(upper, lower)
    reverse!(upper) # counter-clockwise
    append!(indices, upper)
    return nothing
end

function _add_to_hull!(hull, i, points, point_order, p, q, u, v)
    w = point_order[i]
    r = get_point(points, w)
    push!(hull, w)
    turn_cert = point_position_relative_to_line(p, q, r)
    n = length(hull)
    while n > 2 && is_left(turn_cert)
        popat!(hull, n - 1)
        n -= 1
        if n > 2
            v = u
            q = p
            u = hull[n - 2]
            p = get_point(points, u)
            turn_cert = point_position_relative_to_line(p, q, r)
        end
    end
    u = v
    p = q
    v = w
    q = r
    return p, q, u, v
end

Base.sizehint!(ch::ConvexHull, n) = Base.sizehint!(get_indices(ch), n)
