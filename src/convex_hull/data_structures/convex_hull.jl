struct ConvexHull{P, I}
    points::P 
    vertices::Vector{I}
end

function Base.:(==)(ch::ConvexHull, ch2::ConvexHull)
    p1 = get_points(ch)
    p2 = get_points(ch2)
    v1 = get_vertices(ch)
    v2 = get_vertices(ch2)
    p1 ≠ p2 && return false
    length(v1) ≠ length(v2) && return false
    return circular_equality(v1, v2)
end

@inline Base.sizehint!(ch::ConvexHull, n) = Base.sizehint!(get_vertices(ch), n)

@inline Base.copy(ch::ConvexHull) = ConvexHull(copy(get_points(ch)), copy(get_vertices(ch)))

@inline get_points(convex_hull::ConvexHull) = convex_hull.points

@inline get_vertices(convex_hull::ConvexHull) = convex_hull.vertices

function Base.empty!(ch::ConvexHull)
    vertices = get_vertices(ch)
    empty!(vertices)
    return ch
end