@doc """
    ConvexHull{PointsType, IntegerType}
 
Struct for representing a convex hull. See also [`convex_hull`](@ref).

# Fields 
- `points::PointsType`: The underlying point set.
- `vertices::Vector{IntegerType}`: The vertices of the convex hull, in counter-clockwise order. Defined so that `vertices[begin] == vertices[end]`.

# Constructors 
    ConvexHull(points, vertices)
    convex_hull(points; IntegerType=Int)
"""
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
function Base.show(io::IO, m::MIME"text/plain", ch::ConvexHull)
    println(io, "Convex hull.")
    println(io, "   Vertices:")
    v = get_vertices(ch)
    show(io, m, v)
end
Base.sizehint!(ch::ConvexHull, n) = Base.sizehint!(get_vertices(ch), n)


@doc """
    get_points(convex_hull::ConvexHull) -> Points

Returns the underlying point set of `convex_hull`.
"""
get_points(convex_hull::ConvexHull) = convex_hull.points


@doc """
    get_vertices(convex_hull::ConvexHull) -> Vector{Vertices}

Returns the vertices of `convex_hull`. These are given in counter-clockwise order, and 
are defined so that the first and last vertices and equal.
"""
get_vertices(convex_hull::ConvexHull) = convex_hull.vertices


function Base.empty!(ch::ConvexHull)
    vertices = get_vertices(ch)
    empty!(vertices)
    return ch
end