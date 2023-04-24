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
    length(i1) ≠ length(i2) && return false
    return circular_equality(i1, i2)
end
function Base.show(io::IO, m::MIME"text/plain", ch::ConvexHull)
    println(io, "Convex hull.")
    println(io, "    Indices:")
    show(io, m, get_indices(ch))
end

"""
    convex_hull(points; IntegerType::Type{I}=Int64) where {I}

Computes the convex hull of `points` using Graham's scan. Returns a [`ConvexHull`](@ref) object.

Note that if there are a trio of points on the convex hull that are collinear, they will 
all be included, instead of only taking the endpoints of the collinear points.
"""
function convex_hull(points; IntegerType::Type{I}=Int64) where {I}
    ch = ConvexHull(points, I[])
    sizehint!(ch, num_points(points))
    convex_hull!(ch)
    return ch
end
function convex_hull!(ch::ConvexHull{P,A}) where {P,I,A<:AbstractVector{I}}
    indices = get_indices(ch)
    points = get_points(ch)
    empty!(indices)
    sizehint!(indices, num_points(ch))
    point_order = lexicographic_order(points)
    n = num_points(ch)
    if n == 1
        push!(indices, point_order[begin])
        return nothing
    elseif n == 2
        push!(indices, point_order[begin], point_order[begin+1], point_order[begin])
        return nothing
    elseif n == 3
        i, j, k = construct_positively_oriented_triangle(NTuple{3,I}, point_order[begin], point_order[begin+1], point_order[begin+2], points)
        push!(indices, i, j, k, i)
        return nothing
    end
    lower = I[]
    upper = I[]
    sizehint!(lower, num_points(points) ÷ 2)
    sizehint!(upper, num_points(points) ÷ 2)
    for i in eachindex(point_order)
        while length(upper) ≥ 2 && is_left(point_position_relative_to_line(get_point(points, upper[end-1]), get_point(points, upper[end]), get_point(points, point_order[i])))
            pop!(upper)
        end
        push!(upper, point_order[i])
    end
    for i in lastindex(point_order):-1:firstindex(point_order)
        while length(lower) ≥ 2 && is_left(point_position_relative_to_line(get_point(points, lower[end-1]), get_point(points, lower[end]), get_point(points, point_order[i])))
            pop!(lower)
        end
        push!(lower, point_order[i])
    end
    popfirst!(lower)
    append!(upper, lower)
    reverse!(upper) # counter-clockwise
    append!(indices, upper)
    unique!(indices)
    push!(indices,indices[begin])
    return nothing
end

Base.sizehint!(ch::ConvexHull, n) = Base.sizehint!(get_indices(ch), n)