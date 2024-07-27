"""
    convex_hull(points; predicates::AbstractPredicateKernel=AdaptiveKernel(), IntegerType::Type{I}=Int) where {I} -> ConvexHull 

Computes the convex hull of `points`. The monotone chain algorithm is used.

# Arguments 
- `points`: The set of points.

# Keyword Arguments 
- `IntegerType=Int`: The integer type to use for the vertices.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
- `ch`: The [`ConvexHull`](@ref). 
"""
function convex_hull(points; predicates::AbstractPredicateKernel=AdaptiveKernel(), IntegerType::Type{I}=Int) where {I}
    ch = ConvexHull(points, I[])
    sizehint!(ch, num_points(points))
    return convex_hull!(ch; predicates)
end

"""
    convex_hull!(ch::ConvexHull{P,I}; predicates::AbstractPredicateKernel=AdaptiveKernel()) where {P,I}

Using the points in `ch`, computes the convex hull in-place. 

The `predicates` keyword argument determines how predicates are computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.

See also [`convex_hull`](@ref).
"""
function convex_hull!(ch::ConvexHull{P,I}; predicates::AbstractPredicateKernel=AdaptiveKernel()) where {P,I}
    indices = get_vertices(ch)
    points = get_points(ch)
    empty!(indices)
    n = num_points(points)
    sizehint!(indices, n)
    insertion_order = lexicographic_order(points)
    unique!(i -> get_point(points, i), insertion_order)
    if n == 1
        push!(indices, insertion_order[begin])
        return ch
    elseif n == 2
        push!(indices, insertion_order[begin], insertion_order[begin+1], insertion_order[begin])
        return ch
    elseif n == 3
        i, j, k = construct_positively_oriented_triangle(NTuple{3,I}, insertion_order[begin], insertion_order[begin+1], insertion_order[begin+2], points, predicates)
        push!(indices, i, j, k, i)
        return ch
    end
    lower = I[]
    upper = I[]
    sizehint!(lower, max(4, floor(I, cbrt(n))))
    sizehint!(upper, max(4, floor(I, cbrt(n))))
    for i in eachindex(insertion_order)
        while length(upper) ≥ 2 && is_left(point_position_relative_to_line(predicates, get_point(points, upper[end-1]), get_point(points, upper[end]), get_point(points, insertion_order[i])))
            pop!(upper)
        end
        push!(upper, insertion_order[i])
    end
    for i in lastindex(insertion_order):-1:firstindex(insertion_order)
        while length(lower) ≥ 2 && is_left(point_position_relative_to_line(predicates, get_point(points, lower[end-1]), get_point(points, lower[end]), get_point(points, insertion_order[i])))
            pop!(lower)
        end
        push!(lower, insertion_order[i])
    end
    popfirst!(lower)
    append!(upper, lower)
    reverse!(upper) # counter-clockwise
    append!(indices, upper)
    unique!(indices)
    push!(indices, indices[begin])
    return ch
end