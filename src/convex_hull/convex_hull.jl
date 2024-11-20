function convex_hull(points; predicates::AbstractPredicateKernel = AdaptiveKernel(), IntegerType::Type{I}=Int) where {I}
    ch = ConvexHull(points, I[])
    return convex_hull!(ch; predicates)
end

function convex_hull!(ch::ConvexHull{P, I}; predicates::AbstractPredicateKernel = AdaptiveKernel())
    empty!(ch)
    n = num_points(get_points(ch))
    sizehint!(get_vertices(ch), n)
    vertices, points = get_vertices(ch), get_points(ch)
    insertion_order = lexicographic_order(points)
    unique!(i -> get_point(points, i), insertion_order)
    if n == 1 
        push!(vertices, insertion_order[begin])
        return ch
    elseif n == 2
        push!(vertices, insertion_order[begin], insertion_order[begin + 1], insertion_order[begin])
        return ch
    elseif n == 3
        i, j, k = construct_positively_oriented_triangle(insertion_order[begin], insertion_order[begin + 1], insertion_order[begin + 2], points, predicates)
        push!(vertices, i, j, k, i)
        return ch
    end
    lower, upper = I[], I[]
    sizehint!(lower, max(4, floor(I, cbrt(n))))
    sizehint!(upper, max(4, floor(I, cbrt(n))))
    compute_upper_hull!(upper, insertion_order, points, predicates)
    compute_lower_hull!(lower, insertion_order, points, predicates)
    merge_hulls!(vertices, upper, lower)
    return ch
end

function compute_upper_hull!(upper, insertion_order, points, predicates)
    for i in eachindex(insertion_order)
        while length(upper) ≥ 2 && is_left(point_position_relative_to_line(predicates, get_point(points, upper[end - 1]), get_point(points, upper[end]), get_point(points, insertion_order[i])))
            pop!(upper)
        end
        push!(upper, insertion_order[i])
    end
    return upper
end

function compute_lower_hull!(lower, insertion_order, points, predicates)
    for i in lastindex(insertion_order):-1:firstindex(insertion_order)
        while length(lower) ≥ 2 && is_left(point_position_relative_to_line(predicates, get_point(points, lower[end - 1]), get_point(points, lower[end]), get_point(points, insertion_order[i])))
            pop!(lower)
        end
        push!(lower, insertion_order[i])
    end
    return lower
end

function merge_hulls!(vertices, upper, lower)
    popfirst!(lower)
    append!(upper, lower)
    reverse!(upper) # counter-clockwise
    append!(vertices, upper)
    unique!(vertices)
    push!(vertices, vertices[begin])
    return vertices
end