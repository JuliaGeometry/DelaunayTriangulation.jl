# Interface
@inline geti(T) = T[1]
@inline getj(T) = T[2]
@inline getk(T) = T[3]

# Derived
@inline triangle_vertices(T) = (geti(T), getj(T), getk(T))

@inline triangle_edges(i, j, k) = ((i, j), (j, k), (k, i))
@inline triangle_edges(T) = triangle_edges(geti(T), getj(T), getk(T))

@inline function sort_triangle(T)
    i, j, k = triangle_vertices(T)
    minijk = min(i, j, k)
    return minijk == i ? (j, k, i) : minijk == j ? (k, i, j) : (i, j, k)
end
@inline sort_triangle(i::Integer, j::Integer, k::Integer) = sort_triangle((i, j, k))
@inline function sort_triangle(i, j, k) # in case one of the vertices is a point rather than an integer 
    # Note that the return type is a 3-Union 
    return is_ghost_vertex(i) ? (j, k, i) : is_ghost_vertex(j) ? (k, i, j) : (i, j, k)
end

@inline rotate_triangle(T, ::Val{1}) = (getj(T), getk(T), geti(T))
@inline rotate_triangle(T, ::Val{2}) = (getk(T), geti(T), getj(T))

@inline function compare_triangles(T, V)
    return T == V || T == rotate_triangle(V, Val(1)) || T == rotate_triangle(V, Val(2))
end

@inline function construct_positively_oriented_triangle(i, j, k, points, predicates::AbstractPredicateKernel=AdaptiveKernel())
    p, q, r = get_point(points, i, j, k)
    orientation = triangle_orientation(predicates, p, q, r)
    if is_negatively_oriented(orientation)
        return (j, i, k)
    else
        return (i, j, k)
    end
end
