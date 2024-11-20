@inline function get_ghost_vertex(i, j, k)
    is_ghost_vertex(i) && return i
    is_ghost_vertex(j) && return j
    return k
end
@inline function get_ghost_vertex(i, j)
    is_ghost_vertex(i) && return i
    return j
end

@inline function replace_boundary_triangle_with_ghost_triangle(tri::Triangulation, V)
    u, v, w = triangle_vertices(V)
    is_boundary_edge(tri, u, v) && return (v, u, get_adjacent(tri, v, u))
    is_boundary_edge(tri, v, w) && return (w, v, get_adjacent(tri, w, v))
    return (u, w, get_adjacent(tri, u, w))
end

@inline function replace_ghost_triangle_with_boundary_triangle(tri::Triangulation, V)
    T = sort_triangle(V)
    u, v, w = triangle_vertices(T) # w is the ghost vertex 
    return (v, u, get_adjacent(tri, v, u))
end

@inline is_ghost_vertex(i::Integer) = i ≤ 𝒢
@inline is_ghost_vertex(i) = false

@inline function is_boundary_edge(tri::Triangulation, ij)
    e = reverse_edge(ij)
    k = get_adjacent(tri, e)
    return is_ghost_vertex(k)
end
@inline function is_boundary_edge(tri::Triangulation, i, j)
    ij = construct_edge(edge_type(tri), i, j)
    return is_boundary_edge(tri, ij)
end

@inline function is_boundary_triangle(tri::Triangulation, i, j, k)
    is_boundary_edge(tri, i, j) && return true
    is_boundary_edge(tri, j, k) && return true
    is_boundary_edge(tri, k, i) && return true
    return false
end
@inline function is_boundary_triangle(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    return is_boundary_triangle(tri, i, j, k)
end

@inline is_ghost_edge(i, j) = is_ghost_vertex(i) || is_ghost_vertex(j)
@inline is_ghost_edge(ij) = is_ghost_edge(initial(ij), terminal(ij))

@inline function is_ghost_triangle(i, j, k)
    (is_ghost_vertex(i) || is_ghost_vertex(j) || is_ghost_vertex(k)) && return true
    return false
end
@inline function is_ghost_triangle(T)
    i, j, k = triangle_vertices(T)
    return is_ghost_triangle(i, j, k)
end

@inline function is_exterior_ghost_triangle(tri::Triangulation, i, j, k)
    !is_ghost_triangle(i, j, k) && return false
    ℓ = get_ghost_vertex(i, j, k)
    return is_exterior_ghost_vertex(tri::Triangulation, ℓ)
end

@inline function is_exterior_ghost_edge(tri::Triangulation, i, j)
    !is_ghost_edge(i, j) && return false
    ℓ = get_ghost_vertex(i, j)
    return is_exterior_ghost_vertex(tri, ℓ)
end

@inline function is_exterior_boundary_node(tri::Triangulation, i)
    for g in each_ghost_vertex(tri)
        curve_index = get_curve_index(tri, g)
        is_exterior_curve(tri, curve_index) && i ∈ get_neighbours(tri, g) && return true
    end
    return false
end

@inline function is_boundary_node(tri::Triangulation, i)
    for g in each_ghost_vertex(tri)
        i ∈ get_neighbours(tri, g) && return (true, g)
    end
    I = integer_type(tri)
    return (false, I(∅))
end
