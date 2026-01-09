"""
    is_ghost_vertex(i) -> Bool 

Tests if `i` is a ghost vertex, meaning `i ≤ $𝒢`.
"""
is_ghost_vertex(i::I) where {I <: Integer} = i ≤ I(𝒢)
is_ghost_vertex(i) = false # in case we provide a point instead of an integer 

"""
    is_boundary_edge(tri::Triangulation, ij) -> Bool
    is_boundary_edge(tri::Triangulation, i, j) -> Bool 

Tests if the edge `(i, j)` is a boundary edge of `tri`, meaning `(j, i)` adjoins a ghost vertex.
"""
function is_boundary_edge(tri::Triangulation, ij)
    e = reverse_edge(ij)
    k = get_adjacent(tri, e)
    return is_ghost_vertex(k)
end
function is_boundary_edge(tri::Triangulation, i, j)
    ij = construct_edge(edge_type(tri), i, j)
    return is_boundary_edge(tri, ij)
end

"""
    is_boundary_triangle(tri::Triangulation, T) -> Bool
    is_boundary_triangle(tri::Triangulation, i, j, k) -> Bool

Returns `true` if the triangle `T = (i, j, k)` of `tri` has an edge on the boundary, and `false` otherwise.
"""
function is_boundary_triangle(tri::Triangulation, i, j, k)
    is_boundary_edge(tri, i, j) && return true
    is_boundary_edge(tri, j, k) && return true
    is_boundary_edge(tri, k, i) && return true
    return false
end
function is_boundary_triangle(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    return is_boundary_triangle(tri, i, j, k)
end

"""
    is_ghost_edge(ij) -> Bool 
    is_ghost_edge(i, j) -> Bool

Tests if the edge `(i, j)` is a ghost edge, meaning `i` or `j` is a ghost vertex.
"""
is_ghost_edge(i, j) = is_ghost_vertex(i) || is_ghost_vertex(j)
is_ghost_edge(ij) = is_ghost_edge(initial(ij), terminal(ij))

"""
    is_ghost_triangle(T) -> Bool
    is_ghost_triangle(i, j, k) -> Bool 

Tests if `T = (i, j, k)` is a ghost triangle, meaning `i`, `j` or `k` is a ghost vertex.
"""
function is_ghost_triangle(i, j, k)
    (is_ghost_vertex(i) || is_ghost_vertex(j) || is_ghost_vertex(k)) && return true
    return false
end
function is_ghost_triangle(T)
    i, j, k = triangle_vertices(T)
    return is_ghost_triangle(i, j, k)
end

"""
    is_exterior_ghost_triangle(tri::Triangulation, i, j, k) -> Bool 

Tests if the triangle `(i, j, k)` is an exterior ghost triangle of `tri`.

See also [`is_exterior_ghost_vertex`](@ref).
"""
function is_exterior_ghost_triangle(tri::Triangulation, i, j, k)
    !is_ghost_triangle(i, j, k) && return false
    ℓ = get_ghost_vertex(i, j, k)
    return is_exterior_ghost_vertex(tri::Triangulation, ℓ)
end

"""
    is_exterior_ghost_edge(tri::Triangulation, i, j) -> Bool

Tests if the edge `(i, j)` is an exterior ghost edge of `tri`.

See also [`is_exterior_ghost_vertex`](@ref).
"""
function is_exterior_ghost_edge(tri::Triangulation, i, j)
    !is_ghost_edge(i, j) && return false
    ℓ = get_ghost_vertex(i, j)
    return is_exterior_ghost_vertex(tri, ℓ)
end

"""
    is_exterior_boundary_node(tri::Triangulation, i) -> Bool

Tests if the vertex `i` is an exterior boundary node of `tri`.
"""
function is_exterior_boundary_node(tri::Triangulation, i)
    for g in each_ghost_vertex(tri)
        curve_index = get_curve_index(tri, g)
        is_exterior_curve(tri, curve_index) && i ∈ get_neighbours(tri, g) && return true
    end
    return false
end

"""
    is_boundary_node(tri::Triangulation, i) -> (Bool, Vertex)

Tests if the vertex `i` is a boundary node of `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The vertex to test. 

# Outputs 
- `flag`: `true` if `i` is a boundary node, and `false` otherwise.
- `g`: Either the ghost vertex corresponding with the section that `i` lives on if `flag` is true, or $∅ otherwise.
"""
function is_boundary_node(tri::Triangulation, i)
    for g in each_ghost_vertex(tri)
        i ∈ get_neighbours(tri, g) && return (true, g)
    end
    I = integer_type(tri)
    return (false, I(∅))
end
