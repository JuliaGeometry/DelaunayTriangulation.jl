"""
    is_boundary_index(i::I) where {I}

Given an index `i`, returns `i ≤ I(BoundaryIndex)`. 
"""
is_boundary_index(i::I) where {I<:Integer} = i ≤ I(BoundaryIndex)
is_boundary_index(i) = false # In case we provide a point instead of an integer

"""
    is_boundary_edge(ij, adj::Adjacent)
    is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}

Tests if the edge `(i, j)` is a boundary edge, meaning `get_adjacent(adj, i, j)` 
is a boundary index.

!!! notes 

    The orientation of `(i, j)` is important: even if `(i, j)` is an edge 
    on the boundary, if there is a triangle `(i, j, k)` in the triangulation then `(i, j)`
    is not a boundary edge but `(j, i)` would be.
"""
is_boundary_edge(ij, adj::Adjacent) = is_boundary_index(get_adjacent(adj, ij))
function is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj)
end

"""
    is_boundary_triangle(i, j, k, adj)
    is_boundary_triangle(T, adj)

Given a triangle `T = (i, j, k)` and an adjacent map `adj`, returns `true` if `T` is a boundary triangle. 

!!! notes 

    A boundary triangle is still part of the triangulation, but it has at least one edge that
    forms part of the boundary (so that at least one of the edges `(u, v)` satisfies 
    `is_boundary_edge(v, u, adj)`). This is similar to, but not the same as, a ghost triangle.
"""
function is_boundary_triangle(i, j, k, adj)
    for (u, v) in triangle_edges(i, j, k)
        is_boundary_edge(v, u, adj) && return true
    end
    return false
end
function is_boundary_triangle(T, adj)
    for (u, v) in triangle_edges(T)
        is_boundary_edge(v, u, adj) && return true
    end
    return false
end

"""
    is_ghost_edge(i, j)

Given an edge `(i, j)`, returns `true` if `(i, j)` is a ghost edge. 

!!! notes

    A ghost edge is an edge in which either `is_boundary_index(i)` or `is_boundary_index(j)` is true.
"""
is_ghost_edge(i, j) = is_boundary_index(i) || is_boundary_index(j)
is_ghost_edge(ij) = is_ghost_edge(initial(ij), terminal(ij))

"""
    is_ghost_triangle(i, j, k)
    is_ghost_triangle(T)

Given a triangle `T = (i, j, k)`, returns `true` if `T` is a 
ghost triangle and `false` otherwise. 

!!! notes 

    A ghost triangle is one in which any of the vertices `(i, j, k)`
    are a boundary index, as tested via [`is_boundary_index`](@ref).
"""
is_ghost_triangle(i, j, k) = any(is_boundary_index, (i, j, k))
is_ghost_triangle(T) = is_ghost_triangle(geti(T), getj(T), getk(T))

"""
    is_interior_curve(i)
    is_interior_curve(i, boundary_map)

Given an index `i`, tests if the curve is an interior curve, i.e. if `i > 1`. Alternatively, 
if a `boundary_map` is provided from [`construct_boundary_map`](@ref), `i` should be a boundary map so that `is_interior_curve(j)`
is tested, where `j = get_curve_index(boundary_map, i)`.
"""
is_interior_curve(i) = i > 1 # Note that we can't just test e.g. if the boundary index is BoundaryIndex, since the outer-most boundary could have multiple segments
is_interior_curve(i, boundary_map) = is_interior_curve(get_curve_index(boundary_map, i))

"""
    is_outer_boundary_index(i, boundary_map)

Given an index `i` and a `boundary_map` from [`construct_boundary_map`](@ref), tests if the index is a boundary index referring to the outermost boundary.
"""
is_outer_boundary_index(i, boundary_map) = is_boundary_index(i) && !is_interior_curve(i, boundary_map)

"""
    is_outer_ghost_triangle(i, j, k, boundary_map)

Given a ghost triangle `(i, j, k)` and a boundary map `boundary_map` from [`construct_boundary_map`](@ref),
tests if the ghost triangle is on the outermost boundary (`true`) or on an interior boundary (`false`).
"""
function is_outer_ghost_triangle(i, j, k, boundary_map)
    if is_ghost_triangle(i, j, k)
        ℓ = get_boundary_index(i, j, k)
        is_interior = is_interior_curve(ℓ, boundary_map)
        return !is_interior
    else
        return false
    end
end

"""
    is_outer_ghost_edge(i, j, boundary_map)

Given a ghost edge `(i, j)` and a boundary map `boundary_map` from [`construct_boundary_map`](@ref),
tests if the ghost edge is attached to the outermost boundary (`true`) or on an interior boundary (`false`).
"""
function is_outer_ghost_edge(i, j, boundary_map)
    if is_ghost_edge(i, j)
        ℓ = get_boundary_index(i, j)
        is_interior = is_interior_curve(ℓ, boundary_map)
        return !is_interior
    else
        return false
    end
end

"""
    is_outer_boundary_node(i, graph::Graph{I}, boundary_index_ranges) where {I}

Tests if `i` is a node appearing on the outermost boundary. 

# Arguments 
- `i`: The node to test. 
- `graph::Graph`: The graph. 
- `boundary_index_ranges`: A dictionary from [`construct_boundary_index_ranges`](@ref).

# Outputs 
- `is_outer_boundary_node`: A Boolean indicating whether `i` is a node on the outermost boundary.

See also [`is_boundary_node`](@ref).
"""
function is_outer_boundary_node(i, graph::Graph{I}, boundary_index_ranges) where {I}
    outer_range = map_boundary_index(boundary_index_ranges, I(BoundaryIndex))
    for boundary_index in outer_range
        i ∈ get_neighbours(graph, boundary_index) && return true
    end
    return false
end

"""
    is_boundary_node(i, graph::Graph{I}, boundary_index_ranges) where {I}

Tests if `i` is a node appearing on the boundary. 

# Arguments 
- `i`: The node to test. 
- `graph::Graph`: The graph. 
- `boundary_index_ranges`: A dictionary from [`construct_boundary_index_ranges`](@ref).

# Outputs 
- `is_outer_boundary_node`: A Boolean indicating whether `i` is a node on the boundary.
- `boundary_index`: The boundary index of the boundary to which `i` belongs. If there is no such boundary, `boundary_index = I(DefaultAdjacentValue)`.

See also [`is_outer_boundary_node`](@ref).
"""
function is_boundary_node(i, graph::Graph{I}, boundary_index_ranges) where {I}
    for boundary_index_range in values(boundary_index_ranges)
        for boundary_index in boundary_index_range
            i ∈ get_neighbours(graph, boundary_index) && return (true, boundary_index)
        end
    end
    return (false, I(DefaultAdjacentValue))
end

"""
    edge_exists(i::I) where {I}

Returns `i ≠ I(DefaultAdjacentValue)`.
"""
edge_exists(i::I) where {I} = i ≠ I(DefaultAdjacentValue)

"""
    edge_exists(i, j, adj::Adjacent{I,E}) where {I,E}
    edge_exists(ij, adj::Adjacent{I,E}) where {I,E}

Tests if the edge `ij = (i, j)` exists in the triangulation corresponding to the [`Adjacent`](@ref) map `adj`.

See also [`edge_exists(::I) where I`](@ref).
"""
function edge_exists(ij, adj::Adjacent{I,E}) where {I,E}
    k = get_adjacent(adj, ij)
    return edge_exists(k)
end
function edge_exists(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return edge_exists(ij, adj)
end

"""
    has_ghost_triangles(adj::Adjacent{I,E}, adj2v) where {I,E}

Tests if the triangle represented by the [`Adjacent`](@ref) map `adj` and the [`Adjacent2Vertex`](@ref) map `adj2v` contains ghost triangles.
"""
function has_ghost_triangles(adj::Adjacent{I,E}, adj2v) where {I,E}
    I(BoundaryIndex) ∉ keys(get_adjacent2vertex(adj2v)) && return false
    outer_boundary_edges = get_adjacent2vertex(adj2v, I(BoundaryIndex))
    is_empty(outer_boundary_edges) && return false
    e = first(each_edge(outer_boundary_edges))
    return edge_exists(terminal(e), I(BoundaryIndex), adj)
end
