"""
    is_boundary_index(i::I) where {I}

Given an index `i`, returns `i ≤ I(BoundaryIndex)`.
"""
is_boundary_index(i::I) where {I<:Integer} = i ≤ I(BoundaryIndex)
is_boundary_index(i) = false # In case we provide a point instead of an integer

"""
    is_boundary_edge(ij, adj::Adjacent)
    is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}

Given an edge `(i, j)` and an adjacent map `adj`, returns `true` if 
`(i, j)` is a boundary edge and `false` otherwise.

Note that the orientation of `(i, j)` is important: even if `(i, j)` is an edge
on the boundary, if there is a triangle `(i, j, k)` in the triangulation then 
`(i, j)` is not a boundary edge but `(j, i)` would be.
"""
is_boundary_edge(ij, adj::Adjacent) = is_boundary_index(get_adjacent(adj, ij))
function is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj)
end

"""
    is_boundary_triangle(i, j, k, adj)
    is_boundary_triangle(T, adj)

Given a triangle `T = (i, j, k)` and an adjacent map `adj`, 
returns `true` if `T` is a boundary triangle. A boundary triangle 
is still part of the triangulation, but it has at least one edge that 
forms part of the boundary (so that at least one of `is_boundary_edge(j, i)`,
`is_boundary_edge(k, j)`, and `is_boundary_edge(i, k)` is `true`).
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

Given an edge `(i, j)`, returns `true` if `(i, j)` is a 
ghost edge. A ghost edge is an edge in which either `is_boundary_index(i)`
or `is_boundary_index(j)` is true.
"""
is_ghost_edge(i, j) = is_boundary_index(i) || is_boundary_index(j)
is_ghost_edge(ij) = is_ghost_edge(initial(ij), terminal(ij))

"""
    is_ghost_triangle(i, j, k)
    is_ghost_triangle(T)

Given a triangle `T = (i, j, k)`, returns `true` 
if `T` is a ghost triangle and `false` otherwise. A ghost 
triangle is one in which any of the vertices `(i, j, k)` 
are a boundary index, as tested via [`is_boundary_index`](@ref).
"""
function is_ghost_triangle(i, j, k)
    return any(is_boundary_index, (i, j, k))
end
function is_ghost_triangle(T)
    i, j, k = indices(T)
    return is_ghost_triangle(i, j, k)
end

"""
    is_interior_curve(i)
    is_interior_curve(i, boundary_map)

Given an index `i`, tests if the curve is an interior curve, i.e. if `i > 1`. If 
a map `boundary_map` is provided, `i` should be a boundary map so that `is_interior_curve(j)`
is tested, where `j = get_curve_index(boundary_map, i)`.
"""
is_interior_curve(i) = i > 1 # Note that we can't just test e.g. if the boundary index is BoundaryIndex, since the outer-most boundary could have multiple segments
is_interior_curve(i, boundary_map) = is_interior_curve(get_curve_index(boundary_map, i))

"""
    is_outer_boundary_index(i, boundary_map)

Given an index `i`, tests if the index is a boundary index referring to the outermost boundary, making 
use of the `boundary_map` from [`construct_boundary_map`](@ref).
"""
function is_outer_boundary_index(i, boundary_map)
    return is_boundary_index(i) && !is_interior_curve(i, boundary_map)
end

"""
    is_outer_ghost_triangle(i, j, k, boundary_map) 

Given a ghost triangle `(i, j, k)` and a boundary map `boundary_map`
taking boundary indices to their location in the boundary node array, 
tests if `(i, j, k)` is a ghost triangle on the outermost boundary 
(`true`) or on an interior boundary (`false`)`.
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

Given a ghost edge `(i, j)` and a boundary map `boundary_map`
taking boundary indices to their location in the boundary node array, 
tests if `(i, j)` is a ghost edge on the outermost boundary 
(`true`) or on an interior boundary (`false`)`.
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

Given a node index `i`, a `graph::Graph`, a `Dict` from [`construct_boundary_index_ranges`](@ref),
returns true if `i` corresponds to a node on the outermost boundary, 
and `false` otherwise.
"""
function is_outer_boundary_node(i, graph::Graph{I}, boundary_index_ranges) where {I}
    outer_range = map_boundary_index(boundary_index_ranges, I(BoundaryIndex))
    for boundary_index in outer_range
        i ∈ get_neighbours(graph, boundary_index) && return true
    end
    return false
end

"""
    edge_exists(i::I) where {I}

Returns `i ≠ I(DefaultAdjacentValue)`.
"""
edge_exists(i::I) where {I} = i ≠ I(DefaultAdjacentValue)

"""
    edge_exists(ij, adj::Adjacent{I,E}) where {I,E}

Given an edge `ij` and an [`Adjacent`](@ref) map `adj`, 
tests if the edge exists in the corresponding triangulation.
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

Given an [`Adjacent`](@ref) map `adj` and an [`Adjacent2Vertex`](@ref) map `adj2v`,
tests if the corresponding triangulation contains ghost triangles.
"""
function has_ghost_triangles(adj::Adjacent{I,E}, adj2v) where {I,E}
    outer_boundary_edges = get_adjacent2vertex(adj2v, I(BoundaryIndex))
    e = first(each_edge(outer_boundary_edges))
    return edge_exists(terminal(e), I(BoundaryIndex), adj)
end
