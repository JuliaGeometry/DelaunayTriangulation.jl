"""
is_boundary_edge(i, j, adj::Adjacent{I, E}) where {I, E}
is_boundary_edge(i, j, adj2v::Adjacent{I, E}) where {I, E}

Returns `true` if the edge `(i, j)` is a boundary edge of the triangulation, 
and `false` otherwise.
"""
function is_boundary_edge(ij, adj::Adjacent{I,E}) where {I,E}
    return get_edge(adj, ij) == I(BoundaryIndex)
end
function is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj)
end
function is_boundary_edge(ij, adj2v::Adjacent2Vertex{I,Es,E}) where {I,Es,E}
    return ij ∈ get_edge(adj2v, I(BoundaryIndex))
end
function is_boundary_edge(i, j, adj2v::Adjacent2Vertex{I,Es,E}) where {I,Es,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj2v)
end

"""
    edge_exists(i, j, adj::Adjacent{I, E}) where {I, E}

Returns `true` if the edge `(i, j)` is an edge in the triangulation, 
and `false` otherwise.
"""
function edge_exists(ij, adj::Adjacent{I,E}) where {I,E}
    return get_edge(adj, ij) ≠ I(DefaultAdjacentValue)
end
function edge_exists(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return edge_exists(ij, adj)
end

"""
    choose_uvw(e1, e2, e3, i, j, k)

Chooses values for `(u, v, w)` based on the Booleans 
`(e1, e2, e3)`, assuming only one is true. 

- If `e1`, returns `(i, j, k)`
- If `e2`, returns `(j, k, i)`
- If `e3`, returns `(k, i, j)`
"""
function choose_uvw(e1, e2, e3, i, j, k)
    if e1
        return i, j, k
    elseif e2
        return j, k, i
    elseif e3
        return k, i, j
    end
end

function is_delaunay(adj, pts)
    tri_edges = edges(adj)
    for (i, j) in tri_edges
        edge_exists(i, j, adj) && !is_boundary_edge(i, j, adj) && !is_boundary_edge(j, i, adj) && !islegal(i, j, adj, pts) && return false
    end
    return true
end
function validate_triangulation(T, adj::Adjacent{I,E}, adj2v, DG, pts) where {I,E}
    ## Check triangle orientation 
    for T in T
        isoriented(T, pts) == -1 && return false
    end
    ## Test that all edges are legal 
    !is_delaunay(adj, pts) && return false
    ## Check that the adjacent map and the adjacent-to-vertex map edge 
    for (w, S) in adjacent2vertex(adj2v)
        for (i, j) in S
            get_edge(adj, i, j) ≠ w && return false
        end
    end
    ## Check the graph 
    for (i, j) in graph(DG).E
        (i, j) ∉ edges(adj) && return false
        (j, i) ∉ edges(adj) && return false
    end
    ## Done 
    return true
end

"""
    clear_empty_keys!(adj::Adjacent, DG::DelaunayGraph)

Deletes all the keys `(i, j)` in `adj` such that `get_edge(adj, i, j) = $(DefaultAdjacentValue)`.
"""
function clear_empty_keys!(adj::Adjacent{I,E}, DG::DelaunayGraph) where {I,E}
    num_edges_unoriented = length(graph(DG).E)
    all_edges = Vector{E}(undef, 2num_edges_unoriented)
    k = 1
    for (i, j) in graph(DG).E
        all_edges[k] = construct_edge(E, i, j)
        all_edges[num_edges_unoriented+k] = construct_edge(E, j, i)
        k += 1
    end
    invalid_edges = setdiff(edges(adj), all_edges)
    for (i, j) in invalid_edges 
        delete_edge!(adj, i, j)
    end
    return nothing
end

"""
    compare_triangle_sets(T, V)

Compares the collections of triangles `T` and `V`, checking for equality. This function is 
needed to check for equality even with circularly-shifted indices in the triangles. Returns `true` 
if the two collections are equivalent, and false otherwise.
"""
function compare_triangle_sets(T, V)
    length(T) ≠ length(V) && return false 
    for T in T 
        (T ∉ V && shift_triangle_1(T) ∉ V && shift_triangle_2(T) ∉ V) && return false  
    end
    return true
end