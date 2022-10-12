"""
    is_boundary_edge(i, j, adj::Adjacent{I, E}) where {I, E}

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

"""
    edge_exists(i, j, adj::Adjacent{I, E}) where {I, E}

Returns `true` if the edge `(i, j)` is an edge in the triangulation, 
and `false` otherwise.
"""
function edge_exists(ij, adj::Adjacent{I, E}) where {I, E}
    return get_edge(adj, ij) ≠ I(DefaultAdjacentValue)
end 
function edge_exists(i, j, adj::Adjacent{I, E}) where {I, E}
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
function validate_triangulation(T, adj::Adjacent{I,E},adj2v,DG,pts) where {I, E}
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