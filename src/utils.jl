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
    clear_empty_keys!(adj2v::Adjacent2Vertex)

Deletes any keys from `adj2v` that map to empty sets.
"""
function clear_empty_keys!(adj2v::Adjacent2Vertex)
    for (w, S) in adjacent2vertex(adj2v)
        if isempty(S)
            delete_point!(adj2v, w)
        end
    end
    return nothing
end

"""
    clear_empty_points!(DG::DelaunayGraph)

Deletes points from the graph `DG` that have degree 0. 
"""
function clear_empty_points!(DG::DelaunayGraph)
    for (w, S) in graph(DG).N
        if isempty(S)
            delete_point!(DG, w)
        end
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

"""
    compare_unconstrained_triangulations(T1, adj1, adj2v1, DG1, T2, adj2, adj2v2, DG2)

Compares the triangulation `(T1, adj1, adj2v1, DG1)` to the triangulation `(T2, adj2, adj2v2, DG2)`.
Could cause issues with triangulations of collinear points depending on the choices made during 
the algorithm.
"""
function compare_unconstrained_triangulations(T1, adj1, adj2v1, DG1, T2, adj2, adj2v2, DG2)
    adj = deepcopy(adj1)
    _adj = deepcopy(adj2)
    clear_empty_keys!(adj, DG1)
    clear_empty_keys!(_adj, DG2)
    adj2v = deepcopy(adj2v1)
    _adj2v = deepcopy(adj2v2)
    clear_empty_keys!(adj2v)
    clear_empty_keys!(_adj2v)
    return compare_triangle_sets(T1, T2) &&
           adjacent(adj) == adjacent(_adj) &&
           adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
           graph(DG1) == graph(DG2)
end

"""
    check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)

Checks if `adj` and `adj2v` are related so that, if `get_edge(adj, i, j) = k`,
then `(i, j) ∈ get_edge(adj2v, k)`. Returns `true` if so, and `false` otherwise.
"""
function check_adjacent_is_adjacent2vertex_inverse(adj::Adjacent{I,E}, adj2v) where {I,E}
    # Check adj2v 
    for (k, S) in adjacent2vertex(adj2v)
        for ij in S
            get_edge(adj, ij) ≠ k && return false
        end
    end
    # Check adj
    for (ij, k) in adjacent(adj)
        if k ≠ I(DefaultAdjacentValue)
            ij ∉ get_edge(adj2v, k) && return false
        end
    end
    return true
end

"""
    clear_centroid_coordinates!(::Type{T} = Float64) where {T}

Resets the centroid coordinates in `CentroidCoordinates` to be 
`(zero(T), zero(T))`.
"""
function clear_centroid_coordinates!(::Type{T}=Float64) where {T}
    CentroidCoordinates.x = zero(T)
    CentroidCoordinates.y = zero(T)
    CentroidCoordinates.n = 0
    return nothing
end

"""
    update_centroid_after_new_point!(pts, i)

Updates the centroid coordinates in `CentroidCoordinates` after a new point 
`get_point(pts, i)` was added.
"""
function update_centroid_after_new_point!(pts, i)
    new_pt = get_point(pts, i)
    n = CentroidCoordinates.n
    CentroidCoordinates.x = 1 / (n + 1) * (n * CentroidCoordinates.x + getx(new_pt))
    CentroidCoordinates.y = 1 / (n + 1) * (n * CentroidCoordinates.y + gety(new_pt))
    CentroidCoordinates.n += 1
    return nothing
end

"""
    update_centroid_after_deleted_point!(pts, i)

Updates the centroid coordinates in `CentroidCoordinates` after the point
`get_point(pts, i)` was deleted.
"""
function update_centroid_after_deleted_point!(pts, i)
    new_pt = get_point(pts, i)
    n = CentroidCoordinates.n
    CentroidCoordinates.x = 1 / (n - 1) * (n * CentroidCoordinates.x - getx(new_pt))
    CentroidCoordinates.y = 1 / (n - 1) * (n * CentroidCoordinates.y - gety(new_pt))
    CentroidCoordinates.n -= 1
    return nothing
end