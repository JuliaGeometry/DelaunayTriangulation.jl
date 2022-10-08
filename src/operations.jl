"""
    add_triangle!(T, adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours = true)

Adds the triangle `(i, j, r)` into the triangulation, assuming the `r`th point is to the left 
of the directed edge `(i, j)`. and `k = get_edge(adj, i, j)`. If `(i, j, k)` should all remain 
neighbours, set `delete_adjacent_neighbours = false`.
"""
function add_triangle!(T::Ts, adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, HG::HistoryGraph,
    i, j, r, k; delete_adjacent_neighbours=true) where {Ts}
    add_triangle!(T, i, j, r, k)
    add_triangle!(adj, i, j, r)
    add_triangle!(adj2v, i, j, r, k)
    add_triangle!(DG, i, j, r, k; delete_adjacent_neighbours)
    V = triangle_type(Ts)
    add_triangle!(HG, V, i, j, r, k)
    return nothing
end
function add_triangle!(T::Ts, i, j, r, k) where {Ts}
    V = triangle_type(Ts)
    newT = construct_triangle(V, i, j, r)
    oldT = construct_triangle(V, i, j, k)
    add_triangle!(T, newT)
    delete_triangle!(T, oldT)
    return nothing
end
function add_triangle!(adj::Adjacent, i, j, r)
    add_edge!(adj, i, j, r)
    add_edge!(adj, j, r, i)
    add_edge!(adj, r, i, j)
    return nothing
end
function add_triangle!(adj2v::Adjacent2Vertex, i, j, r, k)
    add_edge!(adj2v, i, j, r)
    add_edge!(adj2v, j, r, i)
    add_edge!(adj2v, r, i, j)
    delete_triangle!(adj2v, i, j, k) # By adding the triangle inside (i, j, k), we need to delete (i, j, k)
    return nothing
end
function add_triangle!(DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours=true)
    add_neighbour!(DG, r, i, j)
    delete_adjacent_neighbours && delete_neighbour!(DG, k, i, j)
    return nothing
end
function add_triangle!(HG::HistoryGraph, V, i, j, r, k)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    add_triangle!(HG, Tᵢⱼᵣ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢⱼᵣ)
    return nothing
end

"""
    delete_triangle!(T, adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, k)

Deletes the triangle `(i, j, k)` from the triangulation, assuming `(i, j, k)` is a positively oriented triangle. 
"""
function delete_triangle!(T, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, i, j, k) where {I,E}
    delete_triangle!(T, i, j, k)
    delete_triangle!(adj, i, j, k)
    delete_triangle!(adj2v, i, j, k)
    ij_is_valid = is_valid_edge(j, i, adj)
    jk_is_valid = is_valid_edge(k, j, adj)
    ki_is_valid = is_valid_edge(i, k, adj)
    ij_is_boundary = ij_is_valid && is_boundary_edge(j, i, adj) # check valid to avoid keyerrors
    jk_is_boundary = jk_is_valid && is_boundary_edge(k, j, adj)
    ki_is_boundary = ki_is_valid && is_boundary_edge(i, k, adj)
    ij_is_boundary && delete_boundary_edge!(i, j, k, adj, adj2v, DG; protect_jk=jk_is_boundary, protect_ki=ki_is_boundary)
    jk_is_boundary && delete_boundary_edge!(j, k, i, adj, adj2v, DG; protect_jk=ki_is_boundary, protect_ki=ij_is_boundary)
    ki_is_boundary && delete_boundary_edge!(k, i, j, adj, adj2v, DG; protect_jk=ij_is_boundary, protect_ki=jk_is_boundary)
    !ij_is_valid && delete_neighbour!(DG, i, j)
    !jk_is_valid && delete_neighbour!(DG, j, k)
    !ki_is_valid && delete_neighbour!(DG, k, i)
    return nothing
end
function delete_triangle!(T::Ts, i, j, k) where {Ts}
    V = triangle_type(Ts)
    delT = construct_triangle(V, i, j, k)
    delete_triangle!(T, delT)
    return nothing
end
function delete_triangle!(adj::Adjacent, i, j, k)
    delete_edge!(adj, i, j; protect_boundary=false, delete_uv_only=true)
    delete_edge!(adj, j, k; protect_boundary=false, delete_uv_only=true)
    delete_edge!(adj, k, i; protect_boundary=false, delete_uv_only=true)
    return nothing
end
function delete_triangle!(adj2v::Adjacent2Vertex, i, j, k)
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    return nothing
end
"""
    delete_boundary_edge!(i, j, k, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex, DG::DelaunayGraph;
        protect_jk=false, protect_ki=false) where {I,E}

Deletes the edge `(i, j)`, assuming it is a boundary edge. The vertex `k` is 
such that `(i, j, k)` was originally a positively oriented triangle. If the 
edge `(j, k)` should be retained in `get_edge(adj2v, I(BoundaryIndex))`, set 
`protect_jk = true`, and similarly for `(k, i)`.
"""
function delete_boundary_edge!(i, j, k, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex, DG::DelaunayGraph;
    protect_jk=false, protect_ki=false) where {I,E}
    delete_neighbour!(DG, i, j)
    delete_edge!(adj, j, i; protect_boundary=false, delete_uv_only=true)
    delete_edge!(adj2v, I(BoundaryIndex), j, i)
    if !protect_jk
        add_edge!(adj2v, I(BoundaryIndex), j, k)
        add_edge!(adj, j, k, I(BoundaryIndex))
    end
    if !protect_ki
        add_edge!(adj2v, I(BoundaryIndex), k, i)
        add_edge!(adj, k, i, I(BoundaryIndex))
    end
    return nothing
end







"""
    split_triangle!(T, HG::HistoryGraph, adj, adj2v, DG, Tᵢⱼₖ::V, r) where V

Given a triangulation `T`, adds the `r`th point of the point set into the triangulation, assuming 
that `r` is in the interior of the triangle `Tᵢⱼₖ`.

# Arguments 
- `T`: The current set of triangles defining the triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function split_triangle!(T, HG::HistoryGraph, adj, adj2v, DG, Tᵢⱼₖ::V, r) where {V}
    i, j, k = indices(Tᵢⱼₖ)
    add_triangle!(T, adj, adj2v, DG, HG, i, j, r, k; delete_adjacent_neighbours=false)
    add_triangle!(T, adj, adj2v, DG, HG, j, k, r, i; delete_adjacent_neighbours=false)
    add_triangle!(T, adj, adj2v, DG, HG, k, i, r, j; delete_adjacent_neighbours=false)
    return nothing
end