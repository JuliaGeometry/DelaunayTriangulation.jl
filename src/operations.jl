"""
    add_triangle!(T, adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours = true)

Adds the triangle `(i, j, r)` into the triangulation, assuming the `r`th point is to the left 
of the directed edge `(i, j)` and `k = get_edge(adj, i, j)`. If `(i, j, k)` should all remain 
neighbours, set `delete_adjacent_neighbours = false`.
"""
function add_triangle!(T::Ts, adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, pts, i, j, k) where {Ts}
    V = triangle_type(Ts)
    add_triangle!(T, V, i, j, k)
    add_triangle!(adj, i, j, k)
    add_triangle!(adj2v, i, j, k)
    add_triangle!(DG, i, j, k)
    ji_is_valid, kj_is_valid, ik_is_valid,
    ji_is_boundary, kj_is_boundary, ik_is_boundary =
        check_triangle_edge_validity_and_boundary(i, j, k, adj)
    if ji_is_boundary && !kj_is_boundary && !ik_is_boundary
        add_edge!
    end
    return nothing
end
function add_triangle!(T, V, i, j, k)
    newT = construct_triangle(V, i, j, k)
    add_triangle!(T, newT)
    return nothing
end
function add_triangle!(adj::Adjacent, i, j, k)
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    return nothing
end
function add_triangle!(adj2v::Adjacent2Vertex, i, j, k)
    add_edge!(adj2v, i, j, k)
    add_edge!(adj2v, j, k, i)
    add_edge!(adj2v, k, i, j)
    return nothing
end
function add_triangle!(DG::DelaunayGraph, i, j, k)
    add_neighbour!(DG, k, i, j)
    add_neighbour!(DG, i, j)
    return nothing
end
function add_boundary_edge!(i, j, adj, adj2v)

end

"""
    delete_triangle!(T, adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, k)

Deletes the triangle `(i, j, k)` from the triangulation, assuming `(i, j, k)` is a positively oriented triangle. 
"""
function delete_triangle!(T::Ts, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, i, j, k) where {Ts,I,E}
    V = triangle_type(Ts)
    delete_triangle!(T, V, i, j, k)
    delete_triangle!(adj, i, j, k)
    delete_triangle!(adj2v, i, j, k)
    ij_is_valid, jk_is_valid, ki_is_valid,
    ij_is_boundary, jk_is_boundary, ki_is_boundary = 
    check_triangle_edge_validity_and_boundary(i, j, k, adj)
    ij_is_valid = is_valid_edge(j, i, adj)
    jk_is_valid = is_valid_edge(k, j, adj)
    ki_is_valid = is_valid_edge(i, k, adj)
    ij_is_boundary = ij_is_valid && is_boundary_edge(j, i, adj)
    jk_is_boundary = jk_is_valid && is_boundary_edge(k, j, adj)
    ki_is_boundary = ki_is_valid && is_boundary_edge(i, k, adj)
    if ij_is_boundary && jk_is_boundary && ki_is_boundary # This only occurs if the triangulation is made up of a single triangle 
        delete_neighbour!(DG, i, j, k)
        delete_neighbour!(DG, j, k)
        return nothing
    end
    ij_is_boundary && delete_boundary_edge!(i, j, k, adj, adj2v, DG; protect_jk=jk_is_boundary, protect_ki=ki_is_boundary)
    jk_is_boundary && delete_boundary_edge!(j, k, i, adj, adj2v, DG; protect_jk=ki_is_boundary, protect_ki=ij_is_boundary)
    ki_is_boundary && delete_boundary_edge!(k, i, j, adj, adj2v, DG; protect_jk=ij_is_boundary, protect_ki=jk_is_boundary)
    !ij_is_valid && delete_neighbour!(DG, i, j)
    !jk_is_valid && delete_neighbour!(DG, j, k)
    !ki_is_valid && delete_neighbour!(DG, k, i)
    return nothing
end
function delete_triangle!(T, V, i, j, k)
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
    flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `T`: The current set of triangles defining the triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The current edge.
- `k, r`: Indices for the points the edge is flipped onto.

It is assumed that `(i, k, j)` and `(i, j, r)` are positively oriented triangles.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function flip_edge!(T::Ts, HG::HistoryGraph, adj::Adjacent,
    adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, k, r) where {Ts}
    V = triangle_type(Ts)
    Tᵢₖⱼ = construct_triangle(V, i, k, j)
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    Tᵣₖⱼ = construct_triangle(V, r, k, j)
    Tᵣᵢₖ = construct_triangle(V, r, i, k)
    flip_edge!(T, Tᵢₖⱼ, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    flip_edge!(adj, i, j, k, r)
    flip_edge!(adj2v, i, j, k, r)
    flip_edge!(DG, i, j, k, r)
    flip_edge!(HG, Tᵢₖⱼ, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    return nothing
end
function flip_edge!(T, Tᵢₖⱼ::V, Tᵢⱼᵣ::V, Tᵣₖⱼ::V, Tᵣᵢₖ::V) where {V}
    delete_triangle!(T, Tᵢₖⱼ, Tᵢⱼᵣ)
    add_triangle!(T, Tᵣₖⱼ, Tᵣᵢₖ)
    return nothing
end
function flip_edge!(adj::Adjacent{I,E}, i::I, j::I, k::I, r::I) where {I,E}
    delete_edge!(adj, i, j)
    add_triangle!(adj, r, k, j)
    add_triangle!(adj, r, i, k)
    return nothing
end
function flip_edge!(adj2v::Adjacent2Vertex{I,Es,E}, i::I, j::I, k::I, r::I) where {I,Es,E}
    delete_edge!(adj2v, i, k, j)
    delete_edge!(adj2v, i, j, r)
    delete_edge!(adj2v, j, r, i)
    delete_edge!(adj2v, j, i, k)
    delete_edge!(adj2v, k, j, i)
    delete_edge!(adj2v, r, i, j)
    add_edge!(adj2v, i, k, r)
    add_edge!(adj2v, j, r, k)
    add_edge!(adj2v, k, j, r)
    add_edge!(adj2v, k, r, i)
    add_edge!(adj2v, r, k, j)
    add_edge!(adj2v, r, i, k)
    return nothing
end
function flip_edge!(DG::DelaunayGraph{I}, i::I, j::I, k::I, r::I) where {I}
    delete_neighbour!(DG, i, j)
    add_neighbour!(DG, r, k)
    return nothing
end
function flip_edge!(HG::HistoryGraph{V}, Tᵢₖⱼ::V, Tᵢⱼᵣ::V, Tᵣₖⱼ::V, Tᵣᵢₖ::V) where {V}
    add_triangle!(HG, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(HG, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(HG, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
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
function split_triangle!(T, adj, adj2v, DG, Tᵢⱼₖ::V, r) where {V}
    i, j, k = indices(Tᵢⱼₖ)
    delete_triangle!(T, adj, adj2v, DG, i, j, k)
    add_triangle!(T, adj, adj2v, DG, i, j, r)
    add_triangle!(T, adj, adj2v, DG, j, k, r)
    add_triangle!(T, adj, adj2v, DG, k, i, r)
    return nothing
end
function split_triangle!(T, HG::HistoryGraph, adj, adj2v, DG, Tᵢⱼₖ::V, r) where {V}
    i, j, k = indices(Tᵢⱼₖ)
    split_triangle!(T, adj, adj2v, DG, Tᵢⱼₖ, r)
    split_triangle!(HG, V, i, j, k, r)
    return nothing
end
function split_triangle!(HG::HistoryGraph, V, i, j, k, r)
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    Tⱼₖᵣ = construct_triangle(V, j, k, r)
    Tₖᵢᵣ = construct_triangle(V, k, i, r)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(HG, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
end

"""
    legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `T`: The current set of triangles defining the triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The edge to make legal. Nothing happens if `is_legal(i, j, 𝒜, pts)`.
- `r`: The point being added into the triangulation. 
- `pts`: The point set of the triangulation.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(T, HG, adj, adj2v, DG, i, j, e, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, e, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, e, j, r, pts)
    end
    return nothing
end

function locate_triangle(T, pts, r::I) where {I}
    for V in T
        if isintriangle(V, pts, r) ≠ I(-1)
            return V
        end
    end
end
function add_point_bowyer!(T, adj, adj2v, DG, pts, r)
    V = locate_triangle(T, pts, r)
    i, j, k = indices(V)
    delete_triangle!(T, adj, adj2v, DG, i, j, k)
    dig_cavity!(T, adj, adj2v, DG, pts, r, i, j, k)
    dig_cavity!(T, adj, adj2v, DG, pts, r, j, k, i)
    dig_cavity!(T, adj, adj2v, DG, pts, r, k, i, j)
    return nothing
end
function dig_cavity!(T, adj::Adjacent{I,E}, adj2v, DG, pts, r, i, j, k) where {I,E}
    ℓ = get_edge(adj, j, i) # (j, i, ℓ) is the triangle on the other side of the edge (i, j) from r 
    @show ℓ
    if ℓ == DefaultAdjacentValue || ℓ == BoundaryIndex
        return nothing # The triangle has already been deleted in this case 
    end
    δ = isincircle(pts, r, i, j, ℓ)
    if δ == I(1)
        delete_triangle!(T, adj, adj2v, DG, j, i, ℓ)
        dig_cavity!(T, adj, adj2v, DG, pts, r, i, ℓ, j)
        dig_cavity!(T, adj, adj2v, DG, pts, r, ℓ, j, i)
    else
        add_triangle!(T, adj, adj2v, DG, r, i, j, k)
    end
end