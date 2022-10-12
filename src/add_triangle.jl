function add_triangle!(i, j, k, T::Ts,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E},
    DG::DelaunayGraph{I}) where {I,E,Es,Ts}
    V = triangle_type(Ts)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(T, Tᵢⱼₖ)
    ij_bnd = is_boundary_edge(i, j, adj)
    jk_bnd = is_boundary_edge(j, k, adj)
    ki_bnd = is_boundary_edge(k, i, adj)
    num_bnd_edges = count((ij_bnd, jk_bnd, ki_bnd))
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    add_edge!(adj2v, i, j, k)
    add_edge!(adj2v, j, k, i)
    add_edge!(adj2v, k, i, j)
    add_neighbour!(DG, k, i, j)
    add_neighbour!(DG, i, j)
    if num_bnd_edges == 1
        add_boundary_edges_single!(i, j, k, ij_bnd, jk_bnd, ki_bnd, adj, adj2v)
    elseif num_bnd_edges == 2
        add_boundary_edges_double!(i, j, k, ij_bnd, jk_bnd, ki_bnd, adj, adj2v)
    elseif length(T) == 1 # If this is true, then we just have a single triangle, and thus num_bnd_edges should be 3 
        add_boundary_edges_triple!(i, j, k, adj, adj2v)
    end
    return nothing
end
function add_triangle!(i, j, k, T, adj, adj2v, DG, HG::HistoryGraph{V}) where {V}
    add_triangle!(i, j, k, T, adj, adj2v, DG)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(HG, Tᵢⱼₖ)
    return nothing
end
function add_boundary_edges_single!(i, j, k, ij_bnd, jk_bnd, ki_bnd,
    adj::Adjacent{I,E}, adj2v) where {I,E}
    u, v, w = choose_uvw(ij_bnd, jk_bnd, ki_bnd, i, j, k)
    add_edge!(adj, u, w, I(BoundaryIndex))
    add_edge!(adj, w, v, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), u, w)
    add_edge!(adj2v, I(BoundaryIndex), w, v)
    delete_edge!(adj2v, I(BoundaryIndex), u, v)
    return nothing
end
function add_boundary_edges_double!(i, j, k, ij_bnd, jk_bnd, ki_bnd,
    adj::Adjacent{I,E}, adj2v) where {I,E}
    u, v, w = choose_uvw(!ij_bnd, !jk_bnd, !ki_bnd, i, j, k)
    add_edge!(adj, v, u, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), v, u)
    delete_edge!(adj2v, I(BoundaryIndex), v, w)
    delete_edge!(adj2v, I(BoundaryIndex), w, u)
    return nothing
end
function add_boundary_edges_triple!(i, j, k, adj::Adjacent{I,E}, adj2v) where {I,E}
    add_edge!(adj, j, i, I(BoundaryIndex))
    add_edge!(adj, i, k, I(BoundaryIndex))
    add_edge!(adj, k, j, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), j, i)
    add_edge!(adj2v, I(BoundaryIndex), i, k)
    add_edge!(adj2v, I(BoundaryIndex), k, j)
    return nothing
end