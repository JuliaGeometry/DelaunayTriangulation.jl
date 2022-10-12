function delete_triangle!(i, j, k, T::Ts,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E},
    DG::DelaunayGraph{I}; protect_boundary = false) where {I,E,Es,Ts}
    V = triangle_type(Ts)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    delete_triangle!(T, Tᵢⱼₖ)
    delete_edge!(adj, i, j)
    delete_edge!(adj, j, k)
    delete_edge!(adj, k, i)
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    ji_bnd = is_boundary_edge(j, i, adj)
    ik_bnd = is_boundary_edge(i, k, adj)
    kj_bnd = is_boundary_edge(k, j, adj)
    num_bnd_edges = !protect_boundary ? count((ji_bnd, ik_bnd, kj_bnd)) : 0 # If we are protecting the boundary, setting this count to zero will do that 
    ji_exists = edge_exists(j, i, adj) # Could also write this as ji_is_bnd || edge_exists(j, i, adj),
    ik_exists = edge_exists(i, k, adj) #    but this way makes the code below simpler.
    kj_exists = edge_exists(k, j, adj)
    protect_ji_neighbour = ji_exists && !ji_bnd
    protect_ik_neighbour = ik_exists && !ik_bnd
    protect_kj_neighbour = kj_exists && !kj_bnd
    !protect_ji_neighbour && delete_neighbour!(DG, i, j)
    !protect_ik_neighbour && delete_neighbour!(DG, k, i)
    !protect_kj_neighbour && delete_neighbour!(DG, j, k)
    if num_bnd_edges == 1
        delete_boundary_edges_single!(i, j, k, ji_bnd, ik_bnd, kj_bnd, adj, adj2v)
    elseif num_bnd_edges == 2
        delete_boundary_edges_double!(i, j, k, ji_bnd, ik_bnd, kj_bnd, adj, adj2v)
    elseif num_bnd_edges == 3 # length(T) == 0
        delete_boundary_edges_triple!(i, j, k, adj, adj2v)
    end
    return nothing
end
function delete_boundary_edges_single!(i, j, k, ji_bnd, ik_bnd, kj_bnd,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex) where {I,E}
    u, v, w = choose_uvw(ji_bnd, kj_bnd, ik_bnd, i, j, k)
    delete_edge!(adj, v, u)
    delete_edge!(adj2v, I(BoundaryIndex), v, u)
    add_edge!(adj, v, w, I(BoundaryIndex))
    add_edge!(adj, w, u, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), v, w)
    add_edge!(adj2v, I(BoundaryIndex), w, u)
    return nothing
end
function delete_boundary_edges_double!(i, j, k, ji_bnd, ik_bnd, kj_bnd,
    adj::Adjacent{I,E}, adj2v::Adjacent2Vertex) where {I,E}
    u, v, w = choose_uvw(!ji_bnd, !kj_bnd, !ik_bnd, i, j, k)
    delete_edge!(adj, u, w)
    delete_edge!(adj, w, v)
    delete_edge!(adj2v, I(BoundaryIndex), u, w)
    delete_edge!(adj2v, I(BoundaryIndex), w, v)
    add_edge!(adj, u, v, I(BoundaryIndex))
    add_edge!(adj2v, I(BoundaryIndex), u, v)
    return nothing
end
function delete_boundary_edges_triple!(i, j, k, adj::Adjacent{I,E}, adj2v) where {I,E}
    delete_edge!(adj, k, j)
    delete_edge!(adj, j, i)
    delete_edge!(adj, i, k)
    delete_edge!(adj2v, I(BoundaryIndex), k, j)
    delete_edge!(adj2v, I(BoundaryIndex), j, i)
    delete_edge!(adj2v, I(BoundaryIndex), i, k)
    return nothing
end