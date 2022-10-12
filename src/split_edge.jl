function split_edge!(i, j, r, T, adj, adj2v, DG)
    k = get_edge(adj, i, j)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=!is_boundary_edge(j, i, adj))
    add_triangle!(i, r, k, T, adj, adj2v, DG)
    add_triangle!(r, j, k, T, adj, adj2v, DG)
    return nothing
end
function split_edge!(i, j, r, T, adj, adj2v, DG, HG::HistoryGraph{V}) where {V}
    split_edge!(i, j, r, T, adj, adj2v, DG)
    k = get_edge(adj, i, j)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    Tᵢᵣₖ = construct_triangle(V, i, r, k)
    Tᵣⱼₖ = construct_triangle(V, r, j, k)
    add_triangle!(HG, Tᵢᵣₖ, Tᵣⱼₖ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢᵣₖ, Tᵣⱼₖ)
    return nothing
end