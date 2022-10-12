function split_triangle!(i, j, k, r, T, adj, adj2v, DG)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true)
    add_triangle!(i, j, r, T, adj, adj2v, DG)
    add_triangle!(j, k, r, T, adj, adj2v, DG)
    add_triangle!(k, i, r, T, adj, adj2v, DG)
    return nothing
end
function split_triangle!(i, j, k, r, T, adj, adj2v, DG, HG::HistoryGraph) 
    split_triangle!(i, j, k, r, T, adj, adj2v, DG)
    split_triangle!(i, j, k, r, HG)
    return nothing
end
function split_triangle!(i, j, k, r, HG::HistoryGraph{V}) where {V}
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    Tⱼₖᵣ = construct_triangle(V, j, k, r)
    Tₖᵢᵣ = construct_triangle(V, k, i, r)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    add_triangle!(HG, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢⱼᵣ, Tⱼₖᵣ, Tₖᵢᵣ)
end