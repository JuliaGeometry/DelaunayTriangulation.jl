function flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    delete_triangle!(i, k, j, T, adj, adj2v, DG; protect_boundary=true)
    delete_triangle!(i, j, ℓ, T, adj, adj2v, DG; protect_boundary=true)
    add_triangle!(ℓ, k, j, T, adj, adj2v, DG)
    add_triangle!(ℓ, i, k, T, adj, adj2v, DG)
    return nothing
end
function flip_edge!(i, j, T, adj, adj2v, DG)
    ℓ = get_edge(adj, i, j)
    k = get_edge(adj, j, i)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    return nothing
end
function flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG, HG::HistoryGraph)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG)
    flip_edge!(i, j, k, ℓ, HG)
    return nothing
end
function flip_edge!(i, j, T, adj, adj2v, DG, HG::HistoryGraph)
    ℓ = get_edge(adj, i, j)
    k = get_edge(adj, j, i)
    flip_edge!(i, j, k, ℓ, T, adj, adj2v, DG, HG)
    return nothing
end
function flip_edge!(i, j, k, ℓ, HG::HistoryGraph{V}) where {V}
    Tᵢₖⱼ = construct_triangle(V, i, k, j)
    Tᵢⱼₗ = construct_triangle(V, i, j, ℓ)
    Tₗₖⱼ = construct_triangle(V, ℓ, k, j)
    Tₗᵢₖ = construct_triangle(V, ℓ, i, k)
    add_triangle!(HG, Tₗₖⱼ, Tₗᵢₖ)
    add_edge!(HG, Tᵢₖⱼ, Tₗₖⱼ, Tₗᵢₖ)
    add_edge!(HG, Tᵢⱼₗ, Tₗₖⱼ, Tₗᵢₖ)
    return nothing
end
