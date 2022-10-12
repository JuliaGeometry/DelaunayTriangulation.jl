function legalise_edge!(i, j, r, T, HG, adj, adj2v, DG, pts)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(i, j, e, r, T, adj, adj2v, DG, HG)
        legalise_edge!(i, e, r, T, HG, adj, adj2v, DG, pts)
        legalise_edge!(e, j, r, T, HG, adj, adj2v, DG, pts)
    end
    return nothing
end

