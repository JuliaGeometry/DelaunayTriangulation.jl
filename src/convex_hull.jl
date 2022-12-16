function convex_hull(DG::DelaunayGraph{I}, pts) where {I}
    ch_pts_idx = get_neighbour(DG, I(BoundaryIndex))
    compute_centroid!(pts)
    cx, cy = get_point(pts, I(BoundaryIndex))
    θ = zeros(length(ch_pts_idx))
    j = 1
    for i in ch_pts_idx
        x, y = get_point(pts, i)
        θ[j] = atan(y - cy, x - cx)
        j += 1
    end
    sort_idx = sortperm(θ)
    ch_pts_idx = collect(ch_pts_idx)
    permute!(ch_pts_idx, sort_idx)
    return ch_pts_idx
end

function convex_hull(adj::Adjacent{I,E}, adj2v::Adjacent2Vertex) where {I,E}
    num_bnd_pts = length(get_edge(adj2v, I(BoundaryIndex)))
    bnd_pts = zeros(I, num_bnd_pts)
    i, _ = iterate(get_edge(adj2v, I(BoundaryIndex)))[1] # need to start somewhere 
    bnd_pts[1] = i
    for j in 2:num_bnd_pts
        i = get_edge(adj, i, I(BoundaryIndex))
        bnd_pts[j] = i
    end
    return bnd_pts
end

#TODO: Add method that just uses ghost triangles. There'd be no need to compute 
# anything.