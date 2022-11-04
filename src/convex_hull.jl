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