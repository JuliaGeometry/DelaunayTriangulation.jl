function add_point_bowyer!(T::Ts, adj, adj2v, DG, pts, r;
    pt_idx=graph(DG).V, m=ceil(Int64, length(pt_idx)^(1 / 3)),
    initial_search_point=select_initial_point(pts, r; pt_idx, m)) where {Ts} # only works for inside points currently. also fails for collinear points
    tri_type = triangle_type(Ts)
    V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)
    i, j, k = indices(V)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true)
    dig_cavity!(r, i, j, T, adj, adj2v, DG, pts)
    dig_cavity!(r, j, k, T, adj, adj2v, DG, pts)
    dig_cavity!(r, k, i, T, adj, adj2v, DG, pts)
    return nothing
end
function dig_cavity!(r, i, j, T, adj::Adjacent{I,E}, adj2v, DG, pts) where {I,E}
    ℓ = get_edge(adj, j, i) # (j, i, ℓ) is the triangle on the other side of the edge (i, j) from r 
    if ℓ == I(DefaultAdjacentValue)
        return nothing # The triangle has already been deleted in this case 
    end
    δ = ℓ ≠ I(BoundaryIndex) && isincircle(pts, r, i, j, ℓ)# Boundary edges form part of the boundary on the polygonal cavity if we have reached ℓ = BoundaryIndex
    if δ == I(1)
        delete_triangle!(j, i, ℓ, T, adj, adj2v, DG; protect_boundary=true)
        dig_cavity!(r, i, ℓ, T, adj, adj2v, DG, pts)
        dig_cavity!(r, ℓ, j, T, adj, adj2v, DG, pts)
    else # Must be an edge of the polygonal cavity. Note that this also covers the ℓ ≠ I(BoundaryIndex) case
        add_triangle!(r, i, j, T, adj, adj2v, DG)
    end
    return nothing
end
