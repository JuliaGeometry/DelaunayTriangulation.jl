"""
    structured_triangulation(a, b, c, d, Nˣ, Nʸ; return_boundary_types = false, single_boundary = false, randomly_flip = false)

Computes a structured Delaunay triangulation over the rectangle `[a, b] × [c, d]`, 
using `Nˣ` points to define `[a, b]` and `Nʸ` points to define `[c, d]`. 

If `return_boundary_types = true`, also returns a vector of vectors such `BN` such 
that `BN[i]` gives the node indices that are on the `i`th boundary wall, where `y = 0` is 
`1`, `x = b` is `2`, `y = d` is `3`, and `x = a` is `4`. If `single_boundary = true`, then 
the boundary is instead represented as a contiguous edge, giving only a single set of indices 
rather than 4.

Diagonal edges can be randomly flipped by setting `randomly_flip = true`. This can be useful for avoiding bias in certain simulations 
involving this type of triangulation. Note that only diagonal edges are flipped so that all triangles keep the same area. Edges 
are flipped with probability 1/2.
"""
function triangulate_structured(a, b, c, d, Nˣ, Nʸ; return_boundary_types=false, single_boundary=false, randomly_flip=false)
    # Connectivity matrix
    T = Set{NTuple{3,Int64}}()
    Sub2Ind = LinearIndices((1:Nˣ, 1:Nʸ))
    idx = 1
    for j in 1:Nʸ-1
        for i in 1:Nˣ-1
            u = Sub2Ind[CartesianIndex(i, j)]
            v = Sub2Ind[CartesianIndex(i + 1, j)]
            w = Sub2Ind[CartesianIndex(i, j + 1)]
            add_triangle!(T, (u, v, w))

            idx += 1

            u = Sub2Ind[CartesianIndex(i, j + 1)]
            v = Sub2Ind[CartesianIndex(i + 1, j)]
            w = Sub2Ind[CartesianIndex(i + 1, j + 1)]
            add_triangle!(T, (u, v, w))

            idx += 1
        end
    end
    # Define points 
    points = Matrix{Float64}(undef, 2, Nˣ * Nʸ)
    Δx = (b - a) / (Nˣ - 1)
    Δy = (d - c) / (Nʸ - 1)
    for j in 1:Nʸ
        y = c + (j - 1) * Δy
        for i in 1:Nˣ
            x = a + (i - 1) * Δx
            points[:, Sub2Ind[CartesianIndex(i, j)]] .= [x, y]
        end
    end
    ## Get the boundary nodes so that we can triangulate 
    bnd_nodes = Int64[]
    for i in 1:Nˣ
        push!(bnd_nodes, Sub2Ind[CartesianIndex(i, 1)])
    end
    for j in 1:Nʸ
        push!(bnd_nodes, Sub2Ind[CartesianIndex(Nˣ, j)])
    end
    for i in Nˣ:-1:1
        push!(bnd_nodes, Sub2Ind[CartesianIndex(i, Nʸ)])
    end
    for j in Nʸ:-1:1
        push!(bnd_nodes, Sub2Ind[CartesianIndex(1, j)])
    end
    unique!(bnd_nodes)
    ## Triangulate 
    T, adj, adj2v, DG = triangulate(T, points, bnd_nodes)
    if randomly_flip 
        for (i, j) in edges(DG)
            if !is_boundary_edge(i,j,adj) && !is_boundary_edge(j,i,adj) && !is_ghost_edge(i,j)
                p, q = get_point(points, i, j)
                if (getx(p) ≠ getx(q)) && (gety(p) ≠ gety(q))
                    r = rand()
                    if r > 0.5 
                        flip_edge!(i,j,T,adj,adj2v,DG)
                    end
                end
            end
        end
        clear_empty_keys!(adj)
    end
    # Return boundary_nodes?
    if return_boundary_types
        # Extract boundary and identify points
        ch_idx = convex_hull(DG, points)
        BN = Vector{Vector{Int64}}([[], [], [], []])
        BNx = Vector{Vector{Float64}}([[], [], [], []])
        BNy = Vector{Vector{Float64}}([[], [], [], []])
        if single_boundary
            BN = ch_idx
            return Triangulation(T, adj, adj2v, DG, points), BN
        end
        for i in ch_idx
            if points[1, i] == a
                push!(BN[4], i)
                push!(BNx[4], points[1, i])
                push!(BNy[4], points[2, i])
            end
            if points[1, i] == b
                push!(BN[2], i)
                push!(BNx[2], points[1, i])
                push!(BNy[2], points[2, i])
            end
            if points[2, i] == c
                push!(BN[1], i)
                push!(BNx[1], points[1, i])
                push!(BNy[1], points[2, i])
            end
            if points[2, i] == d
                push!(BN[3], i)
                push!(BNx[3], points[1, i])
                push!(BNy[3], points[2, i])
            end
        end
        # Now put the points in order. 
        # BN[1] should be increasing x, BN[2] increasing y, 
        # BN[3] decreasing x, BN[4] decreasing y.
        idx₁ = sortperm(BNx[1])
        idx₂ = sortperm(BNy[2])
        idx₃ = sortperm(BNx[3]; rev=true)
        idx₄ = sortperm(BNy[4]; rev=true)
        BN[1] .= BN[1][idx₁]
        BN[2] .= BN[2][idx₂]
        BN[3] .= BN[3][idx₃]
        BN[4] .= BN[4][idx₄]
        return Triangulation(T, adj, adj2v, DG, points), BN
    else
        return Triangulation(T, adj, adj2v, DG, points)
    end
end