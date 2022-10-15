"""
    locate_triangle(HG::HistoryGraph, pts, r::I, init=find_root(HG; method=:rng)) where {I}

Finds the triangle `T` containing the query point `r` using the history graph `HG`.
"""
function locate_triangle(HG::HistoryGraph, pts, r::I, init=find_root(HG; method=:rng)) where {I}
    if out_deg(HG, init) == 0
        return init, isintriangle(init, pts, r)
    end
    out = out_neighbors(HG, init)
    for T in out
        isin = isintriangle(T, pts, r)
        (isin == I(1) || isin == I(0)) && return locate_triangle(HG, pts, r, T)
    end
end

"""
    select_initial_point(pts, q; m = ceil(Int64, length(pts)^(1/3)))

Selects an initial point for the jump-and-march algorithm.
"""
function select_initial_point(pts, q; m=ceil(Int64, length(pts)^(1 / 3)))
    current_dist = typemax(eltype(q))
    current_idx = firstindex(pts) - 1
    for _ in 1:m # Not using replacement, but probability of duplicates is approximately 0.5n^(-1/3)
        i = rand(eachindex(pts))
        pᵢ = get_point(pts, i)
        sq_dist = (getx(pᵢ) - getx(q))^2 + (gety(pᵢ) - gety(q))^2
        if sq_dist < current_dist
            current_dist = sq_dist
            current_idx = i
        end
    end
    return current_idx
end
"""
    jump_and_march(q, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E}, pts;
        k = select_initial_point(pts, q),
        TriangleType::Type{V}=NTuple{3,Int64}) where {I,E,Es,V}

Uses the jump and march algorithm to locate the triangle `T` in the triangulation that 
contains the query point, starting at the vertex `k`.
"""
function jump_and_march(q, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E}, pts;
    k=select_initial_point(pts, q),
    TriangleType::Type{V}=NTuple{3,Int64}) where {I,E,Es,V}
    p = get_point(pts, k)
    i, j = rand(get_edge(adj2v, k))
    pᵢ = get_point(pts, i)
    pⱼ = get_point(pts, j)
    # Find the initial triangle to start in
    if orient(p, q, pⱼ) == I(1)
        while orient(p, q, pᵢ) == I(1)
            j = i
            pⱼ = pᵢ
            i = get_edge(adj, i, k)
            i == I(BoundaryIndex) && return jump_and_march(q, adj, adj2v, pts; TriangleType, k=j)
            pᵢ = get_point(pts, i)
        end
    else
        while orient(p, q, pⱼ) == I(-1)
            i = j
            pᵢ = pⱼ
            j = get_edge(adj, k, j)
            j == I(BoundaryIndex) && return jump_and_march(q, adj, adj2v, pts; TriangleType, k=i)
            pⱼ = get_point(pts, j)
        end
    end
    # Now do the straight line search 
    i, j = j, i
    pᵢ, pⱼ = pⱼ, pᵢ # pᵢ is left of pq, pⱼ is right of pq 
    while orient(pᵢ, pⱼ, q) == I(1)
        k = get_edge(adj, i, j)
        k == I(BoundaryIndex) && return jump_and_march(q, adj, adj2v, pts; TriangleType, k=i)
        pₖ = get_point(pts, k)
        if orient(p, q, pₖ) == I(-1)
            j = k
            pⱼ = pₖ
        else
            i = k
            pᵢ = pₖ
        end
    end
    # Swap the orientation to get a positively oriented triangle 
    k = get_edge(adj, j, i)
    return construct_triangle(V, j, i, k)
end