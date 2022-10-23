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

function locate_triangle(T, pts, r::I) where {I}# brute force
    for V in T
        if isintriangle(V, pts, r) ≠ I(-1)
            return V
        end
    end
end

"""
    select_initial_point(pts, q; m = ceil(Int64, length(pts)^(1/3)), pt_idx = eachindex(pts))

Selects an initial point for the jump-and-march algorithm for the query point `q`.
"""
function select_initial_point(pts, q; m=ceil(Int64, length(pts)^(1 / 3)), pt_idx=eachindex(pts))
    current_dist = typemax(eltype(q))
    current_idx = firstindex(pts) - 1 # Just some index not in eachindex(pts)
    for _ in 1:m # Not using replacement, but probability of duplicates is approximately 0.5n^(-1/3)
        i = rand(pt_idx)
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
    select_initial_point(pts, q::Integer; m=ceil(Int64, length(pts)^(1 / 3)))

Selects an initial point for the jump-and-march algorithm for the query point `q`. It is assumed 
that the query point `q` is given by `get_point(pts, q)`, and that `q` is a point being added into 
a triangulation. That is, the initial point will not be `q`.
"""
function select_initial_point(pts, q::Integer; pt_idx=eachindex(pts), m=ceil(Int64, length(pt_idx)^(1 / 3)))
    return select_initial_point(pts, get_point(pts, q); m, pt_idx)
end

"""
    select_initial_triangle(q, adj::Adjacent{I, E}, adj2v, DG, k, pts) where {I, E}

Given a query point `q` and an initial point `k` (from e.g. [`select_initial_point`](@ref)),
finds the triangle `(i, j, k)` to start the jump and march algorithm (see [`jump_and_march`](@ref)) in.
Returns `p, i, j, pᵢ, pⱼ` such that `pᵢ = get_point(pts, i)`, `pⱼ = get_point(pts, j)`, 
`p = get_point(pts, k)`, and `pᵢ` is to the left of `pq` and `pⱼ` is to the right of `pq`.
"""
function select_initial_triangle(q, adj::Adjacent{I,E}, adj2v, DG, k, pts) where {I,E}
    #if !is_boundary_point(k, adj, DG)
    return select_initial_triangle_interior_start(q, adj, adj2v, DG, k, pts)
    #end
end
"""
    select_initial_triangle_interior_start(q, adj::Adjacent{I,E}, adj2v, DG, k, pts) where {I,E}

Selects an initial triangle, assuming that the point that we are starting at is in the interior of the triangulation.
"""
function select_initial_triangle_interior_start(q, adj::Adjacent{I,E}, adj2v, DG, k, pts) where {I,E}
    p = get_point(pts, k)
    i, j = rand(get_edge(adj2v, k))
    pᵢ = get_point(pts, i)
    pⱼ = get_point(pts, j)
    # Find the initial triangle to start in
    if orient(p, q, pⱼ) == 1
        while orient(p, q, pᵢ) == 1
            j = i
            pⱼ = pᵢ
            i = get_edge(adj, i, k)
            i == I(BoundaryIndex) && return select_initial_triangle(q, adj, adj2v, DG, j, pts)
            pᵢ = get_point(pts, i)
        end
    else
        while orient(p, q, pⱼ) == -1
            i = j
            pᵢ = pⱼ
            j = get_edge(adj, k, j)
            j == I(BoundaryIndex) && return select_initial_triangle(q, adj, adj2v, DG, i, pts)
            pⱼ = get_point(pts, j)
        end
    end
    i, j = j, i
    pᵢ, pⱼ = pⱼ, pᵢ # pᵢ is left of pq, pⱼ is right of pq 
    return p, i, j, pᵢ, pⱼ
end
"""
    check_interior_edge_intersections(q, adj::Adjacent{I, E}, DG, k, pts) where {I, E}

Checks if any of the interior edges starting from the boundary point `get_point(pts, k)` 
intersect the line `pq`, where `p = get_point(pts, k)`. Returns:
- `(i, j, true, false)`: If `pq` intersects the edge `(i, j)`.
- `(i, j, false, true)`: If `pq` does not intersect `(i, j)`, but instead `(i, j, k)` contains `q`.
- `(0, 0, false, false)`: If `pq` does not intersect `(i, j)`.
"""
function check_interior_edge_intersections(q, adj::Adjacent{I,E}, DG, k, pts) where {I,E}
    p = get_point(pts, k)
    i = get_edge(adj, k, I(BoundaryIndex)) # Note that this is at the left
    pᵢ = get_point(pts, i)
    o1 = orient(p, q, pᵢ) # Is pᵢ to the left of pq?
    for _ in 1:(deg(graph(DG), k)-2) # There are deg(graph(DG), k) neighbours, but one of those is ∂, and we've already picked one, hence why we subtract 2.
        j = get_edge(adj, k, i)
        pⱼ = get_point(pts, j)
        o2 = orient(p, q, pⱼ) # Is pⱼ to the left of pq?
        if ExactPredicates.opposite_signs(o1, o2) # If they're opposite signs, we have an intersection 
            if meet(p, q, pᵢ, pⱼ) == 1 # Does pq intersect pᵢpⱼ? It's possible that we have opposite signs but q is on the other side, so we do need to check this. Note that we don't need this in the interior edge case because of the ability to completely rotate around.
                return j, i, true, false # Switch i, j so that pi is left of pq and pj is right of pq
            elseif orient(pⱼ, p, q) == 1 && orient(p, pᵢ, q) == 1 # It may not intersect, but it could be inside the triangle. 
                return i, j, false, true
            end
        end
        o1, i, pᵢ = o2, j, pⱼ # Step onto the next triangle
    end
    return I(0), I(0), false, false
end
"""
    straight_line_search_ghost_triangles(q, adj::Adjacent{I, E}, k, pts) where {I, E}

Finds the ghost triangle containing `q`, starting at the `k`th point `get_point(pts, k)`, assuming 
that `get_point(pts, k)` is on the boundary of the triangulation and `q` is outside of the triangulation. 
The returned value is the boundary edge `(i, j)` so that the ghost triangle `(i, j, $BoundaryIndex)` is the 
one that contains `q`. If `q` is inside the triangulation, the returned value is meaningless.
"""
function straight_line_search_ghost_triangles(q, adj::Adjacent{I,E}, k, pts) where {I,E}
    pc = get_point(pts, I(BoundaryIndex))
    pᵢ = get_point(pts, k)
    i = k
    if orient(pc, pᵢ, q) == 1 # This means that q is left of the ghost edge through pᵢ = pₖ, so we'll rotate left around the boundary
        j = get_edge(adj, i, I(BoundaryIndex))
        pⱼ = get_point(pts, j)
        while orient(pc, pⱼ, q) == 1
            i = j
            pᵢ = pⱼ
            j = get_edge(adj, i, I(BoundaryIndex))
            pⱼ = get_point(pts, j)
        end
        return j, i # In this case, we need to swap the orientation so that i, j is a boundary edge 
    else # if orient(pc, pᵢ, q) == -1, which means that q is right of the ghost edge through pᵢ = pₖ, so we'll rotate right around the boundary 
        j = get_edge(adj, I(BoundaryIndex), i)
        pⱼ = get_point(pts, j)
        while orient(pc, pⱼ, q) == -1
            i = j
            pᵢ = pⱼ
            j = get_edge(adj, I(BoundaryIndex), i)
            pⱼ = get_point(pts, j)
        end
        return i, j
    end
end

"""
    jump_and_march(q, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E}, DG, pts;
        pt_idx = eachindex(pts),
        k = select_initial_point(pts, q), 
        TriangleType::Type{V}=NTuple{3,Int64}) where {I,E,Es,V}

Uses the jump and march algorithm to locate the triangle `T` in the triangulation that 
contains the query point, starting at the vertex `k`.
"""
function jump_and_march(q, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E}, DG, pts;
    pt_idx=eachindex(pts), m=ceil(Int64, length(pt_idx)^(1 / 3)),
    k=select_initial_point(pts, q; m, pt_idx),
    TriangleType::Type{V}=NTuple{3,Int64}) where {I,E,Es,V}
    if !is_boundary_point(k, adj, DG) || !triangulation_has_ghost_triangles(adj, adj2v) # If the triangulation does not have ghost triangles, we cannot use the methods below.
        p, i, j, pᵢ, pⱼ = select_initial_triangle_interior_start(q, adj, adj2v, DG, k, pts)
    else
        i, j, intersects_edge, inside_triangle = check_interior_edge_intersections(q, adj, DG, k, pts)
        if inside_triangle
            return construct_triangle(V, i, j, k)
        elseif !intersects_edge
            i, j = straight_line_search_ghost_triangles(q, adj, k, pts)
            return construct_triangle(V, i, j, I(BoundaryIndex))
        end
        p, pᵢ, pⱼ = get_point(pts, k, i, j)
    end
    # Now do the straight line search 
    while orient(pᵢ, pⱼ, q) == I(1)
        k = get_edge(adj, i, j)
        if k == I(BoundaryIndex)
            if triangulation_has_ghost_triangles(adj, adj2v)
                _i, _j = straight_line_search_ghost_triangles(q, adj, i, pts)
                return construct_triangle(V, _i, _j, I(BoundaryIndex))
            else
                return jump_and_march(q, adj, adj2v, DG, pts; TriangleType, k=i)
            end
        end
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
function jump_and_march(q::I, adj::Adjacent{I,E}, adj2v::Adjacent2Vertex{I,Es,E}, DG, pts;
    pt_idx=eachindex(pts), m=ceil(Int64, length(pt_idx)^(1 / 3)),
    k=select_initial_point(pts, q; m, pt_idx),
    TriangleType::Type{V}=NTuple{3,Int64}) where {I,E,Es,V}
    return jump_and_march(get_point(pts, q), adj, adj2v, DG, pts; k, TriangleType)
end