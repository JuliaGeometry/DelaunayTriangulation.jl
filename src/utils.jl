"""
    is_boundary_edge(i, j, adj::Adjacent{I, E}) where {I, E}
    is_boundary_edge(i, j, adj2v::Adjacent{I, E}) where {I, E}

Returns `true` if the edge `(i, j)` is a boundary edge of the triangulation, 
and `false` otherwise.
"""
function is_boundary_edge(ij, adj::Adjacent{I,E}) where {I,E}
    return get_edge(adj, ij) == I(BoundaryIndex)
end
function is_boundary_edge(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj)
end
function is_boundary_edge(ij, adj2v::Adjacent2Vertex{I,Es,E}) where {I,Es,E}
    return ij ∈ get_edge(adj2v, I(BoundaryIndex))
end
function is_boundary_edge(i, j, adj2v::Adjacent2Vertex{I,Es,E}) where {I,Es,E}
    ij = construct_edge(E, i, j)
    return is_boundary_edge(ij, adj2v)
end

"""
    is_boundary_triangle(T, adj)

Returns `true` if `T` has any edges on the boundary, and `false` otherwise. Note that this is 
different from testing if `T` is a ghost triangle - it is assumed that `T` is a solid triangle.
"""
function is_boundary_triangle(T, adj)
    i, j, k = indices(T)
    return is_boundary_edge(k, j, adj) ||
           is_boundary_edge(j, i, adj) ||
           is_boundary_edge(i, k, adj)

end

"""
    is_ghost_edge(i::I, j::I) where {I, E}

Returns `true` if the edge `(i, j)` is a ghost edge of the triangulation, 
and `false` otherwise.
"""
function is_ghost_edge(i::I, j::I) where {I}
    return i == I(BoundaryIndex) || j == I(BoundaryIndex)
end

"""
    is_ghost_triangle(T) 

Returns `true` if the triangle `T` is a ghost triangle.
"""
function is_ghost_triangle(i, j, k)
    return is_ghost_edge(i, j) || is_ghost_edge(j, k) || is_ghost_edge(k, i)
end
function is_ghost_triangle(T)
    i, j, k = indices(T)
    return is_ghost_triangle(i, j, k)
end

"""
    edge_exists(i, j, adj::Adjacent{I, E}) where {I, E}

Returns `true` if the edge `(i, j)` is an edge in the triangulation, 
and `false` otherwise.
"""
function edge_exists(ij, adj::Adjacent{I,E}) where {I,E}
    return get_edge(adj, ij) ≠ I(DefaultAdjacentValue)
end
function edge_exists(i, j, adj::Adjacent{I,E}) where {I,E}
    ij = construct_edge(E, i, j)
    return edge_exists(ij, adj)
end

"""
    is_boundary_point(u, adj, DG::DelaunayGraph{I}) where {I}

Tests if the point `u` is on the boundary of the triangulation with graph `DG` and
adjacent map `adj`.
"""
function is_boundary_point(u, adj, DG::DelaunayGraph{I}) where {I}
    if I(BoundaryIndex) ∈ graph(DG).V
        return u ∈ get_neighbour(DG, I(BoundaryIndex))
    else
        for v in get_neighbour(DG, u)
            is_boundary_edge(u, v, adj) && return true
        end
        return false
    end
end

"""
    choose_uvw(e1, e2, e3, i, j, k)

Chooses values for `(u, v, w)` based on the Booleans 
`(e1, e2, e3)`, assuming only one is true. 

- If `e1`, returns `(i, j, k)`
- If `e2`, returns `(j, k, i)`
- If `e3`, returns `(k, i, j)`
"""
function choose_uvw(e1, e2, e3, i, j, k)
    if e1
        return i, j, k
    elseif e2
        return j, k, i
    elseif e3
        return k, i, j
    end
end

"""
    rotate_triangle_to_boundary_form(T)

This rotates the ghost triangle `T = (i, j, k)` so that it takes the form 
`T = (u, v, $BoundaryIndex)`. This change is done out-of-place.
"""
function rotate_ghost_triangle_to_boundary_form(i, j, k)
    e1 = is_ghost_edge(i, j)
    e2 = is_ghost_edge(j, k)
    e3 = is_ghost_edge(k, i)
    u, v, w = choose_uvw(!e1, !e2, !e3, i, j, k)
    return u, v, w
end
function rotate_ghost_triangle_to_boundary_form(T::V) where {V}
    i, j, k = indices(T)
    u, v, w = rotate_ghost_triangle_to_boundary_form(i, j, k)
    rotated_T = construct_triangle(V, u, v, w)
    return rotated_T
end

"""
    is_delaunay(T::Ts, pts) where {Ts}

Tests if the triangulation `T` is Delaunay by checking that the open circumdisk of each triangle contains no 
points inside it. This check is slow since we intentionally check every point instead of only 
checking the neighbouring triangles, ensuring that no mistakes have been made e.g. with overlapping 
triangles. Returns `true` if the triangulation is indeed Delaunay, and `false` otherwise.
"""
function is_delaunay(T::Ts, pts) where {Ts}
    V = triangle_type(Ts)
    I = integer_type(V)
    for T in T
        for r in _eachindex(pts)
            r = I(r)
            isincircle(T, pts, r) == -1
        end
    end
    return true
end

"""
    validate_triangulation(T, adj::Adjacent{I,E}, adj2v, DG, pts) where {I,E}

Tests that all the triangles in `T` are positively oriented, if all the triangles 
in `T` are Delaunay, if `adj` and `adj2v` are inverses of each other, and if 
all the edges in `DG` are in `adj`. Returns `true` if so, and `false` otherwise.
"""
function validate_triangulation(T, adj::Adjacent{I,E}, adj2v, DG, pts) where {I,E}
    _adj = deepcopy(adj)
    clear_empty_keys!(_adj)
    ## Check triangle orientation 
    for T in T
        isoriented(T, pts) == -1 && (@show "Not all triangles are positively oriented."; return false)
    end
    ## Test that all edges are legal 
    !is_delaunay(T, pts) && (@show "Not all edges are legal."; return false)
    ## Check that the adjacent map and the adjacent-to-vertex map are inverses
    !check_adjacent_is_adjacent2vertex_inverse(_adj, adj2v) && (@show "The Adjacent and Adjacent2Vertex maps are not inverses of each other."; return false)
    ## Check the graph 
    for (i, j) in graph(DG).E
        if triangulation_has_ghost_triangles(_adj, adj2v)
            (i, j) ∉ edges(_adj) && (@show "The adjacent map is missing an edge."; return false)
            (j, i) ∉ edges(_adj) && (@show "The adjacent map is missing an edge."; return false)
        else
            if i ≠ I(BoundaryIndex) && j ≠ I(BoundaryIndex)
                (i, j) ∉ edges(_adj) && (@show "The adjacent map is missing an edge."; return false)
                (j, i) ∉ edges(_adj) && (@show "The adjacent map is missing an edge."; return false)
            end
        end
    end
    ## Check the edges
    for T in T
        i, j, k = indices(T)
        get_edge(_adj, i, j) ≠ k && (@show "An adjacent key is incorrect."; return false)
        get_edge(_adj, j, k) ≠ i && (@show "An adjacent key is incorrect."; return false)
        get_edge(_adj, k, i) ≠ j && (@show "An adjacent key is incorrect."; return false)
    end
    ## Check that every edge is incident to two triangles if it is not a boundary edge, and one otherwise
    for (i, j) in edges(_adj)
        vij = get_edge(_adj, i, j)
        vji = get_edge(_adj, j, i)
        if is_boundary_edge(i, j, adj)
            vij ≠ I(BoundaryIndex) && (@show "An adjacent key is incorrect."; return false)
        elseif is_boundary_edge(j, i, adj)
            vji ≠ I(BoundaryIndex) && (@show "An adjacent key is incorrect."; return false)
        elseif i ≠ I(BoundaryIndex) && j ≠ I(BoundaryIndex)
            (vij < I(LowerLeftBoundingIndex) || vji < I(LowerLeftBoundingIndex)) && (@show "An adjacent key is incorrect."; return false)
        end
    end
    ## Done 
    return true
end

"""
    clear_empty_keys!(adj::Adjacent)

Deletes all the keys `(i, j)` in `adj` such that `get_edge(adj, i, j) = $(DefaultAdjacentValue)`.
"""
function clear_empty_keys!(adj::Adjacent{I,E}) where {I,E}
    for ((i, j), k) in adjacent(adj)
        if k == I(DefaultAdjacentValue)
            delete_edge!(adj, i, j)
        end
    end
    return nothing
end

"""
    clear_empty_keys!(adj2v::Adjacent2Vertex)

Deletes any keys from `adj2v` that map to empty sets.
"""
function clear_empty_keys!(adj2v::Adjacent2Vertex)
    for (w, S) in adjacent2vertex(adj2v)
        if isempty(S)
            delete_point!(adj2v, w)
        end
    end
    return nothing
end

"""
    clear_empty_points!(DG::DelaunayGraph)

Deletes points from the graph `DG` that have degree 0. 
"""
function clear_empty_points!(DG::DelaunayGraph)
    for (w, S) in graph(DG).N
        if isempty(S)
            delete_point!(DG, w)
        end
    end
    return nothing
end

"""
    compare_triangle_sets(T, V)

Compares the collections of triangles `T` and `V`, checking for equality. This function is 
needed to check for equality even with circularly-shifted indices in the triangles. Returns `true` 
if the two collections are equivalent, and false otherwise.
"""
function compare_triangle_sets(T::Ts, V::Vs) where {Ts,Vs}
    VsT = triangle_type(Vs)
    length(T) ≠ length(V) && return false
    for T in T
        i, j, k = indices(T)
        _T = construct_triangle(VsT, i, j, k)
        (_T ∉ V && shift_triangle_1(_T) ∉ V && shift_triangle_2(_T) ∉ V) && return false
    end
    return true
end

"""
    compare_unconstrained_triangulations(T1, adj1, adj2v1, DG1, T2, adj2, adj2v2, DG2)

Compares the triangulation `(T1, adj1, adj2v1, DG1)` to the triangulation `(T2, adj2, adj2v2, DG2)`.
Could cause issues with triangulations of collinear points depending on the choices made during 
the algorithm.
"""
function compare_unconstrained_triangulations(T1, adj1, adj2v1, DG1, T2, adj2, adj2v2, DG2)
    adj = deepcopy(adj1)
    _adj = deepcopy(adj2)
    clear_empty_keys!(adj)
    clear_empty_keys!(_adj)
    adj2v = deepcopy(adj2v1)
    _adj2v = deepcopy(adj2v2)
    clear_empty_keys!(adj2v)
    clear_empty_keys!(_adj2v)
    e1 = compare_triangle_sets(T1, T2)
    e2 = adjacent(adj) == adjacent(_adj)
    e3 = adjacent2vertex(adj2v) == adjacent2vertex(_adj2v)
    e4 = graph(DG1) == graph(DG2)
    !e1 && (@show "Unequal triangle sets."; return false)
    !e2 && (@show "Unequal adjacent maps."; return false)
    !e3 && (@show "Unequal adjacent2vertex maps."; return false)
    !e4 && (@show "Unequal graphs."; return false)
    return true
end

"""
    setdiff_triangles(T::Ts, _T::Ts) where {Ts}

Allows `setdiff` to be used for triangle sets, checking for equalities under 
circular shifts.
"""
function setdiff_triangles(T::Ts, _T::Ts) where {Ts}
    V = sort_triangles(T)
    _V = sort_triangles(_T)
    setdiff!(V, _V)
    return V
end

"""
    check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)

Checks if `adj` and `adj2v` are related so that, if `get_edge(adj, i, j) = k`,
then `(i, j) ∈ get_edge(adj2v, k)`. Returns `true` if so, and `false` otherwise.
"""
function check_adjacent_is_adjacent2vertex_inverse(adj::Adjacent{I,E}, adj2v) where {I,E}
    # Check adj2v 
    for (k, S) in adjacent2vertex(adj2v)
        for ij in S
            get_edge(adj, ij) ≠ k && (@show (ij, k); return false)
        end
    end
    # Check adj
    for (ij, k) in adjacent(adj)
        if k ≠ I(DefaultAdjacentValue)
            ij ∉ get_edge(adj2v, k) && (@show (ij, k); return false)
        end
    end
    return true
end

"""
    clear_centroid_coordinates!(::Type{T} = Float64) where {T}

Resets the centroid coordinates in `CentroidCoordinates` to be 
`(zero(T), zero(T))`.
"""
function clear_centroid_coordinates!(::Type{T}=Float64) where {T}
    CentroidCoordinates.x = zero(T)
    CentroidCoordinates.y = zero(T)
    CentroidCoordinates.n = 0
    return nothing
end

"""
    update_centroid_after_new_point!(pts, i)

Updates the centroid coordinates in `CentroidCoordinates` after a new point 
`get_point(pts, i)` was added.
"""
function update_centroid_after_new_point!(pts, i)
    new_pt = get_point(pts, i)
    n = CentroidCoordinates.n
    CentroidCoordinates.x = 1 / (n + 1) * (n * CentroidCoordinates.x + getx(new_pt))
    CentroidCoordinates.y = 1 / (n + 1) * (n * CentroidCoordinates.y + gety(new_pt))
    CentroidCoordinates.n += 1
    return nothing
end

"""
    update_centroid_after_deleted_point!(pts, i)

Updates the centroid coordinates in `CentroidCoordinates` after the point
`get_point(pts, i)` was deleted.
"""
function update_centroid_after_deleted_point!(pts, i)
    new_pt = get_point(pts, i)
    n = CentroidCoordinates.n
    CentroidCoordinates.x = 1 / (n - 1) * (n * CentroidCoordinates.x - getx(new_pt))
    CentroidCoordinates.y = 1 / (n - 1) * (n * CentroidCoordinates.y - gety(new_pt))
    CentroidCoordinates.n -= 1
    return nothing
end

"""
    compute_centroid!(pts)

Updates `CentroidCoordinates` with the centroid of `pts`.
"""
function compute_centroid!(pts)
    n = num_points(pts)
    clear_centroid_coordinates!()
    for i in _eachindex(pts)
        p = get_point(pts, i)
        x = getx(p)
        y = gety(p)
        CentroidCoordinates.x += x
        CentroidCoordinates.y += y
    end
    CentroidCoordinates.x /= n
    CentroidCoordinates.y /= n
    CentroidCoordinates.n = n
    return nothing
end

"""
    triangulation_has_ghost_triangles(adj::Adjacent{I, E}, adj2v) where {I, E}

Tests if the triangulation with adjacent map `adj` and adjacent-to-vertex map `adj2v` contains 
ghost triangles, returning `true` if so and `false` otherwise. This test is done by testing some 
element of `get_edge(adj2v, $BoundaryIndex)`, say `(u, v)`, and then seeing if `(v, $BoundaryIndex)`
is a valid key in `adj`.
"""
function triangulation_has_ghost_triangles(adj::Adjacent{I,E}, adj2v) where {I,E}
    u, v = iterate(get_edge(adj2v, I(BoundaryIndex)))[1]
    return edge_exists(v, I(BoundaryIndex), adj)
end

"""
    add_ghost_triangles!(T::Ts, adj::Adjacent{I, E}, adj2v, DG) where {Ts, I, E}

Adds ghost triangles to the triangulation `(T, adj, adj2v, DG)`. These structures 
are updated in-place.
"""
function add_ghost_triangles!(T::Ts, adj::Adjacent{I,E}, adj2v, DG) where {Ts,I,E}
    V = triangle_type(Ts)
    for (u, v) in get_edge(adj2v, I(BoundaryIndex))
        add_edge!(adj, v, I(BoundaryIndex), u)
        add_edge!(adj, I(BoundaryIndex), u, v)
        add_edge!(adj2v, u, v, I(BoundaryIndex))
        add_edge!(adj2v, v, I(BoundaryIndex), u)
        #add_neighbour!(DG, I(BoundaryIndex), u, v)
        Tg = construct_triangle(V, u, v, I(BoundaryIndex))
        add_triangle!(T, Tg)
    end
    return nothing
end

"""
    remove_ghost_triangles!(T::Ts, adj::Adjacent{I, E}, adj2v, DG) where {Ts, I, E}

Removes ghost triangles from the triangulation `(T, adj, adj2v, DG)`. These structures 
are updated in-place.
"""
function remove_ghost_triangles!(T::Ts, adj::Adjacent{I,E}, adj2v, DG) where {Ts,I,E}
    V = triangle_type(Ts)
    for (u, v) in get_edge(adj2v, I(BoundaryIndex))
        delete_edge!(adj, v, I(BoundaryIndex))
        delete_edge!(adj, I(BoundaryIndex), u)
        delete_edge!(adj2v, u, v, I(BoundaryIndex))
        delete_edge!(adj2v, v, I(BoundaryIndex), u)
        Tg = construct_triangle(V, u, v, I(BoundaryIndex))
        delete_triangle!(T, Tg)
    end
    #delete_point!(DG, I(BoundaryIndex))
    return nothing
end

"""
    compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)

Compares the results from a triangulation using de Berg's method, `(_T, _adj, _adj2v, _DG)`, to those 
with the Bowyer-Watson method, `(T, adj, adj2v, DG)`. Ghosts triangles are added to the data structures 
from de Berg's method and then compared to the Bowyer-Watson results, and then we remove the ghost triangles 
from both and compare. If both are identical, return `true`, and `false` otherwise. Any changes are made 
out-of-place.
"""
function compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    Tbw, adjbw, adj2vbw, DGbw = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
    Tdb, adjdb, adj2vdb, DGdb = deepcopy(_T), deepcopy(_adj), deepcopy(_adj2v), deepcopy(_DG)
    add_ghost_triangles!(Tdb, adjdb, adj2vdb, DGdb)
    e1 = compare_unconstrained_triangulations(Tbw, adjbw, adj2vbw, DGbw, Tdb, adjdb, adj2vdb, DGdb)
    remove_ghost_triangles!(Tdb, adjdb, adj2vdb, DGdb)
    remove_ghost_triangles!(Tbw, adjbw, adj2vbw, DGbw)
    e2 = compare_unconstrained_triangulations(Tbw, adjbw, adj2vbw, DGbw, Tdb, adjdb, adj2vdb, DGdb)
    return e1 && e2
end
"""
    compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, pts)

Compares the results from a triangulation using the Bowyer-Watson method, `(T, adj, adj2v, DG)`, using the point 
set `pts` to those from de Berg's method. If both are identical, return `true`, and `false` otherwise. Any changes are made 
out-of-place.
"""
function compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, pts)
    (_T, _adj, _adj2v, _DG), _ = triangulate_berg(pts)
    return compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
end

"""
    select_valid_start_of_vector!(pt_order, skip_pts; r=3)

Circularly shifts parts the vector `pt_order` in place until `pt_order[begin:(begin+r-1)]` contains none 
of the elements in `skip_pts`. In particular, we shift until the first element contains no element 
in `skip_pts`, and then we shift `pt_order[(begin+1):end]` until the second element of the original vector 
contains no element in `skip_pts`, and so on.
"""
function select_valid_start_of_vector!(pt_order, skip_pts; r=3)
    for i in 0:(r-1)
        j = 0
        while pt_order[begin+i] ∈ skip_pts
            @views circshift!(pt_order[(begin+i):end], -1)
            j += 1
        end
    end
    return nothing
end

"""
    sort_triangle(T::V) where {V}

Given a triangle `T = (i, j, k)`, sorts it so that the first index is the smallest, maintaining 
the orientation of `T`.
"""
function sort_triangle(T::V) where {V}
    i, j, k = indices(T)
    minijk = min(i, j, k)
    if minijk == i
        return construct_triangle(V, i, j, k)
    elseif minijk == j
        return construct_triangle(V, j, k, i)
    else
        return construct_triangle(V, k, i, j)
    end
end

"""
    sort_triangles(T::Ts) where {Ts}

Sorts the triangles in the collection `T` so that each triangle's first vertex 
has the smallest value. The orientation of each triangle is preserved.
"""
function sort_triangles(T::Ts) where {Ts}
    tris = construct_triangles(Ts)
    for τ in T
        σ = sort_triangle(τ)
        add_triangle!(tris, σ)
    end
    return tris
end

"""
    remove_duplicate_triangles(T) 

Removes duplicate triangles from `T`. This procedure also sorts the triangles 
so that the first index of each triangle is the smallest. Orientations are 
preserved.
"""
function remove_duplicate_triangles(T)
    V = sort_triangles(T)
    if T isa Set
        return V
    else
        unique!(V)
        return V
    end
end

"""
    sort_boundary!(θ, idx_vec, all_pts, boundary_idx, ref)
    sort_boundary!(all_pts, boundary_idx, ref)

Sorts the boundary defined by `all_pts` at the indices `boundary_idx`, 
using the reference point `ref` to sort. The sorting is counter-clockwise.
The vectors `θ`, `idx_vec`, and `boundary_idx` are all modified in-place. The 
second method can be used to avoid creating these cache vectors manually, noting 
that `boundary_idx` will still be modified in-place.

# Arguments 
- `θ`: The vector which will store the angles from `ref` to the points specified by `boundary_idx`.
- `idx_vec`: This will be updated in-place with the sorting order. 
- `all_pts`: The complete point set. 
- `boundary_idx`: Indices in `all_pts` corresponding to the boundary to be sorted. 
- `ref`: The reference point to compute angles about. 

# Outputs 
There are no outputs.
"""
function sort_boundary!(θ, idx_vec, all_pts, boundary_idx, ref)
    cx, cy = ref
    for (j, i) in pairs(boundary_idx)
        x, y = get_point(all_pts, i)
        θ[j] = atan(y - cy, x - cx) + π
        if i == BoundaryIndex
            θ[j] = (θ[j] + π) % 2π # Remember that this point is the centroid if we are at the boundary index, so reflect it across the boundary edge (add +π) to get the correct angle
        end
    end
    sortperm!(idx_vec, θ)
    permute!(boundary_idx, idx_vec)
    return nothing
end
function sort_boundary!(all_pts, boundary_idx, ref)
    θ = zeros(number_type(all_pts), length(boundary_idx))
    idx_vec = zeros(Int64, length(boundary_idx))
    sort_boundary!(θ, idx_vec, all_pts, boundary_idx, ref)
    return nothing
end

"""
    area(pts)

Computes the area of the polygon defined by `pts`, assuming the points are given in 
counter-clockwise order, using the shoelace formula.
"""
function area(pts)
    F = number_type(pts)
    A = zero(F)
    _indices = _eachindex(pts)
    for i in _indices
        pᵢ_current = get_point(pts, i)
        pᵢ_prev = get_point(pts, i == first(_indices) ? last(_indices) : i - 1)
        pᵢ_next = get_point(pts, i == last(_indices) ? first(_indices) : i + 1)
        xᵢ = getx(pᵢ_current)
        yᵢ₊₁ = gety(pᵢ_next)
        yᵢ₋₁ = gety(pᵢ_prev)
        A += xᵢ * (yᵢ₊₁ - yᵢ₋₁)
    end
    return 0.5A
end

"""
    find_first_boundary_index(v)

Given a vector `v`, find the first `i` such that `v[i] == $BoundaryIndex`, treating the 
vector as circular. It is assumed that `v` only has two elements such that `v[i] == $BoundaryIndex`, and 
that these elements always appear next to each other. If no such index exists, returns `nothing`.
"""
function find_first_boundary_index(v)
    idx = findfirst(v .== BoundaryIndex)
    isnothing(idx) && return idx
    #prev_idx = idx == firstindex(v) ? lastindex(v) : idx - 1
    prev_idx = previndex_circular(v, idx)
    if v[prev_idx] == BoundaryIndex
        return prev_idx
    else
        return idx
    end
end

"""
    corss_2d(u, v)

Computes the product `2d_cross(u, v) = ux * vy - uy * vx`.
"""
function cross_2d(u, v)
    ux, uy = u
    vx, vy = v
    return ux * vy - uy * vx
end

"""
    intersection_of_two_line_segments(p, q, r, s)

Finds the intersection point of the two line segments `pq` and `rs`, assuming 
that an intersection exists. 

https://stackoverflow.com/a/565282
"""
function intersection_of_two_line_segments(p, q, r, s)
    rpsr = cross_2d(r .- p, s .- r)
    qpsr = cross_2d(q .- p, s .- r)
    t = rpsr / qpsr
    return p .+ (q .- p) .* t
end

"""
    nextindex_circular(v, i)

Given a vector `v` and an index `i`, finds the index `i+1`, where `lastindex(v) + 1` gets mapped to `firstindex(v)`.
"""
function nextindex_circular(v, i)
    return i == lastindex(v) ? firstindex(v) : i + 1
end

"""
    previndex_circular(v, i)

Given a vector `v` and an index `i`, finds the index `i-1`, where `firstindex(v) - 1` gets mapped to `lastindex(v)`.
"""
function previndex_circular(v, i)
    return i == firstindex(v) ? lastindex(v) : i - 1
end

"""
    is_vertex_of(T, i)

Returns `true` if `i` is a vertex of `T`, and `false` otherwise.
"""
function is_vertex_of(T, i)
    return i == geti(T) || i == getj(T) || i == getk(T)
end

"""
    compare_triangles(T, V)

Let `T = (u, v, w)` and `V = (i, j, k)`. This function sees if `T == V`, 
returning `true` if so, where we say that `T == V` if any of the following 
hold:

- `(u, v, w) == (i, j, k)`,
- `(v, w, u) == (i, j, k)`,
- `(w, u, v) == (i, j, k)`.

Returns false otherwise.
"""
function compare_triangles(T, V)
    if (geti(T), getj(T), getk(T)) == (geti(V), getj(V), getk(V))
        return true
    elseif (getj(T), getk(T), geti(T)) == (geti(V), getj(V), getk(V))
        return true
    elseif (getk(T), geti(T), getj(T)) == (geti(V), getj(V), getk(V))
        return true
    end
    return false
end

"""
    all_points_are_unique(pts)

Returns `true` if all the points in `pts` are unique, and `false` otherwise. 
See also `Base.allunique`.
"""
function all_points_are_unique(pts) # see also allunique
    # Need to get all the points as tuples in a vector so that allunique can be used 
    vec_pts = [get_point(pts, i) for i in _eachindex(pts)]
    return allunique(vec_pts)
end