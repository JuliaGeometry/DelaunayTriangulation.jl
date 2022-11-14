function circumcenter(a, b, c)
    ax, ay = a
    bx, by = b
    cx, cy = c
    D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    norma_sq = ax^2 + ay^2
    normb_sq = bx^2 + by^2
    normc_sq = cx^2 + cy^2
    cent_x = 1 / D * (norma_sq * (by - cy) + normb_sq * (cy - ay) + normc_sq * (ay - by))
    cent_y = 1 / D * (norma_sq * (cx - bx) + normb_sq * (ax - cx) + normc_sq * (bx - ax))
    return (cent_x, cent_y)
end
function circumcenter(T, pts)
    i, j, k = indices(T)
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    cent_xy = circumcenter(pᵢ, pⱼ, pₖ)
    F = number_type(pts)
    return NTuple{2,F}(cent_xy)
end
function circumcenters(T::Ts, pts) where {Ts}
    V = triangle_type(Ts)
    F = number_type(pts)
    I = integer_type(V)
    cents = Vector{NTuple{2,F}}(undef, length(T))
    triangle_to_idx = Dict{V,I}(T .=> 1:length(T))
    idx_to_triangle = Dict{I,V}(1:length(T) .=> T)
    for (τ, i) in triangle_to_idx
        cents[i] = circumcenter(τ, pts)
    end
    return cents, triangle_to_idx, idx_to_triangle
end

"""
    get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)

Given a generator `i` and a set of triangles `T`, an adjacent map `adj`, 
an adjacent-to-vertex map `adj2v`, and a Delaunay graph `DG`
(all from a Delaunay triangulation), computes the indices corresponding to the 
vertices of the associated Voronoi cell, in counter-clockwise order. The 
indices of these vertices are defined by each triangle's circumcenter, whose index 
is obtained from the mapping `tri_to_idx` which maps a triangle to this index.
"""
function get_voronoi_vertices(T::Ts, adj::Adjacent{I,E}, adj2v, DG, i, tri_to_idx) where {I,E,Ts}
    has_ghosts_flag = true
    if !triangulation_has_ghost_triangles(adj, adj2v)
        has_ghosts_flag = false
        add_ghost_triangles!(T, adj, adj2v, DG)
    end
    V = triangle_type(Ts)
    num_ngh = num_neighbours(DG, i)
    cell_idx = zeros(I, num_ngh)
    j = rand(get_neighbour(DG, i))
    for r in 1:num_ngh
        k = get_edge(adj, i, j)
        if k ≠ I(BoundaryIndex)
            if !is_ghost_triangle(i, j, k)
                τ = construct_triangle(V, i, j, k)
                if τ ∉ T
                    τ = shift_triangle_1(τ)
                end
                if τ ∉ T
                    τ = shift_triangle_1(τ)
                end
                cent_idx = tri_to_idx[τ]
                cell_idx[r] = cent_idx
            else
                cell_idx[r] = BoundaryIndex
            end
        else
            cell_idx[r] = BoundaryIndex
        end
        j = k
    end
    if !has_ghosts_flag
        remove_ghost_triangles!(T, adj, adj2v, DG)
    end
    return cell_idx
end

struct VoronoiTessellation{CCs,PGs,TTI,ITT}
    circumcenters::CCs
    polygons::PGs
    triangle_to_idx::TTI
    idx_to_triangle::ITT
end

function voronoi(T, adj::Adjacent{I,E}, adj2v, DG, pts) where {I,E}
    has_ghosts_flag = triangulation_has_ghost_triangles(adj, adj2v)
    if has_ghosts_flag
        remove_ghost_triangles!(T, adj, adj2v, DG)
    end
    cents, tri_to_idx, idx_to_tri = circumcenters(T, pts)
    polys = Dict{I,Vector{Int64}}()
    sizehint!(polys, num_points(pts))
    for i in _eachindex(pts)
        polys[i] = get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
    end
    return VoronoiTessellation(cents, polys, tri_to_idx, idx_to_tri)
end

function area(vorn::VoronoiTessellation, i)
    idx = vorn.polygons[i]
    F = number_type(vorn.circumcenters)
    if BoundaryIndex ∈ idx
        return typemax(F)
    else
        pts = @views vorn.circumcenters[idx]
        return area(pts)
    end
end
function area(vorn::VoronoiTessellation)
    F = number_type(vorn.circumcenters)
    areas = zeros(F, length(vorn.polygons))
    for i in keys(vorn.polygons)
        areas[i] = area(vorn, i)
    end
    return areas
end