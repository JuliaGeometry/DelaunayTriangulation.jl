struct VoronoiTessellation{CCs,PGs,TTI,ITT,GN,CTP}
    circumcenters::CCs
    polygons::PGs
    triangle_to_idx::TTI
    idx_to_triangle::ITT
    generators::GN
    circumcenter_to_polygons::CTP
end

triangle_type(::VoronoiTessellation{CCs,PGs,TTI,ITT,GN,CTP}) where {CCs,PGs,TTI,ITT,GN,CTP} = triangle_type(TTI)
triangle_type(::Type{Dict{V,I}}) where {V,I<:Integer} = V
triangle_type(::Type{Dict{I,V}}) where {V,I<:Integer} = V
number_type(vorn::VoronoiTessellation) = number_type(vorn.generators)

function add_polygon!(cent_to_poly::Dict{I,E}, r, i) where {I,E}
    existing_polys = get!(E, cent_to_poly, r)
    push!(existing_polys, i)
    return nothing
end# r is the circumcenter, i is the polygon

"""
    get_triangle_idx(vorn::VoronoiTessellation, T)

Returns the index `i` such that `vorn.triangle_to_idx[T] = i`, 
where `T` may need to be shifted (i.e. `(i, j, k)` might need to become 
`(j, k, i)`).
"""
function get_triangle_idx(vorn::VoronoiTessellation, T)
    if T ∈ keys(vorn.triangle_to_idx)
        return vorn.triangle_to_idx[T]
    end
    T = shift_triangle_1(T)
    if T ∈ keys(vorn.triangle_to_idx)
        return vorn.triangle_to_idx[T]
    end
    T = shift_triangle_1(T)
    return vorn.triangle_to_idx[T]
end
function get_triangle_idx(vorn::VoronoiTessellation, i, j, k)
    V = triangle_type(vorn)
    T = construct_triangle(V, i, j, k)
    return get_triangle_idx(vorn, T)
end

"""
    get_circumcenter(vorn::VoronoiTessellation, T)

Returns the circumcenter for the triangle `T`, or the `T`th if `T` is an niteger.
"""
function get_circumcenter(vorn::VoronoiTessellation, T)
    i = get_triangle_idx(vorn, T)
    return vorn.circumcenters[i]
end
function get_circumcenter(vorn::VoronoiTessellation, i, j, k)
    _i = get_triangle_idx(vorn, i, j, k)
    return vorn.circumcenters[_i]
end
function get_circumcenter(vorn::VoronoiTessellation, i::Integer)
    return vorn.circumcenters[i]
end

function circumcenter(a, b, c)#wikipedia
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
    if !is_ghost_triangle(T)
        i, j, k = indices(T)
        pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
        cent_xy = circumcenter(pᵢ, pⱼ, pₖ)
        F = number_type(pts)
        return NTuple{2,F}(cent_xy)
    else
        V = rotate_ghost_triangle_to_boundary_form(T)
        u, v, _ = indices(V)
        pᵤ, pᵥ = get_point(pts, u, v)
        cent_x = 0.5 * (getx(pᵤ) + getx(pᵥ))
        cent_y = 0.5 * (gety(pᵤ) + gety(pᵥ))
        F = number_type(pts)
        return NTuple{2,F}((cent_x, cent_y))
    end
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
function get_voronoi_vertices(T::Ts, adj::Adjacent{I,E}, DG, i, tri_to_idx, cent_to_polygons) where {I,E,Ts}
    V = triangle_type(Ts)
    num_ngh = num_neighbours(DG, i)
    cell_idx = zeros(I, num_ngh)
    j = rand(get_neighbour(DG, i)) # just picking a random neighbour
    for r in 1:num_ngh
        k = get_edge(adj, i, j) # note that we are going counter-clockwise 
        if k ≠ I(BoundaryIndex) # need to take care of the unbounded triangles
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
                add_polygon!(cent_to_polygons, cent_idx, i)
            else
                cell_idx[r] = BoundaryIndex
                add_polygon!(cent_to_polygons, BoundaryIndex, i)
            end
        else
            cell_idx[r] = BoundaryIndex
            add_polygon!(cent_to_polygons, BoundaryIndex, i)

        end
        j = k
    end
    return cell_idx
end

function voronoi(T, adj::Adjacent{I,E}, adj2v, DG, pts; trim=false) where {I,E}
    has_ghosts_flag = triangulation_has_ghost_triangles(adj, adj2v)
    if !has_ghosts_flag
        add_ghost_triangles!(T, adj, adj2v, DG)
    end
    cents, tri_to_idx, idx_to_tri = circumcenters(T, pts)
    polys = Dict{I,Vector{Int64}}()
    cent_to_polygons = Dict{I,Set{Int64}}()
    sizehint!(polys, num_points(pts))
    for i in _eachindex(pts)
        polys[i] = get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_polygons)
    end
    vorn = VoronoiTessellation(cents, polys, tri_to_idx, idx_to_tri, pts, cent_to_polygons)
    if trim
        truncate_bounded_polygon_outside_convex_hull(vorn, T, adj, adj2v, DG)
        trim_voronoi_cell!(vorn, adj, DG)
    end
    if !has_ghosts_flag
        remove_ghost_triangles!(T, adj, adj2v, DG)
    end
    return vorn
end

function trim_voronoi_cell!(vorn::VoronoiTessellation, adj, DG)
    for i in get_neighbour(DG, BoundaryIndex)
        trim_voronoi_cell!(vorn, adj, i)
    end
    return nothing
end
function trim_voronoi_cell!(vorn::VoronoiTessellation, adj, i::Integer)
    delete!(vorn.circumcenter_to_polygons[BoundaryIndex], i)
    # Define the polygon and find where the unbounded part starts
    poly = vorn.polygons[i]
    idx = find_first_boundary_index(poly)
    if idx === nothing
        return nothing
    end
    next_idx = nextindex_circular(poly, idx)
    # Truncate the first unbounded line at the midpoint 
    boundary_point = get_edge(adj, BoundaryIndex, i)
    midpt_idx = get_triangle_idx(vorn, i, boundary_point, BoundaryIndex)
    poly[idx] = midpt_idx
    add_polygon!(vorn.circumcenter_to_polygons, midpt_idx, i)
    # Now truncate the second bounded line at the midpoint 
    boundary_point = get_edge(adj, i, BoundaryIndex)
    midpt_idx = get_triangle_idx(vorn, i, BoundaryIndex, boundary_point)
    poly[next_idx] = midpt_idx
    add_polygon!(vorn.circumcenter_to_polygons, midpt_idx, i)
    # Now connect the two lines by going through the generator 
    push!(vorn.circumcenters, get_point(vorn.generators, i))
    insert!(poly, next_idx, length(vorn.circumcenters))
    add_polygon!(vorn.circumcenter_to_polygons, length(vorn.circumcenters), i)
    return nothing
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

"""
    find_circumcenters_outside_convex_hull(vorn::VoronoiTessellation, adj::Adjacent{I,E}, adj2v, DG) where {I,E}

Given a Voronoi tessellation `vorn`, an adjacent map `adj`, an adjacent-to-vertex map `adj2v`, and the
Delaunay graph `DG`, finds all the circumcenters that lie outside of the convex hull of the generators. 
The result returned is a `Dict` which maps exterior circumcenters to the edge that they are next to. 
"""
function find_circumcenters_outside_convex_hull(vorn::VoronoiTessellation, T, adj::Adjacent{I,E}, adj2v, DG) where {I,E}
    outside_circ_idx = Dict{I,E}()
    for τ in T
        is_ghost_triangle(τ) && continue # Ghost triangles put the circumcenter at the edge
        i = get_triangle_idx(vorn, τ)
        p = vorn.circumcenters[i]
        u, v, w = indices(τ)
        pu, pv, pw = get_point(vorn.generators, u, v, w)
        o = isintriangle(pu, pv, pw, p) # First check if the circumference is inside the triangle it came from 
        if o ≥ 0
            continue
        end
        k = select_initial_point(vorn.generators, p; try_points=(u, v, w))
        V = jump_and_march(p, adj, adj2v, DG, vorn.generators; k=k)
        if is_ghost_triangle(V)
            V_rotated = rotate_ghost_triangle_to_boundary_form(V)
            u, v = geti(V_rotated), getj(V_rotated)
            outside_circ_idx[i] = construct_edge(E, u, v)
        end
    end
    return outside_circ_idx
end

"""
    get_neighbouring_circumcenters(vorn::VoronoiTessellation, i, adj)

Given a circumcenter index `i`, a Voronoi tessellation `vorn`, and an adjacent map `adj` for
the dual triangulation, finds the neighbouring circumcenters. The results are listed in 
counter-clockwise order.

Only works on the untrimmed version of the tessellation.
"""
function get_neighbouring_circumcenters(vorn::VoronoiTessellation, i, adj)
    τ = vorn.idx_to_triangle[i]
    i, j, k = indices(τ)
    u = get_edge(adj, j, i) # Remember that the polygons are defined by connecting neighbouring circumcenters
    v = get_edge(adj, k, j)
    w = get_edge(adj, i, k)
    e1 = u == BoundaryIndex
    e2 = v == BoundaryIndex
    e3 = w == BoundaryIndex # Don't want to include ghost triangles with their actual index, they should go to BoundaryIndex
    if any((e1, e2, e3))
        w, u, v = choose_uvw(e1, e2, e3, u, v, w) # This puts the BoundaryIndex into w
        k, i, j = choose_uvw(e1, e2, e3, i, j, k) # Need to rotate i, j, k in the same way as above
    end
    a = get_triangle_idx(vorn, j, i, u)
    b = get_triangle_idx(vorn, k, j, v)
    if any((e1, e2, e3))
        c = BoundaryIndex
    else
        c = get_triangle_idx(vorn, i, k, w)
    end
    return (a, b, c)
end

"""
    find_intersections_of_exterior_circumcenters_with_convex_hull(vorn::VoronoiTessellation, exterior_circumcenters::Dict{I,E}, adj) where {I,E}

Given some circumcenters that are outside of the triangulation's convex hull, finds the coordinates of the edges connecting that circumcenter to its 
neighbours with the convex hull. If no such intersection exists, stores `(NaN, NaN)`. The returned result is a `Dict` mapping the circumcenter indices 
(those in `exterior_circumcenters`) to the two coordinates.
"""
function find_intersections_of_exterior_circumcenters_with_convex_hull(vorn::VoronoiTessellation, exterior_circumcenters::Dict{I,E}, adj) where {I,E}
    F = number_type(vorn)
    intersection_coordinates = Dict{I,NTuple{3,NTuple{2,F}}}()
    for (j, (u, v)) in exterior_circumcenters
        ## Get the coordinates of the relevant polygonal edges
        a, b, c = get_neighbouring_circumcenters(vorn, j, adj)
        pa, pb, pc = get_point(vorn.circumcenters, a, b, c)
        pu, pv = get_point(vorn.generators, u, v)
        pj = get_point(vorn.circumcenters, j)
        ## We only need to have the coordinates if the edges intersect the convex hull 
        o1 = meet(pu, pv, pj, pa)
        o2 = meet(pu, pv, pj, pb)
        if c ≠ BoundaryIndex
            o3 = meet(pu, pv, pj, pc)
        end
        ## Compute the intersection coordinates 
        intersection_uv_a = o1 == 1 ? intersection_of_two_line_segments(pu, pv, pj, pa) : pa .* NaN
        intersection_uv_b = o2 == 1 ? intersection_of_two_line_segments(pu, pv, pj, pb) : pb .* NaN
        if c ≠ BoundaryIndex
            intersection_uv_c = o3 == 1 ? intersection_of_two_line_segments(pu, pv, pj, pc) : pc .* NaN
        else
            intersection_uv_c = pc .* NaN
        end
        ## Store the coordinates 
        intersection_coordinates[j] = (intersection_uv_a, intersection_uv_b, intersection_uv_c)
    end
    return intersection_coordinates
end

"""
    move_exterior_circumcenter_to_convex_hull_intersection!(vorn::VoronoiTessellation, exterior_circumcenters, 
        intersection_coordinates) 

Given a Voronoi tessellation `vorn`, a `Dict` of `exterior_circumcenters` from [`find_circumcenters_outside_convex_hull`](@ref), a `Dict` of 
coordinates `intersection_coordinates` from `find_intersections_of_exterior_circumcenters_with_convex_hull`](@ref), and a circumcenter index `j` 
(assumed to be in `keys(exterior_circumcenters)`), updates the polygons so that the edge going to this exterior circumcenter are cut off at the 
intersection. This update is done in-place. A value is returned, giving a `Dict` that maps polygons to the number of intersections of its edges 
with the convex hull.
"""
function move_exterior_circumcenter_to_convex_hull_intersection!(vorn::VoronoiTessellation, exterior_circumcenters,
    intersection_coordinates::Dict{I,NTuple{3,NTuple{2,F}}}, num_exterior_verts, j::I) where {I,F}
    ## Start by moving the circumenters to the intersection coordinates
    u, v = exterior_circumcenters[j]
    P = vorn.polygons[u]
    Q = vorn.polygons[v]
    pa, pb, pc = intersection_coordinates[j]
    if !isnan(first(pa))
        num_exterior_verts[u] = u ∈ keys(num_exterior_verts) ? num_exterior_verts[u] + 1 : 1
    end
    if !isnan(first(pb))
        num_exterior_verts[v] = v ∈ keys(num_exterior_verts) ? num_exterior_verts[v] + 1 : 1
    end
    pa_index = length(vorn.circumcenters) + 1
    pb_index = length(vorn.circumcenters) + 2
    push!(vorn.circumcenters, pa)
    push!(vorn.circumcenters, pb)
    P[findfirst(==(j), P)] = pa_index
    Q[findfirst(==(j), Q)] = pb_index
    ## The bounded polygon then needs to be closed
    return nothing
end
function move_exterior_circumcenter_to_convex_hull_intersection!(vorn::VoronoiTessellation, exterior_circumcenters,
    intersection_coordinates::Dict{I,NTuple{3,NTuple{2,F}}}) where {I,F}
    num_exterior_verts = Dict{I,I}()
    for j in keys(exterior_circumcenters)
        move_exterior_circumcenter_to_convex_hull_intersection!(vorn, exterior_circumcenters, intersection_coordinates, num_exterior_verts, j)
    end
    return num_exterior_verts
end

"""
    get_modified_polygons(exterior_circumcenters)

Given a list of `exterior_circumcenters`, returns a `Set` of polygons to be modified.
"""
function get_modified_polygons(exterior_circumcenters)
    modified_polygons = Set{Int64}()
    for (u, v) in values(exterior_circumcenters)
        push!(modified_polygons, u, v)
    end
    return modified_polygons
end

"""
    close_polygon!(vorn::VoronoiTessellation, added_generator_indices, num_exterior_vertices, u)

Given a Voronoi tessellation `vorn` and a generator `u` corresponding to an unbounded polygon `vorn.polygons[u]`, closes it by 
replacing the two unbounded edges with a segment going through the generator. The argument `added_generator_indices` maps a generator 
to its position in `vorn.circumcenters`, modifying in-place if the generator is not already in `vorn.circumcenters`. The 
argument `num_exterior_vertices` is the result from [`get_number_exterior_vertices`](@ref). Note that if a polygon only 
has one intersection with an exterior vertex, rather than with two edges, an intersection is introduced at the midpoint of an edge 
on the convex hull to close through.
"""
function close_polygon!(vorn::VoronoiTessellation, added_generator_indices, num_exterior_vertices, u)
    # Start by adding the generators into the circumcenter list if needed
    P = vorn.polygons[u]
    if u ∉ keys(added_generator_indices)
        pu = get_point(vorn.generators, u)
        push!(vorn.circumcenters, pu)
        added_generator_indices[u] = length(vorn.circumcenters)
    end
    # Now close the polygon
    if num_exterior_vertices[u] == 2
        # If there are two intersections, then the only way to close the polygon is to go through the generator
        pu_index = added_generator_indices[u]
        P_boundary_idx = find_first_boundary_index(P)
        P_second_boundary_idx = nextindex_circular(P, P_boundary_idx)
        P[P_boundary_idx] = pu_index
        deleteat!(P, P_second_boundary_idx)
    else

    end
end

"""
    polygon_is_bounded(vorn, i)

Given a Voronoi tessellation `vorn` and a polygon `i`, 
"""
function polygon_is_bounded(vorn, i)
    return BoundaryIndex ∉ vorn.polygons[i]
end

"""
    find_bounded_polygons(vorn::VoronoiTessellation, i)

Given a Voronoi tessellation `vorn` and a circumcenter `i`, finds 
the bounded triangles that the circumcenter is a vertex of.
"""
function find_bounded_polygons(vorn::VoronoiTessellation, i::I; find_first=false) where {I}
    polys = vorn.circumcenter_to_polygons[i]
    if !find_first
        bounded_polys = I[]
    end
    for p in polys
        if !find_first
            polygon_is_bounded(vorn, p) && push!(bounded_polys, p)
        else
            polygon_is_bounded(vorn, p) && return p
        end
    end
    return bounded_polys
end

"""
    truncate_bounded_polygon_outside_convex_hull(vorn::VoronoiTessellation, adj, adj2v, DG)

Given the Voronoi tessellation `vorn`, modifies all bounded polygons that have a circumcenter outside of the domain 
such that the lines intersecting through the convex hull are chopped off at the intersection point, defining a new 
vertex for the neighbouring polygons. This is done for each circumcenter outside of the convex hull, as computed 
using [`find_circumcenters_outside_convex_hull`](@ref).
"""
function truncate_bounded_polygon_outside_convex_hull(vorn::VoronoiTessellation, T, adj, adj2v, DG)
    exterior_circumcenters = find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, DG)
    added_points = Dict{Int64,Int64}() # to avoid adding the same generator twice when cleaning up the unbounded polygons
    modified_polygons = Dict{Int64,Int64}() # Store all the polygons to be modified after the loop. The value stores the number of times the polygon was encountered.
    orig_num_circumcenters = length(vorn.circumcenters)
    for _j in exterior_circumcenters
        pⱼ = vorn.circumcenters[_j]
        # Start by finding the polygon that is the bounded one. Note that every vertex of the Voronoi tessellation has degree 3
        @show exterior_circumcenters
        bounded_polygon_idx = find_bounded_polygons(vorn, _j; find_first=true)
        bounded_polygon = vorn.polygons[bounded_polygon_idx]
        _j_idx = findfirst(bounded_polygon .== _j) # Finds where _j is in the list of vertices for the bounded polygon
        # Find the circumcenters that connect with _j 
        next_j_idx = nextindex_circular(bounded_polygon, _j_idx)
        prev_j_idx = previndex_circular(bounded_polygon, _j_idx)
        next_idx = bounded_polygon[next_j_idx]
        prev_idx = bounded_polygon[prev_j_idx]
        # Now get the triangle that is on the boundary, and rotate it so that its first two vertices give the boundary edge 
        boundary_τ = vorn.idx_to_triangle[_j]
        u, v, w = indices(boundary_τ)
        ub, vb, wb = is_boundary_point.((u, v, w), Ref(adj), Ref(DG))
        if ub && vb
            u, v, w = u, v, w
        elseif ub && wb
            u, v, w = w, u, v
        elseif vb && wb
            u, v, w = v, w, u
        elseif ub && vb && wb
            u, v, w = u, v, w
        end
        # Get the coordinates for the boundary vertices 
        pᵤ, pᵥ = get_point(vorn.generators, u, v)
        # Get the coordinates for the neighbouring circumcenters 
        pᵣ = vorn.circumcenters[next_idx]
        pₛ = vorn.circumcenters[prev_idx]
        if u ∉ keys(added_points)
            push!(vorn.circumcenters, pᵤ)
            added_points[u] = length(vorn.circumcenters)
        end
        if v ∉ keys(added_points)
            push!(vorn.circumcenters, pᵥ)
            added_points[v] = length(vorn.circumcenters)
        end
        # The lines connecting the neighbouring vertices to pⱼ. Let us find these coordinates 
        coord_1 = intersection_of_two_line_segments(pᵤ, pᵥ, pⱼ, pᵣ)
        coord_2 = intersection_of_two_line_segments(pᵤ, pᵥ, pⱼ, pₛ)
        coord_1_idx = length(vorn.circumcenters) + 1
        coord_2_idx = length(vorn.circumcenters) + 2
        push!(vorn.circumcenters, coord_1)
        push!(vorn.circumcenters, coord_2)
        # Now, pⱼ is connected to other unbounded polygons P and Q. We need to update these so that the circumcenter pⱼ is replaced by the intersections we found above 
        other_polys = setdiff(vorn.circumcenter_to_polygons[_j], bounded_polygon_idx)
        P, Q = other_polys
        modified_polygons[P] = P ∈ keys(modified_polygons) ? modified_polygons[P] + 1 : 1
        modified_polygons[Q] = Q ∈ keys(modified_polygons) ? modified_polygons[Q] + 1 : 1
        P_verts = vorn.polygons[P]
        Q_verts = vorn.polygons[Q]
        delete!(vorn.circumcenter_to_polygons[_j], P)
        delete!(vorn.circumcenter_to_polygons[_j], Q)
        P_j = findfirst(P_verts .== _j)
        Q_j = findfirst(Q_verts .== _j)
        if next_idx ∈ P_verts   # ⟹ prev_idx ∈ Q.
            P_verts[P_j] = coord_1_idx
            add_polygon!(vorn.circumcenter_to_polygons, coord_1_idx, P)
            Q_verts[Q_j] = coord_2_idx
            add_polygon!(vorn.circumcenter_to_polygons, coord_2_idx, Q)
        else                    # prev_idx ∈ P_verts ⟹ next_idx ∈ Q. 
            P_verts[P_j] = coord_2_idx
            @show P, coord_2_idx
            add_polygon!(vorn.circumcenter_to_polygons, coord_2_idx, P)
            Q_verts[Q_j] = coord_1_idx
            add_polygon!(vorn.circumcenter_to_polygons, coord_1_idx, Q)
        end
        # The bounded polygon is now modified so that the circumcenter is replaced by the two intersection points, taking care to preserve the order of the points 
        bounded_polygon[_j_idx] = coord_2_idx
        insert!(bounded_polygon, next_j_idx, coord_1_idx)
    end
    # Now we need to close up the bounded polygons, closing them at the intersection points. The lines that connect the intersection points across edges are connected across the generator in between.
    # We need to be careful, not all of the edges of P and Q need to be chopped yet - if the other edge connects with a circumcenter 
    # that is inside the convex hull (meaning not in exterior_circumcenters), then we should do nother.
    @show modified_polygons, added_points, orig_num_circumcenters
    for (P, num_exterior_vertices) in modified_polygons
        verts = vorn.polygons[P]
        idx = find_first_boundary_index(verts)
        idx_2 = nextindex_circular(verts, idx)
        if num_exterior_vertices == 2
            verts[idx] = added_points[P] # This is the generator 
            deleteat!(verts, idx_2)
            add_polygon!(vorn.circumcenter_to_polygons, added_points[P], P)
        else # Chop the edge that goes to the exterior circumcentre, but the other one gets taken to the midpoint
            j = findfirst(verts .> orig_num_circumcenters) # This will be the one that was taken to a new midpoint
            prev_j_idx = previndex_circular(verts, j)
            if verts[prev_j_idx] == BoundaryIndex
                boundary_point = get_edge(adj, BoundaryIndex, P)
                midpt_idx = get_triangle_idx(vorn, P, boundary_point, BoundaryIndex)
                verts[idx] = midpt_idx
                verts[idx_2] = added_points[P] # This is the generator that connects the truncated points 
                add_polygon!(vorn.circumcenter_to_polygons, midpt_idx, P)
                add_polygon!(vorn.circumcenter_to_polygons, added_points[P], P)
            else
                boundary_point = get_edge(adj, P, BoundaryIndex)
                midpt_idx = get_triangle_idx(vorn, P, BoundaryIndex, boundary_point)
                verts[idx_2] = midpt_idx
                verts[idx] = added_points[P] # This is the generator that connects the truncated points 
                add_polygon!(vorn.circumcenter_to_polygons, midpt_idx, P)
                add_polygon!(vorn.circumcenter_to_polygons, added_points[P], P)
            end
        end
        delete!(vorn.circumcenter_to_polygons[BoundaryIndex], P)
    end
    return nothing
end
