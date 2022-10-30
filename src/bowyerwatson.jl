function add_point_bowyer!(T::Ts, adj::Adjacent{I,E}, adj2v, DG, pts, r;
    pt_idx=graph(DG).V, m=ceil(Int64, length(pt_idx)^(1 / 3)),
    initial_search_point=select_initial_point(pts, r; pt_idx, m)) where {Ts,I,E} # only works for inside points currently. also fails for collinear points
    r = I(r)
    tri_type = triangle_type(Ts)
    V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point, TriangleType=tri_type)
    i, j, k = indices(V)
    delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true, update_ghost_edges=false)
    dig_cavity!(r, i, j, T, adj, adj2v, DG, pts)
    dig_cavity!(r, j, k, T, adj, adj2v, DG, pts)
    dig_cavity!(r, k, i, T, adj, adj2v, DG, pts)
    return nothing
end
function dig_cavity!(r, i, j, T::Ts, adj::Adjacent{I,E}, adj2v, DG, pts) where {Ts,I,E}
    ℓ = get_edge(adj, j, i) # (j, i, ℓ) is the triangle on the other side of the edge (i, j) from r 
    if ℓ == I(DefaultAdjacentValue)
        return nothing # The triangle has already been deleted in this case 
    end
    δ = ℓ ≠ I(BoundaryIndex) && isincircle(pts, r, i, j, ℓ)# Boundary edges form part of the boundary on the polygonal cavity if we have reached ℓ = BoundaryIndex. This will not cause issues with the ghost triangles since, if these ghost triangles are relevant, we will have started in one (since this means the point must be in the exterior), and so one of the dig cavity calls will have found it.  
    if δ == I(1)
        delete_triangle!(j, i, ℓ, T, adj, adj2v, DG; protect_boundary=true, update_ghost_edges=false) # Notice that we don't need the update_ghost_edges argument, since we will step over ghost edges already and add them into the triangulation if necessary
        dig_cavity!(r, i, ℓ, T, adj, adj2v, DG, pts)
        dig_cavity!(r, ℓ, j, T, adj, adj2v, DG, pts)
    else # Must be an edge of the polygonal cavity. Note that this also covers the ℓ ≠ I(BoundaryIndex) case
        #if is_ghost_triangle(i, j, ℓ) && isoriented(construct_triangle(triangle_type(Ts), r, i, j), pts) == 0
        #    split_edge!(j, i, r, T, adj, adj2v, DG) # i, j isn't there anymore
        #else
        add_triangle!(r, i, j, T, adj, adj2v, DG, update_ghost_edges=false)
        #end
    end
    return nothing
end

"""

"""
function triangulate_bowyer(pts;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{V}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    trim=true,
    try_last_inserted_point=true) where {I,E,V,Es,Ts}
    pt_order = randomise ? shuffle(_eachindex(pts)) : _eachindex(pts)
    T = Ts()
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    initial_triangle = construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
    while isoriented(initial_triangle, pts) == 0 # We cannot start with a degenerate triangle
        pt_order = collect(pt_order) # We need to collect so that we can change the order of the points 
        reverse!(@view pt_order[begin:(begin+2)]) # This makes it counter-clockwise
        initial_triangle = construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
    end
    u, v, w = indices(initial_triangle)
    add_triangle!(u, v, w, T, adj, adj2v, DG; update_ghost_edges=true)
    # compute_centroid!(@view pts[pt_order[begin:(begin+2)]]) # Can't use this, else ArrayPartitions support fails
    compute_centroid!((
        get_point(pts, pt_order[firstindex(pt_order)]),
        get_point(pts, pt_order[firstindex(pt_order)+1]),
        get_point(pts, pt_order[firstindex(pt_order)+2])
    ))
    for (num_points, new_point) in enumerate(@view pt_order[(begin+3):end])
        last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
        last_inserted_point = pt_order[last_inserted_point_number]
        pt_idx = @view pt_order[begin:last_inserted_point_number] # + 3 for the first three points already inserted
        m = ceil(Int64, length(pt_idx)^(1 / 3))
        initial_search_point = I(select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
        add_point_bowyer!(T, adj, adj2v, DG, pts, new_point;
            pt_idx, m, initial_search_point)
        update_centroid_after_new_point!(pts, new_point)
    end
    trim && remove_ghost_triangles!(T, adj, adj2v, DG)
    return T, adj, adj2v, DG
end

function lazy_triangulate_bowyer(pts;
    IntegerType::Type{I}=Int64,
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TriangleType::Type{VV}=NTuple{3,IntegerType},
    EdgesType::Type{Es}=Set{EdgeType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    randomise=true,
    trim=true,
    try_last_inserted_point=true) where {I,E,VV,Es,Ts}
    pt_order = randomise ? shuffle(_eachindex(pts)) : _eachindex(pts)
    T = Ts()
    adj = Adjacent{I,E}()
    adj2v = Adjacent2Vertex{I,Es,E}()
    DG = DelaunayGraph{I}()
    initial_triangle = construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
    while isoriented(initial_triangle, pts) == 0 # We cannot start with a degenerate triangle
        pt_order = collect(pt_order) # We need to collect so that we can change the order of the points 
        reverse!(@view pt_order[begin:(begin+2)]) # This makes it counter-clockwise
        initial_triangle = construct_positively_oriented_triangle(TriangleType, I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]), pts)
    end
    u, v, w = indices(initial_triangle)
    add_triangle!(u, v, w, T, adj, adj2v, DG; update_ghost_edges=true)
    compute_centroid!((
        get_point(pts, pt_order[firstindex(pt_order)]),
        get_point(pts, pt_order[firstindex(pt_order)+1]),
        get_point(pts, pt_order[firstindex(pt_order)+2])
    ))
    ## Loop over all the points 
    for (num_points, new_point) in enumerate(@view pt_order[(begin+3):end])
        pr = get_point(pts, new_point)
        ## Find all triangles that contain the new point in its circumdisk 
        bad_triangle_edge_list = Es()
        for V in T
            u, v, w = indices(V)
            pu, pv, pw = get_point(pts, u, v, w)
            ## Is the new point in the circumdisk of V?
            # First consider if V is a ghost triangle
            if u == I(BoundaryIndex)
                if orient(pv, pw, pr) ≥ 0 # Is pr left of (v, w)?
                    push!(bad_triangle_edge_list, (u, v), (v, w), (w, u))
                    delete_triangle!(u, v, w, T, adj, adj2v, DG; protect_boundary=true)
                end
            elseif v == I(BoundaryIndex)
                if orient(pw, pu, pr) ≥ 0 # Is pr left of (w, u)?
                    push!(bad_triangle_edge_list, (u, v), (v, w), (w, u))
                    delete_triangle!(u, v, w, T, adj, adj2v, DG; protect_boundary=true)
                end
            elseif w == I(BoundaryIndex)
                if orient(pu, pv, pr) ≥ 0 # Is pr left of (u, v)?
                    push!(bad_triangle_edge_list, (u, v), (v, w), (w, u))
                    delete_triangle!(u, v, w, T, adj, adj2v, DG; protect_boundary=true)
                end
            else # Otherwise, V is a solid triangle
                if incircle(pu, pv, pw, pr) == 1
                    push!(bad_triangle_edge_list, (u, v), (v, w), (w, u))
                    delete_triangle!(u, v, w, T, adj, adj2v, DG; protect_boundary=true)
                end
            end
        end
        ## Now we need to find the boundary of the cavity 
        polygon_edges = Es()
        for (i, j) in bad_triangle_edge_list
            if (j, i) ∉ bad_triangle_edge_list
                push!(polygon_edges, (i, j))
            end
        end
        ## Add the new triangles 
        for (u, v) in polygon_edges
            add_triangle!(u, v, new_point, T, adj, adj2v, DG; update_ghost_edges=false)
        end
        ## Update the centroid
        update_centroid_after_new_point!(pts, new_point)
    end
    ## Remove the ghost triangles and return 
    trim && remove_ghost_triangles!(T, adj, adj2v, DG)
    return T, adj, adj2v, DG
end