"""
    delete_point!(tri::Triangulation, point; rng::AbstractRNG=Random.default_rng())

Deletes `point` from the triangulation `tri` using Chew's algorithm.
"""
function delete_point!(tri::Triangulation, point; rng::AbstractRNG=Random.default_rng())
    ## Start by determining if the node is an outer boundary node or not 
    ## We need to do this here because the information we delete immediately 
    ## after makes it impossible to classify the node.
    is_outer_bnd = is_outer_boundary_node(tri, point) # Interior ghost vertices have a physical location, outer ghost vertices are at infinity.
    nn = num_neighbours(tri, point)
    if nn == 3 && is_outer_bnd
        _delete_point_2!(tri, point)
        return nothing
    end
    ## In case we do have a boundary node, we need to replace the RepresentativeCoordinates now.
    ## Again, we can't defer this due to the information we delete afterwards. 
    if is_outer_bnd
        cx, cy = _replace_representative_coordinates_for_deletion!(tri, point)
    end
    ## Now begin deleting the neighbouring triangles 
    S = get_surrounding_polygon(tri, point) # Need to get this first before we delete anything
    neighbouring_edges = get_adjacent2vertex(tri, point)
    for uv in each_edge(neighbouring_edges)
        u = initial(uv)
        v = terminal(uv)
        delete_triangle!(tri, point, u, v; protect_boundary=true)
    end
    ## Now fill in the cavity
    if is_outer_bnd
        _delete_boundary_point!(tri, S, rng)
    else
        _delete_interior_point!(tri, S, rng)
    end
    ## Delete the point information entirely
    delete_adjacent2vertex!(tri, point)
    delete_vertex!(tri, point)
    ## Cleanup
    push!(S, S[begin])
    is_outer_bnd && _replace_representative_coordinates_after_deletion!(tri, cx, cy)
    fix_edges_after_deletion!(tri, S)
    return nothing
end

function _delete_interior_point!(tri, S, rng)
    triangulate_convex!(tri, S; rng, update_ghost_edges=false) # Works even for boundary nodes, due to the replacement we make above
    return nothing
end

function _delete_boundary_point!(tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR}, S, rng) where {P,Ts,I,E,Es,BN,BNM,B,BIR}
    local boundary_index
    local boundary_index_in_S
    for (i, s) in pairs(S)
        if is_boundary_index(s)
            S[i] = I(BoundaryIndex) # This function is only used for outer boundary indices, and since we use a replacement below, we need 
            boundary_index = s # to represent all boundary indices on the same curve as a single index. This does this for us.
            boundary_index_in_S = i
            break
        end
    end
    temp_tri = triangulate_convex(get_points(tri), S;
        IntegerType=I,
        EdgeType=E,
        TriangleType=triangle_type(Ts),
        EdgesType=Es,
        TrianglesType=Ts,
        representative_point_list = get_representative_point_list(tri),
        rng,
        add_ghost_triangles=false,
        add_convex_hull=false,
        compute_centers=false,
        delete_empty_features=false)
    for T in each_triangle(temp_tri)
        if is_ghost_triangle(T)
            V = rotate_ghost_triangle_to_standard_form(T)
            u, v, _ = indices(V) # ghost index is last
            T = construct_triangle(triangle_type(Ts), u, v, boundary_index)
        end
        add_triangle!(tri, T; protect_boundary=true)
    end
    # Need to also fix S since we use it later 
    S[boundary_index_in_S] = boundary_index
    return Set(S)
end

function _replace_representative_coordinates_for_deletion!(tri, point)
    I = integer_type(tri)
    boundary_index = I(BoundaryIndex)
    cx, cy = get_point(tri, boundary_index) # same thing as e.g. get_representative_point_coordinates(1, Float64)
    left_bnd = get_left_boundary_node(tri, point, boundary_index)
    right_bnd = get_right_boundary_node(tri, point, boundary_index)
    pℓx, pℓy = get_point(tri, left_bnd)
    prx, pry = get_point(tri, right_bnd)
    mx, my = get_point(tri, point)
    mc_length = sqrt((mx - cx)^2 + (my - cy)^2)
    pℓpr_length = sqrt((pℓx - prx)^2 + (pℓy - pry)^2)
    scaled_extent = (mc_length + 25pℓpr_length) / mc_length
    ĉx = cx + scaled_extent * (mx - cx)
    ĉy = cy + scaled_extent * (my - cy)
    representative_coords = get_representative_point_list(tri)[I(1)] # the outer boundary also has a curve index of 1 
    representative_coords.x = ĉx
    representative_coords.y = ĉy
    return cx, cy
end

function _replace_representative_coordinates_after_deletion!(tri, cx, cy)
    I = integer_type(tri)
    representative_coords = get_representative_point_list(tri)[I(1)]
    representative_coords.x = cx
    representative_coords.y = cy
    return nothing
end

function fix_edges_after_deletion!(tri, S)
    for i in firstindex(S):(lastindex(S)-1)
        u = S[i]
        v = S[i+1]
        w = get_adjacent(tri, v, u)
        if edge_exists(w)
            add_adjacent!(tri, u, w, v)
            add_adjacent!(tri, w, v, u)
            add_neighbour!(tri, u, w)
        end
    end
end

function _delete_point_2!(tri::Triangulation, point)
    I = integer_type(tri)
    boundary_index = I(BoundaryIndex)
    left_bnd = get_left_boundary_node(tri, point, boundary_index)
    right_bnd = get_right_boundary_node(tri, point, boundary_index)
    delete_triangle!(tri, point, right_bnd, left_bnd; protect_boundary=false, update_ghost_edges=true)
    delete_adjacent2vertex!(tri, point)
    delete_vertex!(tri, point)
    return nothing
end

"""
    get_surronding_polygon(tri, u; skip_boundary_indices=false)

Returns the set of neighbours of `u` in counter-clockwise order. If `skip_boundary_indices` is `true`, then boundary indices are not included in the set.
"""
@inline get_surrounding_polygon(tri::Triangulation, u; skip_boundary_indices=false) = get_surrounding_polygon(get_adjacent(tri), get_graph(tri), u, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)); skip_boundary_indices)
