"""
    exterior_jump_and_march(
        pts, 
        adj::Adjacent{I,E}, 
        boundary_index_ranges, 
        boundary_map, 
        k, 
        q, 
        check_existence::V=Val(has_multiple_segments(boundary_map))) where {I,E,V}

Given a collection of points `pts`, an [`Adjacent`](@ref) map `adj`, a `Dict` mapping boundary indices to ranges from 
[`construct_boundary_index_ranges`](@ref), a boundary map 
from [`construct_boundary_map`](@ref), a vertex `k`, and a point `q` outside of 
the triangulation, returns the edge `(i, j)` such that the ghost triangle 
`(i, j, $BoundaryIndex)` contains `q`. 

The `check_existence` argument is used to check that the edge exists when using [`get_adjacent`](@ref), 
helping to correct for incorrect boundary indices in the presence of multiple boundary segments. See [`get_adjacent`](@ref).

The result is meaningless if `q` is inside of the triangulation.
"""
function exterior_jump_and_march(
    pts,
    adj::Adjacent{I,E},
    boundary_index_ranges,
    boundary_map,
    k,
    q,
    check_existence::V=Val(has_multiple_segments(boundary_map))
) where {I,E,V}
    pₘ, pᵢ = get_point(pts, boundary_map, I(BoundaryIndex), k)
    i = k
    q_position = point_position_relative_to_line(pₘ, pᵢ, q)
    if is_left(q_position) # q is left of the ghost edge through pᵢ, so rotate left 
        j = get_right_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges,
            check_existence)
        pⱼ = get_point(pts, boundary_map, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        while is_left(new_q_pos)
            i = j
            pᵢ = pⱼ
            j = get_right_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges,
                check_existence)
            pⱼ = get_point(pts, boundary_map, j)
            new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        end
        return j, i # Swap the orientation so that i, j is a boundary edge 
    else # rotate right 
        j = get_left_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges,
            check_existence)
        pⱼ = get_point(pts, boundary_map, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        while is_right(new_q_pos)
            i = j
            pᵢ = pⱼ
            j = get_left_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges,
                check_existence)
            pⱼ = get_point(pts, boundary_map, j)
            new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        end
        return i, j
    end
end

"""
    jump_and_march(pts, adj, adj2v, graph::Graph{I}, boundary_index_ranges, boundary_map, q;
        m=default_num_samples(num_points(pts)),
        point_indices=each_point_index(pts),
        try_points=(),
        k=select_initial_point(pts, q; m, point_indices, try_points),
        TriangleType::Type{V}=NTuple{3,I},
        check_existence::C=Val(has_multiple_segments(boundary_map)),
        store_visited_triangles::F=Val(false),
        visited_triangles=nothing,
        rng::AbstractRNG = Random.default_rng()) where {I,V,C,F}

Using the jump and march algorithm, finds the triangle in the triangulation containing the 
query point `q`.

If your triangulation does not have ghost triangles, and the point `q` is outside of the triangulation, 
this function may fail to terminate. You may like to add ghost triangles in this case (using 
[`add_ghost_triangles!`](@ref)), noting that there is no actual triangle that `q` is inside 
when it is outside of the triangulation unless ghost triangles are present. 

# Arguments 
- `pts`: The collection of points.
- `adj`: The [`Adjacent`](@ref) map.
- `graph::Graph{I}`: The [`Graph`](@ref).
- `boundary_index_ranges`: The `Dict` mapping boundary indices to ranges from [`construct_boundary_index_ranges`](@ref).
- `boundary_map`: The boundary map from [`construct_boundary_map`](@ref) handling the mapping of boundary indices. 
- `q`: The query point.

# Keyword Arguments 
- `m=default_num_samples(num_points(pts))`: The number of samples to use when sampling an initial point from [`select_initial_point`](@ref). Only relevant if `k` is not specified. 
- `point_indices`: The indices for the points. Only relevant if `k` is not specified. 
- `try_points=()`: Extra points to try when sampling an initial point from [`select_initial_point`](@ref). Only relevant if `k` is not specified. 
- `k=select_initial_point(pts, q; m, point_indices, try_points)`: Where to start the algorithm.
- `TriangleType::Type{V}=NTuple{3,I}`: The type used for representing the triangles. 
- `check_existence::C=Val(has_multiple_segments(boundary_map))`: Whether to check that the edge exists when using [`get_adjacent`](@ref), helping to correct for incorrect boundary indices in the presence of multiple boundary segments. See [`get_adjacent`](@ref).
- `stored_visited_triangles::F=Val(false)`: Whether to store visited triangles. Exterior ghost triangles will not be stored.
- `visited_triangles=nothing`: Object to push visited triangles into.
- `rng::AbstractRNG`: The RNG to use.

# Output 
Returns the triangle `V` containing the query point `q`.
"""
function jump_and_march(pts, adj, adj2v, graph::Graph{I}, boundary_index_ranges,
    boundary_map, q;
    m=default_num_samples(num_points(pts)),
    point_indices=each_point_index(pts),
    try_points=(),
    k=select_initial_point(pts, q; m, point_indices, try_points),
    TriangleType::Type{V}=NTuple{3,I},
    check_existence::C=Val(has_multiple_segments(boundary_map)),
    store_visited_triangles::F=Val(false),
    visited_triangles=nothing,
    rng::AbstractRNG=Random.default_rng()) where {I,V,C,F}
    return _jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges, boundary_map, q,
        k, V, check_existence, store_visited_triangles, visited_triangles, rng)
end
function _jump_and_march(pts, adj, adj2v, graph::Graph{I}, boundary_index_ranges,
    boundary_map, q, k, TriangleType::Type{V},
    check_existence::C=Val(has_multiple_segments(boundary_map)),
    store_visited_triangles::F=Val(false),
    visited_triangles=nothing,
    rng::AbstractRNG=Random.default_rng()) where {I,V,C,F}
    if !is_outer_boundary_node(k, graph, boundary_index_ranges) || !has_ghost_triangles(adj, adj2v)
        # If k is not a boundary node, then we can rotate around the point k to find an initial triangle 
        # to start the search. Additionally, if there are no ghost triangles, we can only hope that q 
        # is inside the interior, meaning we should only search for the initial triangle here anyway.
        p, i, j, pᵢ, pⱼ = select_initial_triangle_interior_node(pts, adj, adj2v,
            boundary_map, k, q,
            boundary_index_ranges,
            check_existence, rng)
        is_true(store_visited_triangles) && add_triangle!(visited_triangles, j, i, k)
    else
        # We have an outer boundary node. First, let us check the neighbouring boundary edges. 
        direction, q_pos, next_vertex, right_cert, left_cert = check_for_intersections_with_adjacent_boundary_edges(pts,
            adj,
            boundary_index_ranges,
            boundary_map,
            k,
            q,
            check_existence)
        if !is_outside(direction)
            # q is collinear with one of the edges, so let's jump down these edges and try to find q
            q_pos, u, v, w = search_down_adjacent_boundary_edges(pts, adj,
                boundary_index_ranges,
                boundary_map, k, q,
                direction, q_pos,
                next_vertex,
                check_existence)
            if is_on(q_pos)
                is_true(store_visited_triangles) && add_triangle!(visited_triangles, u, v, w)
                return construct_triangle(V, u, v, w)
            else
                u, v = exterior_jump_and_march(pts, adj, boundary_index_ranges,
                    boundary_map, u, q, check_existence)
                return construct_triangle(V, u, v,
                    get_adjacent(adj, u, v; check_existence,
                        boundary_index_ranges)) # Can't just use I(BoundaryIndex) here since there could be multiple - just use get_adjacent
            end
        end
        # If we did not find anything from the neighbouring boundary edges, we can search the neighbouring interior edges
        i, j, edge_cert, triangle_cert = check_for_intersections_with_interior_edges_adjacent_to_boundary_node(pts,
            adj,
            graph,
            boundary_index_ranges,
            boundary_map,
            k,
            q,
            right_cert,
            left_cert,
            check_existence)
        if is_inside(triangle_cert)
            is_true(store_visited_triangles) && add_triangle!(visited_triangles, i, j, k)
            return construct_triangle(V, i, j, k)
        elseif is_none(edge_cert)
            u, v = exterior_jump_and_march(pts, adj, boundary_index_ranges, boundary_map, k,
                q, check_existence)
            return construct_triangle(V, u, v, get_adjacent(adj, u, v))
        else
            is_true(store_visited_triangles) && add_triangle!(visited_triangles, j, i, k)
            p, pᵢ, pⱼ = get_point(pts, boundary_map, k, i, j)
        end
    end
    if q == p || q == pᵢ || q == pⱼ
        orientation = triangle_orientation(pᵢ, pⱼ, p)
        if is_positively_oriented(orientation)
            return construct_triangle(V, i, j, k)
        else
            return construct_triangle(V, j, i, k)
        end
    end
    ## Now let us do the straight line search 
    # The idea is to keep jumping until pᵢpⱼ goes past q, meaning pᵢpⱼq is no longer a positively oriented triangle
    arrangement = triangle_orientation(pᵢ, pⱼ, q)
    local last_changed # Need this for deciding which variable to use when we hit a collinear point 
    last_changed = I(DefaultAdjacentValue) # just an initial value
    while is_positively_oriented(arrangement)
        # We need to step forward. To do this, we need to be careful of boundary indices. 
        # Since a boundary curve may be represented by different boundary indices 
        # at different locations, we need to check for this. The first check 
        # is for an outer boundary index.
        if is_outer_boundary_index(k, boundary_map)
            # If this happens, it means we have hit an outer boundary edge, and so we need to go into the exterior. If there are no 
            # ghost triangles, though, we just restart the search. Note that interior boundary indices do not matter since the ghost 
            # triangles there have the same orientation, so we can find them as normal.
            if has_ghost_triangles(adj, adj2v)
                i, j = exterior_jump_and_march(pts, adj, boundary_index_ranges,
                    boundary_map, last_changed == i ? j : i, q,
                    check_existence) # use last_changed to ensure we get the boundary point
                return construct_triangle(V, i, j, k)
            else
                return _jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges,
                    boundary_map, q, k, V, check_existence, store_visited_triangles, visited_triangles, rng)
            end
        end
        # Now we can finally move forward. We use check_existence to protect against the issue mentioned above.
        k = get_adjacent(adj, i, j; check_existence, boundary_index_ranges)
        pₖ = get_point(pts, boundary_map, k)
        pₖ_pos = point_position_relative_to_line(p, q, pₖ)
        is_true(store_visited_triangles) && add_triangle!(visited_triangles, i, j, k)
        if is_right(pₖ_pos)
            j, pⱼ = k, pₖ
            last_changed = j
        elseif is_left(pₖ_pos)
            i, pᵢ = k, pₖ
            last_changed = i
        else
            # pₖ is collinear with pq. We will first check if q is already in the current triangle
            in_cert = point_position_relative_to_triangle(pᵢ, pⱼ, pₖ, q) # ... Maybe there is a better way to compute this, reusing previous certificates? Not sure. ...
            if !is_outside(in_cert)
                return construct_triangle(V, i, j, k)
            end
            # To decide which direction this collinear point is away from the line, we can jsut use the last changed variable:
            # If last_changed = i, this means that the left direction was what caused the collinear point, so make k go on the left. 
            # If not, make it go to the right.
            if last_changed == i
                i, pᵢ = k, pₖ
                last_changed = i
            else
                j, pⱼ = k, pₖ
                last_changed = j
            end
            is_true(store_visited_triangles) && add_triangle!(visited_triangles, i, j, k)
        end
        arrangement = triangle_orientation(pᵢ, pⱼ, q)
    end
    # We can finish the above loop even if q is not in the triangle, in which case pᵢpⱼq was a straight line. 
    # To clear this up, let us just restart. 
    if is_degenerate(arrangement)
        pₖ = get_point(pts, boundary_map, k) # Need to have this here in case we skip the entire loop, above, meaning pₖ won't exist
        in_cert = point_position_relative_to_triangle(pᵢ, pⱼ, pₖ, q) # ... Maybe there is a better way to compute this, reusing previous certificates? Not sure. ...
        is_true(store_visited_triangles) && add_triangle!(visited_triangles, i, j, k)
        if is_outside(in_cert)
            return _jump_and_march(pts, adj, adj2v, graph, boundary_index_ranges,
                boundary_map, q,
                last_changed == I(DefaultAdjacentValue) ? i :
                last_changed, V, check_existence, store_visited_triangles, visited_triangles, rng)
        end
    end
    # Swap the orientation to get a positively oriented triangle, remembering that we kept pᵢ on the left of pq and pⱼ on the right 
    k = get_adjacent(adj, j, i; check_existence, boundary_index_ranges)
    is_true(store_visited_triangles) && add_triangle!(visited_triangles, j, i, k)
    return construct_triangle(V, j, i, k)
end
