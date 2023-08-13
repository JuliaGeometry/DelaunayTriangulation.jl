"""
    exterior_jump_and_march(
        pts, 
        adj, 
        boundary_index_ranges, 
        representative_point_list,
        boundary_map, 
        k, 
        q, 
        check_existence=Val(has_multiple_segments(boundary_map)),
        bnd_idx=I(BoundaryIndex)) 

Given a point `q` outside of the triangulation, finds the ghost triangle containing it.

# Arguments 
- `pts`: The collection of points.
- `adj`: The [`Adjacent`](@ref) map.
- `boundary_index_ranges`: A `Dict` mapping boundary indices to ranges from [`construct_boundary_index_ranges`](@ref).
- `representative_point_list`: A `Dict` mapping curve indices to representative points.
- `boundary_map`: A boundary map from [`construct_boundary_map`](@ref).
- `k`: The vertex to start from.
- `q`: The point outside of the triangulation.
- `check_existence=Val(has_multiple_segments(boundary_map)))`: Used to check that the edge exists when using [`get_adjacent`](@ref), in case there are multiple segments.
- `bnd_idx=I(BoundaryIndex)`: A boundary index corresponding to the boundary that `k` is on.

# Ouptut 
- The edge `(i, j)` such that the ghost triangle `(i, j, g)` contains `q`, and `g = get_adjacent(adj, i, j)`.

!!! warning 

    The result is meaningless if `q` is inside of the triangulation.
"""
function exterior_jump_and_march(tri::Triangulation{P,Ts,I,E}, k, q, bnd_idx=I(BoundaryIndex)) where {P,Ts,I,E}
    pₘ, pᵢ = get_point(tri, bnd_idx, k)
    i = k
    q_position = point_position_relative_to_line(pₘ, pᵢ, q)
    if is_left(q_position) # q is left of the ghost edge through pᵢ, so rotate left 
        j = get_right_boundary_node(tri, i, bnd_idx)
        pⱼ = get_point(tri, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        while is_left(new_q_pos)
            i = j
            pᵢ = pⱼ
            j = get_right_boundary_node(tri, i, bnd_idx)
            pⱼ = get_point(tri, j)
            new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        end
        return j, i # Swap the orientation so that i, j is a boundary edge 
    else # rotate right 
        j = get_left_boundary_node(tri, i, bnd_idx)
        pⱼ = get_point(tri, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        while is_right(new_q_pos)
            i = j
            pᵢ = pⱼ
            j = get_left_boundary_node(tri, i, bnd_idx)
            pⱼ = get_point(tri, j)
            new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
        end
        return i, j
    end
end

"""
    jump_and_march(tri::Triangulation, q;
        point_indices=each_point_index(tri),
        m=default_num_samples(length(point_indices)),
        try_points=(),
        k=select_initial_point(get_points(tri), q; m, point_indices, try_points),
        check_existence::C=Val(has_multiple_segments(tri)),
        store_history::F=Val(false),
        history=nothing,
        rng::AbstractRNG=Random.default_rng(),
        exterior_curve_index=1,
        maxiters=2 + length(exterior_curve_index) - num_solid_vertices(tri) + num_solid_edges(tri)) where {C,F}

Returns the triangle containing `q` using the jump-and-march algorithm.

# Arguments 
- `tri::Triangulation`: The triangulation.
- `q`: The query point.

# Keyword Arguments
- `point_indices=each_point_index(tri)`: The indices of the points in the triangulation.
- `m=default_num_samples(length(point_indices))`: The number of samples to use when sampling the point to start the algorithm at.
- `try_points=()`: Additional points to try when determining which point to start at.
- `k=select_initial_point(get_points(tri), q; m, point_indices, try_points)`: The index of the point to start the algorithm at. 
- `store_history::F=Val(false)`: Whether to record the history of the algorithm. See also [`PointLocationHistory`](@ref).
- `history=nothing`: The object to store the history into, if `is_true(store_history)`.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.
- `exterior_curve_index=1`: The curve (or curves) corresponding to the outermost boundary.
- `maxiters = num_triangles(tri)`: Maximum number of iterations to perform before restarting the algorithm at a new initial point. 
- `concavity_protection=false`: When your triangulation has concave boundaries, it is possible that a ghost triangle is incorrectly classified as containing the point. By setting this to `true`, this will be protected against.

!!! note 

    You shouldn't ever need `maxiters` if your triangulation is convex everywhere, as Delaunay triangulations 
    have no problems with jump-and-march, as the sequence of triangles visited is acyclic (H. Edelsbrunner, An acyclicity theorem for cell complexes in d dimensions, Combinatorica 10 (1990) 251–
    260.) However, if the triangulation is not convex, e.g. if you have a constrained triangulation with boundaries 
    and excavations, then an infinite loop can be found where we just keep walking in circles. In this case, 
    you can use the `maxiters` keyword argument to specify the maximum number of iterations to perform before
    reinitialising the algorithm at a random vertex. When reinitialising, the value of `m` is doubled each time.

# Outputs 
Returns `V`, the triangle in `tri` containing `q`.

!!! warning 

    If your triangulation does not have ghost triangles, and the point `q` is outside of the triangulation, 
    this function may fail to terminate. You may like to add ghost triangles in this case (using 
    [`add_ghost_triangles!`](@ref)), noting that there is no actual triangle that `q` is inside of
    when it is outside of the triangulation unless ghost triangles are present. 
"""
function jump_and_march(tri::Triangulation{P,Ts,I}, q;
    point_indices=each_solid_vertex(tri),
    m=default_num_samples(num_vertices(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng(),
    k=select_initial_point(tri, q; point_indices, m, try_points, rng),
    store_history::F=Val(false),
    history=nothing,
    exterior_curve_index=1,
    maxiters=2 + length(exterior_curve_index) - num_solid_vertices(tri) + num_solid_edges(tri),
    concavity_protection=false) where {P,Ts,I,F}
    return _jump_and_march(
        tri,
        q,
        k,
        store_history,
        history,
        rng,
        exterior_curve_index,
        maxiters,
        zero(maxiters),
        concavity_protection,
        zero(maxiters)
    )
end

function _jump_and_march(
    tri::Triangulation{P,Ts,I,E},
    q,
    k,
    store_history::F=Val(false),
    history=nothing,
    rng::AbstractRNG=Random.default_rng(),
    exterior_curve_index=1,
    maxiters=2 + length(exterior_curve_index) - num_vertices(graph) + num_edges(graph),
    cur_iter=0,
    concavity_protection=false,
    num_restarts=0) where {P,Ts,I,E,F}
    is_bnd, bnd_idx = is_boundary_node(tri, k)
    if !(is_bnd && get_curve_index(tri, bnd_idx) ∈ exterior_curve_index) || !has_ghost_triangles(tri)
        # If k is not a boundary node, then we can rotate around the point k to find an initial triangle 
        # to start the search. Additionally, if there are no ghost triangles, we can only hope that q 
        # is inside the interior, meaning we should only search for the initial triangle here anyway.
        p, i, j, pᵢ, pⱼ = select_initial_triangle_interior_node(tri, k, q, store_history, history, rng)
        if !edge_exists(i) && !edge_exists(j) # When we find a possible infinite loop, we use i==j==DefaultAdjacentValue. Let's reinitialise. 
            return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
        end
        if is_true(store_history)
            add_triangle!(history, j, i, k)
        end
    else
        # We have an outer boundary node. First, let us check the neighbouring boundary edges. 
        direction, q_pos, next_vertex, right_cert, left_cert =
            check_for_intersections_with_adjacent_boundary_edges(tri, k, q, bnd_idx)
        if !is_outside(direction)
            # q is collinear with one of the edges, so let's jump down these edges and try to find q
            q_pos, u, v, w = search_down_adjacent_boundary_edges(tri, k, q, direction, q_pos, next_vertex, store_history, history, bnd_idx)
            if is_on(q_pos)
                return construct_triangle(triangle_type(Ts), u, v, w)
            else
                u, v = exterior_jump_and_march(tri, u, q, bnd_idx)
                V = construct_triangle(triangle_type(Ts), u, v, get_adjacent(tri, u, v)) # Can't just use I(BoundaryIndex) here since there could be multiple - just use get_adjacent
                if _concavity_protection_check(concavity_protection, tri, V, q)
                    return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
                else
                    return V
                end
            end
        end
        # If we did not find anything from the neighbouring boundary edges, we can search the neighbouring interior edges
        i, j, edge_cert, triangle_cert =
            check_for_intersections_with_interior_edges_adjacent_to_boundary_node(tri, k, q, right_cert, left_cert, store_history, history, bnd_idx)
        if is_inside(triangle_cert)
            return construct_triangle(triangle_type(Ts), i, j, k)
        elseif is_none(edge_cert)
            u, v = exterior_jump_and_march(tri, k, q, bnd_idx)
            V = construct_triangle(triangle_type(Ts), u, v, get_adjacent(tri, u, v))
            if _concavity_protection_check(concavity_protection, tri, V, q)
                return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
            else
                return V
            end
        else
            p, pᵢ, pⱼ = get_point(tri, k, i, j)
        end
    end
    if q == p || q == pᵢ || q == pⱼ
        # Just return where we currently are. We do need to be careful, though: 
        # If k was a boundary index, then one of pᵢ or pⱼ could come from the 
        # representative point list, which could mean that q is equal to one of the 
        # vertices, but without meaning that it is actually in that triangle. So, 
        # we need to also check for the type of indices we have. 
        safety_check = (q == p && !is_boundary_index(k)) ||
                       (q == pᵢ && !is_boundary_index(i)) ||
                       (q == pⱼ && !is_boundary_index(j))
        if safety_check
            orientation = triangle_orientation(pᵢ, pⱼ, p)
            if is_positively_oriented(orientation)
                return construct_triangle(triangle_type(Ts), i, j, k)
            else
                return construct_triangle(triangle_type(Ts), j, i, k)
            end
        else
            return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
        end
    end
    original_k = k
    ## Now let us do the straight line search 
    # The idea is to keep jumping until pᵢpⱼ goes past q, meaning pᵢpⱼq is no longer a positively oriented triangle
    arrangement = triangle_orientation(pᵢ, pⱼ, q)
    local last_changed # Need this for deciding which variable to use when we hit a collinear point 
    last_changed = I(DefaultAdjacentValue) # just an initial value
    if is_true(store_history)
        add_left_vertex!(history, i)
        add_right_vertex!(history, j)
    end
    while is_positively_oriented(arrangement)
        cur_iter += 1
        if is_true(store_history)
            if last_changed == i
                add_left_vertex!(history, i)
            elseif last_changed == j
                add_right_vertex!(history, j)
            end
        end
        # We need to step forward. To do this, we need to be careful of boundary indices. 
        # Since a boundary curve may be represented by different boundary indices 
        # at different locations, we need to check for this. The first check 
        # is for an outer boundary index.
        if is_boundary_index(k) && get_curve_index(tri, k) ∈ exterior_curve_index
            # If this happens, it means we have hit an outer boundary edge, and so we need to go into the exterior. If there are no 
            # ghost triangles, though, we just restart the search. Note that interior boundary indices do not matter since the ghost 
            # triangles there have the same orientation, so we can find them as normal.
            if has_ghost_triangles(tri)
                i, j = exterior_jump_and_march(tri, last_changed == i ? j : i, q, k) # use last_changed to ensure we get the boundary point
                V = construct_triangle(triangle_type(Ts), i, j, get_adjacent(tri, i, j))
                if _concavity_protection_check(concavity_protection, tri, V, q)
                    return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
                else
                    return V
                end
            else
                return _jump_and_march(tri, q, k, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
            end
        end
        # Now we can finally move forward. We use check_existence to protect against the issue mentioned above.
        k = get_adjacent(tri, i, j)
        pₖ = get_point(tri, k)
        pₖ_pos = point_position_relative_to_line(p, q, pₖ)
        if is_right(pₖ_pos)
            if is_true(store_history)
                add_triangle!(history, i, j, k)
            end
            j, pⱼ = k, pₖ
            last_changed = j
        elseif is_left(pₖ_pos)
            if is_true(store_history)
                add_triangle!(history, i, j, k)
            end
            i, pᵢ = k, pₖ
            last_changed = i
        else
            # pₖ is collinear with pq. We will first check if q is already in the current triangle
            in_cert = point_position_relative_to_triangle(pᵢ, pⱼ, pₖ, q) # ... Maybe there is a better way to compute this, reusing previous certificates? Not sure. ...
            if is_true(store_history)
                # We need to be careful about whether or not we want to add another 
                # collinear segment in this case, since e.g. q might just be on the triangle 
                # on a separate edge (if !is_outside(in_cert)),
                # but not collinear with pq. We just test it directly, 
                # but I'm sure there's a better way. (This is why we have the 
                # original_k variable.) I'm also sure there's a way to know the 
                # actual edge, rather than trying both as we do below. Usually it's 
                # the edge last_changed == i ? j : i, but not always... We split this 
                # into two cases, one where last_changed ≠ 0 and otherwise. This is because 
                # last_changed ≠ 0 path is a lot more common (> 99%).
                if is_collinear(point_position_relative_to_line(tri, original_k, q, last_changed == i ? j : i))
                    add_edge!(history, last_changed == i ? j : i, k)
                elseif edge_exists(last_changed) && is_collinear(point_position_relative_to_line(tri, original_k, q, last_changed))
                    add_edge!(history, last_changed, k)
                elseif k ≠ original_k && pₖ ≠ q && is_collinear(point_position_relative_to_line(tri, k, original_k, q))
                    # This case is needed in case we have a collinearity, but only because we pass through a point 
                    # rather than because we pass through a collinear segment. We will sort through these later using 
                    # fix_segments!.
                    add_edge!(history, original_k, k)
                    add_index!(history, num_edges(history))
                else
                    # This case here is a lot less likely. 
                    if pₖ == q
                        prev = get_adjacent(tri, j, i)
                        if is_collinear(point_position_relative_to_line(tri, i, p, q))
                            add_edge!(history, i, k)
                        elseif is_collinear(point_position_relative_to_line(tri, j, p, q))
                            add_edge!(history, j, k)
                        elseif prev ≠ original_k && is_collinear(point_position_relative_to_line(tri, prev, p, q)) # prev ≠ original_k because otherwise we'll say the line segment we started with is collinear with itself, which is true but we don't need that case.
                            add_edge!(history, prev, k)
                        end
                    end
                end
                add_triangle!(history, i, j, k)
            end
            if !is_outside(in_cert)
                V = construct_triangle(triangle_type(Ts), i, j, k)
                if _concavity_protection_check(concavity_protection, tri, V, q)
                    return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
                else
                    return V
                end
            end
            # To decide which direction this collinear point is away from the line, we can just use the last changed variable:
            # If last_changed = i, this means that the left direction was what caused the collinear point, so make k go on the left. 
            # If not, make it go to the right.
            if last_changed == i
                i, pᵢ = k, pₖ
                last_changed = i
            else
                j, pⱼ = k, pₖ
                last_changed = j
            end
        end
        arrangement = triangle_orientation(pᵢ, pⱼ, q)
        if cur_iter ≥ maxiters
            return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
        end
    end
    # We can finish the above loop even if q is not in the triangle, in which case pᵢpⱼq was a straight line. 
    # To clear this up, let us just restart. 
    if is_degenerate(arrangement)
        pₖ = get_point(tri, k) # Need to have this here in case we skip the entire loop, above, meaning pₖ won't exist
        in_cert = point_position_relative_to_triangle(pᵢ, pⱼ, pₖ, q) # ... Maybe there is a better way to compute this, reusing previous certificates? Not sure. ...
        if is_true(store_history)
            k′ = get_adjacent(tri, i, j)
            add_triangle!(history, i, j, k′)
            add_left_vertex!(history, i)
            add_right_vertex!(history, j)
        end
        if is_outside(in_cert)
            return _jump_and_march(tri, q, last_changed == I(DefaultAdjacentValue) ? i : last_changed, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
        end
    end
    # Swap the orientation to get a positively oriented triangle, remembering that we kept pᵢ on the left of pq and pⱼ on the right 
    k = get_adjacent(tri, j, i)
    V = construct_triangle(triangle_type(Ts), j, i, k)
    if _concavity_protection_check(concavity_protection, tri, V, q)
        return restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts + 1)
    end
    return V
end

function _concavity_protection_check(concavity_protection, tri, V, q)
    !concavity_protection && return false
    points = get_points(tri)
    bn = get_boundary_nodes(tri)
    δ = distance_to_polygon(q, points, bn)
    cert = point_position_relative_to_triangle(tri, V, q)
    if is_outside(cert)
        return true # need to restart
    end
    is_ghost = is_ghost_triangle(V)
    need_to_restart = (is_ghost && δ > 0.0) || (!is_ghost && δ < 0.0)
    return need_to_restart
end

const RESTART_LIMIT = 25
function restart_jump_and_march(tri, q, store_history, history, rng, exterior_curve_index, maxiters, cur_iter, concavity_protection, num_restarts)
    if num_restarts < RESTART_LIMIT
        m = num_solid_vertices(tri)
        point_indices = each_solid_vertex(tri)
        k = select_initial_point(tri, q; m=(m ÷ 2) + 1, point_indices, rng) # don't want to try all points, still want to give the algorithm a chance
        return _jump_and_march(tri, q, k, store_history, history, rng, exterior_curve_index, maxiters, zero(cur_iter), concavity_protection, num_restarts)
    else
        V₁ = brute_force_search(tri, q)
        V₁_is_bad = _concavity_protection_check(concavity_protection, tri, V₁, q)
        if V₁_is_bad
            if is_ghost_triangle(V₁)
                V₁ = brute_force_search(tri, q; itr=each_solid_triangle(tri))
            else
                V₁ = brute_force_search(tri, q; itr=each_ghost_triangle(tri))
            end
        end
        return V₁
    end
end
