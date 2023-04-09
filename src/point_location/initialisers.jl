"""
    default_num_samples(num_points)

Returns `ceil(cbrt(num_points))`. This is the default number of samples 
to use when sampling in the jump-and-march algorithm.
"""
default_num_samples(num_points::I) where {I} = ceil(I, cbrt(num_points))

function compare_distance(current_dist, current_idx, pts, i, qx, qy)
    p = get_point(pts, i)
    px, py = getxy(p)
    sq_dist = (px - qx)^2 + (py - qy)^2
    if sq_dist < current_dist
        current_dist = sq_dist
        current_idx = i
    end
    return current_dist, current_idx
end

"""
    select_initial_point(pts, q;
        point_indices=each_point_index(pts),
        m=default_num_samples(length(point_indices)),
        try_points=(),
        rng::AbstractRNG=Random.default_rng())

Given a collection of points and a point `q`, select `m` random points from `pts` and return 
the one that is closet to `q`.

# Arguments 
- `pts`: The set of points.
- `q`: The coordinates of the query point, or its index in `pts`.

# Keyword Arguments
- `point_indices=each_point_index(pts)`: The indices of the points to consider.
- `m=default_num_samples(length(point_indices))`: The number of points to sample.
- `try_points=()`: The indices of points to try in addition to the sampled points.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.

# Outputs 
- `k`: The index of the point in `pts` that is closest to `q` out of the points sampled.
"""
function select_initial_point(pts, q;
    point_indices=each_point_index(pts),
    m=default_num_samples(length(point_indices)),
    try_points=(),
    rng::AbstractRNG=Random.default_rng())
    F = number_type(pts)
    current_dist = typemax(F)
    current_idx = first(point_indices) # Just some index, doesn't matter 
    qx, qy = getxy(q)
    for _ in 1:m  # Not using replacement, but probability of duplicates is approximately 0.5length(point_indices)^(-1/3)
        i = rand(rng, point_indices)
        is_boundary_index(i) && continue
        current_dist, current_idx = compare_distance(current_dist, current_idx, pts, i, qx,
            qy)
    end
    for i in try_points
        current_dist, current_idx = compare_distance(current_dist, current_idx, pts, i, qx,
            qy)
    end
    return current_idx
end
function select_initial_point(pts, q::Integer;
    m=default_num_samples(num_points(pts)),
    point_indices=each_point_index(pts),
    try_points=(),
    rng::AbstractRNG=Random.default_rng())
    return select_initial_point(pts, get_point(pts, q); m, point_indices, try_points, rng)
end

function select_random_edge(pts, representative_point_list, boundary_map, edges,
    rng::AbstractRNG=Random.default_rng())
    edge = random_edge(edges, rng)
    i, j = edge_indices(edge)
    pᵢ, pⱼ = get_point(pts, representative_point_list, boundary_map, i, j)
    return i, j, pᵢ, pⱼ
end

function prepare_initial_edge(pts, representative_point_list, boundary_map, edges, p, q,
    rng::AbstractRNG=Random.default_rng())
    i, j, pᵢ, pⱼ = select_random_edge(pts, representative_point_list, boundary_map, edges, rng)
    line_cert_i = point_position_relative_to_line(p, q, pᵢ)
    line_cert_j = point_position_relative_to_line(p, q, pⱼ)
    return i, j, pᵢ, pⱼ, line_cert_i, line_cert_j
end

function select_initial_triangle_clockwise(pts, adj::Adjacent{I,E}, representative_point_list, boundary_map, p, q, pᵢ, pⱼ, i, j, k,
    boundary_index_ranges,
    check_existence::V=Val(has_multiple_segments(boundary_map)),
    store_history::F=Val(false),
    history=nothing) where {V,F,I,E}
    line_cert_i = point_position_relative_to_line(p, q, pᵢ)
    if is_true(store_history) && is_collinear(line_cert_i)
        add_edge!(history, k, i)
    end
    while is_left(line_cert_i)
        j = i
        pⱼ = pᵢ
        i = get_adjacent(adj, i, k; check_existence=check_existence,
            boundary_index_ranges)
        pᵢ = get_point(pts, representative_point_list, boundary_map, i)
        line_cert_i = point_position_relative_to_line(p, q, pᵢ)
        if is_true(store_history) && is_collinear(line_cert_i)
            add_edge!(history, k, i)
        end
    end
    return i, j, pᵢ, pⱼ
end

function select_initial_triangle_counterclockwise(pts, adj::Adjacent{I,E}, representative_point_list, boundary_map, line_cert_j, p, q,
    pᵢ, pⱼ, i, j, k, boundary_index_ranges,
    check_existence::V=Val(has_multiple_segments(boundary_map)),
    store_history::F=Val(false),
    history=nothing) where {V,F,I,E}
    if is_true(store_history) && is_collinear(line_cert_j)
        add_edge!(history, k, j)
    end
    while is_right(line_cert_j)
        i = j
        pᵢ = pⱼ
        j = get_adjacent(adj, k, j; check_existence=check_existence,
            boundary_index_ranges)
        pⱼ = get_point(pts, representative_point_list, boundary_map, j)
        line_cert_j = point_position_relative_to_line(p, q, pⱼ)
        if is_true(store_history) && is_collinear(line_cert_j)
            add_edge!(history, k, j)
        end
    end
    return i, j, pᵢ, pⱼ
end

"""
    select_initial_triangle_interior_node(
        pts, 
        adj, 
        adj2v, 
        representative_point_list, 
        boundary_map, 
        k, 
        q, 
        boundary_index_ranges, 
        check_existence=Val(has_multiple_segments(boundary_map)),
        store_history=Val(false),
        history=nothing,
        rng::AbstractRNG = Random.default_rng())

Selects an initial triangle for the jump-and-march algorithm, starting from a point with index `k` 
and finding the triangle such that the line from the `k`th point to `q` intersects it. It is assumed 
that `k` is a point that is not on the boundary. 

# Arguments 
- `pts`: The set of points. 
- `adj`: The [`Adjacent`](@ref) map.
- `adj2v`: The [`Adjacent2Vertex`](@ref) map.
- `representative_point_list`: The list of representative points, giving the coordinates of points corresponding to boundary indices.
- `boundary_map`: The map taking boundary indices to their corresponding segment in the boundary nodes. See [`construct_boundary_map`](@ref).
- `k`: The index of the point in `pts` that we are starting at.
- `q`: The point being searched for.
- `boundary_index_ranges`: The `Dict` handling the mapping of boundary indices to the range of boundary indices belonging to the same curve. See [`construct_boundary_index_ranges`](@ref).
- `check_existence=Val(has_multiple_segments(boundary_map))`: Checks for different possible boundary indices when there are multiple segments. See [`get_adjacent`](@ref).
- `store_history=Val(false)`: Whether to store visited triangles. Exterior ghost triangles will not be stored.
- `history=nothing`: The history. This should be a [`PointLocationHistory`](@ref) type if `store_history` is `true`.
- `rng::AbstractRNG = Random.default_rng()`: The random number generator.

# Outputs 
- `p`: The `k`th point in `pts`.
- `i, j`: These are indices defining the edge of a triangle including the point `p`, such that `i` is to the left of the line `pq` and `j` is to the right of `pq`.
- `pᵢ, pⱼ`: The points in `pts` corresponding to the indices in `i` and `j`, respectively.
"""
function select_initial_triangle_interior_node(
    pts,
    adj::Adjacent{I,E},
    adj2v::Adjacent2Vertex,
    representative_point_list,
    boundary_map,
    k,
    q,
    boundary_index_ranges,
    check_existence::V=Val(has_multiple_segments(boundary_map)),
    store_history::F=Val(false),
    history=nothing,
    rng::AbstractRNG=Random.default_rng()) where {V,F,I,E}
    p = get_point(pts, representative_point_list, boundary_map, k)
    ## Select the initial edge to rotate about
    neighbouring_edges = get_adjacent2vertex(adj2v, k)
    i, j, pᵢ, pⱼ, line_cert_i, line_cert_j = prepare_initial_edge(pts, representative_point_list, boundary_map,
        neighbouring_edges, p, q,
        rng)
    p == q && return p, j, i, pⱼ, pᵢ
    ## Take care for collinear edges
    while is_collinear(line_cert_j) || is_collinear(line_cert_i)
        if is_collinear(line_cert_j)
            opposing_idx = j
            on_edge_cert = point_position_on_line_segment(p, pⱼ, q)
        else
            opposing_idx = i
            on_edge_cert = point_position_on_line_segment(p, pᵢ, q)
        end
        # We want q to be in the direction of of ppᵢ or ppⱼ, which just means not left
        if !is_left(on_edge_cert)
            if is_on(on_edge_cert) || is_degenerate(on_edge_cert)
                is_true(store_history) && add_edge!(history, k, opposing_idx)
                return p, j, i, pⱼ, pᵢ
            else
                pos_cert = point_position_relative_to_line(p, q, pᵢ)
                if is_left(pos_cert)
                    is_true(store_history) && add_edge!(history, k, opposing_idx)
                    return p, i, j, pᵢ, pⱼ
                else
                    is_true(store_history) && add_edge!(history, k, opposing_idx)
                    return p, j, i, pⱼ, pᵢ
                end
            end
        else
            # If we haven't yet found the edge, let's just pick another one
            i, j, pᵢ, pⱼ, line_cert_i, line_cert_j = prepare_initial_edge(pts, representative_point_list, boundary_map,
                neighbouring_edges,
                p, q, rng)
        end
    end

    ## Now rotate around to find a triangle that the line pq intersects through 
    if is_left(line_cert_j)
        i, j, pᵢ, pⱼ = select_initial_triangle_clockwise(pts, adj, representative_point_list, boundary_map, p, q, pᵢ,
            pⱼ, i, j, k, boundary_index_ranges,
            check_existence,
            store_history,
            history)
    else
        i, j, pᵢ, pⱼ = select_initial_triangle_counterclockwise(pts, adj, representative_point_list, boundary_map,
            line_cert_j, p, q, pᵢ, pⱼ,
            i, j, k,
            boundary_index_ranges,
            check_existence,
            store_history,
            history)
    end

    ## Swap values and return
    i, j = j, i
    pᵢ, pⱼ = pⱼ, pᵢ # pᵢ is left of pq, pⱼ is right of pq 
    return p, i, j, pᵢ, pⱼ
end

"""
    check_for_intersections_with_adjacent_boundary_edges(
        pts, 
        adj, 
        boundary_index_ranges, 
        boundary_map, 
        k, 
        q, 
        check_existence = Val(has_multiple_segments(boundary_map)))

Assuming that `k` is on the outer boundary, this function searches down the boundary edges adjacent to `k` to try and locate a triangle 
or edge containing `q`.

See also [`search_down_adjacent_boundary_edges`](@ref), which uses this function to determine an initial direction to search along a 
straight boundary in case `q` is collinear with it.

# Arguments 
- `pts`: The set of points. 
- `adj`: The [`Adjacent`](@ref) map.
- `boundary_index_ranges`: The list of boundary index ranges from [`construct_boundary_index_ranges`](@ref), giving a range of boundary indices for a given boundary index all belonging to the same curve. 
- `boundary_map`: The boundary map from [`construct_boundary_map`](@ref), which maps a boundary index to its corresponding segment. 
- `k`: The outer boundary point `k` to start searching from.
- `q`: The query point `q`.
- `check_existence`: The same keyword argument from [`get_adjacent`](@ref), and is needed when you use multiple segments in your boundary.

# Outputs 
The output is a 5-tuple, with the first three elements having several possible forms:

- `(Certificate.Outside, Certificate.Outside, k)`: The point `q` is not collinear with either of the adjacent boundary edges.
- `(Certificate.Right, C, r)`, where `C` is either `Certificate.On` or `Certificate.Right` and `r` is the vertex right of `k`: The point `q` is collinear with the edge to the right of `k`. If `C` is `Certificate.On`, then `q` is on the edge, whereas `C` being `Certificate.Right` means it is right of the edge. 
- `(Certificate.Left, C, ℓ)`, where `C` is either `Certificate.On` or `Certificate.Left` and `ℓ` is the vertex left of `k`: The point `q` is collinear with the edge to the left of `k`. If `C` is `Certificate.On`, then `q` is on the edge, whereas `C` being `Certificate.Left` means it is left of the edge.

In the latter two outputs above, `C` could also be `Certificate.Degenerate`, which means that `q` is `get_point(pts, r)` or `get_point(pts, ℓ)`, respectively.

The latter two elements of the tuple are: 

- `right_cert`: The position of `q` relative to the edge `(p, p_right)`, where `p_right` is the point on the boundary to the right of `p`.
- `left_cert`: The position of `q` relative to the edge `(p, p_left)`, where `p_left` is the point on the boundary to the left of `p`.

These returned values are useful in case we need to go to [`check_for_intersections_with_interior_edges_adjacent_to_boundary_node`](@ref), since we can reuse these certificates.
"""
function check_for_intersections_with_adjacent_boundary_edges(
    pts,
    adj::Adjacent{I,E},
    boundary_index_ranges,
    representative_point_list,
    boundary_map,
    k,
    q,
    check_existence::V=Val(has_multiple_segments(boundary_map))) where {I,E,V}
    p = get_point(pts, representative_point_list, boundary_map, k)
    right = get_right_boundary_node(adj, k, I(BoundaryIndex), boundary_index_ranges,
        check_existence)
    left = get_left_boundary_node(adj, k, I(BoundaryIndex), boundary_index_ranges,
        check_existence)
    pright, pleft = get_point(pts, representative_point_list, boundary_map, right, left)
    right_cert = point_position_relative_to_line(p, pright, q)
    left_cert = point_position_relative_to_line(pleft, p, q)
    flipped_left_cert = point_position_relative_to_line(p, pleft, q) # Useful for the returned value, but we use left_cert below to get a good value for direction_cert
    !is_collinear(right_cert) && !is_collinear(left_cert) &&
        return (Cert.Outside, Cert.Outside, k, right_cert, flipped_left_cert)
    if is_collinear(right_cert)
        direction_cert = point_position_on_line_segment(p, pright, q)
        if !is_left(direction_cert)
            return (Cert.Right, direction_cert, right, right_cert, flipped_left_cert)
        elseif is_left(direction_cert) && !is_collinear(left_cert)
            return (Cert.Outside, Cert.Outside, k, right_cert, flipped_left_cert)
        end
    end
    direction_cert = point_position_on_line_segment(pleft, p, q)
    if !is_right(direction_cert)
        return (Cert.Left, direction_cert, left, right_cert, flipped_left_cert)
    else
        return (Cert.Outside, Cert.Outside, k, right_cert, flipped_left_cert)
    end
end

"""
    search_down_adjacent_boundary_edges(
        pts, 
        adj,
        boundary_index_ranges, 
        representative_point_list, 
        boundary_map, 
        k, 
        q,
        direction, 
        q_pos, 
        next_vertex,
        check_existence=Val(has_multiple_segments(boundary_map)),
        store_history=Val(false),
        history=nothing) 

Starting at the outer boundary node `k`, walks down the boundary in the direction of `q` until finding `q` or finding that it is outside of the triangulation. 

# Arguments 
- `pts`: The collection of points.
- `adj`: The [`Adjacent`](@ref) map.
- `boundary_index_ranges`: The list of boundary index ranges from [`construct_boundary_index_ranges`](@ref).
- `representative_point_list`: The list of representative points. 
- `boundary_map`: The boundary map from [`construct_boundary_map`](@ref).
- `k`: The outer boundary index.
- `q`: The point being searched for.
- `direction`: The direction of `q` from `get_point(pts, k)`.
- `q_pos`: The certificate for the position of `q` from this point `k`.
- `next_vertex`: The next vertex in the direction of `q` (this argument comes from [`check_for_intersections_with_adjacent_boundary_edges`](@ref)).
- `check_existence=Val(has_multiple_segments(boundary_map))`: The same keyword argument from [`get_adjacent`](@ref), and is needed when you use multiple segments in your boundary.
- `store_history=Val(false)`: Whether to store the history of the search.
- `history=nothing`: The history of the search.

# Outputs 
The output takes the form `(cert, u, v, w)`, where: 

- `cert`: This is `Certificate.On` if `q` is on the edge `(u, v)`, and `Certificate.Outside` if `q` is outside of the triangulation.
- `(u, v, w)`: If `is_on(cert)`, then this is a positively oriented triangle with `q` on the edge `(u, v)`. Otherwise, `(u, v, w)` is a ghost triangle close to `q`.

!!! warning 

    This function relies on the assumption that the geometry is convex.
"""
function search_down_adjacent_boundary_edges(
    pts,
    adj::Adjacent{I,E},
    boundary_index_ranges,
    representative_point_list,
    boundary_map,
    k,
    q,
    direction,
    q_pos,
    next_vertex,
    check_existence::V=Val(has_multiple_segments(boundary_map)),
    store_history::F=Val(false),
    history=nothing) where {I,E,V,F}
    i = k
    j = next_vertex
    pⱼ = get_point(pts, representative_point_list, boundary_map, j)
    if is_right(direction)
        if is_true(store_history) # in case we don't enter the loop
            k′ = get_adjacent(adj, i, j; check_existence, boundary_index_ranges)
            add_triangle!(history, i, j, k′)
            add_edge!(history, i, j)
        end
        while is_right(q_pos)
            i, pᵢ = j, pⱼ
            j = get_right_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges, check_existence)
            if is_true(store_history)
                k′ = get_adjacent(adj, i, j; check_existence, boundary_index_ranges)
                add_triangle!(history, i, j, k′)
                add_edge!(history, i, j)
            end
            pⱼ = get_point(pts, representative_point_list, boundary_map, j)
            right_cert = point_position_relative_to_line(pᵢ, pⱼ, q)
            q_pos = !is_collinear(right_cert) ? Cert.Outside : point_position_on_line_segment(pᵢ, pⱼ, q)
        end
        if is_outside(q_pos)
            return (Cert.Outside, j, i, get_adjacent(adj, j, i; check_existence, boundary_index_ranges))
        elseif is_on(q_pos) || is_degenerate(q_pos)
            k = get_adjacent(adj, i, j; check_existence, boundary_index_ranges)
            return (Cert.On, i, j, k)
        end
    else
        if is_true(store_history) # in case we don't enter the loop
            k′ = get_adjacent(adj, j, i; check_existence, boundary_index_ranges)
            add_triangle!(history, j, i, k′)
            add_edge!(history, i, j) # i,j → j, i to get the ordering of the segments
        end
        while is_left(q_pos)
            i, pᵢ = j, pⱼ
            j = get_left_boundary_node(adj, i, I(BoundaryIndex), boundary_index_ranges, check_existence)
            if is_true(store_history)
                k′ = get_adjacent(adj, j, i; check_existence, boundary_index_ranges)
                add_triangle!(history, j, i, k′)
                add_edge!(history, i, j)
            end
            pⱼ = get_point(pts, representative_point_list, boundary_map, j)
            left_cert = point_position_relative_to_line(pᵢ, pⱼ, q)
            q_pos = !is_collinear(left_cert) ? Cert.Outside : point_position_on_line_segment(pⱼ, pᵢ, q)
        end
        if is_outside(q_pos)
            return (Cert.Outside, i, j, get_adjacent(adj, i, j; check_existence, boundary_index_ranges))
        elseif is_on(q_pos) || is_degenerate(q_pos)
            k = get_adjacent(adj, j, i; check_existence, boundary_index_ranges)
            return (Cert.On, j, i, k)
        end
    end
    throw("Failed to identify the location of $q.")
end

"""
    check_for_intersections_with_interior_edges_adjacent_to_boundary_node(
        pts,
        adj, 
        graph, 
        boundary_index_ranges, 
        representative_point_list,
        boundary_map, 
        k, 
        q, 
        right_cert, 
        left_cert, 
        check_existence=Val(has_multiple_segments(boundary_map)), 
        store_history=Val(false), 
        history=nothing)

Checks if the line connecting the `k`th point of `pts` to `q` intersects any of the edges neighbouring the boundary node `k`.

# Arguments 
- `pts`: The collection of points. 
- `adj`: The [`Adjacent`](@ref) map.
- `graph`: The [`Graph`](@ref).
- `boundary_index_ranges`: The boundary index range mapping from [`construct_boundary_index_ranges`](@ref).
- `representative_point_list`: The list of representative points.
- `boundary_map`: The map that handles the mapping of boundary indices to boundary segments. Sse [`construct_boundary_map`](@ref).
- `k`: The boundary node.
- `q`: The point we are searching for. 
- `right_cert`: A certificate giving the position of `q` to the right of the `k`th point. This comes from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `left_cert`: A certificate giving the position of `q` to the left of the `k`th point. This comes from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `check_existence=Val(has_multiple_segments(boundary_nodes)))`: Checks for different possible boundary indices when there are multiple segments. See [`get_adjacent`](@ref).
- `store_history=Val(false)`: Whether to store visited triangles. Exterior ghost triangles will not be stored.
- `history=nothing`: The history. This should be a [`PointLocationHistory`](@ref) type if `store_history` is `true`.

# Outputs 
There are several possible forms for the returned values. These are listed below, letting `p` be the `k`th point, `pᵢ` the point corresponding to 
the index `i`, and `pⱼ` the point corresponding to the index `j`:

- `(i, j, Certificate.Single, Certificate.Outside)`

The line `pq` intersects the edge `pᵢpⱼ`, and `(j, i, k)` is a positively oriented triangle. In particular, pᵢ is left of `pq` and `pⱼ` is right of `pq`.
- `(i, j, Certificate.None, Certificate.Inside)`

The point `q` is inside the positively oriented triangle `(i, j, k)`.
- `(zero(I), zero(I), Cert.None, Cert.Outside)`

The point `q` is outside of the triangulation. Note that `I` is the integer type.
- `(i, j, Cert.On, Cert.Inside)`

The point `q` is on the edge `pᵢpⱼ`, and so is inside the positively oriented triangle `(i, j, k)`.
- `(i, j, Cert.Right, Cert.Outside)`

The point `q` is collinear with the edge `pᵢpⱼ`, but is off of it and further into the triangulation. 

!!! warning 

    This function should only be used after [`check_for_intersections_with_adjacent_boundary_edges`](@ref), and currently is only guaranteed to work on convex geometries. 
"""
function check_for_intersections_with_interior_edges_adjacent_to_boundary_node(
    pts,
    adj::Adjacent{I,E},
    graph::Graph{I},
    boundary_index_ranges,
    representative_point_list,
    boundary_map,
    k, 
    q,
    right_cert,
    left_cert,
    check_existence::V=Val(has_multiple_segments(boundary_map)),
    store_history::F=Val(false),
    history=nothing) where {I,E,V,F}
    p = get_point(pts, representative_point_list, boundary_map, k)
    other_boundary_node = get_left_boundary_node(adj, k, I(BoundaryIndex),
        boundary_index_ranges, check_existence)
    num_interior_neighbours = num_neighbours(graph, k) - 3 # - 3 = - the two boundary neighbours - BoundaryIndex
    i = get_right_boundary_node(adj, k, I(BoundaryIndex), boundary_index_ranges,
        check_existence)
    pᵢ = get_point(pts, representative_point_list, boundary_map, i)
    certᵢ = right_cert # Do not need to check if is_collinear(certᵢ) - this function should only be used after checking check_for_intersections_with_adjacent_boundary_edges anyway
    for _ in 1:max(num_interior_neighbours, 1)
        # In the above, we need max(num_interior_neighbours, 1) in case the point 
        # being considered is a corner point on the boundary with no interior edge neighbours. 
        # If we just had 1:num_interior_neighbours, then nothing would happen in this loop 
        # since num_interior_neighbours = 0, but we still need to make sure we check the edge 
        # connecting the two boundary neighbours
        j = get_adjacent(adj, k, i; check_existence=check_existence,
            boundary_index_ranges)
        pⱼ = get_point(pts, representative_point_list, boundary_map, j)
        certⱼ = point_position_relative_to_line(p, pⱼ, q)
        # Start by checking if pq intersects the edge pᵢpⱼ, meaning q is left or right of pᵢ, but the opposite of pⱼ
        if (is_left(certᵢ) && is_right(certⱼ)) || (is_right(certᵢ) && is_left(certⱼ))
            # If we have this intersection, there are three possibilities. The first is that there is indeed 
            # an intersection, which we check by counting the intersections of pq with pᵢpⱼ. 
            # Note that we have to check this in this case since q could be behind the edge - this is not 
            # checked for interior nodes since we can completely rotate around an interior node.
            intersection_cert = line_segment_intersection_type(p, q, pᵢ, pⱼ)
            if has_one_intersection(intersection_cert)
                # q is not in the triangle, but we still have an intersection. 
                # In the returned value, we switch i and j so that pᵢ is left of pq and pⱼ is right of pq
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return j, i, Cert.Single, Cert.Outside
            elseif is_touching(intersection_cert)
                # q is on the edge 
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                    add_edge!(history, i, j)
                end
                return i, j, Cert.On, Cert.Inside
            end
            # There may be no intersection, but it could be inside the triangle. 
            pⱼp_edge = point_position_relative_to_line(pⱼ, p, q)
            ppᵢ_edge = point_position_relative_to_line(p, pᵢ, q)
            if is_left(pⱼp_edge) && is_left(ppᵢ_edge)
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return i, j, Cert.None, Cert.Inside
            end
            # If none of the above lead to a return, then q must be on the other side of p away from the interior, meaning 
            # it is in the exterior. Here, we are making use of the fact that the domain is 
            # convex. When the domain is not convex, this won't necessarily work.
            # TODO: See how we can change these rules to allow for non-convex geometries, providing a better 
            # method for point location in general constrained triangulations.
            return zero(I), zero(I), Cert.None, Cert.Outside
        end
        # If we do not have the condition above, it is still possible that q is on the edge.
        if is_collinear(certⱼ)
            # First, let us check if q is on the edge. 
            q_cert = point_position_on_line_segment(p, pⱼ, q)
            if is_on(q_cert) || is_degenerate(q_cert)
                # Here, q is on the edge, so we consider it as being inside the triangle.
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                    add_edge!(history, i, j)
                end
                return i, j, Cert.On, Cert.Inside
            elseif is_right(q_cert)
                # Here, q is to the right of ppⱼ, which just means that it is inside the triangulation. 
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                    add_edge!(history, k, j)
                end
                return j, i, Cert.Right, Cert.Outside # flip i and j to get the correct orientation
            elseif is_left(q_cert)
                # This means that q is left of ppⱼ, but this means that q is away from k, i.e. away from the 
                # boundary node - it is outside of the triangulation. 
                return zero(I), zero(I), Cert.None, Cert.Outside
            end
        end
        # If we still have not returned, then let us rotate onto the next triangle.
        i, pᵢ, certᵢ = j, pⱼ, certⱼ
    end
    # Now, we've made it to this point, but we still need to check the triangle that contains the boundary vertex to the left of p. 
    # This case is kept separate so that we can reuse left_cert
    if num_interior_neighbours > 0 # It is possible that there are actually no interior nodes around the point k 
        j = other_boundary_node
        pⱼ = get_point(pts, representative_point_list, boundary_map, j)
        certⱼ = left_cert
        if (is_left(certᵢ) && is_right(certⱼ)) || (is_right(certᵢ) && is_left(certⱼ))
            intersection_cert = line_segment_intersection_type(p, q, pᵢ, pⱼ)
            if has_one_intersection(intersection_cert)
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return j, i, Cert.Single, Cert.Outside
            end
            pⱼp_edge = point_position_relative_to_line(pⱼ, p, q)
            ppᵢ_edge = point_position_relative_to_line(p, pᵢ, q)
            if is_left(pⱼp_edge) && is_left(ppᵢ_edge)
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return i, j, Cert.None, Cert.Inside
            end
            return zero(I), zero(I), Cert.None, Cert.Outside
        end
        # We do not need to check for collinearities here - this was already done in check_for_intersections_with_adjacent_boundary_edges 
    end
    # If we've made it to this point in the algorithm, then the point is outside of the triangulation. Again, 
    # this is using the assumption that the geometry is convex.
    return zero(I), zero(I), Cert.None, Cert.Outside
end
