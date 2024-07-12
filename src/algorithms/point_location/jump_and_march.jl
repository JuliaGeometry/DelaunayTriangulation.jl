"""
    default_num_samples(n) -> Integer

Given a number of points `n`, returns `∛n` rounded up to the nearest integer. This is the default number of samples used in the jump-and-march algorithm.
"""
default_num_samples(num_points::I) where {I} = ceil(I, cbrt(num_points))


"""
    compare_distance(current_dist, current_idx, pts, i, qx, qy) -> (Number, Vertex)

Computes the minimum of the distance between the `i`th point of `pts` and `(qx, qy)` and `current_dist`.

# Arguments 
- `current_dist`: The current value for the distance to the point `(qx, qy)`.
- `current_idx`: The point of `pts` corresponding to the distance `current_dist`.
- `pts`: The point set. 
- `i`: The vertex to compare with `current_idx`.
- `qx`: The x-coordinate of the query point.
- `qy`: The y-coordinate of the query point.

# Outputs 
- `current_dist`: The minimum of the distance between the `i`th point of `pts` and `(qx, qy)` and `current_dist`.
- `current_idx`: The point of `pts` corresponding to the distance `current_dist`, which will be either `i` or `current_idx`.
"""
function compare_distance(current_dist, current_idx::I, pts, i, qx, qy) where {I}
    p = get_point(pts, i)
    px, py = getxy(p)
    _dist = dist_sqr((px, py), (qx, qy))
    if _dist < current_dist
        current_dist = _dist
        current_idx = i
    end
    return current_dist, I(current_idx)
end


@doc """
    select_initial_point(tri::Triangulation, q; kwargs...) -> Vertex

Selects the initial point for [`find_triangle`](@ref) to start from.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point. Can be either a point or a vertex - if it is a vertex, the corresponding point `get_point(tri, q)` will be used.

# Keyword Arguments 
- `point_indices=each_solid_vertex(tri)`: The indices to sample from. 
- `m=default_num_samples(num_vertices(point_indices))`: The number of samples to take. Replacement is not used, so there may be duplicates.
- `try_points=()`: A list of points to try in addition to those randomly sampled.
- `allow_boundary_points=!is_disjoint(tri)`: Whether to allow boundary points to be selected.
- `rng=Random.default_rng()`: The random number generator to use.

# Outputs
- `i`: The index of the point closest to `q` out of those queried.
"""
select_initial_point
function select_initial_point(
        tri::Triangulation, _q;
        point_indices = each_solid_vertex(tri),
        m = default_num_samples(num_vertices(point_indices)),
        try_points = (),
        allow_boundary_points = !is_disjoint(tri),
        rng::AbstractRNG = Random.default_rng(),
    )
    F = number_type(tri)
    I = integer_type(tri)
    current_dist = typemax(F)
    current_idx = I(first(point_indices))
    _qx, _qy = getxy(_q)
    qx, qy = F(_qx), F(_qy)
    for _ in 1:m # Not using replacement, but probability of duplicates is approximately 0.5num_solid_vertices(tri)^(-1/3)
        i = I(rand(rng, point_indices))
        (is_ghost_vertex(i) || (!allow_boundary_points && is_boundary_node(tri, i)[1])) && continue
        current_dist, current_idx = compare_distance(current_dist, current_idx, tri, i, qx, qy)
    end
    for i in try_points
        if i ≠ ∅
            current_dist, current_idx = compare_distance(current_dist, current_idx, tri, i, qx, qy)
        end
    end
    return I(current_idx)
end
function select_initial_point(
        tri::Triangulation, q::Integer;
        point_indices = each_solid_vertex(tri),
        m = default_num_samples(num_vertices(tri)),
        allow_boundary_points = !is_disjoint(tri),
        try_points = (),
        rng::AbstractRNG = Random.default_rng(),
    )
    return select_initial_point(tri, get_point(tri, q); point_indices, m, try_points, allow_boundary_points, rng)
end


"""
    select_random_edge(tri::Triangulation, edges, rng::AbstractRNG=Random.default_rng()) -> (Vertex, Vertex, Point, Point)

Selects a random edge from the set of edges `edges`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `edges`: The set of edges to sample from.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.

# Outputs
- `i`: The initial vertex of the edge.
- `j`: The terminal vertex of the edge.
- `pᵢ`: The point corresponding to `i`.
- `pⱼ`: The point corresponding to `j`.
"""
function select_random_edge(tri::Triangulation, edges, rng::AbstractRNG = Random.default_rng())
    edge = random_edge(rng, edges)
    i, j = edge_vertices(edge)
    pᵢ, pⱼ = get_point(tri, i, j)
    return i, j, pᵢ, pⱼ
end


"""
    prepare_initial_edge(tri::Triangulation, edges, p, q, rng::AbstractRNG=Random.default_rng()) -> (Vertex, Vertex, Point, Point, Certificate, Certificate)

Selects a random edge from the set of edges `edges` and computes the certificates for the points corresponding to the initial and terminal vertices of the edge.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `edges`: The set of edges to sample from.
- `p`: The initial point.
- `q`: The query point.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.

# Outputs
- `i`: The initial vertex of the edge.
- `j`: The terminal vertex of the edge.
- `pᵢ`: The point corresponding to `i`.
- `pⱼ`: The point corresponding to `j`.
- `line_cert_i`: The [`Certificate`](@ref) for `pᵢ`'s position relative to the oriented line `pq`.
- `line_cert_j`: The [`Certificate`](@ref) for `pⱼ`'s position relative to the oriented line `pq`.
"""
function prepare_initial_edge(tri::Triangulation, edges, p, q, rng::AbstractRNG = Random.default_rng())
    i, j, pᵢ, pⱼ = select_random_edge(tri, edges, rng)
    line_cert_i = point_position_relative_to_line(p, q, pᵢ)
    line_cert_j = point_position_relative_to_line(p, q, pⱼ)
    return i, j, pᵢ, pⱼ, line_cert_i, line_cert_j
end


"""
    select_initial_triangle_interior_vertex(tri::Triangulation, k, q, store_history=Val(false), history=nothing, rng::AbstractRNG=Random.default_rng()) -> (Point, Vertex, Vertex, Point, Point)

Selects the initial triangle for [`find_triangle`](@ref) to start from, for the case where `k` is an interior vertex.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The vertex to start from.
- `q`: The query point.
- `store_history=Val(false)`: Whether to store the history of the algorithm.
- `history=nothing`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use.

# Outputs
- `p`: The point corresponding to `k`.
- `i`: The initial vertex of the triangle. 
- `j`: The terminal vertex of the triangle.
- `pᵢ`: The point corresponding to `i`.
- `pⱼ`: The point corresponding to `j`.

!!! warning "Non-convex geometries"

    This function assumes that the geometry is convex. To deal with this, when an infinite loop is detected 
    we return `∅` for both `i` and `j`, and then let [`find_triangle`](@ref) handle how to 
    correct the algorithm from there.    

# Extended help
This part of the algorithm works by rotating around the vertex `k`, looking for a triangle whose edges adjoining `k` 
are to the left and to the right of `k`. By choosing the initial edge at random via [`prepare_initial_edge`](@ref), 
and computing the position of `q` relative to this initial edge, the rotation will be either clockwise or counter-clockwise, 
and the triangle is then found using either `select_initial_triangle_clockwise` or `select_initial_triangle_counterclockwise`,
respectively. 

In case the initial edge is collinear with the line `pq`, where `p = get_point(tri, q)`, then `fix_initial_collinear_edge_for_interior_vertex` to find a 
non-collinear edge resample more edges from [`prepare_initial_edge`](@ref) if necessary.
"""
function select_initial_triangle_interior_vertex(tri::Triangulation, k, q, store_history::F = Val(false), history = nothing, rng::AbstractRNG = Random.default_rng()) where {F}
    p = get_point(tri, k)
    ## Select the initial edge to rotate about
    neighbouring_edges = get_adjacent2vertex(tri, k)
    i, j, pᵢ, pⱼ, line_cert_i, line_cert_j = prepare_initial_edge(tri, neighbouring_edges, p, q, rng)
    p == q && return p, j, i, pⱼ, pᵢ
    ## Take care for collinear edges
    return_flag, p, i, j, pᵢ, pⱼ, line_cert_i, line_cert_j = fix_initial_collinear_edge_for_interior_vertex(tri, k, q, store_history, history, rng, p, neighbouring_edges, i, j, pᵢ, pⱼ, line_cert_i, line_cert_j)
    return_flag && return p, i, j, pᵢ, pⱼ
    ## Now rotate around to find a triangle that the line pq intersects through 
    if is_left(line_cert_j)
        i, j, pᵢ, pⱼ = select_initial_triangle_clockwise(tri, p, q, pᵢ, pⱼ, i, j, k, store_history, history)
    else
        i, j, pᵢ, pⱼ = select_initial_triangle_counterclockwise(tri, line_cert_j, p, q, pᵢ, pⱼ, i, j, k, store_history, history)
    end
    ## Swap values and return
    i, j = j, i
    pᵢ, pⱼ = pⱼ, pᵢ # pᵢ is left of pq, pⱼ is right of pq 
    return p, i, j, pᵢ, pⱼ
end
function select_initial_triangle_clockwise(tri::Triangulation, p, q, pᵢ, pⱼ, i, j, k, store_history::F = Val(false), history = nothing) where {F}
    line_cert_i = point_position_relative_to_line(p, q, pᵢ)
    if is_true(store_history) && is_collinear(line_cert_i)
        add_edge!(history, k, i)
    end
    nn = num_neighbours(tri, k)
    iter = 0 # when we have Non-convex geometries, sometimes we get stuck in an infinite loop
    while is_left(line_cert_i) && iter ≤ nn + 1
        iter += 1
        j = i
        pⱼ = pᵢ
        i = get_adjacent(tri, i, k)
        pᵢ = get_point(tri, i)
        line_cert_i = point_position_relative_to_line(p, q, pᵢ)
        if is_true(store_history) && is_collinear(line_cert_i)
            add_edge!(history, k, i)
        end
    end
    I = integer_type(tri)
    if iter > nn + 1
        return I(∅), I(∅), pᵢ, pⱼ
    else
        return i, j, pᵢ, pⱼ
    end
end
function select_initial_triangle_counterclockwise(tri::Triangulation, line_cert_j, p, q, pᵢ, pⱼ, i, j, k, store_history::F = Val(false), history = nothing) where {F}
    if is_true(store_history) && is_collinear(line_cert_j)
        add_edge!(history, k, j)
    end
    nn = num_neighbours(tri, k)
    iter = 0 # when we have Non-convex geometries, sometimes we get stuck in an infinite loop
    while is_right(line_cert_j) && iter ≤ nn + 1
        iter += 1
        i = j
        pᵢ = pⱼ
        j = get_adjacent(tri, k, j)
        pⱼ = get_point(tri, j)
        line_cert_j = point_position_relative_to_line(p, q, pⱼ)
        if is_true(store_history) && is_collinear(line_cert_j)
            add_edge!(history, k, j)
        end
    end
    I = integer_type(tri)
    if iter > nn + 1
        return I(∅), I(∅), pᵢ, pⱼ
    else
        return i, j, pᵢ, pⱼ
    end
end
function fix_initial_collinear_edge_for_interior_vertex(tri::Triangulation, k, q, store_history, history, rng, p, neighbouring_edges, i, j, pᵢ, pⱼ, line_cert_i, line_cert_j)
    while is_collinear(line_cert_j) || is_collinear(line_cert_i)
        if is_collinear(line_cert_j)
            opposing_idx = j
            on_edge_cert = point_position_on_line_segment(p, pⱼ, q)
        else
            opposing_idx = i
            on_edge_cert = point_position_on_line_segment(p, pᵢ, q)
        end
        # We want q to be in the direction of of ppᵢ or ppⱼ, which just means not left
        # For example, suppose we have 
        #=
          p -------- pⱼ --------- q 
           \        /
            \      /
             \    /
              \  /
               \pᵢ
        Then we want to know if q is to the right of `ppⱼ` (here, `line_cert_j` is the collinear certificate) 
        or to the left. If it is to the left, then it is on the edge `ppⱼ`. This is not a problem, but we let 
        [`find_triangle`](@ref) handle this case somewhere else. Thus, we should want to find a case where `q` 
        is right of `ppⱼ` and not handle this interior edge case, so we check `!is_left(on_edge_cert)`. 
        =#
        if !is_left(on_edge_cert)
            if is_on(on_edge_cert) || is_degenerate(on_edge_cert)
                is_true(store_history) && add_edge!(history, k, opposing_idx)
                return true, p, j, i, pⱼ, pᵢ, line_cert_i, line_cert_j
            else
                pos_cert = point_position_relative_to_line(p, q, pᵢ)
                if is_left(pos_cert)
                    is_true(store_history) && add_edge!(history, k, opposing_idx)
                    return true, p, i, j, pᵢ, pⱼ, line_cert_i, line_cert_j
                else
                    is_true(store_history) && add_edge!(history, k, opposing_idx)
                    return true, p, j, i, pⱼ, pᵢ, line_cert_i, line_cert_j
                end
            end
        else
            # If we haven't yet found the edge, let's just pick another one
            i, j, pᵢ, pⱼ, line_cert_i, line_cert_j = prepare_initial_edge(tri, neighbouring_edges, p, q, rng)
        end
    end
    return false, p, i, j, pᵢ, pⱼ, line_cert_i, line_cert_j
end


"""
    check_for_intersections_with_adjacent_boundary_edges(tri::Triangulation, k, q, ghost_vertex=𝒢) -> (Certificate, Certificate, Vertex, Certificate, Certificate)

Given a boundary vertex `k`, find a triangle adjacent to `k` to locate a triangle or edge containing `q`.

See also [`search_down_adjacent_boundary_edges`](@ref), which uses this function to determine an initial direction to search along a 
straight boundary in case `q` is collinear with it.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The boundary vertex to start from.
- `q`: The query point.
- `ghost_vertex=𝒢`: The ghost vertex corresponding to the boundary that `k` resides on.

# Outputs 
- `direction_cert`: The direction of `q` relative to the vertex `k` along the boundary, given as a [`Certificate`](@ref) `Left`, `Right`, or `Outside`. If `is_outside(direction_cert)`, then `q` is not collinear with either of the adjacent boundary edges.
- `q_pos_cert`: The position of `q` relative to the vertex `k` along the boundary, given as a [`Certificate`](@ref) `Left`, `Right`, `On`, `Outside`, or `Degenerate`. This is similar to `direction_cert` in that it will be `Outside` whenever `direction_cert` is, but this certificate can also be `On` to indicate that not only is `q` in the direction given by `direction_cert`, but it is directly on the edge in that direction. If `is_degnerate(q_pos_cert)`, then `q = get_point(tri, next_vertex)`.
- `next_vertex`: The next vertex along the boundary in the direction of `q`, or `k` if `q` is not collinear with either of the adjacent boundary edges.
- `right_cert`: The [`Certificate`](@ref) for the position of `q` relative to the boundary edge right of `k`.
- `left_cert`: The [`Certificate`](@ref) for the position of `q` relative to the boundary edge left of `k`.
"""
function check_for_intersections_with_adjacent_boundary_edges(tri::Triangulation{P, T, BN, W, I}, k, q, ghost_vertex = I(𝒢)) where {P, T, BN, W, I}
    p = get_point(tri, k)
    right = get_right_boundary_node(tri, k, ghost_vertex)
    left = get_left_boundary_node(tri, k, ghost_vertex)
    pright, pleft = get_point(tri, right, left)
    right_cert = point_position_relative_to_line(p, pright, q)
    left_cert = point_position_relative_to_line(pleft, p, q)
    flipped_left_cert = point_position_relative_to_line(p, pleft, q) # Useful for the returned value, but we use left_cert below to get a good value for direction_cert
    !is_collinear(right_cert) && !is_collinear(left_cert) && return (Cert.Outside, Cert.Outside, k, right_cert, flipped_left_cert)
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
    search_down_adjacent_boundary_edges(tri::Triangulation, k, q, direction_cert, q_pos_cert, next_vertex, store_history=Val(false), history=nothing, ghost_vertex=𝒢) -> (Bool, Certificate, Vertex, Vertex, Vertex)

Starting at the boundary vertex `k`, walks down the boundary in the direction of `q` until finding `q` or finding that it is outside of the triangulation.

See also [`check_for_intersections_with_adjacent_boundary_edges`](@ref), which uses this function to determine an initial direction to search along.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The boundary vertex to start from.
- `q`: The query point.
- `direction_cert`: The direction of `q` relative to the vertex `k` along the boundary, defined from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `q_pos_cert`: The position of `q` relative to the vertex `k` along the boundary, defined from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `next_vertex`: The next vertex along the boundary in the direction of `q`, defined from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `store_history=Val(false)`: Whether to store the history of the algorithm.
- `history=nothing`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `ghost_vertex=𝒢`: The ghost vertex corresponding to the boundary that `k` resides on.

# Outputs 
- `return_flag`: Whether to return, or throw an exception.
- `q_pos_cert`: A [`Certificate`](@ref) that is `On` if `q` is on the edge `(u, v)`, and `Outside` if `q` is outside of the triangulation.
- `u`: If `is_on(q_pos_cert)`, this is the first vertex of a positively oriented triangle that `q` is on, so that `q` is on the edge `(u, v)`. Otherwise, `(u, v, w)` is a ghost triangle close to `q`.
- `v`: If `is_on(q_pos_cert)`, this is the second vertex of a positively oriented triangle that `q` is on, so that `q` is on the edge `(u, v)`. Otherwise, `(u, v, w)` is a ghost triangle close to `q`.
- `w`: If `is_on(q_pos_cert)`, this is the third vertex of a positively oriented triangle that `q` is on, so that `q` is on the edge `(u, v)` and `w = get_adjacent(tri, u, v)`. Otherwise, `(u, v, w)` is a ghost triangle close to `q`.

!!! warning "Non-convex geometries"

    This function assumes that the geometry is convex. The function will still be able to return, but `is_outside(q_pos_cert)` may not necessarily mean `q` 
    is outside of the triangulation. The main function [`find_triangle`](@ref) will have to restart the algorithm if it is found that `is_outside(q_pos_cert)` 
    was incorrect.

# Extended help 
This function works by stepping along vertices on the boundaries in the direction specified by `direction_cert`, using `search_right_down_adjacent_boundary_edges`
if `is_right(direction_cert)` and `search_left_down_adjacent_boundary_edges` otherwise. In these functions, a `while` loop is used to keep stepping until `q_pos_cert`,
which is updated at each iteration, changes value.
"""
function search_down_adjacent_boundary_edges(tri::Triangulation, k, q, direction, q_pos, next_vertex, store_history::F = Val(false), history = nothing, ghost_vertex = integer_type(tri)(𝒢)) where {F}
    i = k
    j = next_vertex
    pⱼ = get_point(tri, j)
    if is_right(direction)
        return_flag, cert, i, j, k = search_right_down_adjacent_boundary_edges(tri, q, q_pos, store_history, history, ghost_vertex, i, j, pⱼ)
    else
        return_flag, cert, i, j, k = search_left_down_adjacent_boundary_edges(tri, q, q_pos, store_history, history, ghost_vertex, i, j, pⱼ)
    end
    return_flag && return cert, i, j, k
    throw(PointNotFoundError(tri, q))
end
function search_right_down_adjacent_boundary_edges(tri::Triangulation, q, q_pos, store_history::F, history, ghost_vertex, i, j, pⱼ) where {F}
    if is_true(store_history) # in case we don't enter the loop
        k′ = get_adjacent(tri, i, j)
        add_triangle!(history, i, j, k′)
        add_edge!(history, i, j)
    end
    while is_right(q_pos)
        i, pᵢ = j, pⱼ
        j = get_right_boundary_node(tri, i, ghost_vertex)
        if is_true(store_history)
            k′ = get_adjacent(tri, i, j)
            add_triangle!(history, i, j, k′)
            add_edge!(history, i, j)
        end
        pⱼ = get_point(tri, j)
        right_cert = point_position_relative_to_line(pᵢ, pⱼ, q)
        q_pos = !is_collinear(right_cert) ? Cert.Outside : point_position_on_line_segment(pᵢ, pⱼ, q)
    end
    if is_outside(q_pos)
        return (true, Cert.Outside, j, i, get_adjacent(tri, j, i))
    elseif is_on(q_pos) || is_degenerate(q_pos)
        k = get_adjacent(tri, i, j)
        return (true, Cert.On, i, j, k) # true is the return_flag
    end
    return (false, Cert.Outside, i, j, k)
end
function search_left_down_adjacent_boundary_edges(tri::Triangulation, q, q_pos, store_history::F, history, ghost_vertex, i, j, pⱼ) where {F}
    if is_true(store_history) # in case we don't enter the loop
        k′ = get_adjacent(tri, j, i)
        add_triangle!(history, j, i, k′)
        add_edge!(history, i, j) # i,j → j, i to get the ordering of the segments
    end
    while is_left(q_pos)
        i, pᵢ = j, pⱼ
        j = get_left_boundary_node(tri, i, ghost_vertex)
        if is_true(store_history)
            k′ = get_adjacent(tri, j, i)
            add_triangle!(history, j, i, k′)
            add_edge!(history, i, j)
        end
        pⱼ = get_point(tri, j)
        left_cert = point_position_relative_to_line(pᵢ, pⱼ, q)
        q_pos = !is_collinear(left_cert) ? Cert.Outside : point_position_on_line_segment(pⱼ, pᵢ, q)
    end
    if is_outside(q_pos)
        return (true, Cert.Outside, i, j, get_adjacent(tri, i, j))
    elseif is_on(q_pos) || is_degenerate(q_pos)
        k = get_adjacent(tri, j, i)
        return (true, Cert.On, j, i, k)
    end
    return (false, Cert.Outside, i, j, k)
end


"""
    check_for_intersections_with_interior_edges_adjacent_to_boundary_vertex(tri::Triangulation, k, q, right_cert, left_cert, store_history=Val(false), history=nothing, ghost_vertex=𝒢) -> (Bool, Vertex, Vertex, Certificate, Certificate)

Checks for intersections between the line `pq`, where `p = get_point(tri, k)`, and the edges neighbouring `p`, assuming `k` is a boundary node. This function should only be used 
after using [`check_for_intersections_with_adjacent_boundary_edges`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The boundary vertex to start from.
- `q`: The query point.
- `right_cert`: The [`Certificate`](@ref) for the position of `q` relative to the boundary edge right of `k`, coming from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `left_cert`: The [`Certificate`](@ref) for the position of `q` relative to the boundary edge left of `k`, coming from [`check_for_intersections_with_adjacent_boundary_edges`](@ref).
- `store_history=Val(false)`: Whether to store the history of the algorithm.
- `history=nothing`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `ghost_vertex=𝒢`: The ghost vertex corresponding to the boundary that `k` resides on.

# Outputs 
The output takes the form `(i, j, edge_cert, triangle_cert)`. Rather than defining each output individually, here are the possible froms of the output:
- `(i, j, Single, Outside)`: The line `pq` intersects the edge `pᵢpⱼ` and `(j, i, k)` is a positively oriented triangle so that `pᵢ` is left of `pq` and `pⱼ` is right of `pq`.
- `(i, j, None, Inside)`: The point `q` is inside the positively oriented triangle `(i, j, k)`.
- `(0, 0, None, Outside)`: The point `q` is outside of the triangulation.
- `(i, j, On, Inside)`: The point `q` is on the edge `pᵢpⱼ`, and thus inside the positively oriented triangle `(i, j, k)`.
- `(i, j, Right, Outside)`:` The point `q` is collinear with the edge `pᵢpⱼ`, but is off of it and further into the triangulation. 

!!! warning "Non-convex geometries"

    This function assumes that the geometry is convex.

# Extended help
This function works in two stages. Firstly, using `check_for_intersections_with_single_interior_edge_adjacent_to_boundary_vertex`, we check for the intersection of `pq` with the edges neighbouring 
the vertex `k`, rotating counter-clockwise until we find an intersection or reach the other side of the boundary, starting from the first edge counter-clockwise away from the boundary edge right of 
the vertex `k`. By keeping track of the positions of `pq` relative to the current vertex and the previous, we can identify when an intersection is found. If no intersection is found before 
reaching the boundary edge left of `k`, then `check_for_intersections_with_triangle_left_to_boundary_vertex` is used to check the remaining triangle.
"""
function check_for_intersections_with_interior_edges_adjacent_to_boundary_vertex(tri::Triangulation, k, q, right_cert, left_cert, store_history::F = Val(false), history = nothing, ghost_vertex = integer_type(tri)(𝒢)) where {F}
    I = integer_type(tri)
    p = get_point(tri, k)
    other_boundary_node = get_left_boundary_node(tri, k, ghost_vertex)
    num_interior_neighbours = num_neighbours(tri, k) - 3 # 3 = + the two boundary neighbours + the 𝒢
    i = get_right_boundary_node(tri, k, ghost_vertex)
    pᵢ = get_point(tri, i)
    certᵢ = right_cert # Do not need to check if is_collinear(certᵢ) - this function should only be used after checking check_for_intersections_with_adjacent_boundary_edges anyway
    for _ in 1:max(num_interior_neighbours, 1)
        # In the above, we need max(num_interior_neighbours, 1) in case the point 
        # being considered is a corner point on the boundary with no interior edge neighbours. 
        # If we just had 1:num_interior_neighbours, then nothing would happen in this loop 
        # since num_interior_neighbours = 0, but we still need to make sure we check the edge 
        # connecting the two boundary neighbours
        return_flag, _i, _j, edge_cert, triangle_cert, i, pᵢ, certᵢ = check_for_intersections_with_single_interior_edge_adjacent_to_boundary_vertex(tri, k, q, store_history, history, p, i, pᵢ, certᵢ)
        return_flag && return _i, _j, edge_cert, triangle_cert
    end
    # Now, we've made it to this point, but we still need to check the triangle that contains the boundary vertex to the left of p. 
    # This case is kept separate so that we can reuse left_cert
    return_flag, _i, _j, edge_cert, triangle_cert = check_for_intersections_with_triangle_left_to_boundary_vertex(tri, k, q, left_cert, store_history, history, p, other_boundary_node, num_interior_neighbours, i, pᵢ, certᵢ)
    return_flag && return _i, _j, edge_cert, triangle_cert
    # If we've made it to this point in the algorithm, then the point is outside of the triangulation. Again, 
    # this is using the assumption that the geometry is convex.
    return zero(I), zero(I), Cert.None, Cert.Outside
end
function check_for_intersections_with_single_interior_edge_adjacent_to_boundary_vertex(tri::Triangulation, k, q, store_history::F, history, p, i, pᵢ, certᵢ) where {F}
    j = get_adjacent(tri, k, i)
    pⱼ = get_point(tri, j)
    certⱼ = point_position_relative_to_line(p, pⱼ, q)
    I = integer_type(tri)
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
            is_true(store_history) && add_triangle!(history, i, j, k)
            return true, j, i, Cert.Single, Cert.Outside, i, pᵢ, certᵢ # true is the return_flag, and the last three returns are for type stability at the end of the function, in case we need to continue looping in the main function
        elseif is_touching(intersection_cert)
            # q is on the edge 
            if is_true(store_history)
                add_triangle!(history, i, j, k)
                add_edge!(history, i, j)
            end
            return true, i, j, Cert.On, Cert.Inside, i, pᵢ, certᵢ
        end
        # There may be no intersection, but it could be inside the triangle. 
        pⱼp_edge = point_position_relative_to_line(pⱼ, p, q)
        ppᵢ_edge = point_position_relative_to_line(p, pᵢ, q)
        if is_left(pⱼp_edge) && is_left(ppᵢ_edge)
            if is_true(store_history)
                add_triangle!(history, i, j, k)
            end
            return true, i, j, Cert.None, Cert.Inside, i, pᵢ, certᵢ
        end
        # If none of the above lead to a return, then q must be on the other side of p away from the interior, meaning 
        # it is in the exterior. Here, we are making use of the fact that the domain is 
        # convex. When the domain is not convex, this won't necessarily work.
        return true, zero(I), zero(I), Cert.None, Cert.Outside, i, pᵢ, certᵢ
    end
    # If we do not have the condition above, it is still possible that q is on the edge.
    if is_collinear(certⱼ)
        # First, let us check if q is on the edge. 
        q_cert = point_position_on_line_segment(p, pⱼ, q)
        if is_on(q_cert) || is_degenerate(q_cert)
            # Here, q is on the edge, so we consider it as being inside the triangle.
            if is_true(store_history)
                add_triangle!(history, i, j, k)
                pᵢ ≠ q && pⱼ ≠ q && add_edge!(history, i, j) # be careful not to add a collinear edge that is just starting at q
            end
            return true, i, j, Cert.On, Cert.Inside, i, pᵢ, certᵢ
        elseif is_right(q_cert)
            # Here, q is to the right of ppⱼ, which just means that it is inside the triangulation. 
            if is_true(store_history)
                add_triangle!(history, i, j, k)
                p ≠ q && pⱼ ≠ q && add_edge!(history, k, j) # be careful not to add a collinear edge that is just starting at q
            end
            return true, j, i, Cert.Right, Cert.Outside, i, pᵢ, certᵢ # flip i and j to get the correct orientation
        elseif is_left(q_cert)
            # This means that q is left of ppⱼ, but this means that q is away from k, i.e. away from the 
            # boundary node - it is outside of the triangulation. 
            return true, zero(I), zero(I), Cert.None, Cert.Outside, i, pᵢ, certᵢ
        end
    end
    # If we still have not returned, then let us rotate onto the next triangle.
    return false, zero(I), zero(I), Cert.None, Cert.None, j, pⱼ, certⱼ
end
function check_for_intersections_with_triangle_left_to_boundary_vertex(tri::Triangulation, k, q, left_cert, store_history::F, history, p, other_boundary_node, num_interior_neighbours, i, pᵢ, certᵢ) where {F}
    I = integer_type(tri)
    if num_interior_neighbours > 0 # It is possible that there are actually no interior nodes around the point k 
        j = other_boundary_node
        pⱼ = get_point(tri, j)
        certⱼ = left_cert
        if (is_left(certᵢ) && is_right(certⱼ)) || (is_right(certᵢ) && is_left(certⱼ))
            intersection_cert = line_segment_intersection_type(p, q, pᵢ, pⱼ)
            if has_one_intersection(intersection_cert)
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return true, j, i, Cert.Single, Cert.Outside
            end
            pⱼp_edge = point_position_relative_to_line(pⱼ, p, q)
            ppᵢ_edge = point_position_relative_to_line(p, pᵢ, q)
            if is_left(pⱼp_edge) && is_left(ppᵢ_edge)
                if is_true(store_history)
                    add_triangle!(history, i, j, k)
                end
                return true, i, j, Cert.None, Cert.Inside
            end
            return true, zero(I), zero(I), Cert.None, Cert.Outside
        end
        # We do not need to check for collinearities here - this was already done in check_for_intersections_with_adjacent_boundary_edges 
    end
    return false, zero(I), zero(I), Cert.None, Cert.None
end


"""
    exterior_find_triangle(tri::Triangulation, k, q, ghost_vertex=𝒢) -> (Vertex, Vertex)

Given a query point `q` outside of the triangulation, find the ghost triangle containing `q`.

# Arguments
- `tri`: The [`Triangulation`](@ref).
- `k`: The exterior boundary vertex to start from.
- `q`: The query point. 
- `ghost_vertex=𝒢`: The ghost vertex corresponding to the boundary that `k` resides on.

# Outputs
- `i`: The first vertex of the edge on the boundary adjoining the positively oriented ghost triangle. 
- `j`: The second vertex of the edge on the boundary adjoining the positively oriented ghost triangle.

!!! warning "Non-convex geometries and internal query points"

    This function assumes that the geometry is convex. If the geometry is not convex, the returned value may not be correct and should be  
    checked separately. Additionally, if `q` is actually inside the triangulation, then the result is meaningless.

# Extended help
This function works by first finding the position of `q` relative to `pₘp`, where `pₘ` is the representative point for 
the `ghost_vertex` and `p = get_point(tri, k)`. Depending on this position, we rotate counter-clockwise if `q` is 
left of the line using `exterior_find_triangle_rotate_left` and clockwise otherwise using `exterior_find_triangle_rotate_right`.
By keeping track of the current position of `q` and its position relative to the next ghost edge, we can identify when `q` 
resides inside a ghost triangle.
"""
function exterior_find_triangle(tri::Triangulation, k, q, ghost_vertex = integer_type(tri)(𝒢))
    pₘ, pᵢ = get_point(tri, ghost_vertex, k)
    i = k
    q_position = point_position_relative_to_line(pₘ, pᵢ, q)
    if is_left(q_position) # q is left of the ghost edge through pᵢ, so rotate left 
        return exterior_find_triangle_rotate_left(tri, q, i, pₘ, ghost_vertex)
    else # rotate right 
        return exterior_find_triangle_rotate_right(tri, q, i, pₘ, ghost_vertex)
    end
end
function exterior_find_triangle_rotate_left(tri, q, i, pₘ, ghost_vertex)
    j = get_right_boundary_node(tri, i, ghost_vertex)
    pⱼ = get_point(tri, j)
    new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
    while is_left(new_q_pos)
        i = j
        j = get_right_boundary_node(tri, i, ghost_vertex)
        pⱼ = get_point(tri, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
    end
    return j, i  # Swap the orientation so that i, j is a boundary edge 
end
function exterior_find_triangle_rotate_right(tri, q, i, pₘ, ghost_vertex)
    j = get_left_boundary_node(tri, i, ghost_vertex)
    pⱼ = get_point(tri, j)
    new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
    while is_right(new_q_pos)
        i = j
        j = get_left_boundary_node(tri, i, ghost_vertex)
        pⱼ = get_point(tri, j)
        new_q_pos = point_position_relative_to_line(pₘ, pⱼ, q)
    end
    return i, j
end


"""
    find_triangle(tri, q; kwargs...) -> Triangle[, Bool] 

Find the triangle in the triangulation `tri` containing the query point `q` using the jump-and-march algorithm.

# Arguments
- `tri`: The [`Triangulation`](@ref).
- `q`: The query point.

# Keyword Arguments 
- `point_indices=each_solid_vertex(tri)`: The indices of the vertices to consider as possible starting points for the algorithm.
- `m=default_num_samples(num_vertices(point_indices))`: The number of samples to use when selecting the initial point.
- `try_points=()`: A list of points to try as the initial point in addition to the `m` sampled.
- `rng=Random.default_rng()`: The random number generator to use.
- `k=select_initial_point(tri, q; point_indices, m, try_points, rng)`: The initial point to start the algorithm from. See [`select_initial_point`](@ref).
- `store_history=Val(false)`: Whether to store the history of the algorithm.
- `history=nothing`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `maxiters=2 + num_exterior_curves(tri) - num_solid_vertices(tri) + num_solid_edges(tri)`: The maximum number of iterations to perform before restarting the algorithm with [`restart_find_triangle`](@ref).

!!! note "Restarting the algorithm"

    If the algorithm restarts, then the initial point `k` is selected again using [`select_initial_point`](@ref), and the algorithm is restarted from there. 
    This is done if the algorithm gets stuck in a loop, or if the algorithm is not able to find the correct triangle containing `q` after `maxiters` iterations. For a convex 
    geometry, `maxiters` can be safely ignored, as the sequence of triangles visited is acyclic [see H. Edelsbrunner, An acyclicity theorem for cell complexes in d dimensions, Combinatorica 10 (1990) 251-260.)].
    
- `concavity_protection=false`: Whether to use concavity protection. See [`concavity_protection_check`](@ref). This is only needed if your triangulation is not convex. 
- `use_barriers::Val{U}=Val(false)`: Whether to stop searching beyond any segments in the triangulation. 

!!! warning "Found triangles with barriers"

    If you are using barriers, it will be your responsibility to verify that any found triangle from this function actually 
    contains the triangle. This can be verified using the returned `flag` (see below), although the point might still be on the triangle's 
    boundary.

!!! warning "Walking past vertices of barriers"

    If you are using barriers, it is possible that the algorithm can walk past vertices of barriers. One way this can happen is if 
    the initial search line intersects a vertex, in which case the segment might not be considered. Another way this can happen is if you 
    start the algorithm directly on a segment vertex, in which case the algorithm can go past it (e.g. this means that it is possible that a 
    ghost triangle might still be returned if you start the algorithm on a boundary node).

# Output
- `V`: The triangle containing `q`, with type given by `triangle_type(tri)`.

!!! danger "Hitting barriers"

    If a barrier is hit before any initial triangle is properly identified, the returned triangle is 
    `($∅, $∅, $∅)`; this is only possible if `use_barriers == Val(true)`. Moreover, if `use_barriers == Val(true)`, 
    the final triangle may not even be valid if `invisible_flag == true` (defined below).

If you have `use_barriers == Val(true)`, then we also return 

- `invisible_flag`: `false` if the triangle was found without hitting a barrier, and `true` otherwise.

!!! warning "Non-convex geometries"

    While this function does still work for non-convex geometries, it may be significantly slower than for convex geometries, as most of the details 
    of the algorithm assume that the geometry is convex, and so the algorithm may have to restart many times at new initial vertices `k`.

!!! warning "Ghost triangles"

    For this function to work best, the triangulation should have ghost triangles, which you can add using `add_ghost_triangles!` in case `tri` does not already have them. 
    Without ghost triangles, the function may not be able to find the correct triangle containing `q`.

# Extended help 
The algorithm underlying this function is complicated and broken into many parts. Here, we describe a brief overview of the algorithm, but note that the 
    documentation contains a much more detailed description.
    
1. Firstly, the algorithm is initialised depending on whether `k` is a boundary or an interior vertex, using 
     [`initialise_find_triangle_boundary_vertex`](@ref) or [`initialise_find_triangle_interior_vertex`](@ref) respectively.
2. From the initial triangle `(i, j, k)` chosen, we then check if `q` is one of `pᵢ`, `pⱼ`, and `p = pₖ` and then return according to [`find_triangle_return_on_vertex`](@ref) if needed.
3. If we do not return above, we need to step from the initial triangle towards `q`. Since we put `pᵢ` and `pⱼ`
     to the left and right of the line `pq`, respectively, this means that we step until the triangle `pᵢpⱼq` is no longer 
     positively oriented. So, while the triangle is positively oriented, we step according to [`find_triangle_across_triangle`](@ref).
4. If we have not yet returned and the triangle is no longer positively oriented, we check if the triangle is degenerate using [`find_triangle_degenerate_arrangement`](@ref)
     and reinitialise the algorithm if needed. Otherwise, we have found the triangle containing `q` and return the triangle.
"""
function find_triangle(
        tri::Triangulation, _q;
        point_indices = each_solid_vertex(tri),
        m = default_num_samples(num_vertices(point_indices)),
        try_points = (),
        rng::AbstractRNG = Random.default_rng(),
        k = select_initial_point(tri, _q; point_indices, m, try_points, rng),
        store_history::F = Val(false),
        history = nothing,
        maxiters = 2 + num_exterior_curves(tri) - num_solid_vertices(tri) + num_solid_edges(tri),
        concavity_protection = false,
        use_barriers::Val{U} = Val(false),
    ) where {F, U}
    I = integer_type(tri)
    maxiters = Int(maxiters)
    G = number_type(tri)
    _qx, _qy = getxy(_q)
    q = (G(_qx), G(_qy))
    return _find_triangle(tri, q, I(k), store_history, history, rng, maxiters, zero(maxiters), concavity_protection, zero(maxiters), use_barriers)
end


"""
    initialise_find_triangle_interior_vertex(tri::Triangulation, q, k, store_history::F, history, rng) -> (Bool, Point, Vertex, Vertex, Point, Point)

Initialise the jump-and-march algorithm for an interior vertex `k`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `k`: The interior vertex to start from.
- `store_history`: Whether to store the history of the algorithm.
- `history`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `rng`: The random number generator to use.

# Outputs
- `restart_flag`: Whether the algorithm needs to be restarted.
- `p`: The point corresponding to the vertex `k`.
- `i`: The first vertex of the triangle adjoining `k` to start from.
- `j`: The second vertex of the triangle adjoining `k` to start from.
- `pᵢ`: The point corresponding to the vertex `i`.
- `pⱼ`: The point corresponding to the vertex `j`.

# Extended help
This function works by simply using [`select_initial_triangle_interior_vertex`](@ref) to find the initial triangle to start from. A check is made 
to see if the edge `(i, j)` refers to a non-existent edge `($∅, $∅)`, in which case the algorithm needs to be restarted.
"""
function initialise_find_triangle_interior_vertex(tri::Triangulation, q, k, store_history::F, history, rng) where {F}
    # If k is not a boundary node, then we can rotate around the point k to find an initial triangle 
    # to start the search. Additionally, if there are no ghost triangles, we can only hope that q 
    # is inside the interior, meaning we should only search for the initial triangle here anyway.
    p, i, j, pᵢ, pⱼ = select_initial_triangle_interior_vertex(tri, k, q, store_history, history, rng)
    if !edge_exists(i) && !edge_exists(j) # When we find a possible infinite loop, we use i==j==∅. Let's reinitialise. 
        return true, p, i, j, pᵢ, pⱼ # true is the restart_flag
    end
    is_true(store_history) && add_triangle!(history, j, i, k)
    return false, p, i, j, pᵢ, pⱼ
end


"""
    initialise_find_triangle_boundary_vertex(tri::Triangulation, q, k, store_history:, history, ghost_vertex, concavity_protection) -> (Bool, Bool, Triangle, Point, Vertex, Vertex, Point, Point)

Initialise the jump-and-march algorithm for a boundary vertex `k`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `k`: The boundary vertex to start from.
- `store_history`: Whether to store the history of the algorithm.
- `history`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `ghost_vertex`: The ghost vertex corresponding to the boundary that `k` resides on.
- `concavity_protection`: Whether to use concavity protection. See [`concavity_protection_check`](@ref). This is only needed if your triangulation is not convex.

# Outputs
- `restart_flag`: Whether the algorithm needs to be restarted.
- `return_flag`: Whether the algorithm can return immediately, returning `V`.
- `V`: Either a triangle that can returned if `return_flag = true`, or some triangle that is used for type stability for this return value.
- `p`: The point corresponding to the vertex `k`, or it may be `q` if the algorithm is going to be restarted or `return_flag = true`.
- `i`: The first vertex of the triangle adjoining `k` to start from, or `k` if the algorithm is going to be restarted or `return_flag = true`.
- `j`: The second vertex of the triangle adjoining `k` to start from, or `k` if the algorithm is going to be restarted or `return_flag = true`.
- `pᵢ`: The point corresponding to the vertex `i`, or it may be `q` if the algorithm is going to be restarted or `return_flag = true`.
- `pⱼ`: The point corresponding to the vertex `j`, or it may be `q` if the algorithm is going to be restarted or `return_flag = true`.

# Extended help
There are multiple stages to this initialisation, starting from [`check_for_intersections_with_adjacent_boundary_edges`](@ref). 

- If it is found that `q` is not outside of the triangulation, so that `q` is collinear with one of the boundary edges, then we use [`search_down_adjacent_boundary_edges`](@ref) to find where to start, noting 
   that we can return immediately if `q` is found to be on an adjacent boundary edge. Otherwise, [`exterior_find_triangle`](@ref) can then be used to find the ghost triangle containing 
   `q`; if `concavity_protection = true`, then [`concavity_protection_check`](@ref) is used to determine if a restart is needed.

- If is is not found that `q` is outside of the triangulation yet based on information from the adjacent boundary edges, then we need to check the neighbouring 
   interior edges using [`check_for_intersections_with_interior_edges_adjacent_to_boundary_vertex`](@ref), returning early if `q` is found to be inside one of 
   the neighbouring triangles. If the line `pq`, where `p = get_point(tri, k)`, does not intersect any of the neighbouring edges and it is not inside any of 
   the neighbouring triangles, then it must be outside of the triangulation and so we use [`exterior_find_triangle`](@ref) to find the triangle; as before, [`concavity_protection_check`](@ref)
   is used on the found ghost triangle if needed. If there is an intersection, then we return the triangle containing the intersection point that we can start the algorithm from, 
   and its associated vertices and points. 
"""
function initialise_find_triangle_boundary_vertex(tri::Triangulation, _q, k, store_history::F, history, ghost_vertex, concavity_protection) where {F}
    q = getxy(_q) # type stability in case e.g. a user provides a vector into jump and march
    direction, q_pos, next_vertex, right_cert, left_cert =
        check_for_intersections_with_adjacent_boundary_edges(tri, k, q, ghost_vertex)
    Ttype = triangle_type(tri)
    if !is_outside(direction)
        # q is collinear with one of the edges, so let's jump down these edges and try to find q
        q_pos, u, v, w = search_down_adjacent_boundary_edges(tri, k, q, direction, q_pos, next_vertex, store_history, history, ghost_vertex)
        if is_on(q_pos)
            return false, true, construct_triangle(Ttype, u, v, w), q, k, k, q, q # false is the restart_flag, true is the return_flag. We return q, q, q just to get type stability with all returns
        else
            u, v = exterior_find_triangle(tri, u, q, ghost_vertex)
            w = get_adjacent(tri, u, v)
            V = construct_triangle(Ttype, u, v, w) # Can't just use I(𝒢) here since there could be multiple - just use get_adjacent
            if concavity_protection_check(tri, concavity_protection, V, q)
                return true, false, V, q, k, k, q, q # the extra returns are just for type stability
            else
                return false, true, V, q, k, k, q, q
            end
        end
    end
    # If we did not find anything from the neighbouring boundary edges, we can search the neighbouring interior edges
    i, j, edge_cert, triangle_cert =
        check_for_intersections_with_interior_edges_adjacent_to_boundary_vertex(tri, k, q, right_cert, left_cert, store_history, history, ghost_vertex)
    if is_inside(triangle_cert)
        return false, true, construct_triangle(Ttype, i, j, k), q, i, j, q, q
    elseif is_none(edge_cert)
        u, v = exterior_find_triangle(tri, k, q, ghost_vertex)
        w = get_adjacent(tri, u, v)
        V = construct_triangle(Ttype, u, v, w)
        if concavity_protection_check(tri, concavity_protection, V, q)
            return true, false, V, q, i, j, q, q
        else
            return false, true, V, q, i, j, q, q
        end
    else
        p, pᵢ, pⱼ = get_point(tri, k, i, j)
        return false, false, construct_triangle(Ttype, i, j, k), p, i, j, pᵢ, pⱼ # the triangle is just for type stability
    end
end


"""
    find_triangle_return_on_vertex(tri::Triangulation, q, k, p, pᵢ, pⱼ, i, j) -> (Bool, Bool, Triangle)

Check if `q` is one of the vertices of the triangle `(i, j, k)` and return if needed.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `k`: The vertex `k` that the algorithm started from. 
- `p`: The point corresponding to the vertex `k`.
- `pᵢ`: The point corresponding to the vertex `i`.
- `pⱼ`: The point corresponding to the vertex `j`.
- `i`: The first vertex of the triangle adjoining `k` to start from.
- `j`: The second vertex of the triangle adjoining `k` to start from.

# Outputs
- `restart_flag`: Whether the algorithm needs to be restarted.
- `return_flag`: Whether the algorithm can return immediately, returning `V`.
- `V`: The triangle `(i, j, k)`.

# Extended help 
An extra check is made in this algorithm for the case that the point that `q` is equal to is one of the points corresponding to a ghost vertex, 
so it may be for example that `q == pᵢ` but `is_ghost_vertex(i)`, in which case the algorithm would need to restart.
"""
function find_triangle_return_on_vertex(tri::Triangulation, q, k, p, pᵢ, pⱼ, i, j)
    # Just return where we currently are. We do need to be careful, though: 
    # If k was a ghost vertex, then one of pᵢ or pⱼ could come from the 
    # representative point list, which could mean that q is equal to one of the 
    # vertices, but without meaning that it is actually in that triangle. So, 
    # we need to also check for the type of indices we have. 
    safety_check = (q == p && !is_ghost_vertex(k)) ||
        (q == pᵢ && !is_ghost_vertex(i)) ||
        (q == pⱼ && !is_ghost_vertex(j))
    if safety_check
        orientation = triangle_orientation(pᵢ, pⱼ, p)
        if is_positively_oriented(orientation)
            return false, true, construct_triangle(triangle_type(tri), i, j, k)
        else
            return false, true, construct_triangle(triangle_type(tri), j, i, k)
        end
    else
        return true, false, construct_triangle(triangle_type(tri), i, j, k)
    end
end


"""
    find_triangle_across_triangle(tri::Triangulation, q, k, store_history, history, maxiters, cur_iter, concavity_protection, arrangement, original_k, last_changed, p, i, j, pᵢ, pⱼ) -> (Bool, Bool, Bool, Triangle, Integer, Certificate, Vertex, Vertex, Vertex, Point, Point, Integer, Integer)

Walks across the current triangle past the edge `(i, j)`. progressing the [`find_triangle`](@ref) algorithm.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `k`: The vertex that the algorithm started from.
- `store_history`: Whether to store the history of the algorithm.
- `history`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `maxiters`: The maximum number of iterations to perform before restarting the algorithm with [`restart_find_triangle`](@ref).
- `cur_iter`: The current iteration of the algorithm.
- `concavity_protection`: Whether to use concavity protection. See [`concavity_protection_check`](@ref). This is only needed if your triangulation is not convex.
- `arrangement`: A [`Certificate`](@ref) defining the orientation of the triangle `pᵢpⱼq`.
- `original_k`: The original vertex that the algorithm started from.
- `last_changed`: The last vertex that was changed in the algorithm.
- `p`: The point corresponding to the vertex `original_k`.
- `i`: The first vertex of the triangle adjoining `k` to step from.
- `j`: The second vertex of the triangle adjoining `k` to step from.
- `pᵢ`: The point corresponding to the vertex `i`.
- `pⱼ`: The point corresponding to the vertex `j`.

# Outputs 
- `restart_flag`: Whether the algorithm needs to be restarted.
- `return_flag`: Whether the algorithm can return immediately, returning `V`.
- `reinitialise_flag`: Whether the algorithm needs to be reinitialised at a new vertex `k`. (This would only be needed if `!has_ghost_triangles(tri)`.)
- `V`: The triangle stepped into. 
- `cur_iter`: The new number of iterations of the algorithm.
- `arrangement`: A [`Certificate`](@ref) defining the orientation of the triangle `pᵢpⱼq` with the updated values of `i` and `j`.
- `k`: The new vertex that the algorithm is at.
- `last_changed`: The last vertex that was changed in the algorithm.
- `original_k`: The original vertex that the algorithm started from.
- `pᵢ`: The point corresponding to the vertex `i`.
- `pⱼ`: The point corresponding to the vertex `j`.
- `i`: The first vertex of the triangle adjoining `k` to step from.
- `j`: The second vertex of the triangle adjoining `k` to step from.

# Extended help 
This part of the algorithm is relatively complicated because there are many cases that need to be accounted for. Here we give a brief description of how this step works, 
and note that the documentation contains a much more detailed description.

1. Firstly, we need to check whether `k` is an exterior ghost vertex or not. If `k` is an exterior ghost vertex, then this means that we are stepping outside of the 
   triangulation. Thus, we use [`exterior_find_triangle`](@ref) to find where `q` is, starting from the `last_changed` vertex. If `concavity_protection = true`, then 
   [`concavity_protection_check`](@ref) is used to determine if a restart is needed, or if we can return safely. If we reach this step but `!has_ghost_triangles(tri)`,
   then the algorithm should need to be reinitialised since `q` should not be outside of the triangulation, and so we return with `reinitialise_flag = true`.
2. Now we consider the case where `k` is not an exterior ghost vertex. We move forward by updating the value of `k` so that `k = get_adjacent(tri, i, j)`, and then consider where `pₖ` is relative 
   to the line `pq`.

        2a. If `pₖ` is to the right of `pq`, then we should update `j` by `j = k`, ensuring that `j` is always to the right of `pq`.
        2b. If `pₖ` is to the left of `pq`, then we should update `i` by `i = k`, ensuring that `i` is always to the left of `pq`.
        2c. The alternative to 2a and 2b is that `pₖ` is collinear with the edge of `pq`, which could mean that `q` is in the current triangle or it is in a triangle further away. We compute a 
            [`Certificate`](@ref) that determines where `q` is relative to the triangle `pᵢpⱼpₖ`. If `q` is inside or on this triangle, then we return, restarting if necessary according to 
            `concavity_protection` and [`concavity_protection_check`](@ref). If we do not yet need to return, then we need to make a decision as to which of `i` and `j` to update, noting that 
            we want `i` to be left of `pq` and `j` to be right of `pq`, but this is no longer unambiguous since `pₖ` is collinear with `pq`. We make this decision according to `last_changed`:
            If `last_changed = i`, then moving left is what caused us to find this collinear edge, and so we send `k` left by letting `i = k`. Otherwise, we send `k` right by letting `j = k`.
3. Now having stepped forward, we recompute the [`Certificate`](@ref) for arrangement and return, setting `restart_flag = true` if `cur_iters ≥ maxiters`.
"""
function find_triangle_across_triangle(tri::Triangulation, q, k, store_history::F, history, maxiters, cur_iter, concavity_protection, arrangement, original_k, last_changed, p, i, j, pᵢ, pⱼ) where {F}
    cur_iter += 1
    if is_true(store_history)
        if last_changed == i
            add_left_vertex!(history, i)
        elseif last_changed == j
            add_right_vertex!(history, j)
        end
    end
    # We need to step forward. To do this, we need to be careful of ghost vertices.
    if is_exterior_ghost_vertex(tri, k)
        # If this happens, it means we have hit an outer boundary edge, and so we need to go into the exterior. If there are no 
        # ghost triangles, though, we just restart the search. Note that interior boundary indices do not matter since the ghost 
        # triangles there have the same orientation, so we can find them as normal.
        if has_ghost_triangles(tri)
            i, j = exterior_find_triangle(tri, last_changed == i ? j : i, q, k) # use last_changed to ensure we get the boundary point
            V = construct_triangle(triangle_type(tri), i, j, get_adjacent(tri, i, j))
            if concavity_protection_check(tri, concavity_protection, V, q)
                return true, false, false, V, cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j # restart_flag, return_flag, reinitialise_flag, V
            else
                return false, true, false, V, cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j
            end
        else
            return false, false, true, construct_triangle(triangle_type(tri), last_changed, last_changed, last_changed), cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j # just returning a triangle for stability
        end
    end
    # Now we can move forward.
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
            V = construct_triangle(triangle_type(tri), i, j, k)
            if concavity_protection_check(tri, concavity_protection, V, q)
                return true, false, false, V, cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j # restart_flag, return_flag, reinitialise_flag, V
            else
                # return V
                return false, true, false, V, cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j
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
    return cur_iter ≥ maxiters, false, false, construct_triangle(triangle_type(tri), i, j, k), cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j
end


"""
    find_triangle_degenerate_arrangement(tri::Triangulation, q, k, store_history::F, history, pᵢ, pⱼ, i, j) -> Bool

Given a degenerate arrangement `pᵢpⱼq`, reinitialise the jump and march algorithm.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `k`: The vertex that the algorithm started from.
- `store_history`: Whether to store the history of the algorithm.
- `history`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `pᵢ`: The point corresponding to the vertex `i`.
- `pⱼ`: The point corresponding to the vertex `j`.
- `i`: The first vertex of the triangle adjoining `k` to step from.
- `j`: The second vertex of the triangle adjoining `k` to step from.

# Outputs
- `Bool`: Whether the algorithm needs to be restarted.
"""
function find_triangle_degenerate_arrangement(tri::Triangulation, q, k, store_history::F, history, pᵢ, pⱼ, i, j) where {F}
    pₖ = get_point(tri, k) # Need to have this here in case we skip the entire loop, above, meaning pₖ won't exist
    in_cert = point_position_relative_to_triangle(pᵢ, pⱼ, pₖ, q) # ... Maybe there is a better way to compute this, reusing previous certificates? Not sure. ...
    if is_true(store_history)
        k′ = get_adjacent(tri, i, j)
        add_triangle!(history, i, j, k′)
        add_left_vertex!(history, i)
        add_right_vertex!(history, j)
    end
    return is_outside(in_cert) # if outside, we need to reinitialise. Otherwise, we're all good.
end


function _find_triangle(tri::Triangulation, q, k, store_history::F, history, rng::AbstractRNG, maxiters, cur_iter, concavity_protection, num_restarts, use_barriers::Val{U}) where {F, U}
    is_bnd, ghost_vertex = is_boundary_node(tri, k)
    I = integer_type(tri)
    trit = triangle_type(tri)
    if !(is_bnd && is_exterior_ghost_vertex(tri, ghost_vertex)) || !has_ghost_triangles(tri)
        restart_flag, p, i, j, pᵢ, pⱼ = initialise_find_triangle_interior_vertex(tri, q, k, store_history, history, rng)
        if restart_flag && (!(is_bnd && is_interior_ghost_vertex(tri, ghost_vertex)) || !is_true(use_barriers)) # if we are not at an interior boundary node, then we should not have reached a barrier yet. but if we are at such a node, then the only reason to restart is if we have reached a barrier.
            return restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        elseif restart_flag && is_true(use_barriers)
            V = construct_triangle(trit, I(∅), I(∅), I(∅))
            return V, true
        end
    else
        restart_flag, return_flag, V, p, i, j, pᵢ, pⱼ = initialise_find_triangle_boundary_vertex(tri, q, k, store_history, history, ghost_vertex, concavity_protection)
        if restart_flag && !is_true(use_barriers)
            return restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        elseif restart_flag && is_true(use_barriers)
            return V, true
        end
        return_flag && return is_true(use_barriers) ? (V, is_ghost_triangle(V)) : V
    end
    if q == p || q == pᵢ || q == pⱼ
        restart_flag, return_flag, V = find_triangle_return_on_vertex(tri, q, k, p, pᵢ, pⱼ, i, j)
        if restart_flag && !is_true(use_barriers)
            return restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        elseif restart_flag && is_true(use_barriers)
            return V, true
        end
        return_flag && return is_true(use_barriers) ? (V, is_ghost_triangle(V)) : V
    end
    if is_true(use_barriers) && any(is_ghost_vertex, (i, j))
        V = construct_triangle(trit, I(∅), I(∅), I(∅))
        return V, true
    end
    ## Now let us do the straight line search 
    # The idea is to keep jumping until pᵢpⱼ goes past q, meaning pᵢpⱼq is no longer a positively oriented triangle
    original_k = k
    arrangement = triangle_orientation(pᵢ, pⱼ, q)
    local last_changed # Need this for deciding which variable to use when we hit a collinear point 
    I = integer_type(tri)
    last_changed = I(∅) # just an initial value
    reached_barrier = is_true(use_barriers) && is_positively_oriented(arrangement) && contains_segment(tri, i, j)
    if is_true(store_history)
        add_left_vertex!(history, i)
        add_right_vertex!(history, j)
    end
    while is_positively_oriented(arrangement) && !reached_barrier
        restart_flag, return_flag, reinitialise_flag, V, cur_iter, arrangement, k, last_changed, original_k, pᵢ, pⱼ, i, j = find_triangle_across_triangle(tri, q, k, store_history, history, maxiters, cur_iter, concavity_protection, arrangement, original_k, last_changed, p, i, j, pᵢ, pⱼ)
        if restart_flag && !is_true(use_barriers)
            return restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        elseif restart_flag && is_true(use_barriers)
            return V, true
        end
        return_flag && return is_true(use_barriers) ? (V, reached_barrier) : V
        reinitialise_flag && return _find_triangle(tri, q, k, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        reached_barrier = is_true(use_barriers) && is_positively_oriented(arrangement) && contains_segment(tri, i, j) # check the arrangement cert since, if it is now negative, then we don't care about the next edge (i, j) because we have reached the triangle
    end
    # We can finish the above loop even if q is not in the triangle, in which case pᵢpⱼq was a straight line. 
    # To clear this up, let us just restart. 
    if is_degenerate(arrangement)
        reinitialise_flag = find_triangle_degenerate_arrangement(tri, q, k, store_history, history, pᵢ, pⱼ, i, j)
        if reinitialise_flag && !is_true(use_barriers)
            return _find_triangle(tri, q, last_changed == I(∅) ? i : last_changed, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
        elseif reinitialise_flag && is_true(use_barriers)
            return V, true
        end
    end
    # Swap the orientation to get a positively oriented triangle, remembering that we kept pᵢ on the left of pq and pⱼ on the right 
    k = get_adjacent(tri, j, i)
    V = construct_triangle(trit, j, i, k)
    if !reached_barrier && concavity_protection_check(tri, concavity_protection, V, q)
        return restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts + 1, use_barriers)
    end
    if is_true(use_barriers)
        return V, reached_barrier
    else
        return V
    end
end


"""
    concavity_protection_check(tri::Triangulation, concavity_protection, V, q) -> Bool

Check whether the [`find_triangle`](@ref) algorithm needs to restart. This is only needed if `tri` is not convex.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `concavity_protection`: Whether this check is needed.
- `V`: The triangle that the algorithm has found.
- `q`: The query point.

# Outputs 
- `need_to_restart`: Whether the algorithm needs to restart. Will also be `false` if `concavity_protection`.

# Extended help 
This function uses [`dist`](@ref) to determine whether the query point `q` is inside or outside of the polygon defined by the triangulation,
and also checks the position of `q` relative to `V` via [`point_position_relative_to_triangle`](@ref). If `q` is outside of this triangle, 
then `need_to_restart = true`. If `q` is inside this triangle, then issues can still arise due to overlapping ghost triangles from the non-convexity. Thus, 
depending on the result from [`dist`](@ref) and whether `V` is a ghost triangle, `need_to_restart` will be set accordingly.
"""
function concavity_protection_check(tri::Triangulation, concavity_protection, V, q)
    !concavity_protection && return false
    cert = point_position_relative_to_triangle(tri, V, q)
    is_outside(cert) && return true
    δ = dist(tri, q)
    is_ghost = is_ghost_triangle(V)
    need_to_restart = (is_ghost && δ > 0.0) || (!is_ghost && δ < 0.0)
    return need_to_restart
end


const RESTART_LIMIT = 25
"""
    restart_find_triangle(tri::Triangulation, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts, use_barriers) -> Triangle[, Bool]

Restart the [`find_triangle`](@ref) algorithm, or use [`brute_force_search`](@ref) to find `q` if `num_restarts ≥ $RESTART_LIMIT`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `q`: The query point.
- `store_history`: Whether to store the history of the algorithm.
- `history`: The history of the algorithm. If `store_history`, then this should be a [`PointLocationHistory`](@ref) object.
- `rng`: The random number generator to use.
- `maxiters`: The maximum number of iterations to perform before restarting the algorithm with [`restart_find_triangle`](@ref).
- `cur_iter`: The current iteration of the algorithm.
- `concavity_protection`: Whether to use concavity protection. See [`concavity_protection_check`](@ref). This is only needed if your triangulation is not convex.
- `num_restarts`: The number of times the algorithm has been restarted. 
- `use_barriers`: Whether to use barriers, stopping the algorithm at any segment.

# Outputs
- `V`: The triangle containing `q`.

In addition, if `use_barriers = Val(true)`, then a second output is returned, which is a boolean indicating whether the algorithm reached a barrier (`true`) or not (`false`).
"""
function restart_find_triangle(tri, q, store_history, history, rng, maxiters, cur_iter, concavity_protection, num_restarts, use_barriers)
    if num_restarts < RESTART_LIMIT
        m = num_solid_vertices(tri)
        point_indices = each_solid_vertex(tri)
        k = select_initial_point(tri, q; m = (m >> 1) + 1, point_indices, rng) # don't want to try all points, still want to give the algorithm a chance
        return _find_triangle(tri, q, k, store_history, history, rng, maxiters, zero(cur_iter), concavity_protection, num_restarts, use_barriers)
    else
        V = brute_force_search(tri, q)
        V_is_bad = concavity_protection_check(tri, concavity_protection, V, q)
        if V_is_bad
            if is_ghost_triangle(V)
                V = brute_force_search(tri, q; itr = each_solid_triangle(tri))
            else
                V = brute_force_search(tri, q; itr = each_ghost_triangle(tri))
            end
        end
        if is_true(use_barriers)
            return V, false
        else
            return V
        end
    end
end