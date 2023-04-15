"""
    try_circumcenter_insertion!(tri::Triangulation, T, events::InsertionEventHistory, queue::RefinementQueue, exterior_curve_index=1, insertion_strategy=:off, rng::AbstractRNG=Random.default_rng())

Attempt to insert a new point into the triangulation at the circumcenter of the triangle `T`. If the insertion causes edges to be encroached, 
the insertion is undone and the encroached edges are queued for refinement. Otherwise, the insertion is accepted and the new point is added to the 
triangulation.

If `insertion_strategy` is `:off`, then we insert the off-centre rather than the circumcenter, but if `:on` then we use the circumcenter.

The returned value is a [`Certificate`](@ref) type, used for distinguishing between three possible results:

- `Cert.PrecisionFailure`: The circumcenter wasn't inserted due to a failure relating to precision (e.g. the circumcenter was another vertex).
- `Cert.EncroachmentFailure`: The circumcenter wasn't inserted due to encroachment.
- `Cert.SuccessfulInsertion`: The circumcenter was inserted successfully.

We need to have `Cert.PrecisionFailure` be separate to avoid requeueing the triangle later.
"""
function try_circumcenter_insertion!(tri::Triangulation, T, events::InsertionEventHistory, queue::RefinementQueue, exterior_curve_index=1, insertion_strategy=:off, rng::AbstractRNG=Random.default_rng())
    empty!(events)
    ε = sqrt(eps(number_type(tri)))
    ## Find the circumcenter
    i, j, k = indices(T)
    p, q, r = get_point(tri, i, j, k)
    A² = squared_triangle_area(p, q, r)
    A² ≤ ε^2 && return Cert.PrecisionFailure

    cx, cy = triangle_circumcenter(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    # Just as we need when checking for precision issues when splitting a segment, we need to check that the circumcenter is not one of the vertices.
    if ((px == cx) && (py == cy)) || ((qx == cx) && (qy == cy)) || ((rx == cx) && (ry == cy)) || # equality is actually valid to check here, sometimes the precison really does push them to be equal. We do it separately as it's faster.
       (!iszero(cx) && !iszero(cy) && ((abs((px - cx) / cx) ≤ ε) && (abs((py - cy) / cy) ≤ ε)) || ((abs((qx - cx) / cx) ≤ ε) && (abs((qy - cy) / cy) ≤ ε)) || ((abs((rx - cx) / cx) ≤ ε) && (abs((ry - cy) / cy) ≤ ε))) ||
       ((abs(px - cx) ≤ ε) && (abs(py - cy) ≤ ε)) || ((abs(qx - cx) ≤ ε) && (abs(qy - cy) ≤ ε)) || ((abs(rx - cx) ≤ ε) && (abs(ry - cy) ≤ ε))
        return Cert.PrecisionFailure
    end

    ## Where should we start the point location?
    δpc² = (cx - px)^2 + (cy - py)^2
    δqc² = (cx - qx)^2 + (cy - qy)^2
    δrc² = (cx - rx)^2 + (cy - ry)^2
    min_δ = min(δpc², δqc², δrc²)
    if min_δ == δpc²
        point_location_initial_vertex = i
    elseif min_δ == δqc²
        point_location_initial_vertex = j
    else
        point_location_initial_vertex = k
    end

    ## Find the triangle where the vertex is located 
    # First, check if it's in the current triangle
    q = (cx, cy)
    flag = point_position_relative_to_triangle(tri, T, q)
    if !is_outside(flag)
        V = T
    else
        # If not, we need to search for it
        V = jump_and_march(
            tri,
            q;
            m=nothing,
            point_indices=nothing,
            try_points=nothing,
            k=point_location_initial_vertex,
            rng,
            check_existence=Val(has_multiple_segments(tri)),
            exterior_curve_index
        )
    end

    # Check if there are any precision issues
    if is_ghost_triangle(V)
        return Cert.PrecisionFailure
    end
    v1, v2, v3 = indices(V)
    p1, p2, p3 = get_point(tri, v1, v2, v3)
    p1x, p1y = getxy(p1)
    p2x, p2y = getxy(p2)
    p3x, p3y = getxy(p3)
    if ((p1x == cx) && (p1y == cy)) || ((p2x == cx) && (p2y == cy)) || ((p3x == cx) && (p3y == cy)) ||
       (!iszero(cx) && !iszero(cy) && ((abs((p1x - cx) / cx) ≤ ε) && (abs((p1y - cy) / cy) ≤ ε)) || ((abs((p2x - cx) / cx) ≤ ε) && (abs((p2y - cy) / cy) ≤ ε)) || ((abs((p3x - cx) / cx) ≤ ε) && (abs((p3y - cy) / cy) ≤ 1e-14))) ||
       ((abs(p1x - cx) ≤ ε) && (abs(p1y - cy) ≤ ε)) || ((abs(p2x - cx) ≤ ε) && (abs(p2y - cy) ≤ ε)) || ((abs(p3x - cx) ≤ ε) && (abs(p3y - cy) ≤ ε))
        return Cert.PrecisionFailure
    end

    ## Now insert the vertex
    V = add_point!(tri, cx, cy;
        point_indices=nothing,
        m=nothing,
        try_points=nothing,
        rng,
        initial_search_point=nothing,
        update_representative_point=false,
        store_event_history=Val(true),
        event_history=events,
        exterior_curve_index=exterior_curve_index,
        V=V
    )

    ## Now we need to search the edges opposite to the circumcenter to see if we have introduced any new encroached edges 
    any_encroached = false
    circumcenter_index = num_points(tri)
    for edge in get_adjacent2vertex(tri, circumcenter_index)
        if contains_constrained_edge(tri, edge) && is_encroached(tri, edge)
            any_encroached = true
            u, v = edge_indices(edge)
            p, q = get_point(tri, u, v)
            px, py = getxy(p)
            qx, qy = getxy(q)
            ℓ² = (qx - px)^2 + (qy - py)^2
            encroachment_enqueue!(queue, edge, ℓ²)
        end
    end

    ## If we found an encroached edge, we need to completely reverse the insertion! Thankfully, our event storage can be used.
    any_encroached && undo_circumcenter_insertion!(tri, events, circumcenter_index)
    if any_encroached
        return Cert.EncroachmentFailure
    else
        return Cert.SuccessfulInsertion
    end
end

function undo_circumcenter_insertion!(tri::Triangulation, events, circumcenter_index)
    for T in events.added_triangles
        delete_triangle!(tri, T; protect_boundary=true, update_ghost_edges=false)
    end
    for T in events.deleted_triangles
        add_triangle!(tri, T; protect_boundary=true, update_ghost_edges=false)
    end
    undo_constrained_segment_changes!(tri, events)
    undo_boundary_segment_changes!(tri, events, circumcenter_index)
    delete_adjacent2vertex!(tri, circumcenter_index)
    delete_vertex!(tri, circumcenter_index)
    pop_point!(tri)
    return nothing
end

function undo_constrained_segment_changes!(tri::Triangulation, events)
    constrained_edges = get_constrained_edges(tri)
    all_constrained_edges = get_all_constrained_edges(tri)
    E = edge_type(tri)
    for edges in (constrained_edges, all_constrained_edges)
        for added_edges in events.added_segments
            u, v = edge_indices(added_edges)
            delete_edge!(edges, construct_edge(E, u, v))
            delete_edge!(edges, construct_edge(E, v, u))
        end
        for deleted_edges in events.deleted_segments
            u, v = edge_indices(deleted_edges)
            e = construct_edge(E, u, v)
            e′ = reverse_edge(e)
            if !contains_edge(e′, edges)
                add_edge!(edges, construct_edge(E, u, v))
            end
        end
    end
    return nothing
end

function undo_boundary_segment_changes!(tri::Triangulation, events, circumcenter_index)
    # Only one boundary edge will ever be changed. So, just extract it. 
    deleted_boundary_segments = events.deleted_boundary_segments
    isempty(deleted_boundary_segments) && return nothing
    e = pop!(deleted_boundary_segments)
    merge_boundary_edge!(tri, e, circumcenter_index)
    return nothing
end

function compute_split_position(tri::Triangulation, p, q, e, segment_list)
    if e ∈ segment_list
        # Split at the midpoint when splitting for the first time
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = px + (qx - px) / 2, py + (qy - py) / 2
        return mx, my
    else
        num_adjoin = segment_vertices_adjoin_other_segments(tri, e)
        if num_adjoin == 2
            # If both vertices adjoin another segment, then the segments can undergo up to two unbalanced splits 
            # on each other. To overcome this, we instead split the segment to be between 1/4 and 1/2 of the split subsegment.
            t = compute_concentric_shell_quarternary_split_position(p, q)
        elseif num_adjoin == 1
            # If we just have one other adjoining segment, we can split the segment to be between 1/3 and 2/3 of the split subsegment.
            t = compute_concentric_shell_ternary_split_position(p, q)
        else
            # If there are no other adjoining segments, there's no need to worry about ping-pong encroachment.
            t = inv(2.0)
        end
        # Need to place the split on the side that is closer to the midpoint
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = px + (qx - px) / 2, py + (qy - py) / 2
        if abs(t - 1 / 2) < 1e-6 # close enough to the midpoint, so just bisect - every third split on a segment should be a bisection anyway. 
            return mx, my
        else
            t′ = 1.0 - t
            ax, ay = px + t * (qx - px), py + t * (qy - py)
            bx, by = px + t′ * (qx - px), py + t′ * (qy - py)
            if (ax - mx)^2 + (ay - my)^2 > (bx - mx)^2 + (by - my)^2
                return bx, by
            else
                return ax, ay
            end
        end
    end
end

function split_subsegment!(tri::Triangulation, queue, events, targets, segment_list, e)
    empty!(events)
    ε = eps(number_type(tri))
    u, v = edge_indices(e)
    p, q = get_point(tri, u, v)
    mx, my = compute_split_position(tri, p, q, e, segment_list)
    # We can run into precision issues (e.g. I found this was needed when reproducing Figure 6.14 in the Delaunay 
    # Mesh Generation book). Let's ensure that the computed splitting point isn't just one of p or q.
    px, py = getxy(p)
    qx, qy = getxy(q)
    if ((mx == px) && (my == py)) || ((mx == qx) && (my == qy)) ||
       (abs(mx - px) ≤ ε) && (abs(my - py) ≤ ε) || (abs(mx - qx) ≤ ε) && (abs(my - qy) ≤ ε)
        return nothing
    end
    push_point!(tri, mx, my)
    r = num_points(tri)
    complete_split_edge_and_legalise!(tri, u, v, r, Val(true), events)
    assess_added_triangles!(tri, queue, events, targets)
    return nothing
end

function split_all_encroached_segments!(tri::Triangulation, queue, events, targets, segment_list)
    iters = 0
    while !encroachment_queue_is_empty(queue) && !compare_points(targets, num_points(tri))
        iters += 1
        e = encroachment_dequeue!(queue)
        if !edge_exists(tri, e) && !edge_exists(tri, reverse_edge(e))
            continue
        end
        split_subsegment!(tri, queue, events, targets, segment_list, e)
    end
    return nothing
end

function split_triangle!(tri::Triangulation, queue::RefinementQueue, events::InsertionEventHistory, T, exterior_curve_index, insertion_strategy, rng::AbstractRNG=Random.default_rng())
    empty!(events)
    success_certificate = try_circumcenter_insertion!(tri, T, events, queue, exterior_curve_index, insertion_strategy, rng)
    return success_certificate
end