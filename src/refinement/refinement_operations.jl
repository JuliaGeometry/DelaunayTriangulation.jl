"""
    InsertionEventHistory{T, E}

A struct to store the events that occur during the insertion of a new point into a triangulation.

# Fields 
- `added_triangles::Set{T}`

A set of triangles that are added to the triangulation after the insertion.
- `deleted_triangles::Set{T}`

A set of triangles that are deleted from the triangulation after the insertion.
- `added_segments::Set{E}`

A set of segments that are added to the triangulation after the insertion.
- `deleted_segments::Set{E}`

A set of segments that are deleted from the triangulation after the insertion.
- `added_boundary_segments::Set{E}`

A set of boundary segments that are added to the triangulation after the insertion.
- `deleted_boundary_segments::Set{E}`

A set of boundary segments that are deleted from the triangulation after the insertion.
"""
struct InsertionEventHistory{T,E}
    added_triangles::Set{T}
    deleted_triangles::Set{T}
    added_segments::Set{E}
    deleted_segments::Set{E}
    added_boundary_segments::Set{E}
    deleted_boundary_segments::Set{E}
end
add_triangle!(events::InsertionEventHistory, T) = push!(events.added_triangles, T)
delete_triangle!(events::InsertionEventHistory, T) = push!(events.deleted_triangles, T)
add_edge!(events::InsertionEventHistory, e) = push!(events.added_segments, e)
delete_edge!(events::InsertionEventHistory, e) = push!(events.deleted_segments, e)
function split_boundary_edge!(events::InsertionEventHistory{T,E}, u, v, new_point) where {T,E}
    push!(events.deleted_boundary_segments, construct_edge(E, u, v))
    push!(events.added_boundary_segments, construct_edge(E, u, new_point))
    push!(events.added_boundary_segments, construct_edge(E, new_point, v))
    return nothing
end
function has_segment_changes(events::InsertionEventHistory)
    return any(!isempty, (events.added_segments, events.deleted_segments,
        events.added_boundary_segments, events.deleted_boundary_segments))
end

function initialise_event_history(tri::Triangulation)
    T = triangle_type(tri)
    E = edge_type(tri)
    add_set = Set{T}()
    delete_set = Set{T}()
    add_edge_set = Set{E}()
    delete_edge_set = Set{E}()
    add_bnd_set = Set{E}()
    delete_bnd_set = Set{E}()
    sizehint!(add_set, 16)
    sizehint!(delete_set, 16)
    sizehint!(add_edge_set, 8)
    sizehint!(delete_edge_set, 8)
    sizehint!(add_bnd_set, 8)
    sizehint!(delete_bnd_set, 8)
    return InsertionEventHistory{T,E}(Set{T}(), Set{T}(), Set{E}(), Set{E}(), Set{E}(), Set{E}())
end
function Base.empty!(events::InsertionEventHistory)
    empty!(events.added_triangles)
    empty!(events.deleted_triangles)
    empty!(events.added_segments)
    empty!(events.deleted_segments)
    empty!(events.added_boundary_segments)
    empty!(events.deleted_boundary_segments)
    return nothing
end

"""
    try_circumcenter_insertion!(tri::Triangulation, T, events::InsertionEventHistory, queue::RefinementQueue, rng::AbstractRNG=Random.default_rng())

Attempt to insert a new point into the triangulation at the circumcenter of the triangle `T`. If the insertion causes edges to be encroached, 
the insertion is undone and the encroached edges are queued for refinement. Otherwise, the insertion is accepted and the new point is added to the 
triangulation.
"""
function try_circumcenter_insertion!(tri::Triangulation, T, events::InsertionEventHistory, queue::RefinementQueue, rng::AbstractRNG=Random.default_rng())
    empty!(events)
    ## Find the circumcenter
    i, j, k = indices(T)
    p, q, r = get_point(tri, i, j, k)
    cx, cy = triangle_circumcenter(p, q, r)

    ## Where should we start the point location?
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
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

    ## Now insert the vertex
    add_point!(tri, cx, cy;
        point_indices=nothing,
        m=nothing,
        try_points=nothing,
        rng,
        initial_search_point=point_location_initial_vertex,
        update_representative_point=false,
        store_event_history=Val(true),
        event_history=events
    )

    ## Now we need to search the edges opposite to the circumcenter to see if we have introduced any new encroached edges 
    any_encroached = false
    circumcenter_index = num_points(tri)
    for edge in get_adjacent2vertex(tri, circumcenter_index)
        if is_encroached(tri, edge)
            @show edge
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
    return any_encroached
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