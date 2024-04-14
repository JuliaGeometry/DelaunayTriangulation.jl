"""
    InsertionEventHistory{T,E}

A data structure for storing the changes to the triangulation during the insertion of a point.

# Fields
- `added_triangles::Set{T}`: The triangles that were added.
- `deleted_triangles::Set{T}`: The triangles that were deleted.
- `added_segments::Set{E}`: The interior segments that were added.
- `deleted_segments::Set{E}`: The interior segments that were deleted.
- `added_boundary_segments::Set{E}`: The boundary segments that were added.
- `deleted_boundary_segments::Set{E}`: The boundary segments that were deleted.

# Constructor 
The default constructor is available, but we also provide 

    InsertionEventHistory(tri::Triangulation)

which will initialise this struct with empty, appropriately `sizehint!`ed, sets.
"""
struct InsertionEventHistory{T,E}
    added_triangles::Set{T}
    deleted_triangles::Set{T}
    added_segments::Set{E}
    deleted_segments::Set{E}
    added_boundary_segments::Set{E}
    deleted_boundary_segments::Set{E}
end
function Base.show(io::IO, ::MIME"text/plain", events::InsertionEventHistory)
    println(io, "InsertionEventHistory")
    println(io, "   $(length(events.added_triangles)) added_triangles: $(events.added_triangles)")
    println(io, "   $(length(events.deleted_triangles)) deleted_triangles: $(events.deleted_triangles)")
    println(io, "   $(length(events.added_segments)) added_segments: $(events.added_segments)")
    println(io, "   $(length(events.deleted_segments)) deleted_segments: $(events.deleted_segments)")
    println(io, "   $(length(events.added_boundary_segments)) added_boundary_segments: $(events.added_boundary_segments)")
    print(io, "   $(length(events.deleted_boundary_segments)) deleted_boundary_segments: $(events.deleted_boundary_segments)")
end

"""
    InsertionEventHistory(tri::Triangulation) -> InsertionEventHistory

Initialises an `InsertionEventHistory` for the triangulation `tri`.
"""
function InsertionEventHistory(tri::Triangulation)
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
    return InsertionEventHistory{T,E}(add_set, delete_set, add_edge_set, delete_edge_set, add_bnd_set, delete_bnd_set)
end

"""
    add_triangle!(events::InsertionEventHistory, T)

Add the triangle `T` to the `added_triangles` of `events`.
"""
add_triangle!(events::InsertionEventHistory, T) = push!(events.added_triangles, T)

"""
    delete_triangle!(events::InsertionEventHistory, T)

Add the triangle `T` to the `deleted_triangles` of `events`.
"""
delete_triangle!(events::InsertionEventHistory, T) = push!(events.deleted_triangles, T)

"""
    add_edge!(events::InsertionEventHistory, e)

Add the edge `e` to the `added_segments` of `events`.
"""
add_edge!(events::InsertionEventHistory, e) = push!(events.added_segments, e)

"""
    delete_edge!(events::InsertionEventHistory, e)

Add the edge `e` to the `deleted_segments` of `events`.
"""
delete_edge!(events::InsertionEventHistory, e) = push!(events.deleted_segments, e)

"""
    split_boundary_edge!(events::InsertionEventHistory, u, v, new_point)

Add the edge `(u, v)` to the `deleted_boundary_segments` of `events` and add the edges `(u, new_point)` and `(new_point, v)` to the `added_boundary_segments` of `events`.
"""
function split_boundary_edge!(events::InsertionEventHistory{T,E}, u, v, new_point) where {T,E}
    !contains_edge(construct_edge(E, v, u), events.deleted_boundary_segments) && push!(events.deleted_boundary_segments, construct_edge(E, u, v))
    !contains_edge(construct_edge(E, new_point, u), events.added_boundary_segments) && push!(events.added_boundary_segments, construct_edge(E, u, new_point))
    !contains_edge(construct_edge(E, v, new_point), events.added_boundary_segments) && push!(events.added_boundary_segments, construct_edge(E, new_point, v))
    return events
end

"""
    has_triangle_changes(events::InsertionEventHistory) -> Bool 

Returns `true` if there are any changes to the segments in `events`, and `false` otherwise.
"""
function has_segment_changes(events::InsertionEventHistory)
    return any(!isempty, (events.added_segments, events.deleted_segments,
        events.added_boundary_segments, events.deleted_boundary_segments))
end

"""
    each_added_triangle(events::InsertionEventHistory) -> Iterator

Returns an iterator over the triangles that were added to the triangulation during the insertion of a point.
"""
each_added_triangle(events::InsertionEventHistory) = each_triangle(events.added_triangles)

"""
    each_deleted_triangle(events::InsertionEventHistory) -> Iterator

Returns an iterator over the triangles that were deleted from the triangulation during the insertion of a point.
"""
each_added_segment(events::InsertionEventHistory) = each_edge(events.added_segments)

"""
    each_added_boundary_segment(events::InsertionEventHistory) -> Iterator

Returns an iterator over the boundary segments that were added to the triangulation during the insertion of a point.
"""
each_added_boundary_segment(events::InsertionEventHistory) = each_edge(events.added_boundary_segments)

"""
    triangle_type(::InsertionEventHistory{T}) where {T} = T

Returns the type of the triangles in `events`, `T`.
"""
triangle_type(::InsertionEventHistory{T}) where {T} = T

"""
    empty!(events::InsertionEventHistory)

Empties `events` by emptying all of its fields.
"""
function Base.empty!(events::InsertionEventHistory)
    empty!(events.added_triangles)
    empty!(events.deleted_triangles)
    empty!(events.added_segments)
    empty!(events.deleted_segments)
    empty!(events.added_boundary_segments)
    empty!(events.deleted_boundary_segments)
    return events
end

"""
    undo_insertion!(tri::Triangulation, events::InsertionEventHistory, pop=Val(true))

Undoes the insertion of the most recent vertex into `tri`, assuming that its insertion history has been 
recorded into `events` and the vertex is `num_points(tri)`. 

If you do not want to delete the latest vertex from the triangulation, set `pop` to `Val(false)`.
"""
function undo_insertion!(tri::Triangulation, events::InsertionEventHistory, pop=Val(true))
    vertex = num_points(tri)
    for T in events.added_triangles
        delete_triangle!(tri, T; protect_boundary=true, update_ghost_edges=false)
    end
    for T in events.deleted_triangles
        add_triangle!(tri, T; protect_boundary=true, update_ghost_edges=false)
    end
    undo_segment_changes!(tri, events)
    undo_boundary_segment_changes!(tri, events)
    if is_true(pop)
        delete_adjacent2vertex!(tri, vertex)
        delete_vertex!(tri, vertex)
        pop_point!(tri)
    end
    return tri
end

"""
    undo_segment_changes!(tri::Triangulation, events::InsertionEventHistory)

Undoes any changes to the segments in `tri` that were made after an insertion event, as recorded in `events`.
"""
function undo_segment_changes!(tri::Triangulation, events::InsertionEventHistory)
    interior_segments = get_interior_segments(tri)
    all_segments = get_all_segments(tri)
    E = edge_type(tri)
    for edges in (interior_segments, all_segments)
        for added_edges in events.added_segments
            u, v = edge_vertices(added_edges)
            delete_unoriented_edge!(edges, construct_edge(E, u, v))
        end
        for deleted_edges in events.deleted_segments
            u, v = edge_vertices(deleted_edges)
            e = construct_edge(E, u, v)
            e′ = reverse_edge(e)
            if !contains_edge(e′, edges)
                add_edge!(edges, construct_edge(E, u, v))
            end
        end
    end
    return tri
end

"""
    undo_boundary_segment_changes!(tri::Triangulation, events::InsertionEventHistory)

Undoes any changes to the boundary segments in `tri` that were made after an insertion event, as recorded in `events`,
assuming the inserted vertex is `num_points(tri)`.
"""
function undo_boundary_segment_changes!(tri::Triangulation, events::InsertionEventHistory)
    vertex = num_points(tri)
    # Only one boundary edge will ever be changed. So, just extract it. 
    deleted_boundary_segments = events.deleted_boundary_segments
    isempty(deleted_boundary_segments) && return tri
    e = pop!(deleted_boundary_segments)
    merge_boundary_edge!(tri, e, vertex)
    return tri
end