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
each_added_triangle(events::InsertionEventHistory) = each_triangle(events.added_triangles)
each_added_segment(events::InsertionEventHistory) = each_edge(events.added_segments)
each_added_boundary_segment(events::InsertionEventHistory) = each_edge(events.added_boundary_segments)

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