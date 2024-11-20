struct InsertionEventHistory{I, E}
    added_triangles::Set{NTuple{3,I}}
    deleted_triangles::Set{NTuple{3,I}}
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

function InsertionEventHistory(tri::Triangulation)
    I = integer_type(tri)
    E = edge_type(tri)
    add_set = Set{NTuple{3,I}}()
    delete_set = Set{NTuple{3,I}}()
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
    return InsertionEventHistory{I, E}(add_set, delete_set, add_edge_set, delete_edge_set, add_bnd_set, delete_bnd_set)
end

@inline add_triangle!(events::InsertionEventHistory, T) = push!(events.added_triangles, T)
@inline delete_triangle!(events::InsertionEventHistory, T) = push!(events.deleted_triangles, T)
@inline add_edge!(events::InsertionEventHistory, e) = push!(events.added_segments, e)
@inline delete_edge!(events::InsertionEventHistory, e) = push!(events.deleted_segments, e)
@inline function split_boundary_edge!(events::InsertionEventHistory{T, E}, u, v, new_point) where {T, E}
    !contains_edge(construct_edge(E, v, u), events.deleted_boundary_segments) && push!(events.deleted_boundary_segments, construct_edge(E, u, v))
    !contains_edge(construct_edge(E, new_point, u), events.added_boundary_segments) && push!(events.added_boundary_segments, construct_edge(E, u, new_point))
    !contains_edge(construct_edge(E, v, new_point), events.added_boundary_segments) && push!(events.added_boundary_segments, construct_edge(E, new_point, v))
    return events
end
@inline function has_segment_changes(events::InsertionEventHistory)
    return any(
        !isempty, (
            events.added_segments, events.deleted_segments,
            events.added_boundary_segments, events.deleted_boundary_segments,
        ),
    )
end
@inline each_added_triangle(events::InsertionEventHistory) = each_triangle(events.added_triangles)
@inline each_added_segment(events::InsertionEventHistory) = each_edge(events.added_segments)
@inline each_added_boundary_segment(events::InsertionEventHistory) = each_edge(events.added_boundary_segments)
@inline integer_type(::InsertionEventHistory{I}) where {I} = I

@inline function Base.empty!(events::InsertionEventHistory)
    empty!(events.added_triangles)
    empty!(events.deleted_triangles)
    empty!(events.added_segments)
    empty!(events.deleted_segments)
    empty!(events.added_boundary_segments)
    empty!(events.deleted_boundary_segments)
    return events
end


function undo_insertion!(tri::Triangulation, events::InsertionEventHistory, pop = Val(true))
    vertex = num_points(tri)
    for T in events.added_triangles
        delete_triangle!(tri, T)
    end
    for T in events.deleted_triangles
        add_triangle!(tri, T)
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

function undo_boundary_segment_changes!(tri::Triangulation, events::InsertionEventHistory)
    vertex = num_points(tri)
    # Only one boundary edge will ever be changed. So, just extract it. 
    deleted_boundary_segments = events.deleted_boundary_segments
    isempty(deleted_boundary_segments) && return tri
    e = pop!(deleted_boundary_segments)
    merge_boundary_edge!(tri, e, vertex)
    return tri
end
