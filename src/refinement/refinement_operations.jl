"""
    InsertionEventStorage{T, E}

A struct to store the events that occur during the insertion of a new point into a triangulation.

# Fields 
- `added_triangles::Set{T}`

A set of triangles that are added to the triangulation after the insertion.
- `deleted_triangles::Set{T}`

A set of triangles that are deleted from the triangulation after the insertion.

- `new_segments{T}`
A set of segments that were added after the insertion.
"""
struct InsertionEventStorage{T, E}
    added_triangles::Set{T}
    deleted_triangles::Set{T}
    new_segments::Set{T}
end 
add_triangle!(events::InsertionEventStorage, T) = push!(events.added_triangles, T)
delete_triangle!(events::InsertionEventStorage, T) = push!(events.deleted_triangles, T)
add_edge!(events::InsertionEventStorage, e) = push!(events.new_encroached_edges, e)

function initialise_event_storage(tri::Triangulation)
    T = triangle_type(tri)
    E = edge_type(tri)
    add_set = Set{T}()
    delete_set = Set{T}()
    encroached_set = Set{E}()
    sizehint!(add_set, 16)
    sizehint!(delete_set, 16)
    sizehint!(encroached_set, 32)
    return InsertionEventStorage{T,E}(Set{T}(), Set{T}(), Set{E}())
end

function insert_circumcenter!(tri::Triangulation, T, events::InsertionEventStorage, queue::RefinementQueue, targets::RefinementTargets)
    
end
