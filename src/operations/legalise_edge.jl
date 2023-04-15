"""
    legalise_edge!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing)

Given a triangulation `tri`, an edge `(i, j)`, and a 
point `r` that was added into a triangle that `(i, j)` 
belongs to, legalises the edge `(i, j)` and other neighbouring 
edges recursively.

If `store_event_history` is `Val(true)`, then the event history is stored in
`event_history`.

!!! warning 

    Edge flipping can lead to final event_historys that have triangles both in added_triangles and deleted_triangles
"""
function legalise_edge!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing)
    cert = is_legal(tri, i, j)
    if is_illegal(cert) 
        e = get_adjacent(tri, j, i)
        flip_edge!(tri, i, j, e, r, store_event_history, event_history)
        legalise_edge!(tri, i, e, r, store_event_history, event_history)
        legalise_edge!(tri, e, j, r, store_event_history, event_history)
    end
    return nothing
end