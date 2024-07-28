"""
    legalise_edge!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing; predicates::AbstractPredicateKernel=AdaptiveKernel())

Legalises the edge `(i, j)` and other neighbouring edges in `tri` if they are illegal, assuming the vertex `r` 
was just added into a triangle that contains `(i, j)`. [`flip_edge!`](@ref) is used.

See also [`is_legal`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the edge to legalise.
- `j`: The second vertex of the edge to legalise.
- `r`: The vertex that was just added into a triangle that contains `(i, j)`.
- `store_event_history=Val(false)`: Whether to store the event history of the flip.
- `event_history=nothing`: The event history. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
There is no output, as `tri` is updated in-place.

!!! warning "Invalid event histories"

    Edge flipping can lead to `event_history` having triangles both in `event_history.added_triangles` and `event_history.deleted_triangles`.
    To get around this, we only store in these fields the triangles necessary to allow [`undo_insertion!`](@ref) to work,
    so that at a triangle that might have appeared in both will only appear in one.
"""
function legalise_edge!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing; predicates::AbstractPredicateKernel=AdaptiveKernel())
    cert = is_legal(predicates, tri, i, j)
    if is_illegal(cert) 
        e = get_adjacent(tri, j, i)
        flip_edge!(tri, i, j, e, r, store_event_history, event_history)
        legalise_edge!(tri, i, e, r, store_event_history, event_history; predicates)
        legalise_edge!(tri, e, j, r, store_event_history, event_history; predicates)
    end
    return tri
end