"""
    flip_edge!(tri::Triangulation, i, j, store_event_history=Val(false), event_history=nothing)
    flip_edge!(tri::Triangulation, i, j, k, ℓ, store_event_history=Val(false), event_history=nothing)

Flips the edge between vertices `i` and `j` in `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the edge to flip.
- `j`: The second vertex of the edge to flip.
- `k`: The vertex `k = get_adjacent(tri, j, i)`. This is only used in the second method.
- `ℓ`: The vertex `ℓ = get_adjacent(tri, i, j)`. This is only used in the second method.
- `store_event_history=Val(false)`: Whether to store the event history of the flip. 
- `event_history=nothing`: The event history. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object. This storage is done using `store_flip_edge_history!`.

# Outputs 
There is no output, as `tri` is updated in-place.

!!! warning "Invalid flips"

    If `(i, j, k, ℓ)`, where `ℓ = get_adjacent(tri, i, j)` and `k = get_adjacent(tri, j, i)`, is not a convex quadrilateral, then this edge flip will make the triangulation non-planar.
"""
function flip_edge!(tri::Triangulation, i::I, j::I, store_event_history::Val{B} = Val(false), event_history = nothing) where {I <: Integer, B}
    ℓ = get_adjacent(tri, i, j)
    k = get_adjacent(tri, j, i)
    flip_edge!(tri, i, j, k, ℓ, store_event_history, event_history)
    return tri
end
function flip_edge!(tri::Triangulation, i::I, j::I, k::I, ℓ::I, store_event_history::Val{B} = Val(false), event_history = nothing) where {I <: Integer, B}
    delete_triangle!(tri, i, k, j; protect_boundary = true, update_ghost_edges = false)
    delete_triangle!(tri, i, j, ℓ; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, ℓ, k, j; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, ℓ, i, k; protect_boundary = true, update_ghost_edges = false)
    is_true(store_event_history) && store_flip_edge_history!(event_history, i, j, k, ℓ)
    return tri
end


function store_flip_edge_history!(event_history, i, j, k, ℓ)
    trit = triangle_type(event_history)
    Tikj = construct_triangle(trit, i, k, j)
    Tijℓ = construct_triangle(trit, i, j, ℓ)
    Tℓkj = construct_triangle(trit, ℓ, k, j)
    Tℓik = construct_triangle(trit, ℓ, i, k)
    delete_triangle!(event_history, Tikj)
    delete_triangle!(event_history, Tijℓ)
    add_triangle!(event_history, Tℓkj)
    add_triangle!(event_history, Tℓik)
    # Remove duplicate entries 
    delete_triangle!(event_history.added_triangles, Tikj)
    delete_triangle!(event_history.added_triangles, Tijℓ)
    delete_triangle!(event_history.deleted_triangles, Tℓkj)
    delete_triangle!(event_history.deleted_triangles, Tℓik)
    return event_history
end