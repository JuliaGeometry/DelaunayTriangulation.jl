"""
    flip_edge!(tri::Triangulation, i, j, store_event_history=Val(false), event_history=nothing)
    flip_edge!(tri::Triangulation, i, j, k, ℓ)

Given a triangulation `tri` and an edge `(i, j)` appearing in the triangulation, 
flips the edge `(i, j)` so that it becomes `(ℓ, k)`, where `ℓ = get_adjacent(tri, i, j)`
and `k = get_adjacent(tri, j, i)`. 

!!! warning 

    If `(i, j, ℓ, k)` is not a convex quadrilateral, than this edge flip makes the triangulation non-planar.
"""
function flip_edge!(tri::Triangulation, i::I, j::I, store_event_history::Val{B}=Val(false), event_history=nothing) where {I<:Integer,B}
    ℓ = get_adjacent(tri, i, j)
    k = get_adjacent(tri, j, i)
    flip_edge!(tri, i, j, k, ℓ, store_event_history, event_history)
    return nothing
end
function flip_edge!(tri::Triangulation, i::I, j::I, k::I, ℓ::I, store_event_history::Val{B}=Val(false), event_history=nothing) where {I<:Integer,B}
    delete_triangle!(tri, i, k, j; protect_boundary=true, update_ghost_edges=false)
    delete_triangle!(tri, i, j, ℓ; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, ℓ, k, j; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, ℓ, i, k; protect_boundary=true, update_ghost_edges=false)
    if is_true(store_event_history)
        trit = triangle_type(tri)
        delete_triangle!(event_history, construct_triangle(trit, i, k, j))
        delete_triangle!(event_history, construct_triangle(trit, i, j, ℓ))
        add_triangle!(event_history, construct_triangle(trit, ℓ, k, j))
        add_triangle!(event_history, construct_triangle(trit, ℓ, i, k))
    end
    return nothing
end