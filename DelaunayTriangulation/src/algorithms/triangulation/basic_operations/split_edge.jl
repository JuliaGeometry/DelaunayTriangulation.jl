"""
    get_edges_for_split_edge(tri::Triangulation, i, j, r)

Returns the edges `(i, j)`, `(j, i)`, `(i, r)`, `(r, i)`, `(r, j)`, and `(j, r)`.
"""
function get_edges_for_split_edge(tri::Triangulation, i, j, r)
    E = edge_type(tri)
    eᵢⱼ = construct_edge(E, i, j)
    eⱼᵢ = reverse_edge(eᵢⱼ)
    eᵢᵣ = construct_edge(E, i, r)
    eᵣᵢ = reverse_edge(eᵢᵣ)
    eᵣⱼ = construct_edge(E, r, j)
    eⱼᵣ = reverse_edge(eᵣⱼ)
    return eᵢⱼ, eⱼᵢ, eᵢᵣ, eᵣᵢ, eᵣⱼ, eⱼᵣ
end

"""
    split_edge!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing)

Splits the edge `(i, j)` in `tri` at the vertex `r`. For the triangulation to be valid after this splitting, it is assumed that `r` is collinear with, 
or at least very close to collinear with, the edge `(i, j)`.

See also [`legalise_split_edge!`](@ref) and [`complete_split_edge_and_legalise!`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).    
- `i`: The first vertex of the edge to split.
- `j`: The second vertex of the edge to split.
- `r`: The vertex to split the edge at.
- `store_event_history=Val(false)`: Whether to store the event history of the flip.
- `event_history=nothing`: The event history. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.

# Outputs 
There is no output, as `tri` is updated in-place.

!!! warning "Handling unoriented edges"

    The triangulation will only be updated as if `(i, j)` has been split rather than also `(j, i)`. You will need to call `split_edge!` again with `(j, i)` if you want to split that edge as well.
"""
function split_edge!(tri::Triangulation, i, j, r, store_event_history = Val(false), event_history = nothing)
    if is_constrained(tri) && contains_segment(tri, i, j)
        is_bnd = contains_boundary_edge(tri, i, j) || contains_boundary_edge(tri, j, i)
        eᵢⱼ, eⱼᵢ, eᵢᵣ, eᵣᵢ, eᵣⱼ, eⱼᵣ = get_edges_for_split_edge(tri, i, j, r)
        all_segments = get_all_segments(tri)
        interior_segments = get_interior_segments(tri)
        delete_edge!(all_segments, eᵢⱼ)
        delete_edge!(all_segments, eⱼᵢ)
        if is_true(store_event_history) && !is_bnd
            !contains_edge(eᵢⱼ, event_history.deleted_segments) && delete_edge!(event_history, eᵢⱼ)
        end
        if !contains_segment(tri, r, i)
            add_edge!(all_segments, eᵢᵣ)
            is_true(store_event_history) && !is_bnd && !contains_edge(eᵣᵢ, event_history.added_segments) && add_edge!(event_history, eᵢᵣ)
        end
        if !contains_segment(tri, j, r)
            add_edge!(all_segments, eᵣⱼ)
            is_true(store_event_history) && !is_bnd && !contains_edge(eⱼᵣ, event_history.added_segments) && add_edge!(event_history, eᵣⱼ)
        end
        if contains_boundary_edge(tri, i, j)
            split_boundary_edge!(tri, i, j, r)
            is_true(store_event_history) && split_boundary_edge!(event_history, i, j, r)
        elseif contains_boundary_edge(tri, j, i)
            split_boundary_edge!(tri, j, i, r)
            is_true(store_event_history) && split_boundary_edge!(event_history, j, i, r)
        else # The edge isn't a boundary edge, so let's make sure the edge is removed from interior_segments
            delete_edge!(interior_segments, eᵢⱼ)
            delete_edge!(interior_segments, eⱼᵢ)
            if is_true(store_event_history)
                !contains_edge(eᵢⱼ, event_history.deleted_segments) && delete_edge!(event_history, eᵢⱼ)
            end
            if !contains_edge(eᵣᵢ, interior_segments)
                add_edge!(interior_segments, eᵢᵣ)
                is_true(store_event_history) && add_edge!(event_history, eᵢᵣ)
            end
            if !contains_edge(eⱼᵣ, interior_segments)
                add_edge!(interior_segments, eᵣⱼ)
                is_true(store_event_history) && add_edge!(event_history, eᵣⱼ)
            end
        end
    end
    k = get_adjacent(tri, i, j)
    delete_triangle!(tri, i, j, k; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, i, r, k; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, r, j, k; protect_boundary = true, update_ghost_edges = false)
    if is_true(store_event_history)
        trit = triangle_type(tri)
        Tᵢⱼₖ = construct_triangle(trit, i, j, k)
        Tᵢᵣₖ = construct_triangle(trit, i, r, k)
        Tᵣⱼₖ = construct_triangle(trit, r, j, k)
        delete_triangle!(event_history, Tᵢⱼₖ)
        add_triangle!(event_history, Tᵢᵣₖ)
        add_triangle!(event_history, Tᵣⱼₖ)
    end
    return tri
end

"""
    legalise_split_edge!(tri::Triangulation, i, j, k, r, store_event_history=Val(false), event_history=nothing; predicates::AbstractPredicateKernel=AdaptiveKernel())

Legalises the newly added edges in `tri` after the edge `(i, j)` was split using [`split_edge!`](@ref).

See also [`complete_split_edge_and_legalise!`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the edge that was split.
- `j`: The second vertex of the edge that was split.
- `k`: The vertex that was originally adjacent to `(i, j)`.
- `r`: The vertex that `(i, j)` was split at.
- `store_event_history=Val(false)`: Whether to store the event history of the flip.
- `event_history=nothing`: The event history. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.

# Keyword Arguments
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, as `tri` is updated in-place.
"""
function legalise_split_edge!(tri::Triangulation, i, j, k, r, store_event_history = Val(false), event_history = nothing; predicates::AbstractPredicateKernel = AdaptiveKernel())
    legalise_edge!(tri, j, k, r, store_event_history, event_history; predicates)
    legalise_edge!(tri, k, i, r, store_event_history, event_history; predicates)
    return tri
end

"""
    complete_split_edge_and_legalise!(tri::Triangulation, i, j, r, store_event_history=Val(false), event_history=nothing; predicates::AbstractPredicateKernel=AdaptiveKernel())

Given a triangulation `tri`, an edge `(i, j)`, and a point `r`, splits both `(i, j)` and `(j, i)` at `r` using [`split_edge!`](@ref) and then subsequently legalises the new edges with [`legalise_split_edge!`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `i`: The first vertex of the edge to split.
- `j`: The second vertex of the edge to split.
- `r`: The vertex to split the edge at.
- `store_event_history=Val(false)`: Whether to store the event history of the flip.
- `event_history=nothing`: The event history. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, as `tri` is updated in-place.
"""
function complete_split_edge_and_legalise!(tri::Triangulation, i, j, r, store_event_history = Val(false), event_history = nothing; predicates::AbstractPredicateKernel = AdaptiveKernel())
    k = get_adjacent(tri, i, j)
    ℓ = get_adjacent(tri, j, i)
    split_edge!(tri, i, j, r, store_event_history, event_history)
    split_edge!(tri, j, i, r, store_event_history, event_history)
    legalise_split_edge!(tri, i, j, k, r, store_event_history, event_history; predicates)
    legalise_split_edge!(tri, j, i, ℓ, r, store_event_history, event_history; predicates)
    if is_true(store_event_history)
        # The deleted_triangles in event_history will contain triangles with the vertex r, which didn't actually appear 
        # initially. So, we need to deal any triangles with r as a vertex from event_history.deleted_triangles.
        filter!(
            T -> let r = r
                r ∉ triangle_vertices(T)
            end, event_history.deleted_triangles,
        )
    end
    return tri
end
