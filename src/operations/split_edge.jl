"""
    split_edge!(tri::Triangulation, i, j, r)

Given a triangulation `tri` and an edge `(i, j)`, splits the edge at 
the point `r` so that the edges `(i, r)` and `(r, j)` now appear in `tri` 
(with the triangles updated accordingly). It is assumed that `r` is (at least 
very close to) collinear with `(i, j)`.
"""
function split_edge!(tri::Triangulation, i, j, r)
    if is_constrained(tri) && contains_constrained_edge(tri, i, j)
        E = edge_type(tri)
        constrained_edges = get_constrained_edges(tri)
        all_constrained_edges = get_all_constrained_edges(tri)
        delete_edge!(all_constrained_edges, construct_edge(E, i, j))
        delete_edge!(all_constrained_edges, construct_edge(E, j, i))
        if !contains_constrained_edge(tri, r, i)
            add_edge!(all_constrained_edges, construct_edge(E, i, r))
        end
        if !contains_constrained_edge(tri, j, r)
            add_edge!(all_constrained_edges, construct_edge(E, r, j))
        end
        if contains_boundary_edge(tri, i, j)
            split_boundary_edge!(tri, i, j, r)
        elseif contains_boundary_edge(tri, j, i)
            split_boundary_edge!(tri, j, i, r)
        else # The edge isn't a boundary edge, so let's make sure the edge is removed from constrained_edges 
            delete_edge!(constrained_edges, construct_edge(E, i, j))
            delete_edge!(constrained_edges, construct_edge(E, j, i))
            if !contains_edge(r, i, constrained_edges) 
                add_edge!(constrained_edges, construct_edge(E, i, r))
            end
            if !contains_edge(j, r, constrained_edges)
                add_edge!(constrained_edges, construct_edge(E, r, j))
            end
        end
    end
    k = get_adjacent(tri, i, j)
    delete_triangle!(tri, i, j, k; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, i, r, k; protect_boundary=true, update_ghost_edges=false)
    add_triangle!(tri, r, j, k; protect_boundary=true, update_ghost_edges=false)
    return nothing
end

"""
    legalise_split_edge!(tri::Triangulation, i, j, k, r)

Given a triangulation `tri`, an edge `(i, j)` that has already 
been split by [`split_edge!`](@ref) at the point `r`, 
legalises the new edges using [`legalise_edge`](@ref), letting 
`k` be the vertex that was originally adjacent to `(i, j)`.
"""
function legalise_split_edge!(tri::Triangulation, i, j, k, r)
    legalise_edge!(tri, j, k, r)
    legalise_edge!(tri, k, i, r)
    return nothing
end

"""
    complete_split_and_legalise!(tri::Triangulation, i, j, r)

Given a triangulation `tri`, an edge `(i, j)`, and a point `r`,
splits both `(i, j)` and `(j, i)` at `r` using [`split_edge!`](@ref) and then legalises
the new edges using [`legalise_split_edge!`](@ref).
"""
function complete_split_and_legalise!(tri::Triangulation, i, j, r)
    k = get_adjacent(tri, i, j)
    ℓ = get_adjacent(tri, j, i)
    split_edge!(tri, i, j, r)
    split_edge!(tri, j, i, r)
    legalise_split_edge!(tri, i, j, k, r)
    legalise_split_edge!(tri, j, i, ℓ, r)
    return nothing
end