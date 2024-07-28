"""
    edge_exists(tri::Triangulation, ij) -> Bool 
    edge_exists(tri::Triangulation, i, j) -> Bool

Tests if the edge `(i, j)` is in `tri`, returning `true` if so and `false` otherwise. 

See also [`unoriented_edge_exists`](@ref).
"""
edge_exists(i::I) where {I <: Integer} = i ≠ I(∅)
edge_exists(tri::Triangulation, ij) = edge_exists(get_adjacent(tri, ij))
edge_exists(tri::Triangulation, i, j) = edge_exists(get_adjacent(tri, i, j))

"""
    unoriented_edge_exists(tri::Triangulation, ij) -> Bool 
    unoriented_edge_exists(tri::Triangulation, i, j) -> Bool

Tests if the unoriented edge `(i, j)` is in `tri`, returning `true` if so and `false` otherwise.
"""
unoriented_edge_exists(tri::Triangulation, ij) = edge_exists(tri, ij) || edge_exists(tri, reverse_edge(ij))
unoriented_edge_exists(tri::Triangulation, i, j) = unoriented_edge_exists(tri, construct_edge(edge_type(tri), i, j))

"""
    has_ghost_triangles(tri::Triangulation) -> Bool

Returns `true` if `tri` has ghost triangles, and `false` otherwise.
"""
function has_ghost_triangles(tri::Triangulation)
    num_ghost_vertices(tri) == 0 && return false
    I = integer_type(tri)
    some_outer_boundary_edges = get_adjacent2vertex(tri, I(𝒢))
    isempty(some_outer_boundary_edges) && return false
    e = (first ∘ each_edge)(some_outer_boundary_edges)
    return edge_exists(tri, terminal(e), I(𝒢))
end

"""
    is_constrained(tri::Triangulation) -> Bool

Returns `true` if `tri` has constrained edges (segments), and `false` otherwise.
"""
is_constrained(tri::Triangulation) = !isempty(get_all_segments(tri))

"""
    has_multiple_curves(tri::Triangulation) -> Bool

Returns `true` if `tri` has multiple boundary curves, and `false` otherwise.
"""
has_multiple_curves(tri::Triangulation) = has_multiple_curves(get_boundary_nodes(tri))

"""
    has_multiple_sections(tri::Triangulation) -> Bool

Returns `true` if `tri` has multiple boundary sections, and `false` otherwise.
"""
has_multiple_sections(tri::Triangulation) = has_multiple_sections(get_boundary_nodes(tri))

"""
    has_boundary_nodes(tri::Triangulation) -> Bool 

Returns `true` if `tri` has boundary nodes, and `false` otherwise.
"""
function has_boundary_nodes(boundary_nodes)
    return has_multiple_sections(boundary_nodes) || num_boundary_edges(boundary_nodes) ≠ 0 || eltype(boundary_nodes) <: AbstractParametricCurve
end
has_boundary_nodes(tri::Triangulation) = has_boundary_nodes(get_boundary_nodes(tri))
has_boundary_nodes(::Nothing) = false

"""
    is_weighted(tri::Triangulation) -> Bool 

Returns `true` if `tri` is weighted, and `false` otherwise.
"""
is_weighted(tri::Triangulation) = is_weighted(get_weights(tri))
