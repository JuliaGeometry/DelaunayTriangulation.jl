"""
    edge_type(tri::Triangulation)

Returns the type used for representing edges in `tri`.
"""
@inline edge_type(::Triangulation{P,Ts,I,E}) where {P,Ts,I,E} = E

"""
    num_edges(tri::Triangulation)

Returns `num_edges(get_graph(tri))`, the number of edges in `tri`.
"""
@inline num_edges(tri::Triangulation) = num_edges(get_graph(tri))

"""
    each_edge(tri::Triangulation)

Returns `get_edges(get_graph(tri))`. This iterator over the edges could include ghost edges.

See also [`each_solid_edge`](@ref) and [`each_ghost_edge`](@ref)
"""
@inline each_edge(tri::Triangulation) = get_edges(get_graph(tri))

"""
    each_constrained_edge(tri::Triangulation)

Returns `each_edge(get_all_constrained_edges(tri))`, an iterator over all `constrained_edges` 
(including all constrained boundary edges)`.
"""
@inline each_constrained_edge(tri::Triangulation) = each_edge(get_all_constrained_edges(tri))

"""
    contains_constrained_edge(tri::Triangulation, e)
    contains_constrained_edge(tri::Triangulation, i, j)

Tests if `tri` contains the constrained edge `e = (i, j)` (or `reverse_edge(e) = (j, i)`).
"""
@inline function contains_constrained_edge(tri::Triangulation, e)
    constrained_edges = get_all_constrained_edges(tri)
    return contains_edge(e, constrained_edges) || contains_edge(reverse_edge(e), constrained_edges)
end
@inline function contains_constrained_edge(tri::Triangulation, i, j)
    E = edge_type(tri)
    e = construct_edge(E, i, j)
    return contains_constrained_edge(tri, e)
end

"""
    num_ghost_edges(tri::Triangulation)

Returns the number of ghost edges in `tri`.
"""
function num_ghost_edges(tri::Triangulation)
    all_bnd_idx = all_boundary_indices(tri)
    num_ghosts = 0
    for i in all_bnd_idx
        bnd_ngh = get_neighbours(tri, i)
        num_ghosts += length(bnd_ngh)
    end
    return num_ghosts
end

"""
    num_solid_edges(tri::Triangulation)

Returns the number of solid edges in `tri`.
"""
num_solid_edges(tri::Triangulation) = num_edges(tri) - num_ghost_edges(tri)

abstract type AbstractEachEdge{E} end
Base.IteratorSize(::Type{<:AbstractEachEdge}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachEdge{E}}) where {E} = Base.IteratorEltype(E)
Base.eltype(::Type{<:AbstractEachEdge{E}}) where {E} = edge_type(E)
initialise_edges(::Type{<:AbstractEachEdge{E}}) where {E} = initialise_edges(E)
each_edge(itr::AbstractEachEdge) = itr
struct EachSolidEdge{T,E} <: AbstractEachEdge{E}
    tri::T
    edges::E
end
struct EachGhostEdge{T,E} <: AbstractEachEdge{E}
    tri::T
    edges::E
end
Base.length(itr::EachSolidEdge) = num_solid_edges(itr.tri)
Base.length(itr::EachGhostEdge) = num_ghost_edges(itr.tri)
function Base.iterate(itr::EachSolidEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end
function Base.iterate(itr::EachGhostEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while !is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end

"""
    each_solid_edge(tri)

Returns an iterator over all solid edges of the triangulation `tri`, i.e. 
over all edges that are not ghost edges.

See also [`each_edge`](@ref) and [`each_ghost_edge`](@ref).
"""
each_solid_edge(tri) = EachSolidEdge(tri, each_edge(tri))

"""
    each_ghost_edge(tri)

Returns an iterator over all ghost edges of the triangulation `tri`, i.e.
over all edges that have a boundary index as one of its indices.

See also [`each_edge`](@ref) and [`each_solid_edge`](@ref).
"""
each_ghost_edge(tri) = EachGhostEdge(tri, each_edge(tri))

"""
    sort_edge_by_degree(tri::Triangulation, e)

Returns `sort_edge_by_degree(e, get_graph(tri))`.
"""
@inline sort_edge_by_degree(tri::Triangulation, e) = sort_edge_by_degree(e, get_graph(tri))

"""
    split_constrained_edge!(tri::Triangulation, constrained_edge, collinear_segments)

Calls `split_constrained_edge!(get_constrained_edges(tri), constrained_edge, collinear_segments)`, splitting 
the `constrained_edge` which is assumed to represent the union of the `collinear_segments` so that 
we instead store those `collinear_segments`.
"""
@inline split_constrained_edge!(tri::Triangulation, constrained_edge, collinear_segments) = split_constrained_edge!(get_constrained_edges(tri), constrained_edge, collinear_segments)
