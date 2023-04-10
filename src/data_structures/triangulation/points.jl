"""
    get_point(tri::Triangulation, i...)

Returns the `i`th point of `tri. Boundary indices are automatically handled.
"""
@inline function get_point(tri::Triangulation, i)
    return get_point(get_points(tri), get_representative_point_list(tri), get_boundary_map(tri), i)
end
@inline function get_point(tri::Triangulation, i::Vararg{Any,N}) where {N}
    return get_point(get_points(tri), get_representative_point_list(tri), get_boundary_map(tri), i...)
end

"""
    each_point_index(tri::Triangulation)

Returns an iterator over the indices of the points of `tri`.

See also [`each_vertex`](@ref), [`each_solid_vertex](@ref) and [`each_ghost_vertex`](@ref).
"""
@inline each_point_index(tri::Triangulation) = each_point_index(get_points(tri))

"""
    each_point(tri::Triangulation)

Returns an iterator over the points of `tri`. 

See also [`each_vertex`](@ref), [`each_solid_vertex](@ref) and [`each_ghost_vertex`](@ref).
"""
@inline each_point(tri::Triangulation) = each_point(get_points(tri))

"""
    num_points(tri::Triangulation)

Returns the number of points of `tri`.

!!! note 

    Note that this is just the size of `get_points(tri)`,
    but if there are some missing points then this will not match the number in the triangulation itself. 
    Use [`num_vertices`](@ref) to get the number of vertices in the triangulation.
"""
@inline num_points(tri::Triangulation) = num_points(get_points(tri))

"""
    push_point!(tri::Triangulation, x, y) = push_point!(get_points(tri), x, y)
    push_point!(tri::Triangulation, p) = push_point!(get_points(tri), p)

Pushes the point `p = (x, y)` into `get_points(tri0`.

!!! note

    This does not add the point `p` into the triangulation itself. See [`add_point!`](@ref) for this.
"""
@inline push_point!(tri::Triangulation, x, y) = push_point!(get_points(tri), x, y)
@inline push_point!(tri::Triangulation, p) = push_point!(get_points(tri), p)

abstract type AbstractEachVertex{V} end
Base.IteratorSize(::Type{<:AbstractEachVertex}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachVertex{V}}) where {V} = Base.IteratorEltype(V)
Base.eltype(::Type{<:AbstractEachVertex{V}}) where {V} = eltype(V)
each_vertex(itr::AbstractEachVertex) = itr
num_vertices(itr::AbstractEachVertex) = length(itr)
struct EachSolidVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end
Base.length(verts::EachSolidVertex) = num_vertices(verts.tri) - length(all_boundary_indices(verts.tri)) * (integer_type(verts.tri)(BoundaryIndex) ∈ verts.vertices)
struct EachGhostVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end
Base.length(verts::EachGhostVertex) = length(all_boundary_indices(verts.tri))
function Base.iterate(itr::EachSolidVertex, state...)
    vertices_state = iterate(itr.vertices, state...)
    vertices_state === nothing && return nothing
    vertices, state = vertices_state
    while is_boundary_index(vertices)
        vertices_state = iterate(itr.vertices, state...)
        vertices_state === nothing && return nothing
        vertices, state = vertices_state
    end
    return vertices, state
end
Base.iterate(itr::EachGhostVertex, state...) = Base.iterate(itr.vertices, state...)

"""
    each_solid_vertex(tri::Triangulation)

Returns an iterator over the solid vertices in the triangulation `tri`, i.e. 
those that are not ghost vertices.
"""
each_solid_vertex(tri) = EachSolidVertex(each_vertex(tri), tri)

"""
    each_ghost_vertex(tri::Triangulation)

Returns an iterator over the ghost vertices in the triangulation `tri`.
""" 
each_ghost_vertex(tri) = EachGhostVertex(all_boundary_indices(tri), tri)

function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachSolidVertex})
    itr = v[]
    verts = itr.vertices 
    r = rand(rng, verts)
    while is_boundary_index(r)
        r = rand(rng, verts)
    end
    return r
end
function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachGhostVertex})
    itr = v[]
    verts = itr.vertices 
    r = rand(rng, verts)
    while !is_boundary_index(r)
        r = rand(rng, verts)
    end
    return r
end