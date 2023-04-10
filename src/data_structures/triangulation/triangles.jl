"""
    triangle_type(tri::Triangulation)

Returns the type used for representing triangles in `tri`.
"""
@inline triangle_type(::Triangulation{P,Ts}) where {P,Ts} = triangle_type(Ts)

"""
    num_triangles(tri::Triangulation) = num_triangles(get_triangles(tri))

Returns the number of triangles in `tri`, including ghost triangles.
"""
@inline num_triangles(tri::Triangulation) = num_triangles(get_triangles(tri))

"""
    each_triangle(tri::Triangulation)

Returns an iterator over each triangle in `tri`. This iterator could include ghost triangles. 

See also [`each_solid_triangle`](@ref) and [`each_ghost_triangle`](@ref).
"""
@inline each_triangle(tri::Triangulation) = each_triangle(get_triangles(tri))

"""
    contains_triangle(tri::Triangulation, T)
    contains_triangle(tri::Triangulation, i, j, k)

Tests if `tri` contains the triangle `T`.
"""
@inline contains_triangle(tri::Triangulation, T) = contains_triangle(T, get_triangles(tri))
@inline contains_triangle(tri::Triangulation, i, j, k) = contains_triangle(i, j, k, get_triangles(tri))

"""
    construct_positively_oriented_triangle(tri::Triangulation, i, j, k)

Constructs a triangle `T` with vertices `(i, j, k)` in `tri` such that `T` is 
positively oriented.
"""
@inline function construct_positively_oriented_triangle(tri::Triangulation, i, j, k)
    return construct_positively_oriented_triangle(triangle_type(tri), i, j, k, get_points(tri))
end

"""
    num_ghost_triangles(tri::Triangulation)

Returns the number of ghost triangles in `tri`.
"""
@inline function num_ghost_triangles(tri::Triangulation)
    if has_ghost_triangles(tri)
        num_ghosts = 0
        all_bnd_idx = all_boundary_indices(tri)
        for i in all_bnd_idx
            S = get_adjacent2vertex(tri, i)
            num_ghosts += num_edges(S)
        end
    else
        num_ghosts = 0
    end
    return num_ghosts
end

"""
    num_solid_triangles(tri::Triangulation)

Returns the number of solid triangles in `tri`.
"""
@inline num_solid_triangles(tri::Triangulation) = num_triangles(tri) - num_ghost_triangles(tri)

abstract type AbstractEachTriangle{T} end
Base.IteratorSize(::Type{<:AbstractEachTriangle}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachTriangle{T}}) where {T} = Base.IteratorEltype(T)
Base.eltype(::Type{<:AbstractEachTriangle{T}}) where {T} = triangle_type(T)
initialise_triangles(::Type{<:AbstractEachTriangle{T}}) where {T} = initialise_triangles(T)
each_triangle(itr::AbstractEachTriangle) = itr
struct EachSolidTriangle{T,V} <: AbstractEachTriangle{V}
    tri::T
    triangles::V
end
struct EachGhostTriangle{T,V} <: AbstractEachTriangle{V}
    tri::T
    triangles::V
end
Base.length(itr::EachSolidTriangle) = num_solid_triangles(itr.tri)
Base.length(itr::EachGhostTriangle) = num_ghost_triangles(itr.tri)
function Base.iterate(itr::EachSolidTriangle, state...)
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing && return nothing
    tri, state = tri_state
    while is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end
function Base.iterate(itr::EachGhostTriangle, state...)
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing || !has_ghost_triangles(itr.tri) && return nothing
    tri, state = tri_state
    while !is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end

"""
    each_solid_triangle(tri)

Returns an iterator over all triangles of `tri` that are not ghost triangles.
"""
each_solid_triangle(tri) = EachSolidTriangle(tri, each_triangle(tri))

"""
    each_ghost_triangle(tri)

Returns an iterator over all triangles of `tri` that are ghost triangles.
"""
each_ghost_triangle(tri) = EachGhostTriangle(tri, each_triangle(tri))