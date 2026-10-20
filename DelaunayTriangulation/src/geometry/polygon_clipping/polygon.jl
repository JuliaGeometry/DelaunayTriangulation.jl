"""
    Polygon{T,V,P} <: AbstractVector{T}

A struct for representing a polygon. The `vertices` are to be a counter-clockwise list of integers, where 
the integers themselves refer to points in `points`.

# Fields
- `vertices::V`: A list of integers that refer to points in `points`. The last vertex may be the same as the first, in which case it is ignored.
- `points::P`: A list of points.

!!! warning "Aliasing"

    The `vertices` field is exactly the same as the input `vertices`, not a copy. Where `vertices[begin] == vertices[end]`,
    the last element is ignored when indexing into the polygon (this is tracked by the `is_circular` field).
"""
struct Polygon{T, V, P} <: AbstractVector{T}
    vertices::V
    points::P
    is_circular::Bool
    @inline function Polygon(vertices::V, points::P) where {V, P}
        p = get_point(points, vertices[begin])
        T = typeof(p)
        return new{T, V, P}(vertices, points, is_circular(vertices))
    end
end
Base.length(P::Polygon) = length(P.vertices) - P.is_circular
Base.size(P::Polygon) = (length(P),)
Base.getindex(P::Polygon, i::Int) = get_point(P.points, P.vertices[i])
Base.getindex(P::Polygon, i::Vararg{Int, N}) where {N} =
    map(i) do j
    P[j]
end
