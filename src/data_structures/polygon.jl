"""
    Polygon{T,V,P} <: AbstractVector{T}

A struct for representing a polygon. The `vertices` are to be a counter-clockwise list of integers, where 
the integers themselves refer to points in `points`.

# Fields
- `vertices::V`: A list of integers that refer to points in `points`. The last vertex shokuld not be the same as the first.
- `points::P`: A list of points.

!!! warning "Aliasing"

    In the case where `vertices[begin] ≠ vertices[end]`, the `vertices` field is exactly the same as the input `vertices`. Where 
    `vertices[begin] = vertices[end]`, the `vertices` field is a view of `vertices` that excludes the last element.
"""
struct Polygon{T,V,P} <: AbstractVector{T}
    vertices::V
    points::P
    @stable default_union_limit = 2 @inline function Polygon(vertices::V, points::P) where {V,P}
        p = get_point(points, vertices[begin])
        T = typeof(p)
        if vertices[begin] ≠ vertices[end]
            return new{T,V,P}(vertices, points)
        else
            _verts = @views vertices[begin:(end-1)]
            return new{T,typeof(_verts),P}(_verts, points)
        end
    end
end
Base.size(P::Polygon) = (length(P.vertices),)
Base.getindex(P::Polygon, i::Int) = get_point(P.points, P.vertices[i])
Base.getindex(P::Polygon, i::Vararg{Int,N}) where {N} =
    map(i) do j
        P[j]
    end