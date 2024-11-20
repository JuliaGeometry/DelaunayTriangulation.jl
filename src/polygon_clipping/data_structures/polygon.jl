struct Polygon{T,V,P} <: AbstractVector{T}
    vertices::V
    points::P
    is_circular::Bool
    @inline function Polygon(vertices::V, points::P) where {V,P}
        p = get_point(points, vertices[begin])
        T = typeof(p)
        return new{T,V,P}(vertices, points, is_circular(vertices))
    end
end
get_vertices(P::Polygon) = P.vertices

get_points(P::Polygon) = P.points

is_circular(P::Polygon) = P.is_circular

Base.length(P::Polygon) = length(get_vertices(P)) - is_circular(P)

Base.size(P::Polygon) = (length(P),)

Base.getindex(P::Polygon, i::Int) = get_point(get_points(P), get_vertices(P)[i])