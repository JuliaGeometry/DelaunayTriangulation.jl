struct Polygon{T,V,P} <: AbstractVector{T}
    vertices::V
    points::P
    function Polygon(vertices::V, points::P) where {V,P}
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

"""
    clip_polygon(vertices, points, clip_vertices, clip_points)

Clips the counter-clockwise polygon defined by `(vertices, points)` against the clip 
polygon `(clip_vertices, clip_points)`, assumed to be convex and counter-clockwise. The 
Sutherland-Hodgman algorithm is used. It is assumed that 
`vertices[begin] == vertices[end]` and `clip_vertices[begin] == clip_vertices[end]`.

The returned result is another counter-clockwise polygon `P` that also satisfies `P[begin] == P[end]`.
"""
function clip_polygon(vertices, points, clip_vertices, clip_points)
    return clip_polygon(Polygon(vertices, points), Polygon(clip_vertices, clip_points))
end

function clip_polygon(poly::Polygon, clip_poly::Polygon{T}) where {T}
    output_list = poly
    q = clip_poly[end] 
    for p in clip_poly 
        input_list = output_list 
        output_list = T[]
        s = input_list[end] 
        for vertex in input_list 
            if (is_left ∘ point_position_relative_to_line)(q, p, vertex) 
                flag = point_position_relative_to_line(q, p, s)
                if !is_left(flag) # allow for is_on 
                    r = segment_intersection_coordinates(q, p, s, vertex)
                    push!(output_list, r)
                end
                push!(output_list, vertex)
            elseif (is_left ∘ point_position_relative_to_line)(q, p, s)
                r = segment_intersection_coordinates(q, p, s, vertex)
                push!(output_list, r)
            end
            s = vertex 
        end
        q = p 
    end
    push!(output_list, output_list[begin])
    return output_list::Vector{T}
end