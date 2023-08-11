struct Polygon{T,V,P} <: AbstractVector{T}
    vertices::V
    points::P
    @inline function Polygon(vertices::V, points::P) where {V,P}
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

Note that this algorithm may potentially return overlapping edges, but this is fine for rendering.
"""
function clip_polygon(vertices, points, clip_vertices, clip_points)
    return clip_polygon(Polygon(vertices, points), Polygon(clip_vertices, clip_points))
end

"""
    clip_polygon_to_edge(input_list, q, p)

Clips the input polygon `input_list` against the edge `qp`. This is the 
intermediate step used in the Sutherland-Hodgman algorithm.
"""
function clip_polygon_to_edge(input_list, q, p)
    output_list = eltype(input_list)[]
    s = input_list[end]
    for vertex in input_list
        # By considering s = input_list[end], we can consider each edge of the 
        # input polygon one at a time, with the initial edge being s-vertex. 
        # The first step is to check if vertex is to the left or to the right of the 
        # edge qp. If it's to the left, then it is outside of the original polygon  
        # (note the assumption here that the poylgon is counter-clockwise).
        if (is_left ∘ point_position_relative_to_line)(q, p, vertex)
            # Now that we know that vertex is outside of the polygon, we need to know 
            # if there is any intersection to consider. In particular, does s-vertex intersect 
            # the edge qp? If s is also to the left of the line, then no there is no intersection 
            # and we can safely ignore the intersection (note that this check makes the implicit 
            # assumption that the clipping polygon is convex). If s is to the right of the line,
            # then there is an intersection and we need to consider it.
            flag = point_position_relative_to_line(q, p, s)
            if !is_left(flag) # allow for is_on. If left, then the edge s-vertex is not inside the clipping polygon, so there is no intersection to consider
                r = segment_intersection_coordinates(q, p, s, vertex)
                push!(output_list, r)
            end
            # We still need to add back in the vertex to the output list regardless so that 
            # we can keep going with the next edge of the clipping polygon easily.
            push!(output_list, vertex)
        elseif (is_left ∘ point_position_relative_to_line)(q, p, s)
            # In this case, vertex is to the right of qp, but now s is outside. This means there 
            # is an intersection of s-vertex with qp.
            r = segment_intersection_coordinates(q, p, s, vertex)
            push!(output_list, r)
        end
        s = vertex
    end
    return output_list
end

function clip_polygon(poly::Polygon, clip_poly::Polygon{T}) where {T}
    output_list = poly
    q = clip_poly[end]
    for p in clip_poly
        input_list = output_list 
        isempty(output_list) && break
        output_list = clip_polygon_to_edge(input_list, q, p)
        q = p
    end
    !isempty(output_list) && push!(output_list, output_list[begin])
    return output_list::Vector{T}
end