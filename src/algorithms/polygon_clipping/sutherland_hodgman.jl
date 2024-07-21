"""
    clip_polygon(vertices, points, clip_vertices, clip_points) -> Vector

Clip a polygon defined by `(vertices, points)` to a convex clip polygon defined by `(clip_vertices, clip_points)` with the Sutherland-Hodgman algorithm.
The polygons should be defined in counter-clockwise order.

# Arguments 
- `vertices`: The vertices of the polygon to be clipped.
- `points`: The underlying point set that the vertices are defined over. 
- `clip_vertices`: The vertices of the clipping polygon.
- `clip_points`: The underlying point set that the clipping vertices are defined over.

# Output 
- `clipped_polygon`: The coordinates of the clipped polygon, given in counter-clockwise order and `clipped_polygon[begin] == clipped_polygon[end]`.
"""
function clip_polygon(vertices, points, clip_vertices, clip_points)
    return clip_polygon(Polygon(vertices, points), Polygon(clip_vertices, clip_points))
end

function clip_polygon_to_edge(input_list, q, p)
    output_list = eltype(input_list)[]
    s = input_list[end]
    for vertex in input_list
        # By considering s = input_list[end], we can consider each edge of the 
        # input polygon one at a time, with the initial edge being s-vertex. 
        # The first step is to check if vertex is to the left or to the right of the 
        # edge qp. If it's to the left, then it is outside of the original polygon  
        # (note the assumption here that the polygon is counter-clockwise).
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

function clip_polygon(poly::Polygon{T}, clip_poly::Polygon) where {T}
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
