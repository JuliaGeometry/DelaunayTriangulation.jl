"""
    clip_polygon(vertices, points, clip_vertices, clip_points; predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Vector

Clip a polygon defined by `(vertices, points)` to a convex clip polygon defined by `(clip_vertices, clip_points)` with the Sutherland-Hodgman algorithm.
The polygons should be defined in counter-clockwise order.

# Arguments 
- `vertices`: The vertices of the polygon to be clipped.
- `points`: The underlying point set that the vertices are defined over. 
- `clip_vertices`: The vertices of the clipping polygon.
- `clip_points`: The underlying point set that the clipping vertices are defined over.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output 
- `clipped_polygon`: The coordinates of the clipped polygon, given in counter-clockwise order and `clipped_polygon[begin] == clipped_polygon[end]`.
"""
function clip_polygon(vertices, points, clip_vertices, clip_points; predicates::AbstractPredicateKernel = AdaptiveKernel())
    return clip_polygon(Polygon(vertices, points), Polygon(clip_vertices, clip_points); predicates)
end

function clip_polygon_to_edge(input_list, q, p, predicates::AbstractPredicateKernel)
    output_list = eltype(input_list)[]
    s = input_list[end]
    for vertex in input_list
        # By considering s = input_list[end], we can consider each edge of the
        # input polygon one at a time, with the initial edge being s-vertex.
        # The first step is to check if vertex is to the left or to the right of the
        # edge qp. If it's to the left, then it is on the same side of qp as the clipping polygon
        # (note the assumption here that the clipping polygon is counter-clockwise), so we keep it.
        if is_left(point_position_relative_to_line(predicates, q, p, vertex))
            # Now that we know that vertex is kept, we need to know if there is any intersection
            # to consider. In particular, does s-vertex intersect the edge qp? If s is also to the
            # left of the line, then there is no intersection and we can safely ignore it (note that
            # this check makes the implicit assumption that the clipping polygon is convex). If s is
            # not to the left of the line, then there is an intersection and we need to consider it.
            flag = point_position_relative_to_line(predicates, q, p, s)
            if !is_left(flag) # allow for is_on. If left, then both s and vertex are kept, so there is no intersection to consider
                r = segment_intersection_coordinates(q, p, s, vertex)
                push!(output_list, r)
            end
            # We then add the vertex itself to the output list, since it is kept.
            push!(output_list, vertex)
        elseif is_left(point_position_relative_to_line(predicates, q, p, s))
            # In this case, vertex is not to the left of qp and so is clipped away, but s is kept.
            # This means there is an intersection of s-vertex with qp.
            r = segment_intersection_coordinates(q, p, s, vertex)
            push!(output_list, r)
        end
        s = vertex
    end
    return output_list
end

function clip_polygon(poly::Polygon{T}, clip_poly::Polygon; predicates::AbstractPredicateKernel = AdaptiveKernel()) where {T}
    output_list = poly
    q = clip_poly[end]
    for p in clip_poly
        input_list = output_list
        isempty(output_list) && break
        output_list = clip_polygon_to_edge(input_list, q, p, predicates)
        q = p
    end
    !isempty(output_list) && push!(output_list, output_list[begin])
    return output_list::Vector{T}
end
