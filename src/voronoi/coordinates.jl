"""
    get_polygon_coordinates(vorn::VoronoiTessellation, i; bounding_box = nothing)

Returns a vector for the coordinates of the `i`th polygon in `vorn`. If `bounding_box` 
is provided, the polygon will be clipped to the bounding box, assuming that it takes 
the form `(xmin, xmax, ymin, ymax)`. Some specific cases:

- If the polygon is unbounded but `bounding_box` is `nothing`, then an error will be thrown. 
- If the polygon is bounded and `bounding_box` is `nothing`, then the polygon coordinates will be returned as if without any clipping.
- If the polygon is outside of the bounding box entirely, then an empty vector will be returned.
"""
function get_polygon_coordinates(vorn::VoronoiTessellation, i; bounding_box=nothing)
    if !isnothing(bounding_box)
        a, b, c, d = bounding_box
        @assert a < b && c < d "The bounding box must be of the form (xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax."
    end
    if i ∈ get_unbounded_polygons(vorn)
        isnothing(bounding_box) && throw(ArgumentError("The polygon is unbounded, so a bounding box must be provided."))
        return get_unbounded_polygon_coordinates(vorn, i, bounding_box)
    else
        return get_bounded_polygon_coordinates(vorn, i, bounding_box)
    end
end

"""
    polygon_position_relative_to_box(vorn::VoronoiTessellation, bounding_box, i)

Tests the position of the `i`th polygon of `vorn` relative to the bounding box`, ignoring 
boundary indices. Returns:

- `Cert.Inside`: The polygon is contained entirely within the box.
- `Cert.Outside`: The polygon is entirely outside the box.
- `Cert.Touching`: The polygon is contained entirely within the box, but touches the boundary of the box. Note that if a polygon touches the boundary of the box but is entirely outside otherwise, `Cert.Outside` will be returned instead. 
- `Cert.Multiple`: The polygon intersects the box in multiple places.
"""
function polygon_position_relative_to_box(vorn::VoronoiTessellation, bounding_box, i)
    a, b, c, d = bounding_box
    vertices = get_polygon(vorn, i)
    points = get_polygon_points(vorn)
    flag = polygon_position_relative_to_box(a, b, c, d, vertices, points)
    return flag
end