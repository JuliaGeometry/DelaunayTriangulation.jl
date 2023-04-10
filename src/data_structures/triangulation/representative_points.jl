"""
    compute_representative_points!(tri::Triangulation; use_convex_hull=!has_multiple_segments(tri) && num_boundary_edges(get_boundary_nodes(tri)) == 0)

Updates `get_representative_point_list(tri)` to match the current position of the boundaries. If there are no boundary nodes, `use_convex_hull`
instead represents them using the indices of the convex hull.
"""
@inline function compute_representative_points!(tri::Triangulation; use_convex_hull=!has_multiple_segments(tri) && num_boundary_edges(get_boundary_nodes(tri)) == 0)
    if !use_convex_hull
        compute_representative_points!(get_representative_point_list(tri), get_points(tri), get_boundary_nodes(tri))
    else
        compute_representative_points!(get_representative_point_list(tri), get_points(tri), get_convex_hull_indices(tri))
    end
end

"""
    get_representative_point_coordinates(tri::Triangulation, i)

Returns the coordinates of the `i`th representative point in `tri`.
"""
get_representative_point_coordinates(tri::Triangulation, i) = get_representative_point_coordinates(get_representative_point_list(tri), i)

"""
    reset_representative_points!(tri::Triangulation)

Resets `get_representative_point_list(tri)`.
"""
reset_representative_points!(tri::Triangulation) = reset_representative_points!(get_representative_point_list(tri))

"""
    update_centroid_after_addition!(tri::Triangulation, i, p)

After the point `p` has been added into `tri`, this updates the `i`th representative point in 
`get_representative_point_list(tri)`, treating it as if it were a centroid.
"""
update_centroid_after_addition!(tri::Triangulation, i, p) = update_centroid_after_addition!(get_representative_point_list(tri), i, p)

"""
    update_centroid_after_deletion!(tri::Triangulation, i, p)

After the point `p` has been deleted from `tri`, this updates the `i`th representative point in
`get_representative_point_list(tri)`, treating it as if it were a centroid.
"""
update_centroid_after_deletion!(tri::Triangulation, i, p) = update_centroid_after_deletion!(get_representative_point_list(tri), i, p)

"""
    new_representative_point!(tri::Triangulation, i)

Adds a new representative point to `get_representative_point_list(tri)` with index `i`.
"""
new_representative_point!(tri::Triangulation, i) = new_representative_point!(get_representative_point_list(tri), i)
