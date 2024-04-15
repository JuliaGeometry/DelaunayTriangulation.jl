"""
    reset_representative_points!(tri::Triangulation)

Resets each representative point of `tri` to the origin.
"""
function reset_representative_points!(tri::Triangulation)
    representative_point_list = get_representative_point_list(tri)
    for c in values(representative_point_list)
        reset!(c)
    end
    return representative_point_list
end

"""
    empty_representative_points!(tri::Triangulation)

Empties the `Dict` of representative points of `tri`.
"""
function empty_representative_points!(tri::Triangulation)
    representative_point_list = get_representative_point_list(tri)
    empty!(representative_point_list)
    return representative_point_list
end

"""
    update_centroid_after_addition!(tri::Triangulation, curve_index, p)

Updates the centroid of the `curve_index`th curve in `tri` after the addition of the point `p`.
""" 
function update_centroid_after_addition!(tri::Triangulation, curve_index, p)
    I = integer_type(tri)
    F = number_type(tri)
    representative_point_list = get_representative_point_list(tri)
    centroid = get!(RepresentativeCoordinates{I,F}, representative_point_list, curve_index)
    add_point!(centroid, p)
    return representative_point_list
end

"""
    update_centroid_after_deletion!(tri::Triangulation, curve_index, p)

Updates the centroid of the `curve_index`th curve in `tri` after the deletion of the point `p`.
"""
function update_centroid_after_deletion!(tri::Triangulation, curve_index, p)
    I = integer_type(tri)
    F = number_type(tri)
    representative_point_list = get_representative_point_list(tri)
    centroid = get!(RepresentativeCoordinates{I,F}, representative_point_list, curve_index)
    delete_point!(centroid, p)
    return representative_point_list
end

"""
    get_representative_point_coordinates(tri::Triangulation, curve_index) -> NTuple{2, Number}

Returns the coordinates of the representative point of the `curve_index`th curve in `tri`.
"""
function get_representative_point_coordinates(tri::Triangulation, curve_index)
    representative_point_list = get_representative_point_list(tri)
    F = number_type(tri)
    c = representative_point_list[curve_index]
    x, y = getxy(c)
    return (F(x), F(y))
end

"""
    new_representative_point!(tri::Triangulation, curve_index) 

Creates a new representative point for the `curve_index`th curve in `tri`.
"""
function new_representative_point!(tri::Triangulation, curve_index)
    I = integer_type(tri)
    T = number_type(tri)
    representative_point_list = get_representative_point_list(tri)
    representative_point_list[curve_index] = RepresentativeCoordinates(zero(T), zero(T), zero(I))
    return representative_point_list
end

"""
    compute_representative_points!(tri::Triangulation; use_convex_hull=!has_boundary_nodes(tri), precision=one(number_type(tri)))

Computes a new set of representative points for `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref) for which to compute the representative points.

# Keyword Arguments
- `use_convex_hull=!has_boundary_nodes(tri)`: If `true`, then the representative points are computed using the convex hull of the triangulation. Otherwise, the representative points are computed using the boundary nodes of the triangulation.
- `precision=one(number_type(tri))`: The precision to use when computing the representative points via [`pole_of_inaccessibility`](@ref).

# Output 
There are no outputs as `tri` is updated in-place, but for each curve the representative point is computed using [`pole_of_inaccessibility`](@ref).

!!! warning "Exterior curves"

    While `get_exterior_curve_indices(tri)` does store the curves corresponding to exterior curves, this function still treats the first 
    curve as the most important exterior curve, computing the representative point so that it is in no holes. In particular, other exterior curves 
    might have representative points that are in a hole of one of their interior holes. This isn't much of a problem, indeed it wouldn't be a significant 
    problem even if we had the representative point in a hole of the first curve, but it is something to be aware of.
"""
function compute_representative_points!(tri::Triangulation; use_convex_hull=!has_boundary_nodes(tri), precision=one(number_type(tri)))
    reset_representative_points!(tri)
    if use_convex_hull
        return _compute_representative_points!(tri, get_convex_hull_vertices(tri); precision)
    else
        return _compute_representative_points!(tri, get_boundary_nodes(tri); precision)
    end
end
function _compute_representative_points!(tri::Triangulation, boundary_nodes; precision=one(number_type(tri))) # need this to be separate for type stability.
    points = get_points(tri)
    representative_point_list = get_representative_point_list(tri)
    I = integer_type(tri)
    if has_multiple_curves(tri)
        nc = num_curves(tri)
        x, y = pole_of_inaccessibility(points, boundary_nodes; precision) # this is separate so that we make sure the center is not in any of the holes
        representative_point_list[1] = RepresentativeCoordinates(x, y, zero(I))
        for i in 2:nc
            bn = get_boundary_nodes(boundary_nodes, i)
            x, y = pole_of_inaccessibility(points, bn; precision)
            representative_point_list[i] = RepresentativeCoordinates(x, y, zero(I))
        end
    else
        x, y = pole_of_inaccessibility(points, boundary_nodes; precision)
        representative_point_list[1] = RepresentativeCoordinates(x, y, zero(I))
    end
    return representative_point_list
end