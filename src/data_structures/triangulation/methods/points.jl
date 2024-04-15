"""
    get_point(tri::Triangulation, i) -> NTuple{2, Number}
    get_point(tri::Triangulation, i...) -> NTuple{length(i), NTuple{2, Number}}

Returns the coordinates corresponding to the vertices `i...` of `tri`, given as a `Tuple` of the form `(x, y)` for each point.
If `i` is a ghost vertex, then the coordinates of the representative point of the curve associated with `i` are returned instead.
"""
function get_point(tri::Triangulation, i)
    points = get_points(tri)
    if !is_ghost_vertex(i)
        return get_point(points, i)
    elseif i == ∅
        throw(BoundsError(points, i))
    else
        curve_index = get_curve_index(tri, i)
        c = get_representative_point_coordinates(tri, curve_index)
        return c
    end
end
get_point(tri::Triangulation, i::Vararg{Any,N}) where {N} = ntuple(j -> get_point(tri, i[j]), Val(N))

"""
    num_points(tri::Triangulation) -> Integer

Returns the number of points in `tri`.

!!! danger 

    If `tri` has vertices that are not yet present in the triangulation, e.g. if you have deleted vertices or have 
    some submerged vertices in a weighted triangulation, then the corresponding points will still be counted in this 
    function. It is recommended that you instead consider [`num_vertices`](@ref), [`num_solid_vertices`](@ref), or 
    [`num_ghost_vertices`](@ref).
"""
num_points(tri::Triangulation) = num_points(get_points(tri))

"""
    push_point!(tri::Triangulation, x, y)
    push_point!(tri::Triangulation, p)

Pushes the point `p = (x, y)` into the points of `tri`.
"""
push_point!(tri::Triangulation, x, y) = push_point!(get_points(tri), x, y)
push_point!(tri::Triangulation, p) = push_point!(get_points(tri), p)

"""
    pop_point!(tri::Triangulation)

Pops the last point from the points of `tri`.
"""
pop_point!(tri::Triangulation) = pop_point!(get_points(tri))

"""
    set_point!(tri::Triangulation, i, x, y)
    set_point!(tri::Triangulation, i, p)

Sets the `i`th point of `tri` to `p = (x, y)`.
"""
set_point!(tri::Triangulation, i, x, y) = set_point!(get_points(tri), i, x, y)
set_point!(tri::Triangulation, i, p) = set_point!(get_points(tri), i, p)