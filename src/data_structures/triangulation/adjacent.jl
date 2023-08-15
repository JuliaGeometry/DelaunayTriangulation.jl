"""
    get_adjacent(tri, uv; check_existence=has_multiple_segments(tri))
    get_adjacent(tri, u, v; check_existence=has_multiple_segments(tri))

Returns the vertex `w` adjacent to the edge `(u, v)` such that `(u, v, w)` is a positively oriented triangle in `tri`.
If `check_existence` is `true`, then care is taken in case `tri` has a segmented boundary such that multiple boundary indices 
may be present, in case `(u, v)` is on the boundary.
"""
@inline function get_adjacent(tri::Triangulation, uv; check_existence::V=Val(has_multiple_segments(tri))) where {V}
    return get_adjacent(get_adjacent(tri), uv; check_existence, boundary_index_ranges=get_boundary_index_ranges(tri))
end
@inline function get_adjacent(tri::Triangulation, u, v; check_existence::V=Val(has_multiple_segments(tri))) where {V}
    return get_adjacent(get_adjacent(tri), u, v; check_existence, boundary_index_ranges=get_boundary_index_ranges(tri))
end

"""
    add_adjacent!(tri::Triangulation, uv, w)
    add_adjacent!(tri::Triangulation, u, v, w)

Updates the adjacent map of `tri` so that `w` is adjacent to `(u, v)`.
"""
@inline add_adjacent!(tri::Triangulation, uv, w) = add_adjacent!(get_adjacent(tri), uv, w)
@inline add_adjacent!(tri::Triangulation, u, v, w) = add_adjacent!(get_adjacent(tri), u, v, w)

"""
    delete_adjacent!(tri::Triangulation, uv)
    delete_adjacent!(tri::Triangulation, u, v)

Deletes `(u, v)` from the adjacent map of `tri`.
"""
@inline delete_adjacent!(tri::Triangulation, uv) = delete_adjacent!(get_adjacent(tri), uv)
@inline delete_adjacent!(tri::Triangulation, u, v) = delete_adjacent!(get_adjacent(tri), u, v)