"""
    get_adjacent(tri, uv; check_existence=Val(has_multiple_segments(tri)))
    get_adjacent(tri, u, v; check_existence=Val(has_multiple_segments(tri)))

Returns `get_adjacent(tri, uv; check_existence)` or `get_adjacent(tri, u, v; check_existence)`.
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

Calls `add_adjacent!(get_adjacent(tri), uv, w)` or `add_adjacent!(get_adjacent(tri), u, v, w)`,
adding the edge `uv = (u, v)` to the adjacent map of `tri` with corresponding vertex `w`.
"""
@inline add_adjacent!(tri::Triangulation, uv, w) = add_adjacent!(get_adjacent(tri), uv, w)
@inline add_adjacent!(tri::Triangulation, u, v, w) = add_adjacent!(get_adjacent(tri), u, v, w)

"""
    delete_adjacent!(tri::Triangulation, uv)
    delete_adjacent!(tri::Triangulation, u, v)

Calls `delete_adjacent!(get_adjacent(tri), uv)` or `delete_adjacent!(get_adjacent(tri), u, v)`,
removing the edge `uv = (u, v)` from the adjacent map.
"""
@inline delete_adjacent!(tri::Triangulation, uv) = delete_adjacent!(get_adjacent(tri), uv)
@inline delete_adjacent!(tri::Triangulation, u, v) = delete_adjacent!(get_adjacent(tri), u, v)