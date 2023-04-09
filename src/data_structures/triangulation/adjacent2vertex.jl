"""
    get_adjacent2vertex(tri::Triangulation, w)

Returns `get_adjacent2vertex(get_adjacent2vertex(tri), w)`, the edges `(u, v)` in `tri` 
that are adjacent to `w`.
"""
@inline get_adjacent2vertex(tri::Triangulation, w) = get_adjacent2vertex(get_adjacent2vertex(tri), w)

"""
    add_adjacent2vertex!(tri::Triangulation, w, uv)
    add_adjacent2vertex!(tri::Triangulation, w, u, v)

Calls `add_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)` or `add_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)`,
pushing the edge `uv = (u, v)` into `get_adjacent2vertex(tri, w)`.
"""
@inline add_adjacent2vertex!(tri::Triangulation, w, uv) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
@inline add_adjacent2vertex!(tri::Triangulation, w, u, v) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)

"""
    delete_adjacent2vertex!(tri::Triangulation, w, uv)
    delete_adjacent2vertex!(tri::Triangulation, w, u, v)

Calls `delete_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)` or `delete_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)`,
deleting the edge `uv = (u, v)` from `get_adjacent2vertex(tri, w)`.
"""
@inline delete_adjacent2vertex!(tri::Triangulation, w, uv) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
@inline delete_adjacent2vertex!(tri::Triangulation, w, u, v) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)

"""
    delete_adjacent2vertex!(tri::Triangulation, w)

Calls `delete_adjacent2vertex!(get_adjacent2vertex(tri), w)`, deleting the vertex `w` from `get_adjacent2vertex(tri)`.
"""
@inline delete_adjacent2vertex!(tri::Triangulation, w) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w)
