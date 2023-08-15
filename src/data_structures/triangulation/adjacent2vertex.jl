"""
    get_adjacent2vertex(tri::Triangulation, w)

Returns all edges `(u, v)` in `tri` such that `(u, v, w)` is a positively oriented triangle in `tri`.
"""
@inline get_adjacent2vertex(tri::Triangulation, w) = get_adjacent2vertex(get_adjacent2vertex(tri), w)

"""
    add_adjacent2vertex!(tri::Triangulation, w, uv)
    add_adjacent2vertex!(tri::Triangulation, w, u, v)

Adds the edge `(u, v)` to the list of edges adjacent to `w`, meaning that `(u, v, w)` is a positively oriented triangle in `tri`.
"""
@inline add_adjacent2vertex!(tri::Triangulation, w, uv) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
@inline add_adjacent2vertex!(tri::Triangulation, w, u, v) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)

"""
    delete_adjacent2vertex!(tri::Triangulation, w, uv)
    delete_adjacent2vertex!(tri::Triangulation, w, u, v)

Deletes the edge `(u, v)` from the list of edges adjacent to `w`, meaning that `(u, v, w)` is no longer a positively oriented triangle in `tri`.
"""
@inline delete_adjacent2vertex!(tri::Triangulation, w, uv) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
@inline delete_adjacent2vertex!(tri::Triangulation, w, u, v) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)

"""
    delete_adjacent2vertex!(tri::Triangulation, w)

Deletes the vertex `w` from the `adjacent2vertex` map of `tri`.
"""
@inline delete_adjacent2vertex!(tri::Triangulation, w) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w)
