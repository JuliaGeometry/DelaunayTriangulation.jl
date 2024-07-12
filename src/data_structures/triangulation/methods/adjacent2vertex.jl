"""
    get_adjacent2vertex(tri::Triangulation, w) -> Edges

Returns the set of all edges `(u, v)` in `tri` such that `(u, v, w)` is a positively oriented triangle in `tri`.
"""
get_adjacent2vertex(tri::Triangulation, w) = get_adjacent2vertex(get_adjacent2vertex(tri), w)


"""
    add_adjacent2vertex!(tri::Triangulation, w, uv)
    add_adjacent2vertex!(tri::Triangulation, w, u, v)

Adds the edge `(u, v)` into the set of edges returned by `get_adjacent2vertex(tri, w)`.
"""
add_adjacent2vertex!(tri::Triangulation, w, uv) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
add_adjacent2vertex!(tri::Triangulation, w, u, v) = add_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)


"""
    delete_adjacent2vertex!(tri::Triangulation, w, uv)
    delete_adjacent2vertex!(tri::Triangulation, w, u, v)

Deletes the edge `(u, v)` from the set of edges returned by `get_adjacent2vertex(tri, w)`.
"""
delete_adjacent2vertex!(tri::Triangulation, w, uv) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, uv)
delete_adjacent2vertex!(tri::Triangulation, w, u, v) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w, u, v)


"""
    delete_adjacent2vertex!(tri::Triangulation, w)

Deletes the key `w` from the [`Adjacent2Vertex`](@ref) map of `tri`.
"""
delete_adjacent2vertex!(tri::Triangulation, w) = delete_adjacent2vertex!(get_adjacent2vertex(tri), w)