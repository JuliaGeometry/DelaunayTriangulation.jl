"""
    get_adjacent(tri::Triangulation, uv) -> Vertex
    get_adjacent(tri::Triangulation, u, v) -> Vertex

Returns the vertex `w` such that `(u, v, w)` is a positively oriented triangle in the
underlying triangulation, or `∅` if no such triangle exists.
"""
function get_adjacent(tri::Triangulation, uv)
    adj = get_adjacent(tri)
    check_existence = has_multiple_sections(tri)
    if !check_existence || !is_ghost_edge(uv)
        return get_adjacent(adj, uv)
    else
        return _safe_get_adjacent(tri, uv)
    end
end
get_adjacent(tri::Triangulation, u, v) = get_adjacent(tri, construct_edge(edge_type(tri), u, v))

"""
    _safe_get_adjacent(tri::Triangulation, uv) -> Vertex

This is the safe version of `get_adjacent`, which is used when the triangulation has multiple
sections, ensuring that the correct ghost vertex is returned in case `uv` is a ghost edge.
"""
function _safe_get_adjacent(tri::Triangulation, uv)
    adj = get_adjacent(tri)
    u, v = edge_vertices(uv)
    w = get_adjacent(adj, u, v)
    if !edge_exists(w)
        ℓ = get_ghost_vertex(u, v)
        range = get_ghost_vertex_range(tri, ℓ)
        if ℓ == u
            for i in range
                w = get_adjacent(adj, i, v)
                edge_exists(w) && return w
            end
        else
            for j in range
                w = get_adjacent(adj, u, j)
                edge_exists(w) && return w
            end
        end
    end
    return get_adjacent(adj, uv)
end

"""
    add_adjacent!(tri::Triangulation, uv, w)
    add_adjacent!(tri::Triangulation, u, v, w)

Adds the key-value pair `(u, v) ⟹ w` to the adjacency map of `tri`.
"""
add_adjacent!(tri::Triangulation, uv, w) = add_adjacent!(get_adjacent(tri), uv, w)
add_adjacent!(tri::Triangulation, u, v, w) = add_adjacent!(get_adjacent(tri), u, v, w)

"""
    delete_adjacent!(tri::Triangulation, uv)
    delete_adjacent!(tri::Triangulation, u, v)

Deletes the key `(u, v)` from the adjacency map of `tri`.
"""
delete_adjacent!(tri::Triangulation, uv) = delete_adjacent!(get_adjacent(tri), uv)
delete_adjacent!(tri::Triangulation, u, v) = delete_adjacent!(get_adjacent(tri), u, v)