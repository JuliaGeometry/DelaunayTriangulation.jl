@inline function get_adjacent(tri::Triangulation, uv)
    adj = get_adjacent(tri)
    check_existence = has_multiple_sections(tri)
    if !check_existence || !is_ghost_edge(uv)
        return get_adjacent(adj, uv)
    else
        return _safe_get_adjacent(tri, uv)
    end
end
@inline get_adjacent(tri::Triangulation, u, v) = get_adjacent(tri, construct_edge(edge_type(tri), u, v))

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

@inline add_adjacent!(tri::Triangulation, uv, w) = add_adjacent!(get_adjacent(tri), uv, w)
@inline add_adjacent!(tri::Triangulation, u, v, w) = add_adjacent!(get_adjacent(tri), u, v, w)

@inline delete_adjacent!(tri::Triangulation, uv) = delete_adjacent!(get_adjacent(tri), uv)
@inline delete_adjacent!(tri::Triangulation, u, v) = delete_adjacent!(get_adjacent(tri), u, v)

@inline add_candidate!(tri::Triangulation, u, v) = add_candidate!(get_adjacent2vertex_candidates(tri), u, v)

@inline delete_candidate!(tri::Triangulation, u) = delete_candidate!(get_adjacent2vertex_candidates(tri), u)

@inline get_candidate(tri::Triangulation, u) = get_candidate(get_adjacent(tri), u)