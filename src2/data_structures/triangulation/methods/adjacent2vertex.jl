@doc """
    get_adjacent2vertex(tri::Triangulation, w) -> Iterator

Returns the set of all edges `(u, v)` in `tri` such that `(u, v, w)` is a positively oriented triangle in `tri`.
"""
get_adjacent2vertex(::Triangulation, ::Any)

concretize_adjacent2vertex(tri::Triangulation) = concretize_adjacent2vertex(get_adjacent2vertex(tri))

function Adjacent2VertexIterator(tri::Triangulation, u, v)
    I = integer_type(tri)
    E = edge_type(tri)
    return Adjacent2VertexIterator{typeof(tri), E, I}(tri, u, v)
end

add_candidates!(tri::Triangulation, u, v, w) = add_candidates!(get_adjacent2vertex(tri), u, v, w)

get_candidate(tri::Triangulation, u) = get_candidate(get_adjacent2vertex(tri), u)