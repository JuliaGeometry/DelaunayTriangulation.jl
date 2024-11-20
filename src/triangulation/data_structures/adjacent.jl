struct Adjacent{I,E} <: AbstractDict{E,I}
    map::Dict{E,I}
end
@inline Adjacent{I,E}() where {I,E} = Adjacent{I,E}(Dict{E,I}())

@inline get_adjacent(adj::Adjacent) = adj.adjacent

@inline Base.sizehint!(adj::Adjacent, n) = sizehint!(get_adjacent(adj), n)
@inline Base.iterate(adj::Adjacent, state...) = iterate(get_adjacent(adj), state...)
@inline Base.length(adj::Adjacent) = length(get_adjacent(adj))

@inline edge_type(::Adjacent{I,E}) where {I,E} = E
@inline integer_type(::Adjacent{I,E}) where {I,E} = I

@inline function get_adjacent(adj::Adjacent, uv::E)
    I = integer_type(adj)
    map = get_adjacent(adj)
    return Base.get(map, uv, I(∅))
end
@inline function get_adjacent(adj::Adjacent, u, v)
    E = edge_type(adj)
    e = construct_edge(E, u, v)
    return get_adjacent(adj, e)
end

@inline function add_adjacent!(adj::Adjacent, uv, w)
    map = get_adjacent(adj)
    map[uv] = w
    return adj
end
@inline function add_adjacent!(adj::Adjacent, u, v, w)
    E = edge_type(adj)
    e = construct_edge(E, u, v)
    return add_adjacent!(adj, e, w)
end

@inline function delete_adjacent!(adj::Adjacent, uv)
    map = get_adjacent(adj)
    delete!(map, uv)
    return adj
end
@inline function delete_adjacent!(adj::Adjacent, u, v)
    E = edge_type(adj)
    e = construct_edge(E, u, v)
    return delete_adjacent!(adj, e)
end

@inline function add_triangle!(adj::Adjacent, u, v, w)
    add_adjacent!(adj, u, v, w)
    add_adjacent!(adj, v, w, u)
    add_adjacent!(adj, w, u, v)
    return adj
end
@inline add_triangle!(adj::Adjacent, T) = add_triangle!(adj, triangle_vertices(T)...)

@inline function delete_triangle!(adj::Adjacent, u, v, w)
    delete_adjacent!(adj, u, v)
    delete_adjacent!(adj, v, w)
    delete_adjacent!(adj, w, u)
    return adj
end
@inline delete_triangle!(adj::Adjacent, T) = delete_triangle!(adj, triangle_vertices(T)...)

@inline function Base.empty!(adj::Adjacent)
    map = get_adjacent(adj)
    empty!(map)
    return adj
end

@inline edge_exists(w) = w ≠ ∅
@inline edge_exists(adj, uv) = edge_exists(get_adjacent(adj, uv))
@inline edge_exists(adj, u, v) = edge_exists(get_adjacent(adj, u, v))
