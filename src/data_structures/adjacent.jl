"""
    Adjacent{I,E}

Struct for storing adjacency relationships for mapping edges to vertices that 
together form a positively oriented triangle in an associated triangulation. 
The type `I` is the integer type used, while `E` is the edge type.

See the docs for a description of how boundary edges 
are handled.

See also [`Adjacent2Vertex`](@ref).

# Fields 
- `adjacent::DefaultDict{E,I,I}`

The `Dict` used for storing the edges (the keys) and the associated vertices 
(the values). If `(u, v)` is not a valid edge, then `w = adjacent[(u, v)]`
returns `$DefaultAdjacentValue` (this value is defined in `DefaultAdjacentValue`).
Otherwise, `(u, v, w)` is a positively oriented triangle.

# Constructors 
The adjacent map can be constructed in two ways:

- `Adjacent{I, E}() where {I, E}`

Creates an empty map.
- `Adjacent(adj::DefaultDict{E,I,I}) where {E,I,I}`

Creates an adjacent map from an existing `DefaultDict`.

# Extended help 
You should not work with the `adjacent` field directly. We provide the following 
functions for working with `Adjacent`, where `adj` denotes an `Adjacent{I, E}` type.
(Type information in the function signatures is omitted.)

## Accessors

- `get_adjacent(adj)`
- `get_adjacent(adj, uv)` or `get_adjacent(adj, u, v)`

In the latter methods, you can also use the keyword argument `check_existence` to 
declare whether to check that the edge exists. This would be used if you need 
to be careful about different boundary indices on the same boundary curve. The 
default value is `Val(false)`, meaning this isn't checked.

## Mutators

- `add_adjacent!(adj, uv, w)` or `add_adjacent!(adj, u, v, w)`
- `delete_adjacent!(adj, uv)` or `delete_adjacent!(adj, u, v)`
- `add_triangle!(adj, i, j, k)` or `add_triangle!(adj, T)`
- `add_triangle!(adj, T...)`
- `delete_triangle!(adj, i, j, k)` or `delete_triangle!(adj, T)`
- `delete_triangle!(adj, T...)`
- `clear_empty_keys!(adj)`

## Iteration 
You can also iterate over `Adjacent` maps the same way as you would 
with a `Dict`, e.g.

```julia
for (edge, vertex) in adj 
    get_adjacent(adj, edge) == vertex 
end
```
"""
struct Adjacent{I,E}
    adjacent::DefaultDict{E,I,I}
    function Adjacent{I,E}() where {I,E}
        A = DefaultDict{E,I,I}(I(DefaultAdjacentValue))
        adj = new{I,E}(A)
        return adj
    end
    Adjacent(adj::DefaultDict{E,I,I}) where {I,E} = new{I,E}(adj)
end
Base.:(==)(adj::Adjacent, adj2::Adjacent) = get_adjacent(adj) == get_adjacent(adj2)
function Base.show(io::IO, m::MIME"text/plain", adj::Adjacent{I, E}) where {I, E}
    println(io, "Adjacent{$I, $E}, with map:")
    show(io,m,get_adjacent(adj))
end

"""
    get_adjacent(adj::Adjacent)

Given the adjacent map `adj`, returns the `adjacent` field.
"""
get_adjacent(adj::Adjacent) = adj.adjacent

"""
    get_adjacent(adj::Adjacent{I,E}, uv::E; check_existence::V = Val(false), boundary_index_ranges=nothing) where {I,E,V}
    get_adjacent(adj::Adjacent{I,E}, u, v; check_existence::V = Val(false), boundary_index_ranges=nothing) where {I,E,V}

Given the adjacent map `adj` and an edge `(u, v)`, returns the vertex `w`
such that `(u, v, w)` is a positively oriented triangle in the underlying 
triangulation.

In the case of a ghost edge, `check_existence = Val(true)` may be useful in case the 
boundary curve has multiple segments, meaning multiple boundary indices correspond 
to the same curve. Assuming the boundary index came from a neighbouring triangle, 
setting `check_existence = Val(true)` will check neighbouring boundary indices 
in case the found edge does not exist. When `is_true(check_existence)`, you also need
to provide the range of boundary indices to check via the keyword argument `boundary_index_ranges`. This should be 
the `Dict` from [`construct_boundary_index_ranges`](@ref).
"""
function get_adjacent(adj::Adjacent{I,E}, uv::E; check_existence::V=Val(false),
                      boundary_index_ranges=nothing) where {I,E,V}
    return (!is_true(check_existence) || !is_ghost_edge(uv)) ? _get_adjacent(adj, uv) :
           _safe_get_adjacent(adj, uv, boundary_index_ranges)
end
function get_adjacent(adj::Adjacent{I,E}, u, v; check_existence=Val(false),
                      boundary_index_ranges=nothing) where {I,E}
    return get_adjacent(adj, construct_edge(E, u, v); check_existence,
                        boundary_index_ranges)
end

@inline _get_adjacent(adj::Adjacent{I,E}, uv::E) where {I,E} = get_adjacent(adj)[uv]
@inline function _safe_get_adjacent(adj::Adjacent{I,E}, uv::E,
                                    boundary_index_ranges=nothing) where {I,E}
    u = initial(uv)
    v = terminal(uv)
    if !edge_exists(u, v, adj)
        ℓ = get_boundary_index(u, v)
        range = map_boundary_index(boundary_index_ranges, ℓ)
        if ℓ == u
            for i in range
                edge_exists(i, v, adj) && return get_adjacent(adj, i, v)
            end
        else
            for j in range
                edge_exists(u, j, adj) && return get_adjacent(adj, u, j)
            end
        end
    end
    return get_adjacent(adj, u, v)
end

"""
    add_adjacent!(adj::Adjacent, uv, w)
    add_adjacent!(adj::Adjacent{I,E}, u, v, w) where {I,E}

Given the adjacent map `adj`, an edge `(u, v)`, and a vertex `w`, adds 
the edge `(u, v)` with corresponding value `w` into the adjacent map 
so that `(u, v, w)` is a positively oriented triangle in the 
underlying triangulation.
"""
function add_adjacent!(adj::Adjacent, uv, w)
    get_adjacent(adj)[uv] = w
    return nothing
end
function add_adjacent!(adj::Adjacent{I,E}, u, v, w) where {I,E}
    uv = construct_edge(E, u, v)
    add_adjacent!(adj, uv, w)
    return nothing
end

"""
    delete_adjacent!(adj::Adjacent, uv)
    delete_adjacent!(adj::Adjacent{I,E}, u, v) where {I,E}

Given the adjacent map `adj` and an edge `(u, v)`, deletes the 
edge `(u, v)` from the adjacent map.
"""
function delete_adjacent!(adj::Adjacent, uv)
    delete!(get_adjacent(adj), uv)
    return nothing
end
function delete_adjacent!(adj::Adjacent{I,E}, u, v) where {I,E}
    uv = construct_edge(E, u, v)
    delete!(get_adjacent(adj), uv)
    return nothing
end

"""
    add_triangle!(adj::Adjacent, i, j, k)

Given an adjacent map `adj` and indices `(i, j, k)` representing some triangle, 
adds that triangle into the adjacent map. In particular, adds the edges `(i, j)`,
`(j, k)`, and `(k, i)` into `adj` with corresponding values `k`, `i`, and `j`, 
respectively.
"""
function add_triangle!(adj::Ts, i::V, j::V,
                       k::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}} # the signature here is just for resolving a method ambiguity
    for (u, v, w) in ((i, j, k), (j, k, i), (k, i, j))
        add_adjacent!(adj, I(u), I(v), I(w))
    end
    return nothing
end

"""
    add_triangle!(adj::Adjacent, T...)

Given an adjacent map `adj` and triangles `T...`, adds the 
triangles into `adj`. See also [`add_triangle!(::Ts, ::V, ::V, ::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}}`](@ref).
"""
function add_triangle!(adj::Adjacent, T)
    i, j, k = indices(T)
    add_triangle!(adj, i, j, k)
    return nothing
end
function add_triangle!(adj::Adjacent, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        add_triangle!(adj, T[i])
    end
    return nothing
end

"""
    delete_triangle!(adj::Adjacent, i, j, k)

Given an adjacent map `adj` and indices `(i, j, k)` representing some triangle, 
deletes that triangle into the adjacent map. In particular, deletes the edges `(i, j)`,
`(j, k)`, and `(k, i)` from `adj`.
"""
function delete_triangle!(adj::Ts, i::V, j::V,
                          k::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}}
    for (u, v) in triangle_edges(i, j, k)
        delete_adjacent!(adj, u, v)
    end
    return nothing
end

"""
    delete_triangle!(adj::Adjacent, T...)

Given an adjacent map `adj` and triangles `T...`, deletes the 
triangles from `adj`. See also [`delete_triangle!(::Ts, ::V, ::V, ::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}}`](@ref).
"""
function delete_triangle!(adj::Adjacent, T)
    i, j, k = indices(T)
    delete_triangle!(adj, i, j, k)
    return nothing
end
function delete_triangle!(adj::Adjacent, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        delete_triangle!(adj, T[i])
    end
    return nothing
end

Base.iterate(adj::Adjacent, state...) = Base.iterate(get_adjacent(adj), state...)
Base.length(adj::Adjacent) = Base.length(get_adjacent(adj))
Base.eltype(adj::Type{Adjacent{I,E}}) where {I,E} = Pair{E,I}

"""
    clear_empty_keys!(adj::Adjacent) 

Given an [`Adjacent`](@ref) map `adj`, removes any edges that 
map to $DefaultAdjacentValue`.
"""
function clear_empty_keys!(adj::Adjacent)
    for (uv, w) in adj
        !edge_exists(w) && delete_adjacent!(adj, uv)
    end
    return nothing
end

Base.sizehint!(adj::Adjacent, n) = Base.sizehint!(get_adjacent(adj), n)
