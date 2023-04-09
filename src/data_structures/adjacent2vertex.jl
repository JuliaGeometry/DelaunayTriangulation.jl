"""
    Adjacent2Vertex{I, Es, E}

Struct for storing adjacency relationships for mapping vertices to 
all other edges that together form a positively oriented triangle 
in an associated triangulation. The type `I` is the integer type 
used, `Es` is the type used for representing a collection of edges, 
and `E` is the type for a single edge.

See the docs for a description of how boundary edges 
are handled.

See also [`Adjacent`](@ref).

# Fields 
- `adjacent2vertex::Dict{I, Es}`

The `Dict` used for storing the vertices (the keys) and the associated 
collection of edges (the values). In particular, if `w` is a vertex, 
then `(u, v, w)` is a positively oriented triangle for each 
`(u, v)` in `adjacent2vertex[w]`.

# Constructors 
The adjacent2vertex map can be constructed in two ways:

- `Adjacent2Vertex{I,Es,E}()`

Creates an empty map.
- `Adjacent2Vertex(adj2v::Dict{I,E}) where {I,Es}`

Creates an adjacent2vertex map from an existing `Dict`. The type `E` is obtained from 
[`edge_type`](@ref)`(Es)`.

# Extended help 
You should not work with the `adjacent2vertex` field directory. We provide the 
following functions for working with `Adjacent2Vertex`, where `adj2v` 
denotes an `Adjacent2Vertex{I,Es,E}` type. (Type information in the 
function signatures is omitted.)

## Accessors 

- `get_adjacent2vertex(adj2v)`
- `get_adjacent2vertex(adj2v, w)`

## Mutators 

- `add_adjacent2vertex!(adj2v, w, uv)` or `add_adjacent2vertex!(adj2v, w, u, v)`
- `delete_adjacent2vertex!(adj2v, w, uv)` or `delete_adjacent2vertex!(adj2v, w, u, v)`
- `delete_adjacent2vertex!(adj2v, w)`
- `add_triangle!(adj2v, i, j, k)` or `add_triangle!(adj2v, T)`
- `add_triangle!(adj2v, T...)`
- `delete_triangle!(adj2v, i, j, k)` or `delete_triangle!(adj2v, T)`
- `delete_triangle!(adj2v, T...)`
- `clear_empty_keys!(adj2v)`

## Iteration 
You can also iterate over `Adjacent2Vertex` maps the same way as you would 
with a `Dict`, e.g. if `adj` is a corresponding [`Adjacent`](@ref) map,

```julia
for (vertex, edge_list) in adj2v 
    for edge in each_edge(edge_list)
        get_adjacent(adj, edge) == vertex 
    end
end
```
"""
struct Adjacent2Vertex{I,Es,E}
    adjacent2vertex::Dict{I,Es}
    function Adjacent2Vertex{I,Es,E}() where {I,Es,E}
        D = Dict{I,Es}()
        TA2V = new{I,Es,E}(D)
        return TA2V
    end
    Adjacent2Vertex(adj2v::Dict{I,Es}) where {I,Es} = new{I,Es,edge_type(Es)}(adj2v)
end
function Base.:(==)(adj2v::Adjacent2Vertex, adj2v2::Adjacent2Vertex)
    return get_adjacent2vertex(adj2v) == get_adjacent2vertex(adj2v2)
end
function Base.show(io::IO, m::MIME"text/plain", adj2v::Adjacent2Vertex{I,Es,E}) where {I,Es,E}
    println(io, "Adjacent{$I, $Es, $E}, with map:")
    show(io, m, get_adjacent2vertex(adj2v))
end

"""
    get_adjacent2vertex(adj2v::Adjacent2Vertex)

Returns the `adjacent2vertex` field from the adjacent2vertex map `adj2v`.
"""
get_adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    get_adjacent2vertex(adj2v::Adjacent2Vertex, w)

Given the adjacent2vertex map `adj2v` and a vertex `w`, returns the set of 
edges associated with the vertex `w`, i.e. the set of edges `(u, v)` such that 
`(u, v, w)` is a positively oriented triangle in the underlying triangulation.
"""
get_adjacent2vertex(adj2v::Adjacent2Vertex, w) = get_adjacent2vertex(adj2v)[w]

"""
    add_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    add_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}

Given the adjacent2vertex map `adj2v`, a vertex `w`, and an edge `(u, v)`, adds 
the edge `(u, v)` into the set of edges associated with the vertex `w` in the map.
"""
function add_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    existing_edges = get!(adj2v, w)
    add_edge!(existing_edges, uv)
    return nothing
end
function add_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    add_adjacent2vertex!(adj2v, w, uv)
    return nothing
end

"""
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w, uv)
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}

Given the adjacent2vertex map `adj2v`, a vertex `w`, and an edge `(u, v)`, deletes 
the edge `(u, v)` from the set of edges associated with the vertex `w` in the map.
"""
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w, uv)
    existing_edges = get_adjacent2vertex(adj2v, w)
    delete_edge!(existing_edges, uv)
    return nothing
end
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    delete_adjacent2vertex!(adj2v, w, uv)
    return nothing
end

"""
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w)

Given the adjacent2vertex map `adj2v` and a vertex `w`, deletes 
the key `w` from the map.
"""
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w)
    delete!(get_adjacent2vertex(adj2v), w)
    return nothing
end

"""
    add_triangle!(adj2v::Adjacent2Vertex, i, j, k)

Given an adjacent2vertex map `adj2v` and indices `(i, j, k)` 
representing some triangle, adds that triangle from the 
map. In particular, adds the edges `(i, j)`, `(j, k)`, and `(k, i)` 
into the set of edges associated with the vertices `k`, `i`, and `j`,
respectively, in the map.
"""
function add_triangle!(adj2v::Ts, i::V, j::V, k::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}} # the signature here is for resolving a method ambiguity
    for (u, v, w) in ((i, j, k), (j, k, i), (k, i, j))
        add_adjacent2vertex!(adj2v, I(w), I(u), I(v))
    end
    return nothing
end

"""
    add_triangle!(adj::Adjacent2Vertex, T...)

Given an adjacent map `adj2v` and triangles `T...`, adds the 
triangles into `adj2v`. See also [`add_triangle!`](@ref add_triangle!(::Ts, ::V, ::V, ::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}}).
"""
function add_triangle!(adj2v::Adjacent2Vertex, T)
    i, j, k = indices(T)
    add_triangle!(adj2v, i, j, k)
    return nothing
end
function add_triangle!(adj2v::Adjacent2Vertex, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        add_triangle!(adj2v, T[i])
    end
    return nothing
end

"""
    delete_triangle!(adj2v::Adjacent2Vertex, i, j, k)

Given an adjacent2vertex map `adj2v` and indices `(i, j, k)` 
representing some triangle, deletes that triangle from the 
map. In particular, deletes the edges `(i, j)`, `(j, k)`, and `(k, i)` 
from the set of edges associated with the vertices `k`, `i`, and `j`,
respectively, from the map.
"""
function delete_triangle!(adj2v::Ts, i::V, j::V, k::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}} # the signature here is for resolving a method ambiguity
    for (u, v, w) in ((i, j, k), (j, k, i), (k, i, j))
        delete_adjacent2vertex!(adj2v, w, u, v)
    end
    return nothing
end

"""
    delete_triangle!(adj::Adjacent2Vertex, T...)

Given an adjacent2vertex map `adj2v` and triangles `T...`, deletes the 
triangles from `adj2v`. See also [`delete_triangle!`](@ref delete_triangle!(::Ts, ::V, ::V, ::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}}).
"""
function delete_triangle!(adj2v::Adjacent2Vertex, T)
    i, j, k = indices(T)
    delete_triangle!(adj2v, i, j, k)
    return nothing
end
function delete_triangle!(adj2v::Adjacent2Vertex, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        delete_triangle!(adj2v, T[i])
    end
    return nothing
end

function Base.iterate(adj2v::Adjacent2Vertex, state...)
    return Base.iterate(get_adjacent2vertex(adj2v), state...)
end
Base.length(adj2v::Adjacent2Vertex) = Base.length(get_adjacent2vertex(adj2v))
Base.eltype(::Type{Adjacent2Vertex{I,Es,E}}) where {I,Es,E} = Pair{I,Es}

function Base.get!(adj2v::Adjacent2Vertex{I,Es,E}, w) where {I,Es,E}
    dict = get_adjacent2vertex(adj2v)
    index = Base.ht_keyindex2!(dict, w)
    index > 0 && return dict.vals[index]
    age0 = dict.age
    v = initialise_edges(Es) # this is the only line that changed from get!
    if dict.age != age0
        index = Base.ht_keyindex2!(dict, w)
    end
    if index > 0
        dict.age += 1
        @inbounds dict.keys[index] = w
        @inbounds dict.vals[index] = v
    else
        @inbounds Base._setindex!(dict, v, w, -index)
    end
    return v
end

"""
    clear_empty_keys!(adj2v::Adjacent2Vertex)

Given an [`Adjacent2Vertex`](@ref) map `adj2v`, removes 
any keys that map to empty sets.
"""
function clear_empty_keys!(adj2v::Adjacent2Vertex)
    for (w, S) in adj2v
        is_empty(S) && delete_adjacent2vertex!(adj2v, w)
    end
    return nothing
end

Base.sizehint!(adj2v::Adjacent2Vertex, n) = Base.sizehint!(get_adjacent2vertex(adj2v), n)
