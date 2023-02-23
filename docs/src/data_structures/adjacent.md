```@meta
CurrentModule = DelaunayTriangulation
```

# Adjacent 

The `Adjacent` map is used for mapping edges to vertices that together form a positively oriented triangle. The definition of the `Adjacent` map is below:

```julia
struct Adjacent{I,E}
    adjacent::DefaultDict{E,I,I}
    function Adjacent{I,E}() where {I,E}
        A = DefaultDict{E,I,I}(I(DefaultAdjacentValue))
        adj = new{I,E}(A)
        return adj
    end
    Adjacent(adj::DefaultDict{E,I,I}) where {I,E} = new{I,E}(adj)
end
```

We use a `DefaultDict` from DataStructures.jl for storing these relationships. We need to use a `DefaultDict` rather than a `Dict` so that we more efficiently check if an edge exists, simply seeing if `get_adjacent(adj, u, v)` returns a `DefaultAdjacentValue` (see the `edge_exists` function). We list the complete docstring for `Adjacent` below, along with individual docstrings for methods for working with `Adjacent`.

```@docs 
Adjacent
get_adjacent(::Adjacent)
get_adjacent(::Adjacent{I, E}, ::E) where {I, E, V}
add_adjacent!(::Adjacent, ::Any, ::Any)
delete_adjacent!(::Adjacent, ::Any)
add_triangle!(::Ts, ::V, ::V, ::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}}
add_triangle!(::Adjacent, ::Any)
delete_triangle!(::Ts, ::V, ::V, ::V) where {I,E,V<:Integer,Ts<:Adjacent{I,E}}
delete_triangle!(::Adjacent, ::Any)
clear_empty_keys!(::Adjacent)
```