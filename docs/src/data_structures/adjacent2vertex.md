```@meta
CurrentModule = DelaunayTriangulation
```

# Adjacent2Vertex 

The `Adjacent2Vertex` map is closely related to the `Adjacent` map. Instead of mapping edges to vertices that together form positively oriented triangles, we map vertices to all edges that will form a positively oriented triangle with that vertex. The definition is simply via a `Dict`:

```julia
struct Adjacent2Vertex{I,Es,E}
    adjacent2vertex::Dict{I,Es}
    function Adjacent2Vertex{I,Es,E}() where {I,Es,E}
        D = Dict{I,Es}()
        TA2V = new{I,Es,E}(D)
        return TA2V
    end
    Adjacent2Vertex(adj2v::Dict{I,Es}) where {I,Es} = new{I,Es,edge_type(Es)}(adj2v)
end
```

A complete list of methods used for working with this struct, along with the struct's docstring itself, is shown below.

```@docs 
Adjacent2Vertex 
get_adjacent2vertex(::Adjacent2Vertex)
get_adjacent2vertex(::Adjacent2Vertex, ::Any)
add_adjacent2vertex!(::Adjacent2Vertex{I,Es,E}, ::Any, ::Any) where {I,Es,E}
delete_adjacent2vertex!(::Adjacent2Vertex, ::Any, ::Any)
delete_adjacent2vertex!(::Adjacent2Vertex, ::Any)
add_triangle!(::Ts, ::V, ::V, ::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}}
add_triangle!(::Adjacent2Vertex, ::Any)
delete_triangle!(::Ts, ::V, ::V, ::V) where {I,Es,E,V<:Integer,Ts<:Adjacent2Vertex{I,Es,E}}
delete_triangle!(::Adjacent2Vertex, ::Any)
clear_empty_keys!(::Adjacent2Vertex)
```