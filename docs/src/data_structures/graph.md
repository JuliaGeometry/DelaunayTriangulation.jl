```@meta
CurrentModule = DelaunayTriangulation
```

# Graph 

To make it possible to iterate over all points that share an edge with a given vertex, we have a `Graph` struct that we define as follows:

```julia
struct Graph{I}
    graph::UndirectedGraph{I}
    function Graph{I}() where {I}
        G = UndirectedGraph{I}()
        return new{I}(G)
    end
    Graph() = Graph{Int}()
    Graph(G::UndirectedGraph{I}) where {I} = new{I}(G)
end
```

Note that this graph is undirected. 

```@docs 
Graph 
get_graph(::Graph)
get_edges(::Graph)
get_vertices(::Graph)
each_vertex(::Graph)
get_neighbours(::Graph)
get_neighbours(::Graph, ::Any)
num_neighbours(::Graph, ::Any)
num_edges(::Graph)
num_vertices(::Graph)
add_vertex!(::Graph{I}, ::Vararg{I,N}) where {I,N}
add_neighbour!(::Graph{I}, ::I, ::I) where {I}
add_triangle!(::Ts, ::V, ::V, ::V) where {I,V<:Integer,Ts<:Graph{I}}
add_triangle!(::Graph, ::Any)
delete_triangle!(::Ts, ::V, ::V, ::V) where {I,V<:Integer,Ts<:Graph{I}}
delete_triangle!(::Graph, ::Any)
delete_neighbour!(::Graph, ::Any, ::Any)
delete_vertex!(::Graph, ::Any)
delete_boundary_vertices_from_graph!(::Graph{I}) where {I}
clear_empty_points!(::Graph)
```