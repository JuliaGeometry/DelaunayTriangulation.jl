"""
    Graph{I}

Struct for storing neighbourhood relationships that map vertices 
to all other vertices that share an edge with that vertex. 
The type `I` is the integer type.

See the docs for a description of how boundary edges 
are handled.

# Fields 
- `graph::UndirectedGraph{I}`

The `UndirectedGraph` that maps a vertex `u` to a list of edges,
`V`, such that `(u, v)` is an edge of the triangulation for each 
`v` in `V`. 

# Extended help
You should not work with the `graph` field directly. We provide 
the following functions for working with `Graph`, where `G` denotes 
a `Graph{I}` type. (Type information in the function signatures 
is omitted.)

## Accessors 
- `get_graph(G)`
- `get_vertices(G)`
- `each_vertex(G)`
- `get_edges(G)`
- `get_neighbours(G)`
- `get_neighbours(G, u)`

## Mutators
- `add_vertex!(G, u...)`
- `add_neighbour!(G, u, v...)`
- `add_triangle!(G, i, j, k)` or `add_triangle!(G, T)`
- `add_triangle!(G, T...)`
- `delete_triangle!(G, i, j, k)` or `delete_triangle!(G, T)`
- `delete_triangle!(G, T...)`
- `delete_neighbour!(G, u, v...)`
- `delete_vertex!(G, u...)`
- `delete_boundary_vertices_from_graph!(G)`
- `clear_empty_points!(G)`

## Miscellaneous
- `num_edges(G)`
- `num_neighbours(G, u)`
- `num_vertices(G)`
"""
struct Graph{I}
    graph::UndirectedGraph{I}
    function Graph{I}() where {I}
        G = UndirectedGraph{I}()
        return new{I}(G)
    end
    Graph() = Graph{Int64}()
    Graph(G::UndirectedGraph{I}) where {I} = new{I}(G)
end
Base.:(==)(G::Graph, H::Graph) = get_graph(G) == get_graph(H)
function Base.show(io::IO, ::MIME"text/plain", graph::Graph)
    println(io, "Graph")
    println(io, "    Number of edges: $(num_edges(graph))")
    print(io, "    Number of vertices: $(num_vertices(graph))")
end

"""
    get_graph(G::Graph)

Returns the field `graph` of `G`.
"""
get_graph(G::Graph) = G.graph

"""
    get_edges(G::Graph)

Returns all the edges of the graph `G`. Edges 
are unordered.
"""
get_edges(G::Graph) = get_graph(G).E

"""
    get_vertices(G::Graph)

Given a `graph`, returns the current set of vertices.
"""
get_vertices(G::Graph) = get_graph(G).V

"""
    each_vertex(G::Graph)

Given a `graph`, returns the current set of vertices.
"""
each_vertex(G::Graph) = get_vertices(G)

"""
    get_neighbours(G::Graph)

Returns the set of neighbourhoods of the graph `G`.
"""
get_neighbours(G::Graph) = get_graph(G).N

"""
    get_neighbours(G::Graph, u) 

Returns the neighbourhood of the points `u` of 
the graph `G`.
"""
get_neighbours(G::Graph, u) = get_neighbours(G)[u]

"""
    num_neighbours(G::Graph, u) 

Returns the number of neighbours of the point `u` 
of the graph `G`.
"""
num_neighbours(G::Graph, u) = deg(get_graph(G), u)

"""
    num_edges(G::Graph)

Returns the number of edges `G`. The edges 
`(i, j)` and `(j, i)` will only be counted once.
"""
num_edges(G::Graph) = length(get_graph(G).E)

"""
    num_vertices(G::Graph)

Returns the number of vertices in the graph `G`.
"""
num_vertices(G::Graph) = length(get_vertices(G))

"""
    add_vertex!(G::Graph, u...)

Given a graph `G` and vertices `u...`, adds these 
vertices into `G`.
"""
function add_vertex!(G::Graph{I}, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_vertex!(G, u[i])
    end
    return nothing
end
function add_vertex!(G::Graph{I}, u::I) where {I}
    add!(get_graph(G), u)
    return nothing
end

"""
    add_neighbour!(G::Graph, u, v...)

Given a graph `G`, adds `v...` into the neighbourhood of `u`.
"""
function add_neighbour!(G::Graph{I}, u::I, v::I) where {I}
    add!(get_graph(G), u, v)
    return nothing
end
function add_neighbour!(G::Graph{I}, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_neighbour!(G, u, v[i])
    end
    return nothing
end

"""
    add_triangle!(G::Graph, i, j, k)

Given a graph `G`, adds the triangle `(i, j, k)` into `G`. In particular, the 
indices `(i, j, k)` are added into `G`, and the indices are all in each other's 
neighbourhoods.
"""
function add_triangle!(G::Ts, i::V, j::V, k::V) where {I,V<:Integer,Ts<:Graph{I}} # the signature here is for resolving a method ambiguity
    add_vertex!(G, I(i), I(j), I(k))
    add_neighbour!(G, I(i), I(j), I(k))
    add_neighbour!(G, I(j), I(k)) # Don't need to do all the cases, the graph is undirected
    return nothing
end

"""
    add_triangle!(G::Graph, T...)

Adds the triangles `T...` into the graph `G`.
"""
function add_triangle!(G::Graph, T)
    i, j, k = indices(T)
    add_triangle!(G, i, j, k)
    return nothing
end
function add_triangle!(G::Graph, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        add_triangle!(G, T[i])
    end
    return nothing
end

"""
    delete_triangle!(G::Graph, i, j, k)

Given a graph `G`, deletes the triangle `(i, j, k)` deletes `G`. In particular, the 
indices `(i, j, k)` will no longer be neighbours of each other.

!!! note

    Be careful with using this function - you could have a triangle `(j, i, ℓ)`, say, 
    which will also be affected since the graph is undirected. Note also 
    that the vertices `(i, j, k)` will not be removed - only the neighbourhoods 
    are affected.
"""
function delete_triangle!(G::Ts, i::V, j::V, k::V) where {I,V<:Integer,Ts<:Graph{I}} # the signature here is for resolving a method ambiguity
    delete_neighbour!(G, I(i), I(j), I(k))
    delete_neighbour!(G, I(j), I(k)) # Don't need to do all the cases, the graph is undirected
    return nothing
end

"""
    delete_triangle!(G::Graph, T...)

Deletes the triangles `T...` from the graph `G`.
"""
function delete_triangle!(G::Graph, T)
    i, j, k = indices(T)
    delete_triangle!(G, i, j, k)
    return nothing
end
function delete_triangle!(G::Graph, T::Vararg{V,N}) where {V,N}
    for i in 1:N
        delete_triangle!(G, T[i])
    end
    return nothing
end

"""
    delete_neighbour!(G::Graph, u, v...)

Given a graph `G` and a vertex `u`, deletes the vertices `v...`
from the neighbourhood of `u`.
"""
function delete_neighbour!(G::Graph, u, v)
    delete!(get_graph(G), u, v)
    return nothing
end
function delete_neighbour!(G::Graph, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_neighbour!(G, u, v[i])
    end
    return nothing
end

"""
    delete_vertex!(G::Graph, u...)

Given a graph `G` and vertices `u...`, deletes the vertices from `G`.
"""
function delete_vertex!(G::Graph, u)
    delete!(get_graph(G), u)
    return nothing
end
function delete_vertex!(G::Graph, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_vertex!(G, u[i])
    end
    return nothing
end

"""
    delete_boundary_vertices_from_graph!(G::Graph{I}) where {I}

Given a graph `G`, deletes all the boundary indices from `G`, i.e. 
all `u` such that `u ≤ $BoundaryIndex`.
"""
function delete_boundary_vertices_from_graph!(G::Graph{I}) where {I}
    smallest_idx = minimum(get_vertices(G))
    for i in I(BoundaryIndex):-1:smallest_idx
        delete_vertex!(G, i)
    end
    return nothing
end

"""
    clear_empty_points!(G::Graph)

Given a `graph`, deletes any points that have empty neighbourhoods.
"""
function clear_empty_points!(G::Graph)
    for (w, S) in get_graph(G).N
        isempty(S) && delete_vertex!(G, w)
    end
    return nothing
end

function Base.sizehint!(G::Graph, ne, nn, nv)
    E = get_edges(G)
    N = get_neighbours(G)
    V = get_vertices(G)
    sizehint!(E, ne)
    sizehint!(N, nn)
    sizehint!(V, nv)
    return nothing
end
