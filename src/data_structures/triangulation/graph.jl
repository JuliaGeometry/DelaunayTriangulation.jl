abstract type AbstractGraph end

"""
    Graph{IntegerType}

Struct for storing neighbourhood relationships between vertices in a triangulation. This is an undirected graph.

# Fields 
- `vertices::Set{IntegerType}`

The set of vertices in the underlying triangulation.
- `edges::Set{NTuple{2, IntegerType}}`

The set of edges in the underlying triangulation.
- `neighbours::Dict{IntegerType, Set{IntegerType}}`

The map taking vertices `u` to the set of all `v` such that `(u, v)` is an edge in the underlying triangulation.

# Constructors 
    Graph{IntegerType}()
    Graph(vertices::Set{IntegerType}, edges::Set{NTuple{2, IntegerType}}, neighbours::Dict{IntegerType, Set{IntegerType}})
"""
struct Graph{I} <: AbstractGraph
    vertices::Set{I}
    edges::Set{NTuple{2,I}}
    neighbours::Dict{I,Set{I}}
end
Graph{I}() where {I} = Graph(Set{I}(), Set{NTuple{2,I}}(), Dict{I,Set{I}}())
function Base.show(io::IO, ::MIME"text/plain", graph::Graph)
    println(io, "Graph")
    println(io, "    Number of edges: ", num_edges(graph))
    print(io, "    Number of vertices: ", num_vertices(graph))
end
function Base.:(==)(g::Graph, h::Graph)
    g_v = get_vertices(g)
    h_v = get_vertices(h)
    g_e = get_edges(g)
    h_e = get_edges(h)
    g_n = get_neighbours(g)
    h_n = get_neighbours(h)
    length(g_v) == length(h_v) && length(g_e) == length(h_e) && length(g_n) == length(h_n) || return false
    g_v == h_v || return false
    g_n == h_n || return false
    for e in g_e
        (e ∈ h_e || reverse_edge(e) ∈ h_e) || return false
    end
    return true
end
function Base.sizehint!(graph::Graph, n1, n2, n3)
    V = get_vertices(graph)
    E = get_edges(graph)
    N = get_neighbours(graph)
    sizehint!(V, n1)
    sizehint!(E, n2)
    sizehint!(N, n3)
    return graph
end

function Base.copy(graph::Graph)
    V = get_vertices(graph)
    E = get_edges(graph)
    N = get_neighbours(graph)
    return Graph(copy(V), copy(E), copy(N))
end

"""
    get_vertices(graph::Graph) -> Set{Vertex}

Returns the set of vertices in `graph`.
"""
get_vertices(graph::Graph) = graph.vertices

"""
    get_edges(graph::Graph) -> Set{NTuple{2, Vertex}}

Returns the set of edges in `graph`.
"""
get_edges(graph::Graph) = graph.edges

"""
    get_neighbours(graph::Graph) -> Dict{Vertex, Set{Vertex}}

Returns the `neighbours` map of `graph`.
"""
get_neighbours(graph::Graph) = graph.neighbours

"""
    get_neighbours(G::Graph, u) -> Set{Vertex}

Returns the set of neighbours of `u` in `G`.
"""
get_neighbours(G::Graph, u) = get_neighbours(G)[u]

"""
    num_neighbours(G::Graph, u) -> Integer

Returns the number of neighbours of `u` in `G`.
"""
num_neighbours(G::Graph, u) = length(get_neighbours(G, u))

"""
    num_edges(G::Graph) -> Integer

Returns the number of edges in `G`. The  
edges `(i, j)` and `(j, i)` are counted as one edge.
"""
num_edges(G::Graph) = length(get_edges(G))

"""
    num_vertices(G::Graph) -> Integer

Returns the number of vertices in `G`.
"""
num_vertices(G::Graph) = length(get_vertices(G))

"""
    has_vertex(G::Graph, u) -> Bool

Returns `true` if `u` is a vertex in `G`, and `false` otherwise.
"""
has_vertex(G::Graph, u) = u ∈ get_vertices(G)

"""
    add_vertex!(G::Graph, u...)

Adds the vertices `u...` to `G`.
"""
function add_vertex!(G::Graph{I}, v) where {I}
    has_vertex(G, v) && return G
    V = get_vertices(G)
    push!(V, v)
    N = get_neighbours(G)
    get!(N, v) do
        Set{I}() # in case N is empty, let's add it here 
    end
    return G
end
function add_vertex!(G::Graph{I}, u...) where {I}
    foreach(u) do v
        add_vertex!(G, v)
    end
    return G
end

"""
    add_edge!(G::Graph, u, v)

Adds the edge `(u, v)` to `G`.
"""
function add_edge!(G::Graph{I}, u, v) where {I}
    E = get_edges(G)
    (v, u) ∉ E && push!(E, (u, v))
    return G
end

"""
    delete_edge!(G::Graph, u, v)

Deletes the edge `(u, v)` from `G`.
"""
function delete_edge!(G::Graph{I}, u, v) where {I}
    E = get_edges(G)
    delete!(E, (u, v))
    delete!(E, (v, u)) # undirected graph
    return G
end

"""
    add_neighbour!(G::Graph, u, v...)

Adds the neighbours `v...` to `u` in `G`.
"""
function add_neighbour!(G::Graph{I}, u, v) where {I}
    N = get_neighbours(G)
    u_N = get!(Set{I}, N, u)
    v_N = get!(Set{I}, N, v)
    push!(u_N, v)
    push!(v_N, u) # undirected graph
    add_edge!(G, u, v)
    !has_vertex(G, u) && add_vertex!(G, u)
    !has_vertex(G, v) && add_vertex!(G, v)
    return G
end
function add_neighbour!(G::Graph{I}, u, v...) where {I}
    foreach(v) do w
        add_neighbour!(G, u, w)
    end
    return G
end

"""
    _delete!(G::Graph, u, v)

Deletes the neighbour `v` from `u` in `G`. If the neighbours of `u` are empty once `v` is deleted, then `u` is also deleted from the vertex set.

!!! warning "Undirected graph"

    Even though the graph is undirected, this function only deletes `v` from the neighbours of `u`. If you want to delete `u` from the neighbours of `v`, you should call `delete_neighbour!(G, v, u)`.
"""
function _delete!(G::Graph{I}, u, v) where {I}
    N = get_neighbours(G)
    u_N = get!(Set{I}, N, u)
    delete!(u_N, v)
    iszero(num_neighbours(G, u)) && delete_vertex!(G, u) # We need to do this because, for weighted triangulations, this vertex might be forever missing (because it's  submerged) and cause issues with empty sets
    return G
end

"""
    delete_neighbour!(G::Graph, u, v...)

Deletes the neighbours `v...` from `u` in `G`.
"""
function delete_neighbour!(G::Graph{I}, u, v) where {I}
    _delete!(G, u, v)
    _delete!(G, v, u)
    delete_edge!(G, u, v)
    return G
end
function delete_neighbour!(G::Graph{I}, u, v...) where {I}
    foreach(v) do w
        delete_neighbour!(G, u, w)
    end
    return G
end

"""
    add_triangle!(G::Graph, u, v, w)
    add_triangle!(G::Graph, T)

Adds the neighbourhood relationships defined by the triangle `T = (u, v, w)` to the graph `G`.
"""
function add_triangle!(G::Graph, u::Integer, v::Integer, w::Integer)
    add_vertex!(G, u, v, w)
    add_neighbour!(G, u, v, w)
    add_neighbour!(G, v, w) # undirected graph, so don't need to do all the cases
    return G
end
add_triangle!(G::Graph, T) = add_triangle!(G, geti(T), getj(T), getk(T))

"""
    delete_triangle!(G::Graph, u, v, w)
    delete_triangle!(G::Graph, T)

Deletes the neighbourhood relationships defined by the triangle `T = (u, v, w)` from the graph `G`.
"""
function delete_triangle!(G::Graph, u::Integer, v::Integer, w::Integer)
    delete_neighbour!(G, u, v, w)
    delete_neighbour!(G, v, w) # undirected graph, so don't need to do all the cases
    return G
end
delete_triangle!(G::Graph, T) = delete_triangle!(G, geti(T), getj(T), getk(T))

"""
    delete_vertex!(G::Graph, u...)

Deletes the vertices `u...` from `G`.
"""
function delete_vertex!(G::Graph{I}, u) where {I}
    N = get_neighbours(G)
    u_N = get!(Set{I}, N, u)
    foreach(u_N) do v
        v_N = get!(Set{I}, N, v)
        delete!(v_N, u) # don't use delete_neighbour! because it will delete u_N[v] as well
    end
    delete!(N, u)
    u_G = get_vertices(G)
    delete!(u_G, u)
    return G
end
function delete_vertex!(G::Graph, u...)
    foreach(u) do v
        delete_vertex!(G, v)
    end
    return G
end

"""
    delete_ghost_vertices!(G::Graph)

Deletes all ghost vertices from `G`.
"""
function delete_ghost_vertices_from_graph!(G::Graph{I}) where {I}
    V = get_vertices(G)
    min_g = minimum(V)
    for i in min_g:I(𝒢)
        delete_vertex!(G, i)
    end
    return G
end

"""
    has_ghost_vertices(G::Graph) -> Bool

Returns `true` if `G` has ghost vertices, and `false` otherwise.
"""
function has_ghost_vertices(G::Graph)
    V = get_vertices(G)
    return any(is_ghost_vertex, V)
end

"""
    clear_empty_vertices!(G::Graph)

Deletes all empty vertices from `G`.
"""
function clear_empty_vertices!(G::Graph)
    V = get_vertices(G)
    foreach(V) do i
        n = num_neighbours(G, i)
        iszero(n) && delete_vertex!(G, i)
    end
    return G
end

function Base.empty!(G::Graph)
    V = get_vertices(G)
    E = get_edges(G)
    N = get_neighbours(G)
    empty!(V)
    empty!(E)
    empty!(N)
    return G
end
