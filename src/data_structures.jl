###################################################
#/
#/
#/ Adjacent 
#/
#/
###################################################
"""
    Adjacent{I, E}

The adjacent map of the triangulation, mapping an edge `(u, v)` to a vertex `w` 
such that `(u, v, w)` is a positively oriented triangle. The map is stored in the 
field `adjacent` as a `Dict{E, I}`, where `E` is the edge type and `I` is the 
integer type.
"""
struct Adjacent{I,E} # Adjacent{IntegerType, EdgeType}
    adjacent::DefaultDict{E,I,I}
    function Adjacent{I,E}() where {I,E}
        A = DefaultDict{E,I,I}(I(DefaultAdjacentValue))
        adj = new{I,E}(A)
        return adj
    end
    Adjacent(adj::DefaultDict{E,I,I}) where {I,E} = new{I,E}(adj)
end
adjacent(adj::Adjacent) = adj.adjacent
edges(adj::Adjacent) = keys(adjacent(adj))

"""
    get_edge(adj::Adjacent, uv)
    get_edge(adj::Adjacent{I, E}, u, v) where {I, E}

Returns the vertex `w` such that `(u, v, w)` is a positively oriented triangle, indexing 
the adjacent map `adj` at the key `(u, v)`.
"""
@inline get_edge(adj::Adjacent{I,E}, uv::E) where {I,E} = adjacent(adj)[uv]::I
@doc (@doc get_edge(::Adjacent, ::Any))
@inline function get_edge(adj::Adjacent{I,E}, u, v) where {I,E}
    uv = construct_edge(E, u, v)
    w = get_edge(adj, uv)
    return w
end

"""
    add_edge!(adj::Adjacent, uv, w) 
    add_edge!(adj::Adjacent{I, E}, u, v, w) where {I, E}

Adds the edge `(u, v)` into the adjacent map `adj`.
"""
function add_edge!(adj::Adjacent, uv, w)
    adjacent(adj)[uv] = w
    return nothing
end
@doc (@doc add_edge!(::Adjacent, ::Any, ::Any))
function add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E}
    uv = construct_edge(E, u, v)
    add_edge!(adj, uv, w)
    return nothing
end

"""
    delete_edge!(adj::Adjacent{I,E}, u, v) where {I,E}

Deletes the edge `(u, v)` from the map `adj`. 
"""
function delete_edge!(adj::Adjacent{I,E}, u, v) where {I,E}
    Euv = construct_edge(E, u, v)
    delete!(adjacent(adj), Euv)
    return nothing
end

"""
    add_triangle!(adj::Adjacent, T...)

Adds the edges of the triangle `T` to the adjacent map `adj`.`
"""
function add_triangle!(adj::Adjacent, i, j, k)
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    return nothing
end
function add_triangle!(adj::Adjacent, T)
    i, j, k = indices(T)
    add_triangle!(adj, i, j, k)
end
@doc (@doc add_triangle!(::Adjacent, ::Any))
function add_triangle!(adj::Adjacent, T::Vararg{Tri,N}) where {Tri,N}
    for i in 1:N
        add_triangle!(adj, T[i])
    end
    return nothing
end

###################################################
#/
#/
#/ Adjacent2Vertex
#/
#/
###################################################
"""
    Adjacent2Vertex{I,Es,E}

The adjacent-to-vertex map of the triangulation, mapping a vertex `w` 
to the collection of all edges `(u, v)` such that `(u, v, w)` is a positively 
oriented triangle. The map is stored in `adjacent2vertex` as a `Dict{I, Es}`, 
where `I` is the integer type and `Es` is the type representing the collection 
of edges. The extra parameter `E` is the edge type.
"""
struct Adjacent2Vertex{I,Es,E}
    adjacent2vertex::Dict{I,Es}
    function Adjacent2Vertex{I,Es,E}() where {I,Es,E}
        D = Dict{I,Es}()
        TA2V = new{I,Es,E}(D)
        return TA2V
    end
    Adjacent2Vertex(adj2v::Dict{I,Es}) where {I,Es} = new{I,Es,eltype(iterate(adj2v)[1][2])}(adj2v)
end
adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    get_edge(adj2v::Adjacent2Vertex, w)

Gets the set of edges adjacent to `w` from the map `adj2v`.
"""
get_edge(adj2v::Adjacent2Vertex, w) = adjacent2vertex(adj2v)[w]

"""
    add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}

Adds the edge `(u, v)` to the set of edges in `adjacent2vertex[w]`.
"""
function add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    existing_segments = get!(Es, adjacent2vertex(adj2v), w)
    push!(existing_segments, uv)
    return nothing
end
@doc (@doc add_edge!(::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E})
function add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    add_edge!(adj2v, w, uv)
    return nothing
end

function delete_edge!(adj2v::Adjacent2Vertex, w, uv)
    uvs = get_edge(adj2v, w)
    delete!(uvs, uv)
    return nothing
end
@doc (@doc delete_edge!(::Adjacent2Vertex, ::Any, ::Any))
function delete_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    delete_edge!(adj2v, w, uv)
end

"""
    delete_point!(adj2v::Adjacent2Vertex, w) 

Deletes the point `w` from the map `adj2v`.
"""
function delete_point!(adj2v::Adjacent2Vertex, w)
    delete!(adjacent2vertex(adj2v), w)
    return nothing
end

###################################################
#/
#/
#/ DelaunayGraph
#/
#/
###################################################
"""
    struct DelaunayGraph{I}

Graph representation of a Delaunay triangulation, mapping all points to their neighbours, 
where a point `v` is a neighbour of `u` if `(u, v)` is an edge in the triangulation. The 
implementation is through an `UndirectedGraph`. The graph can be accessed using [`graph`](@ref).
"""
struct DelaunayGraph{I}
    graph::UndirectedGraph{I}
    function DelaunayGraph{I}() where {I}
        DG = UndirectedGraph{I}()
        return new{I}(DG)
    end
    DelaunayGraph() = DelaunayGraph{Int64}()
    DelaunayGraph(DG::UndirectedGraph{I}) where {I} = new{I}(DG)
end
graph(G::DelaunayGraph) = G.graph
edges(G::DelaunayGraph) = graph(G).E

"""
    add_point!(DG::DelaunayGraph, u)
    add_point!(DG::DelaunayGraph, u::Vararg{I, N})

Addds the point(s) `u` into the vertex set of `DG`.
"""
function add_point!(DG::DelaunayGraph, u)
    add!(graph(DG), u)
    return nothing
end
@doc (@doc add_point!(::DelaunayGraph, ::Any))
function add_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_point!(DG, u[i])
    end
    return nothing
end

"""
    add_neighbour!(DG::DelaunayGraph{I}, u, v...) 

Adds an edge from `u` to the points in `v` into the graph `DG`. Note that, as 
`DG` is an undirected graph, this also adds an edge from `v` to `u`.
"""
function add_neighbour!(DG::DelaunayGraph, u, v)
    add!(graph(DG), u, v)
    return nothing
end
@doc (@doc add_neighbour!(::DelaunayGraph, ::Any, ::Any))
function add_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_neighbour!(DG, u, v[i])
    end
    return nothing
end

"""
    get_neighbour(DG::DelaunayGraph, u)

Gets the neighbours of `u` from the graph `DG`.
"""
function get_neighbour(DG::DelaunayGraph, u)
    return graph(DG).N[u]
end

"""
    delete_neighbour!(DG::DelaunayGraph, u, v)
    delete_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}

Deletes `v` from the neighbourhood of `u`. Note that, as `DG` is an 
undirected graph, this also deletes `u` from the neighbourhood of `v`.
"""
function delete_neighbour!(DG::DelaunayGraph, u, v)
    delete!(graph(DG), u, v)
    return nothing
end
@doc (@doc delete_neighbour!(::DelaunayGraph, ::Any, ::Any))
function delete_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_neighbour!(DG, u, v[i])
    end
    return nothing
end

"""
    delete_point!(DG::DelaunayGraph, u)
    delete_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}

Deletes the point(s) `u` from the vertex set of `DG`.
"""
function delete_point!(DG::DelaunayGraph, u)
    delete!(graph(DG), u)
    return nothing
end
@doc (@doc delete_point!(::DelaunayGraph, ::Any))
function delete_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_point!(DG, u[i])
    end
    return nothing
end

###################################################
#/
#/
#/ HistoryGraph
#/
#/
###################################################
"""
    struct HistoryGraph{I,T<:AbstractTriangle{I}} 

The point location graph for the Delaunay triangulation. This is a directed graph, 
implemented using `DirectedGraph`, that stores the history of the triangulation as points are 
inserted. The graph can be accessed using [`graph`](@ref).
"""
struct HistoryGraph{T}
    graph::DirectedGraph{T}
    function HistoryGraph{T}() where {T}
        G = DirectedGraph{T}()
        forbid_loops!(G)
        TDAG = new{T}(G)
        return TDAG
    end
    HistoryGraph(HG::DirectedGraph{T}) where {T} = new{T}(HG)
end
graph(G::HistoryGraph) = G.graph

"""
    SimpleGraphs.out_neighbors(G::HistoryGraph, T)

Returns the set of triangles that `T` was replaced by in the triangulation.
"""
SimpleGraphs.out_neighbors(G::HistoryGraph, T) = graph(G).N[T] # The version in SimpleGraphs is allocating...

"""
    SimpleGraphs.in_neighbors(G::HistoryGraph, T)

Returns the triangles that were replaced to form the triangle `T`.
"""
SimpleGraphs.in_neighbors(G::HistoryGraph, T) = SimpleGraphs.in_neighbors(graph(G), T) # Need the allocating version for find_root. Only place it's used anyway.

"""
    SimpleGraphs.in_deg(G::HistoryGraph, T)

Returns the number of triangles that were replaced to form the triangle `T`.
"""
SimpleGraphs.in_deg(G::HistoryGraph, T) = SimpleGraphs.in_deg(graph(G), T)

"""
    SimpleGraphs.out_deg(G::HistoryGraph, T)

Returns the number of triangles that `T` was replaced by.
"""
SimpleGraphs.out_deg(G::HistoryGraph, T) = SimpleGraphs.out_deg(graph(G), T)

"""
    add_triangle!(D::HistoryGraph, T)

Adds the triangles in `T` into the graph `D`. Checks are made for 
circular shifts of the vertices of `T` to avoid duplicates.
"""
function add_triangle!(D::HistoryGraph, T)
    has(graph(D), T) && return nothing
    has(graph(D), shift_triangle_1(T)) && return nothing
    has(graph(D), shift_triangle_2(T)) && return nothing
    add!(graph(D), T)
    return nothing
end
@doc (@doc add_triangle!(::HistoryGraph, ::Any))
function add_triangle!(D::HistoryGraph, T::Vararg{Tri,N}) where {Tri,N}
    for i in 1:N
        add_triangle!(D, T[i])
    end
end

"""
    add_edge!(G::HistoryGraph, T, V)
    add_edge!(G::HistoryGraph, T, V::Vararg{Tri, N}) where {Tri, N}

Adds a directed edge from `T` to the triangle(s) in `V`, making checks for circular shifts 
in the vertices of `T` and `V`.
"""
function add_edge!(G::HistoryGraph, T, V)
    for i in 0:2
        Tshift = shift_triangle(T, i)
        if has(graph(G), Tshift)
            for j in 0:2
                Vshift = shift_triangle(V, j)
                if has(graph(G), Vshift)
                    SimpleGraphs.add!(graph(G), Tshift, Vshift)
                    return nothing
                end
            end
        end
    end
    return nothing
end
@doc (@doc add_edge!(::HistoryGraph, ::Any, ::Any))
function add_edge!(G::HistoryGraph, T, V::Vararg{Tri,N}) where {Tri,N}
    for V in V
        add_edge!(G, T, V)
    end
    return nothing
end