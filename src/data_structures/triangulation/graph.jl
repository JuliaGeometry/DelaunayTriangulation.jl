abstract type AbstractGraph{I,E} end

struct Graph{I,E} <: AbstractGraph{I,E}
    adjacent2vertex::Adjacent2Vertex{I,E}
end

get_adjacent2vertex(graph::Graph) = graph.adjacent2vertex

function Base.show(io::IO, ::MIME"text/plain", graph::AbstractGraph)
    println(io, "Graph")
    println(io, "    Number of edges: ", num_edges(graph))
    print(io, "    Number of vertices: ", num_vertices(graph))
end
Base.copy(graph::Graph) = Graph(copy(get_adjacent2vertex(graph)))

function get_vertices(graph::Graph)
    adj2v = get_adjacent2vertex(graph)
    candidates = get_candidates(adj2v)
    return keys(candidates)
end

function get_edges(graph::Graph)
    adj2v = get_adjacent2vertex(graph)
    adj = get_adjacent(adj2v)
    return EdgeIterator(adj)
end

struct GraphNeighbourhoodIterator{A,E,I}
    adj::A # parametric so that this can be used for both Adjacent and Triangulation
    u::I
    v::I
end
GraphNeighbourhoodIterator(adj::Adjacent{I,E}, u, v) where {I,E} = GraphNeighbourhoodIterator{Adjacent{I,E},E,I}(adj, u, v)
function Base.iterate(itr::GraphNeighbourhoodIterator)
    adj, u, v = itr.adj, itr.u, itr.v
    !edge_exists(v) && return nothing
    w = get_adjacent(adj, u, v)
    return v, w
end
function Base.iterate(itr::GraphNeighbourhoodIterator, w)
    adj, u, v = itr.adj, itr.u, itr.v
    (!edge_exists(v) || w == v) && return nothing
    x = get_adjacent(adj, u, w)
    return w, x
end
Base.IteratorSize(::Type{<:GraphNeighbourhoodIterator}) = Base.SizeUnknown()
Base.eltype(::Type{GraphNeighbourhoodIterator{A,E,I}}) where {A,E,I} = I
function Base.length(itr::GraphNeighbourhoodIterator)
    n = 0
    for _ in itr
        n += 1
    end
    return n
end

function get_neighbours(graph::Graph, u)
    adj2v = get_adjacent2vertex(graph)
    adj = get_adjacent(adj2v)
    v = find_adjacent(adj2v, u)
    return GraphNeighbourhoodIterator(adj, u, v)
end

num_neighbours(graph::AbstractGraph, u) = length(get_neighbours(graph, u))

num_edges(graph::AbstractGraph) = length(get_edges(graph))

num_vertices(graph::AbstractGraph) = length(get_vertices(graph))

has_vertex(graph::AbstractGraph, u) = u ∈ get_vertices(graph)

has_ghost_vertices(graph::AbstractGraph) = any(is_ghost_vertex, get_vertices(graph))
 
struct ConcreteGraph{I,E} <: AbstractGraph{I,E}
    vertices::Set{I}
    edges::Set{E}
    neighbours::Dict{I,Set{I}}
end
function concretize_graph(graph)
    vertices = Set(get_vertices(graph))
    edges = Set(get_edges(graph))
    I = integer_type(graph)
    neighbours = Dict{I,Set{I}}()
    for u in vertices
        neighbours[u] = Set(get_neighbours(graph, u))
    end
    return ConcreteGraph(vertices, edges, neighbours)
end