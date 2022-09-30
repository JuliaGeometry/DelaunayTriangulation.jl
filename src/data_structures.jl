## Adjacent
"""
    Adjacent{A}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `adjacent`, so that 
`(u, v, adjacent[(u, v)])` is a positively oriented triangle. This struct is also callable, e.g. if 
`𝒜::Adjacent`, then `(u, v, 𝒜(u, v))` is a positively oriented triangle.
"""
struct Adjacent{A}
    adjacent::A
    function Adjacent()
        A = Dict{EdgeType,Int64}()
        adj = new{typeof(A)}(A)
        return adj
    end
end
Base.setindex!(adj::Adjacent, w, uv) = Base.setindex!(adj.adjacent, w, uv) # (u, v, w) is a positively oriented triangle
Base.getindex(adj::Adjacent, uv) = Base.getindex(adj.adjacent, uv) # Find the w so that (u, v, w) is a positively oriented triangle
edges(adj::Adjacent) = keys(adj.adjacent) # List of edges - note that there are duplicates, i.e. (u, v) and (v, u) are both stored.
(adj::Adjacent)(u, v) = adj[(u, v)] # Find the w so that (u, v, w) is a positively oriented triangle

## DelaunayGraph
struct DelaunayGraph
    graph::UndirectedGraph{Int64}
    function DelaunayGraph()
        DG = UndirectedGraph{Int64}()
        return new(DG)
    end
    DelaunayGraph(DG) = new(DG)
end
graph(G::DelaunayGraph) = G.graph
function add_neighbour!(DG::DelaunayGraph, u::Integer, v...)
    for w in v
        add!(graph(DG), u, w)
    end
    return nothing
end
neighbours(DG::DelaunayGraph, u) = graph(DG).N[u]
(DG::DelaunayGraph)(u) = neighbours(DG, u)
Base.getindex(DG::DelaunayGraph, u) = neighbours(DG, u)
points(𝒟𝒢::DelaunayGraph) = graph(𝒟𝒢).V

## HistoryDAG 
"""
    HistoryDAG <: AbstractSimpleGraph

Data structure for the Delaunay triangulation. This is a directed acyclic graph 
that is useful for point location.
"""
struct HistoryDAG <: AbstractSimpleGraph
    graph::DirectedGraph{TriangleType}
    function HistoryDAG()
        G = DirectedGraph{TriangleType}()
        forbid_loops!(G)
        TDAG = new(G)
        return TDAG
    end
end

graph(G::HistoryDAG) = G.graph
SimpleGraphs.out_neighbors(G::HistoryDAG, u) = graph(G).N[u] # The version in SimpleGraphs is allocating...
SimpleGraphs.in_neighbors(G::HistoryDAG, u) = SimpleGraphs.in_neighbors(graph(G), u) # Need the allocating version for find_root. Only place it's used anyway.
SimpleGraphs.in_deg(G::HistoryDAG, u) = SimpleGraphs.in_deg(graph(G), u)
SimpleGraphs.out_deg(G::HistoryDAG, u) = SimpleGraphs.out_deg(graph(G), u)

"""
    find_root(G::HistoryDAG; method=:brute)

Finds the root of the graph `G`, assuming only one root exists and the only node with in-degree zero is this root. 
There are two methods:

    - `method=:brute`: In this case, all vertices are searched until one is found with in-degree zero.
    - `method=:rng`: In this case, a random integer is selected, defniing the node to start at. 
    We then go upwards until the root is found; note that this only works as the graph is acyclic.
"""
function find_root(G::HistoryDAG; method=:brute)
    if method == :brute
        return _find_root_brute(G)
    elseif method == :rng
        return _find_root_rng(G)
    end
end
function _find_root_brute(G::HistoryDAG)
    for (k, v) in graph(G).NN
        length(v) == 0 && return k
    end
end
function _find_root_rng(G::HistoryDAG)
    verts = vlist(graph(G))
    num_verts = length(verts)
    starting_node = verts[rand(1:num_verts)]
    for _ in 1:num_verts
        num_choices = in_deg(G, starting_node)
        num_choices == 0 && return starting_node
        starting_node = in_neighbors(G, starting_node)[rand(1:num_choices)]
    end
end

## Adjacent2Vertex 
""" 
    Adjacent2Vertex{A⁻¹}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `adjacent2vertex`, so that, 
for each `(u, v) ∈ adjacent2vertex[w]`, `(u, v, w)` is a positively oriented triangle. Note that 
this is the inverse of `TriangulationAdjacent`.
"""
struct Adjacent2Vertex{A⁻¹}
    adjacent2vertex::A⁻¹
    function Adjacent2Vertex()
        D = Dict{Int64,Adjacent2VertexVector}()
        TA2V = new{typeof(D)}(D)
        return TA2V
    end
end
Base.setindex!(adj::Adjacent2Vertex, uv, w) = Base.setindex!(adj.adjacent2vertex, uv, w) # (u, v, w) is a positively oriented triangle
Base.getindex(adj::Adjacent2Vertex, w) = Base.getindex(adj.adjacent2vertex, w) # Find the (u, v) so that (u, v, w) is a positively oriented triangle
Base.get!(f, adj::Adjacent2Vertex, w) = Base.get!(f, adj.adjacent2vertex, w)

## Triangulation 
"""
    Triangles

Data structure representing a collection of triangles.
"""
struct Triangles
    triangles::Set{TriangleType}
    Triangles() = new(TriangleType[])
    Triangles(tris::Set{TriangleType}) = new(tris)
    Triangles(tris::AbstractVector{TriangleType}) = new(Set(tris))
end
Base.delete!(𝒯::Triangles, T) = Base.delete!(𝒯.triangles, T)

struct Triangulation{A,A2V,DG,H,T,P,R}
    adjacent::A
    adjacent2vertex::A2V
    graph::DG
    history::H
    triangles::T
    points::P
    root::R
end
adjacent(DT::Triangulation) = DT.adjacent
adjacent(DT::Triangulation, u, v) = adjacent(DT)(u, v)
adjacent2vertex(DT::Triangulation) = DT.adjacent2vertex
adjacent2vertex(DT::Triangulation, u) = adjacent2vertex(DT)[u]
graph(DT::Triangulation) = DT.graph
neighbours(DT::Triangulation, u) = graph(DT)(u)
history(DT::Triangulation) = DT.history
triangles(DT::Triangulation) = DT.triangles
points(DT::Triangulation) = DT.points
points(DT::Triangulation, i) = points(DT)[i]
root(DT::Triangulation) = DT.root
num_points(DT::Triangulation) = length(DT.points)
num_triangles(DT::Triangulation) = length(DT.triangles)
