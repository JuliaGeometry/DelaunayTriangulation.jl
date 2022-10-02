############################################
##
## CONSTRUCTORS FOR THE TRIANGULATION DATA STRUCTURES 
##
############################################
"""
    struct Adjacent{I,E<:AbstractEdge{I}}

Data structure for the `adjacent` map of a triangulation. The map can be accessed 
via [`adjacent`](@ref). The map is defined using a `Dict{E, I}`, such that the key 
`E(u, v)` is mapped to an `w::I` such that `Triangle(u, v, w)` is positively oriented.
"""
struct Adjacent{I,E<:AbstractEdge{I}}
    adjacent::Dict{E,I}
    function Adjacent{I,E}() where {I,E}
        A = Dict{E,I}()
        adj = new{I,E}(A)
        return adj
    end
    Adjacent() = Adjacent{Int64,Edge{Int64}}()
end
adjacent(adj::Adjacent) = adj.adjacent
edges(adj::Adjacent) = keys(adj.adjacent)

"""
    struct Adjacent2Vertex{I,E<:AbstractEdge{I}}

Data structure for the `adjacent2vertex` map of a triangulation. The map can be accessed 
via [`adjacent2vertex`](@ref). The map is defined using a `Dict{I,Set{E}}`, such that the key 
`w` is mapped to a `Set` containing all pairs `(u, v)` such that `Triangle(u, v, w)` is a 
positively oriented triangle.
"""
struct Adjacent2Vertex{I,E<:AbstractEdge{I}}
    adjacent2vertex::Dict{I,Set{E}}
    function Adjacent2Vertex{I,E}() where {I,E}
        D = Dict{I,Set{E}}()
        TA2V = new{I,E}(D)
        return TA2V
    end
    Adjacent2Vertex() = Adjacent2Vertex{Int64,Edge{Int64}}()
end
adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    struct DelaunayGraph{I}

Graph representation of a Delaunay triangulation, mapping all points to their neighbours, 
where a point `v` is a neighbour of `u` is `(u, v)` is an edge in the triangulation. The 
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

"""
    struct HistoryGraph{I,T<:AbstractTriangle{I}} 

The point location graph for the Delaunay triangulation. This is a directed graph, 
implemented using `DirectedGraph`, that stores the history of the triangulation as points are 
inserted. The graph can be accessed using [`graph`](@ref).
"""
struct HistoryGraph{I,T<:AbstractTriangle{I}}
    graph::DirectedGraph{T}
    function HistoryGraph{I,T}() where {I,T}
        G = DirectedGraph{T}()
        TDAG = new{I,T}(G)
        return TDAG
    end
    HistoryGraph() = HistoryGraph{Int64,Triangle{Int64}}()
    function HistoryGraph(HG::DirectedGraph{T}) where {I,T<:AbstractTriangle{I}}
        return new{I,T}(HG)
    end
end
graph(G::HistoryGraph) = G.graph

"""
    Triangulation{A,A2V,DG,H,T,P,R}

Struct for a Delaunay triangulation. See also [`triangulate`](@ref).

# Fields 
- `adjacent`: The adjacent map. See [`Adjacent`](@ref).
- `adjacent2vertex`: The adjacent-to-vertex map. See [`Adjacent2Vertex`](@ref).
- `graph`: The graph representation of the triangulation. See also [`DelaunayGraph`](@ref).
- `history`: The history structure for the triangulation. See also [`HistoryGraph`](@ref).
- `triangles`: The triangles that define the triangulation. See also [`Triangles`](@ref) and [`Triangle`](@ref).
- `points`: The point set of the triangulation. See also [`Points`](@ref) and [`Point`](@ref).
- `root`: The root of `history`.
"""
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
adjacent2vertex(DT::Triangulation) = DT.adjacent2vertex
graph(DT::Triangulation) = DT.graph
history(DT::Triangulation) = DT.history
triangles(DT::Triangulation) = DT.triangles
points(DT::Triangulation) = DT.points
root(DT::Triangulation) = DT.root