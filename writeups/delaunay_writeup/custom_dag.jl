const triangle_type = NTuple{3, Int64}
struct TriangulationDAG <: AbstractSimpleGraph
    graph::DirectedGraph{triangle_type}
    function TriangulationDAG()
        G = DirectedGraph{triangle_type}()
        forbid_loops!(G)
        TDAG = new(G)
        return TDAG
    end
end
const TriDAG = TriangulationDAG
graph(G::TriangulationDAG) = G.graph
for op in (:out_deg, :in_deg, :deg,
    :dual_deg, :out_neighbors, :in_neighbors,
    :NE, :has, :add!, :delete!,
    :elist)
    @eval begin
        SimpleGraphs.$op(G::TriDAG) = SimpleGraphs.$op(graph(G))
        SimpleGraphs.$op(G::TriDAG, u) = SimpleGraphs.$op(graph(G), u)
        SimpleGraphs.$op(G::TriDAG, u, v) = SimpleGraphs.$op(graph(G), u, v)
    end
end