############################################
##
## CONSTRUCTORS FOR THE TRIANGULATION DATA STRUCTURES 
##
############################################
@testset "Adjacent" begin
    adj = DT.Adjacent()
    @test adj.adjacent == Dict{Edge{Int64},Int64}()
    @test typeof(adj) == DT.Adjacent{Int64,Edge{Int64}}
    @test isconcretetype(typeof(adj))
    adj = DT.Adjacent{Int16,Edge{Int16}}()
    @test adj.adjacent == Dict{Edge{Int16},Int16}()
    @test typeof(adj) == DT.Adjacent{Int16,Edge{Int16}}
    @test isconcretetype(typeof(adj))
    @test adjacent(adj) == Dict{Edge{Int16},Int16}()
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex()
    @test adj2v.adjacent2vertex == Dict{Edge{Int64},Int64}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int64,Edge{Int64}}
    @test isconcretetype(typeof(adj2v))
    adj2v = DT.Adjacent2Vertex{Int16,Edge{Int16}}()
    @test adj2v.adjacent2vertex == Dict{Edge{Int16},Int16}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int16,Edge{Int16}}
    @test isconcretetype(typeof(adj2v))
    @test adjacent2vertex(adj2v) == Dict{Int64,Set{Edge{Int16}}}()
end

@testset "DelaunayGraph" begin
    dg = DT.DelaunayGraph()
    @test dg.graph == UndirectedGraph{Int64}()
    @test typeof(dg) == DT.DelaunayGraph{Int64}
    @test isconcretetype(typeof(dg))
    dg = DT.DelaunayGraph{Int16}()
    @test dg.graph == UndirectedGraph{Int16}()
    @test typeof(dg) == DT.DelaunayGraph{Int16}
    @test isconcretetype(typeof(dg))
    @test graph(dg) == UndirectedGraph{Int16}()
    tvg = UndirectedGraph{Int64}()
    add!(tvg, 7)
    add!(tvg, 9)
    add!(tvg, 9, 7)
    dg = DT.DelaunayGraph(tvg)
    @test DT.graph(dg) == tvg
end

@testset "HistoryDAG" begin
    hg = DT.HistoryDAG()
    @test hg.graph == DirectedGraph{Triangle{Int16}}()
    @test typeof(hg) == DT.HistoryDAG{Int64,Triangle{Int64}}
    @test isconcretetype(typeof(hg))
    hg = DT.HistoryDAG{Int16,Triangle{Int16}}()
    @test hg.graph == DirectedGraph{Triangle{Int16}}()
    @test typeof(hg) == DT.HistoryDAG{Int16,Triangle{Int16}}
    @test graph(hg) == DirectedGraph{Triangle{Int16}}()
    dg = DirectedGraph{Triangle{Int64}}()
    add!(dg, Triangle(1, 2, 3))
    add!(dg, Triangle(4, 5, 6))
    add!(dg, Triangle(7, 8, 9))
    add!(dg, Triangle(1, 2, 3), Triangle(4, 5, 6))
    add!(dg, Triangle(1, 2, 3), Triangle(7, 8, 9))
    hg = DT.HistoryDAG(dg)
    @test hg.graph == dg
end