############################################
##
## INDEXING AND UPDATING OF THE TRIANGULATION DATA STRUCTURES
##
############################################
@testset "Adjacent" begin
    adj = DT.Adjacent()
    i, j, k, r = 1, 2, 3, 4
    DT.add_edge!(adj, i, j, r)
    DT.add_edge!(adj, j, r, i)
    DT.add_edge!(adj, r, i, j)
    DT.add_edge!(adj, i, k, j)
    DT.add_edge!(adj, k, j, i)
    DT.add_edge!(adj, j, i, k)
    @test DT.get_edge(adj, i, j) == r
    @test DT.get_edge(adj, j, r) == i
    @test DT.get_edge(adj, r, i) == j
    @test DT.get_edge(adj, i, k) == j
    @test DT.get_edge(adj, k, j) == i
    @test DT.get_edge(adj, j, i) == k
    @test DT.get_edge(adj, Edge{Int64}(j, r)) == i
    DT.add_edge!(adj, Edge{Int64}(5, 7), 9)
    DT.add_edge!(adj, Edge{Int64}(7, 5), 3)
    @test DT.get_edge(adj, 5, 7) == DT.get_edge(adj, Edge{Int64}(5, 7)) == 9
    DT.delete_edge!(adj, 5, 7)
    @test_throws KeyError DT.get_edge(adj, 5, 7)
    @test_throws KeyError DT.get_edge(adj, Edge{Int64}(5, 7))
    DT.delete_edge!(adj, i, j)
    @test_throws KeyError DT.get_edge(adj, Edge{Int64}(i, j))
    @test_throws KeyError DT.get_edge(adj, i, j)
    @test_throws KeyError DT.get_edge(adj, Edge{Int64}(j, i))
    @test_throws KeyError DT.get_edge(adj, i, j)
    DT.add_edge!(adj, 16, 13, 5)
    DT.add_edge!(adj, 13, 16, DT.BoundaryIndex)
    DT.delete_edge!(adj, 16, 13)
    @test_throws KeyError DT.get_edge(adj, 16, 13)
    @test DT.get_edge(adj, 13, 16) == DT.BoundaryIndex
    DT.add_edge!(adj, 16, 13, 5)
    DT.delete_edge!(adj, 16, 13; protect_boundary=false)
    @test_throws KeyError DT.get_edge(adj, 16, 13)
    @test_throws KeyError DT.get_edge(adj, 13, 16)
    adj = DT.Adjacent()
    T1 = Triangle(1, 2, 3)
    T2 = Triangle(4, 5, 6)
    T3 = Triangle(7, 8, 9)
    DT.add_triangle!(adj, T1)
    @test DT.get_edge(adj, 1, 2) == 3
    @test DT.get_edge(adj, 2, 3) == 1
    @test DT.get_edge(adj, 3, 1) == 2
    DT.add_triangle!(adj, T2, T3)
    @test DT.get_edge(adj, 4, 5) == 6
    @test DT.get_edge(adj, 5, 6) == 4
    @test DT.get_edge(adj, 6, 4) == 5
    @test DT.get_edge(adj, 7, 8) == 9
    @test DT.get_edge(adj, 8, 9) == 7
    @test DT.get_edge(adj, 9, 7) == 8
    @test length(DT.adjacent(adj)) == 9

    p1 = Point(-1.0, 4.0)
    p2 = Point(4.0, 6.0)
    p3 = Point(2.0, 2.0)
    p4 = Point(6.0, -3.0)
    p5 = Point(7.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)
    adj = DT.Adjacent()
    DT.add_triangle!(adj, Triangle(1, 3, 2))
    DT.add_triangle!(adj, Triangle(1, 4, 3))
    DT.add_triangle!(adj, Triangle(3, 4, 5))
    DT.add_triangle!(adj, Triangle(3, 5, 2))
    @test DT.get_edge(adj, 1, 3) == 2
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, 4) == 3
    @test DT.get_edge(adj, 4, 3) == 1
    @test DT.get_edge(adj, 3, 1) == 4
    @test DT.get_edge(adj, 3, 4) == 5
    @test DT.get_edge(adj, 4, 5) == 3
    @test DT.get_edge(adj, 5, 3) == 4
    @test DT.get_edge(adj, 3, 5) == 2
    @test DT.get_edge(adj, 5, 2) == 3
    @test DT.get_edge(adj, 2, 3) == 5

    adj = DT.Adjacent()
    DT.add_edge!(adj, 2, 3, 7)
    DT.add_edge!(adj, 5, 7, 17)
    @test collect(DT.edges(adj)) == [Edge{Int64}(2, 3), Edge{Int64}(5, 7)] || collect(DT.edges(adj)) == [Edge{Int64}(5, 7), Edge{Int64}(2, 3)]
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex()
    DT.add_edge!(adj2v, 3, Edge{Int64}(1, 2))
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([Edge{Int64}(1, 2)])
    DT.add_edge!(adj2v, 10, 5, 7)
    @test DT.get_edge(adj2v, 10) == Set{Edge{Int64}}([Edge{Int64}(5, 7)])
    DT.add_edge!(adj2v, 3, Edge{Int64}(1, 4))
    DT.add_edge!(adj2v, 3, Edge{Int64}(10, 11))
    DT.add_edge!(adj2v, 3, Edge{Int64}(13, 5))
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([Edge{Int64}(1, 2), Edge{Int64}(1, 4), Edge{Int64}(10, 11), Edge{Int64}(13, 5)])
    DT.delete_edge!(adj2v, 3, Edge{Int64}(1, 4))
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([Edge{Int64}(1, 2), Edge{Int64}(10, 11), Edge{Int64}(13, 5)])
    DT.delete_edge!(adj2v, 3, 1, 2)
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([Edge{Int64}(10, 11), Edge{Int64}(13, 5)])
    DT.delete_point!(adj2v, 10)
    @test_throws KeyError DT.get_edge(adj2v, 10)

    adj2v = DT.Adjacent2Vertex()
    i, j, k, r = 5, 9, 13, 6
    DT.add_edge!(adj2v, i, k, j)
    DT.add_edge!(adj2v, i, j, r)
    DT.add_edge!(adj2v, i, 20, 30) # need to make sure we aren't just removing all the edges
    DT.add_edge!(adj2v, j, i, k)
    DT.add_edge!(adj2v, j, r, i)
    DT.add_edge!(adj2v, j, 50, 70)
    DT.add_edge!(adj2v, k, j, i)
    DT.add_edge!(adj2v, k, 100, 110)
    DT.add_edge!(adj2v, r, i, j)
    DT.add_edge!(adj2v, r, 130, 150)
    DT.add_edge!(adj2v, 17, 19, 32)
    DT.update_after_flip!(adj2v, i, j, k, r)
    @test DT.get_edge(adj2v, i) == Set{Edge{Int64}}([Edge{Int64}(k, r), Edge{Int64}(20, 30)])
    @test DT.get_edge(adj2v, j) == Set{Edge{Int64}}([Edge{Int64}(r, k), Edge{Int64}(50, 70)])
    @test DT.get_edge(adj2v, k) == Set{Edge{Int64}}([Edge{Int64}(r, i), Edge{Int64}(j, r), Edge{Int64}(100, 110)])
    @test DT.get_edge(adj2v, r) == Set{Edge{Int64}}([Edge{Int64}(i, k), Edge{Int64}(k, j), Edge{Int64}(130, 150)])
    @test DT.get_edge(adj2v, 17) == Set{Edge{Int64}}([Edge{Int64}(19, 32)])
    @test length(adjacent2vertex(adj2v)) == 5

    adj2v = DT.Adjacent2Vertex()
    i, j, k, r = 17, 13, 29, 9
    DT.add_edge!(adj2v, i, j, k)
    DT.add_edge!(adj2v, i, 170, 302)
    DT.add_edge!(adj2v, j, k, i)
    DT.add_edge!(adj2v, j, 500, 501)
    DT.add_edge!(adj2v, k, i, j)
    DT.add_edge!(adj2v, k, 177, 111)
    DT.add_edge!(adj2v, 37, 400, 1)
    DT.update_after_insertion!(adj2v, i, j, k, r)
    @test DT.get_edge(adj2v, i) == Set{Edge{Int64}}([Edge{Int64}(j, r), Edge{Int64}(r, k), Edge{Int64}(170, 302)])
    @test DT.get_edge(adj2v, j) == Set{Edge{Int64}}([Edge{Int64}(k, r), Edge{Int64}(r, i), Edge{Int64}(500, 501)])
    @test DT.get_edge(adj2v, k) == Set{Edge{Int64}}([Edge{Int64}(i, r), Edge{Int64}(r, j), Edge{Int64}(177, 111)])
    @test DT.get_edge(adj2v, r) == Set{Edge{Int64}}([Edge{Int64}(j, k), Edge{Int64}(k, i), Edge{Int64}(i, j)])
    @test DT.get_edge(adj2v, 37) == Set{Edge{Int64}}([Edge{Int64}(400, 1)])
    @test length(adjacent2vertex(adj2v)) == 5
end

@testset "DelaunayGraph" begin
    DG = DT.DelaunayGraph()
    DT.add_point!(DG, 5)
    @test 5 ∈ DT.graph(DG).V
    DT.add_point!(DG, 17, 53, 32)
    @test all(v ∈ DT.graph(DG).V for v in [5, 17, 53, 32])
    DT.add_neighbour!(DG, 5, 17)
    @test DT.graph(DG).N[5] == Set(17)
    DT.add_neighbour!(DG, 5, 53, 32)
    @test DT.graph(DG).N[5] == Set([17, 53, 32])
    @test DT.get_neighbour(DG, 5) == Set([17, 53, 32])
    @test DT.get_neighbour(DG, 53) == Set(5)
    @test DT.get_neighbour(DG, 32) == Set(5)
    DT.delete_neighbour!(DG, 32, 5)
    @test DT.get_neighbour(DG, 5) == Set([17, 53])
    @test DT.get_neighbour(DG, 32) == Set()
    DT.delete_point!(DG, 5)
    @test_throws KeyError DT.get_neighbour(DG, 5)
    DT.add_point!(DG, 17, 101, 291)
    @test 17 ∈ DT.graph(DG).V
    @test 101 ∈ DT.graph(DG).V
    @test 291 ∈ DT.graph(DG).V
    DT.delete_point!(DG, 17, 101, 291)
    @test 17 ∉ DT.graph(DG).V
    @test 101 ∉ DT.graph(DG).V
    @test 291 ∉ DT.graph(DG).V

    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    DG = DT.DelaunayGraph(tvn)
    @test DT.graph(DG) == tvn
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3])
    @test DT.get_neighbour(DG, 2) == Set([11, 9, 8, 15, 16])
    @test DT.get_neighbour(DG, 3) == Set([1])
    DT.add_point!(DG, 40)
    DT.add_neighbour!(DG, 40, 50)
    @test DT.get_neighbour(DG, 40) == Set([50])
    DT.add_neighbour!(DG, 60, 70)
    @test DT.get_neighbour(DG, 60) == Set([70])
    @test DT.get_neighbour(DG, 70) == Set([60])
    DT.add_point!(DG, 70, -1, -2, -3, -4, -5, -6, -7)
    @test all(has(DT.graph(DG), u) for u ∈ [70, -1, -2, -3, -4, -5, -6, -7])
    DT.add_neighbour!(DG, 70, 4, 5, 9, 3)
    @test DT.get_neighbour(DG, 70) == Set([4, 60, 5, 9, 3])
    DT.add_neighbour!(DG, 1, 10)
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3, 10])
    DT.add_neighbour!(DG, 1, 10, 11, 12, 13)
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3, 10, 11, 12, 13])
    DG_empty = DT.DelaunayGraph()
    @test isempty(DG_empty.graph.V)
    DT.add_neighbour!.(Ref(DG_empty), 2, [3, 4, 5, 1])
    @test DT.get_neighbour(DG_empty, 2) == Set([3, 4, 5, 1])

    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    DG = DT.DelaunayGraph(tvn)
    DT.delete_neighbour!(DG, 1, 6)
    @test DT.get_neighbour(DG, 1) == Set([3, 4, 5, 9])
    @test DT.get_neighbour(DG, 2) == Set([11, 9, 8, 15, 16])
    @test DT.get_neighbour(DG, 3) == Set([1])
end

@testset "HistoryGraph" begin
    dg = DirectedGraph{Triangle{Int64}}()
    add!(dg, Triangle(1, 2, 3))
    add!(dg, Triangle(4, 5, 6))
    add!(dg, Triangle(7, 8, 9))
    add!(dg, Triangle(1, 2, 3), Triangle(4, 5, 6))
    add!(dg, Triangle(1, 2, 3), Triangle(7, 8, 9))
    hg = DT.HistoryGraph(dg)
    @test out_neighbors(hg, Triangle(1, 2, 3)) == Set([Triangle(4, 5, 6), Triangle(7, 8, 9)])
    @test in_neighbors(hg, Triangle(7, 8, 9)) == [Triangle(1, 2, 3)]
    @test in_deg(hg, Triangle(1, 2, 3)) == 0
    @test in_deg(hg, Triangle(7, 8, 9)) == 1
    @test out_deg(hg, Triangle(1, 2, 3)) == 2
    @test out_deg(hg, Triangle(7, 8, 9)) == 0
    DT.add_triangle!(hg, Triangle(11, 13, 15))
    @test Triangle(11, 13, 15) ∈ hg.graph.V
    DT.add_triangle!(hg, Triangle(17, 18, 19), Triangle(22, 25, 29))
    @test Triangle(17, 18, 19) ∈ hg.graph.V
    @test Triangle(22, 25, 29) ∈ hg.graph.V
    DT.add_triangle!(hg, Triangle(18, 19, 17))
    DT.add_edge!(hg, Triangle(17, 18, 19), Triangle(22, 25, 29))
    @test Triangle(22, 25, 29) ∈ out_neighbors(hg, Triangle(17, 18, 19))
    DT.add_edge!(hg, Triangle(18, 19, 17), Triangle(1, 2, 3))
    @test Triangle(1, 2, 3) ∈ out_neighbors(hg, Triangle(17, 18, 19))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    og = out_deg(hg, Triangle(7, 8, 9))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    @test og == out_deg(hg, Triangle(7, 8, 9))

    H = DT.HistoryGraph()
    DT.add_triangle!(H, Triangle(1, 2, 3))
    DT.add_triangle!(H, Triangle(4, 5, 6))
    DT.add_edge!(H, Triangle(1, 2, 3), Triangle((4, 5, 6)))
    Htrue = deepcopy(H)
    DT.add_edge!(H, Triangle(2, 3, 1), Triangle((4, 5, 6)))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, Triangle(2, 3, 1), Triangle((5, 6, 4)))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, Triangle(1, 2, 3), Triangle((6, 4, 5)))
    @test Htrue.graph == H.graph
    DT.add_triangle!(H, Triangle(3, 1, 2))
    @test Htrue.graph == H.graph

    H = DT.HistoryGraph()
    HH = DT.HistoryGraph()
    T₁ = Triangle((1, 2, 3))
    T₂ = Triangle((4, 5, 6))
    T₃ = Triangle((7, 8, 9))
    DT.add_triangle!(H, T₁, T₂, T₃)
    DT.add_triangle!(HH, T₁)
    DT.add_triangle!(HH, T₂)
    DT.add_triangle!(HH, T₃)
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, Triangle((2, 3, 1))) # same triangle as T₁
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, Triangle((2, 3, 1)), Triangle((9, 7, 8))) # same as T₁ and T₃
    @test DT.graph(H) == DT.graph(HH)
    DT.add_edge!(H, T₁, T₂, T₃)
    DT.add_edge!(HH, T₁, T₂)
    DT.add_edge!(HH, T₁, T₃)
    @test DT.graph(H) == DT.graph(HH)
end