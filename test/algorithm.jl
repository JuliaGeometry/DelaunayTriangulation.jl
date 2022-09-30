############################################
##
## SPECIFIC DATA TYPES 
##
############################################
@testset "Edge" begin
    e = Edge(1, 2)
    @test e.initial == 1
    @test e.terminal == 2
    @test initial(e) == 1
    @test terminal(e) == 2
    e = Edge((3, 10))
    @test initial(e) == 3
    @test terminal(e) == 10
    @test Edge{Int} <: DT.AbstractEdge{Int}
    @test typeof(e) == Edge{Int}
    @test isconcretetype(typeof(e))
    e = Edge{Int64}(2, 3)
    @test e == Edge((2, 3))
    @test typeof(e) == Edge{Int64}
end

@testset "Triangle" begin
    T = Triangle(1, 2, 3)
    T′ = Triangle((1, 2, 3))
    @test Triangle{Int64}((1, 2, 3)) == T
    @test T.indices == T′.indices
    @test T == T′
    @test T.indices == (1, 2, 3)
    @test indices(T) == (1, 2, 3)
    @test geti(T) == 1
    @test getj(T) == 2
    @test getk(T) == 3
    @test Triangle{Int} <: DT.AbstractTriangle{Int}
    @test typeof(T) == Triangle{Int}
    @test isconcretetype(typeof(T))
    @test DT.shift_triangle_indices_1(T) == (2, 3, 1)
    @test DT.shift_triangle_indices_2(T) == (3, 1, 2)
    T = Triangle(4, 10, 2)
    @test DT.shift_triangle_indices_1(T) == (10, 2, 4)
    @test DT.shift_triangle_indices_2(T) == (2, 4, 10)
    @test DT.shift_triangle_indices_0(T) == (4, 10, 2)
    @test DT.shift_triangle_0(T) == T
    @test DT.shift_triangle_1(T) == Triangle(10, 2, 4)
    @test DT.shift_triangle_2(T) == Triangle(2, 4, 10)
    @test DT.shift_triangle_indices(T, 2) == (2, 4, 10)
    @test DT.shift_triangle_indices(T, 1) == (10, 2, 4)
    @test DT.shift_triangle_indices(T, 0) == (4, 10, 2)
    @test DT.shift_triangle(T, 0) == T
    @test DT.shift_triangle(T, 1) == Triangle(10, 2, 4)
    @test DT.shift_triangle(T, 2) == Triangle(2, 4, 10)
end

@testset "Point" begin
    p = Point(2.0, 3.0)
    @test p.coords == Float64[2.0, 3.0]
    @test getx(p) == 2.0
    @test gety(p) == 3.0
    @test coords(p) == Float64[2.0, 3.0]
    @test Point{Float64,Vector{Float64}} <: DT.AbstractPoint{Float64,Vector{Float64}}
    @test typeof(p) == Point{Float64,Vector{Float64}}
    @test isconcretetype(typeof(p))

    p = Point(2, 71)
    @test p.coords == Int[2, 71]
    @test getx(p) == 2
    @test gety(p) == 71
    @test coords(p) == Int[2, 71]
    @test typeof(p) == Point{Int64,Vector{Int64}}

    p = Point{Vector{Float64}}(2, 3)
    q = Point(2.0, 3.0)
    @test p.coords == q.coords
    @test typeof(p) == typeof(q) == Point{Float64,Vector{Float64}}

    p = Point((2, 3))
    @test typeof(p) == Point{Int64,NTuple{2,Int64}}
    @test getx(p) == 2
    @test gety(p) == 3
    @test coords(p) == (2, 3)
    @test p.coords == (2, 3)

    p = [2.0, 5.0]
    p = OffsetVector(p, -3)
    pt = Point(p)
    @test getx(pt) == 2.0
    @test gety(pt) == 5.0
    @test coords(pt) == p
    @test typeof(pt) == Point{Float64,typeof(p)}

    p = Point(5, 4)
    tp = Tuple(p)
    @test tp == (5, 4)

    p = Point(2.7, 13.01)
    tp = Tuple(p)
    @test tp == (2.7, 13.01)
end

############################################
##
## PRIMITIVE COLLECTIONS 
##
############################################
@testset "Triangles" begin
    T₁ = Triangle(1, 2, 3)
    T₂ = Triangle(5, 8, 2)
    T₃ = Triangle(17, 29, 15)
    T₄ = Triangle(-2, 0, 5)
    T₅ = Triangle(17, 30, 72)
    T₆ = Triangle(5, 50, 101)
    T_set = Set{Triangle{Int}}([T₁, T₂, T₃, T₄, T₅, T₆])
    T_set_struct = Triangles(T_set)
    @test triangles(T_set_struct) == T_set
    T_vec = [T₁, T₂, T₃, T₄, T₅, T₆]
    T_vec_struct = Triangles(T_vec)
    @test triangles(T_vec_struct) == T_set
    T_tup = [T₁, T₂, T₃, T₄, T₅, T₆]
    T_tup_struct = Triangles(T_tup)
    @test triangles(T_tup_struct) == T_set
    @test isconcretetype(typeof(T_tup_struct))

    T = Triangle(1, 2, 3)
    Ts = Triangles([T])
    DT.add_triangle!(Ts, Triangle(4, 5, 7))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(1, 2, 3), Triangle(4, 5, 7)])
    DT.delete_triangle!(Ts, Triangle(1, 2, 3))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7)])
    DT.add_triangle!(Ts, Triangle(5, 3, 2), Triangle(10, 11, 12), Triangle(13, 15, 19))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7), Triangle(5, 3, 2), Triangle(10, 11, 12), Triangle(13, 15, 19)])
    DT.delete_triangle!(Ts, Triangle(2, 5, 3), Triangle(19, 13, 15))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(4, 5, 7), Triangle(10, 11, 12)])
    DT.delete_triangle!(Ts, Triangle(4, 5, 7))
    @test Ts.triangles == Set{Triangle{Int64}}([Triangle(10, 11, 12)])
end

@testset "Points" begin
    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    @test p₁ ≠ p₂
    @test Point(2.0, 3.0) == Point(2.0, 3.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts_vec_struct = Points(pts_vec)
    @test points(pts_vec_struct) == pts_vec
    pts_tup = (p₁, p₂, p₃, p₄)
    pts_tup_struct = Points(pts_tup)
    @test points(pts_tup_struct) == pts_vec
    @test isconcretetype(typeof(pts_tup_struct))
    @test DT.get_point(pts_tup_struct, 1) == p₁
    @test DT.get_point(pts_tup_struct, 2) == p₂
    @test DT.get_point(pts_tup_struct, 3) == p₃
    @test DT.get_point(pts_tup_struct, 4) == p₄
    @test length(pts_tup_struct) == 4
    add_point!(pts_tup_struct, Point(5.7, 13.3))
    @test points(pts_tup_struct) == [p₁, p₂, p₃, p₄, Point(5.7, 13.3)]
    @test length(pts_tup_struct) == 5

    p0 = Float64[5, 5]
    p1 = Float64[4.5, 2.5]
    p2 = Float64[2.5, 1.5]
    p3 = Float64[3, 3.5]
    p4 = Float64[0, 2]
    p5 = Float64[1, 5]
    p6 = Float64[1, 3]
    p7 = Float64[4, -1]
    p8 = Float64[-1, 4]
    pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]
    qts= Points(pts)
    @test getx.(points(qts)) == first.(pts)
    @test gety.(points(qts)) == last.(pts)
end

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

@testset "HistoryDAG" begin
    dg = DirectedGraph{Triangle{Int64}}()
    add!(dg, Triangle(1, 2, 3))
    add!(dg, Triangle(4, 5, 6))
    add!(dg, Triangle(7, 8, 9))
    add!(dg, Triangle(1, 2, 3), Triangle(4, 5, 6))
    add!(dg, Triangle(1, 2, 3), Triangle(7, 8, 9))
    hg = DT.HistoryDAG(dg)
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
    @test Triangle(18, 19, 17) ∉ hg.graph.V
    DT.add_edge!(hg, Triangle(17, 18, 19), Triangle(22, 25, 29))
    @test Triangle(22, 25, 29) ∈ out_neighbors(hg, Triangle(17, 18, 19))
    DT.add_edge!(hg, Triangle(18, 19, 17), Triangle(1, 2, 3))
    @test Triangle(1, 2, 3) ∈ out_neighbors(hg, Triangle(17, 18, 19))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    @test Triangle(11, 13, 15) ∈ out_neighbors(hg, Triangle(7, 8, 9))
    og = out_deg(hg, Triangle(7, 8, 9))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    DT.add_edge!(hg, Triangle(9, 7, 8), Triangle(13, 15, 11))
    @test og == out_deg(hg, Triangle(7, 8, 9))

    H = DT.HistoryDAG()
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

    H = DT.HistoryDAG()
    HH = DT.HistoryDAG()
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

############################################
##
## PREDICATES
##
############################################
p0 = Float64[5, 5]
p1 = Float64[4.5, 2.5]
p2 = Float64[2.5, 1.5]
p3 = Float64[3, 3.5]
p4 = Float64[0, 2]
p5 = Float64[1, 5]
p6 = Float64[1, 3]
p7 = Float64[4, -1]
p8 = Float64[-1, 4]
pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]
pts = Points(pts)
@test orient(Triangle(4, 6, 7), pts) == 1
@test orient(Triangle(4, 7, 6), pts) == -1
@test orient(Triangle(4, 2, 3), pts) == -1
@test orient(Triangle(4, 7, 3), pts) == 1
@test orient(Triangle(5, 7, 9), pts) == 1
@test orient(Triangle(5, 9, 7), pts) == -1
@test orient(Triangle(3, 8, 5), pts) == -1
@test orient(Triangle(1, 2, 3), [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]) == 0

############################################
##
## UTILITY FUNCTIONS 
##
############################################
@testset "Can we ocrrectly find the root of a DAG?" begin
    𝒟 = DT.HistoryDAG()
    add!(DT.graph(𝒟), Triangle((1, 2, 3)))
    add!(DT.graph(𝒟), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)))
    add!(DT.graph(𝒟), Triangle((10, 11, 12)))
    add!(DT.graph(𝒟), Triangle((13, 14, 15)))
    add!(DT.graph(𝒟), Triangle((1, 2, 3)), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((1, 2, 3)), Triangle((7, 8, 9)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)), Triangle((10, 11, 12)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((4, 5, 6)), Triangle((13, 14, 15)))
    @test DT.find_root(𝒟; method=:brute) == Triangle((1, 2, 3))
    @test all(DT.find_root(𝒟; method=:rng) == Triangle((1, 2, 3)) for _ in 1:10)
end

@testset "Can we correctly compute whether a point is higher than another?" begin
    p = [2.0, 3.71] |> Point
    q = [2.0, 4.81] |> Point
    @test !DT.is_point_higher(p, q)
    q = [2.0, 3.60] |> Point
    @test DT.is_point_higher(p, q)
    q = [2.381, 3.71] |> Point
    @test DT.is_point_higher(p, q)
    p = [1.2999, 1.0] |> Point
    q = [1.981, 1.71] |> Point
    @test !DT.is_point_higher(p, q)
    @test DT.is_point_higher(q, p)
    p = [57.131, 4.0] |> Point
    q = [2.0, 3.1] |> Point
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
    p = [-2.31, 4.0] |> Point
    q = [5.0, 4.0] |> Point
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
end

@testset "Can we correctly sort points by height?" begin
    v = Vector{Point{Float64, Vector{Float64}}}(undef, 100)
    for _ in 1:500
        v .= [Point(rand(2)) for _ in 1:100]
        DT.partial_highest_point_sort!(v, 1)
        @test all(DT.is_point_higher(v[1], v[j]) for j in 2:lastindex(v))
    end
end

@testset "Can we count the number of negative values?" begi 
    v = [-1, -2, 3, 5, 4.0]
    @test DT.num_less(0, v) == 2 
    @test DT.num_less(7, v) == 5 
    v = [2, 5, 10, 29.0]
    @test DT.num_less(0, v) == 0 
    v = [2, 1, 5, 4]
    @test DT.num_less(10, v) == 4
    @test DT.num_less(5, v) == 2 
end