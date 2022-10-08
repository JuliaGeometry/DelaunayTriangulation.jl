############################################
##
## PRIMITIVES AND COLLECTIONS 
##
############################################
@testset "Edge" begin
    e = (1, 3)
    @test DT.construct_edge(NTuple{2,Int64}, 1, 3) == e
    e = (3, 1)
    @test DT.construct_edge(NTuple{2,Int64}, 3, 1) == e
    e = [1, 10]
    DT.construct_edge(::Type{Vector{Int64}}, i, j) = [i, j]
    @test DT.construct_edge(Vector{Int64}, 1, 10) == [1, 10]
end

@testset "Triangle" begin
    T = (1, 5, 10)
    @test DT.construct_triangle(NTuple{3,Int64}, 1, 5, 10) == T
    T = (5, 23, 15)
    @test DT.construct_triangle(NTuple{3,Int64}, 5, 23, 15) == T
    @test DT.indices(T) == (5, 23, 15)
    @test DT.geti(T) == 5
    @test DT.getj(T) == 23
    @test DT.getk(T) == 15
    @test DT.shift_triangle_1(T) == (23, 15, 5)
    @test DT.shift_triangle_2(T) == (15, 5, 23)
    @test DT.shift_triangle(T, 1) == (23, 15, 5)
    @test DT.shift_triangle(T, 2) == (15, 5, 23)
    @test DT.shift_triangle(T, 0) == (5, 23, 15)
end

@testset "Point" begin
    p = [1.0, 5.0]
    @test DT.getx(p) == 1.0
    @test DT.gety(p) == 5.0
    p = (1.388182, 5.0001)
    @test DT.getx(p) == 1.388182
    @test DT.gety(p) == 5.0001
end

@testset "Triangles" begin
    T = DT.construct_triangles(Set{NTuple{3,Int64}})
    @test T == Set{NTuple{3,Int64}}([])
    DT.add_triangle!(T, (2, 3, 1), (7, 10, 15), (18, 19, 23))
    @test T == Set{NTuple{3,Int64}}([(2, 3, 1), (7, 10, 15), (18, 19, 23)])
    DT.add_triangle!(T, (7, 5, 2))
    @test T == Set{NTuple{3,Int64}}([(2, 3, 1), (7, 10, 15), (18, 19, 23), (7, 5, 2)])
    DT.delete_triangle!(T, (7, 5, 2))
    @test T == Set{NTuple{3,Int64}}([(2, 3, 1), (7, 10, 15), (18, 19, 23)])
    DT.delete_triangle!(T, (2, 3, 1), (7, 10, 15))
    @test T == Set{NTuple{3,Int64}}([(18, 19, 23)])
    @test DT.triangle_type(typeof(T)) == NTuple{3,Int64}
end

@testset "Points" begin
    p = [(1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0)]
    @test DT.get_point(p, 1) == (1.0, 2.0)
    @test DT.get_point(p, 2) == (5.0, 2.0)
    @test DT.get_point(p, 3) == (5.0, 1.0)
    @test DT.get_point(p, 4) == (17.0, 2.0)
    @test DT.get_point(p, 1, 2, 3, 4) == ((1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0))
    DT.add_point!(p, (1.0, 5.0))
    @test p == [(1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0), (1.0, 5.0)]
    DT.add_point!(p, (1.0, 2.0), (5.0, 17.0))
    @test p == [(1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0), (1.0, 5.0), (1.0, 2.0), (5.0, 17.0)]

    p₁ = [2.0, 5.0]
    p₂ = [5.0, 1.7]
    p₃ = [2.2, 2.2]
    p₄ = [-17.0, 5.0]
    pts = [p₁, p₂, p₃, p₄]
    @test DT.get_point(pts, DT.LowerRightBoundingIndex)[1] ≈ -6.0 + DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerRightBoundingIndex)[2] ≈ 3.35 - 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex)[1] ≈ -6.0 - DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex)[2] ≈ 3.35 - 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.UpperBoundingIndex)[1] ≈ -6.0 rtol = 1e-1
    @test DT.get_point(pts, DT.UpperBoundingIndex)[2] ≈ 3.35 + DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == DT.lower_right_bounding_triangle_coords(pts)
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == DT.lower_left_bounding_triangle_coords(pts)
    @test DT.get_point(pts, DT.UpperBoundingIndex) == DT.upper_bounding_triangle_coords(pts)
    @test_throws BoundsError DT.get_point(pts, 0)
    @test_throws BoundsError DT.get_point(pts, -5)
    @test_throws BoundsError DT.get_point(pts, 17)


end

############################################
##
## CONSTRUCTORS FOR THE TRIANGULATION DATA STRUCTURES 
##
############################################
@testset "Adjacent" begin
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    @test adj.adjacent == DefaultDict{NTuple{2,Int64},Int64}(DT.DefaultAdjacentValue)
    @test typeof(adj) == DT.Adjacent{Int64,NTuple{2,Int64}}
    @test isconcretetype(typeof(adj))
    adj = DT.Adjacent{Int16,NTuple{2,Int16}}()
    @test adj.adjacent == DefaultDict{NTuple{2,Int16},Int16}(Int16(DT.DefaultAdjacentValue))
    @test typeof(adj) == DT.Adjacent{Int16,NTuple{2,Int16}}
    @test isconcretetype(typeof(adj))
    @test DT.adjacent(adj) == DefaultDict{NTuple{2,Int16},Int16}(Int16(DT.DefaultAdjacentValue))
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    @test adj2v.adjacent2vertex == Dict{NTuple{2,16},Int64}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}
    @test isconcretetype(typeof(adj2v))
    adj2v = DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}()
    @test adj2v.adjacent2vertex == Dict{NTuple{2,Int16},Int16}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test isconcretetype(typeof(adj2v))
    @test DT.adjacent2vertex(adj2v) == Dict{Int64,Set{NTuple{2,Int16}}}()
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
    @test DT.graph(dg) == UndirectedGraph{Int16}()
    tvg = UndirectedGraph{Int64}()
    add!(tvg, 7)
    add!(tvg, 9)
    add!(tvg, 9, 7)
    dg = DT.DelaunayGraph(tvg)
    @test DT.graph(dg) == tvg
end

@testset "HistoryGraph" begin
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    @test hg.graph == DirectedGraph{NTuple{3,Int64}}()
    @test typeof(hg) == DT.HistoryGraph{NTuple{3,Int64}}
    @test isconcretetype(typeof(hg))
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    @test hg.graph == DirectedGraph{NTuple{3,Int64}}()
    @test typeof(hg) == DT.HistoryGraph{NTuple{3,Int64}}
    @test DT.graph(hg) == DirectedGraph{NTuple{3,Int64}}()
    dg = DirectedGraph{NTuple{3,Int64}}()
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
end

############################################
##
## INDEXING AND UPDATING OF THE TRIANGULATION DATA STRUCTURES
##
############################################
@testset "Adjacent" begin
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
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
    @test DT.get_edge(adj, (j, r)) == i
    DT.add_edge!(adj, (5, 7), 9)
    DT.add_edge!(adj, (7, 5), 3)
    @test DT.get_edge(adj, 5, 7) == DT.get_edge(adj, (5, 7)) == 9
    DT.delete_edge!(adj, 5, 7)
    @test DT.get_edge(adj, 5, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (5, 7)) == DT.DefaultAdjacentValue
    DT.delete_edge!(adj, i, j)
    @test DT.get_edge(adj, i, j) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, j, i) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (i, j)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (j, i)) == DT.DefaultAdjacentValue

    DT.add_edge!(adj, 16, 13, 5)
    DT.add_edge!(adj, 13, 16, DT.BoundaryIndex)
    DT.delete_edge!(adj, 16, 13)
    @test DT.get_edge(adj, (16, 13)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 13, 16) == DT.BoundaryIndex
    DT.add_edge!(adj, 16, 13, 5)
    DT.delete_edge!(adj, 16, 13; protect_boundary=false)
    @test DT.get_edge(adj, (16, 13)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (13, 16)) == DT.DefaultAdjacentValue
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    T1 = DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)
    T2 = DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6)
    T3 = DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)
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

    p1 = (-1.0, 4.0)
    p2 = (4.0, 6.0)
    p3 = (2.0, 2.0)
    p4 = (6.0, -3.0)
    p5 = (7.0, 3.0)
    pts = [p1, p2, p3, p4, p5]
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 3, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 4, 3))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 4, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 5, 2))
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
    DT.delete_edge!(adj, 2, 3; delete_uv_only=true)
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj, (2, 3)) == DT.DefaultAdjacentValue

    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_edge!(adj, 2, 3, 7)
    DT.add_edge!(adj, 5, 7, 17)
    @test collect(DT.edges(adj)) == [(2, 3), (5, 7)] || collect(DT.edges(adj)) == [(5, 7), (2, 3)]
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, 3, (1, 2))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2)])
    DT.add_edge!(adj2v, 10, 5, 7)
    @test DT.get_edge(adj2v, 10) == Set{NTuple{2,Int64}}([(5, 7)])
    DT.add_edge!(adj2v, 3, (1, 4))
    DT.add_edge!(adj2v, 3, (10, 11))
    DT.add_edge!(adj2v, 3, (13, 5))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2), (1, 4), (10, 11), (13, 5)])
    DT.delete_edge!(adj2v, 3, (1, 4))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2), (10, 11), (13, 5)])
    DT.delete_edge!(adj2v, 3, 1, 2)
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(10, 11), (13, 5)])
    DT.delete_point!(adj2v, 10)
    @test_throws KeyError DT.get_edge(adj2v, 10)

    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
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
    @test DT.get_edge(adj2v, i) == Set{NTuple{2,Int64}}([(k, r), (20, 30)])
    @test DT.get_edge(adj2v, j) == Set{NTuple{2,Int64}}([(r, k), (50, 70)])
    @test DT.get_edge(adj2v, k) == Set{NTuple{2,Int64}}([(r, i), (j, r), (100, 110)])
    @test DT.get_edge(adj2v, r) == Set{NTuple{2,Int64}}([(i, k), (k, j), (130, 150)])
    @test DT.get_edge(adj2v, 17) == Set{NTuple{2,Int64}}([(19, 32)])
    @test length(DT.adjacent2vertex(adj2v)) == 5

    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    i, j, k, r = 17, 13, 29, 9
    DT.add_edge!(adj2v, i, j, k)
    DT.add_edge!(adj2v, i, 170, 302)
    DT.add_edge!(adj2v, j, k, i)
    DT.add_edge!(adj2v, j, 500, 501)
    DT.add_edge!(adj2v, k, i, j)
    DT.add_edge!(adj2v, k, 177, 111)
    DT.add_edge!(adj2v, 37, 400, 1)
    DT.update_after_insertion!(adj2v, i, j, k, r)
    @test DT.get_edge(adj2v, i) == Set{NTuple{2,Int64}}([(j, r), (r, k), (170, 302)])
    @test DT.get_edge(adj2v, j) == Set{NTuple{2,Int64}}([(k, r), (r, i), (500, 501)])
    @test DT.get_edge(adj2v, k) == Set{NTuple{2,Int64}}([(i, r), (r, j), (177, 111)])
    @test DT.get_edge(adj2v, r) == Set{NTuple{2,Int64}}([(j, k), (k, i), (i, j)])
    @test DT.get_edge(adj2v, 37) == Set{NTuple{2,Int64}}([(400, 1)])
    @test length(DT.adjacent2vertex(adj2v)) == 5
end

@testset "DelaunayGraph" begin
    DG = DT.DelaunayGraph{Int64}()
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
    dg = DirectedGraph{NTuple{3,Int64}}()
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    hg = DT.HistoryGraph(dg)
    @test out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == Set([DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)])
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == [DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)]
    @test in_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == 0
    @test in_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == 1
    @test out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == 2
    @test out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == 0
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 11, 13, 15))
    @test DT.construct_triangle(NTuple{3,Int64}, 11, 13, 15) ∈ hg.graph.V
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19), DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29))
    @test DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19) ∈ hg.graph.V
    @test DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29) ∈ hg.graph.V
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 18, 19, 17))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19), DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29))
    @test DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29) ∈ out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 18, 19, 17), DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    @test DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3) ∈ out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    og = out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    @test og == out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))

    H = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    Htrue = deepcopy(H)
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 5, 6, 4))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 6, 4, 5))
    @test Htrue.graph == H.graph
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 2))
    @test Htrue.graph == H.graph

    H = DT.HistoryGraph{NTuple{3,Int64}}()
    HH = DT.HistoryGraph{NTuple{3,Int64}}()
    T₁ = DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)
    T₂ = DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6)
    T₃ = DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)
    DT.add_triangle!(H, T₁, T₂, T₃)
    DT.add_triangle!(HH, T₁)
    DT.add_triangle!(HH, T₂)
    DT.add_triangle!(HH, T₃)
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1)) # same triangle as T₁
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8)) # same as T₁ and T₃
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
@testset "Orientation of a triangle's points" begin
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
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 6, 7), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 7, 6), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 2, 3), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 7, 3), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 5, 7, 9), pts) == 1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 5, 9, 7), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 3, 8, 5), pts) == -1
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]) == 0
    pts = [[0, -3.0], [3.0, 0.0], [0.0, 3.0], [-3.0, 0.0]]
    @test DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), pts) == DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 2, 3, 4), pts) == DT.isoriented(DT.construct_triangle(NTuple{3,Int64}, 4, 1, 2), pts) == 1
    p₁ = [5.7044025422189, 1.801603986463]
    p₂ = [8.3797127128527, 5.8924221838871]
    p₃ = [2.8875415689061, 6.2038339497809]
    @test orient(p₁, p₂, p₃) == 1
    @test orient(p₁, p₂, [10.0, 4.0]) == -1
    p₁ = [5.0, 1.0]
    p₂ = [5.0, 6.0]
    @test orient(p₁, p₂, [5.0, 5.0]) == 0
    @test orient(p₁, p₂, [5.0, 2.0]) == 0
    @test orient(p₂, p₁, [5.0, 2.0]) == 0
end

@testset "Testing if points are in a given circle" begin
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
    @test DT.isincircle(pts, 5, 7, 6, 9) == 1
    @test DT.isincircle(pts, 5, 7, 6, 3) == -1
    @test DT.isincircle(pts, 5, 7, 6, 3) == -1
    @test DT.isincircle(pts, 5, 7, 6, 6) == 0
    @test DT.isincircle(pts, 3, 2, 1, 4) == 1
    @test DT.isincircle(pts, 3, 2, 1, 6) == 1
    @test DT.isincircle(pts, 3, 2, 1, 7) == 1
    @test DT.isincircle(pts, 3, 2, 1, 5) == -1
    @test DT.isincircle(pts, 3, 2, 1, 8) == -1
    @test incircle(p4, p6, p5, p8) == 1
end

@testset "Testing if a point is to the left of a line" begin
    A = [4.6, 3.2]
    B = [3.2, 2.2]
    C = [3.4, 3.2]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    C = [5.8, 3.6]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == -1
    A = [1.0, 7.0]
    B = [1.0, 1.0]
    C = [1.0, 5.0]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 0
    @test DT.isleftofline([A, B, C], 3, 1, 2) == 0
    A = [2.123933267613, 7.1892809338214]
    B = [-1.5542939635314, 3.3935384556756]
    C = [2.8172732249214, 5.085758012496]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == -1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == 1
    C = [-2.8172732249214, 5.085758012496]
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    @test DT.isleftofline([A, B, C], 3, 2, 1) == 1
    @test DT.isleftofline([A, B, C], 3, 1, 2) == -1
    D = [-6.3, -2.77]
    E = [-3.46, -2.07]
    F = [-5.22, -1.41]
    @test DT.isleftofline([D, E, F], 1, 2, 3) == 1
    @test DT.isleftofline([D, E, F], 2, 1, 3) == -1
    F = [-4.88, -4.57]
    @test DT.isleftofline([D, E, F], 1, 2, 3) == -1
    @test DT.isleftofline([D, E, F], 2, 1, 3) == 1
end

@testset "Testing if a point is to the left of a line defined by edges of the bounding triangle" begin
    p1 = [2.0, 3.0]
    p2 = [5.7, 2.3]
    p3 = [17.0, -2.0]
    p = [0.0, 0.0]
    pts = [p1, p2, p3, p]
    i, j = DT.LowerRightBoundingIndex, DT.UpperBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
    i, j = DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.UpperBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
    i, j = DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == 1
    i, j = DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex
    @test DT.isleftofline(pts, i, j, 4) == -1
end

@testset "Testing if we can find where a point lies to a segment connecting to the bounding triangle" begin
    pᵢ = [17.0, 2.0]
    p = [2.0, -2.5]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [-8.7999865, 21.89396]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [-12.49835, -10.738]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == -1
    pᵢ = [20.0, 20.0]
    p = [0.0, 20.0]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerRightBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerRightBoundingIndex, 1, 2) == 1
    pᵢ = [20.569333, 8.405]
    p = [8.743, -17.13]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerLeftBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.LowerLeftBoundingIndex, 1, 2) == -1
    p = [-98.23, 26.9]
    pts = [pᵢ, p]
    @test DT.isleftofline(pts, 1, DT.LowerLeftBoundingIndex, 2) == -1
    @test DT.isleftofline(pts, DT.LowerLeftBoundingIndex, 1, 2) == 1
    @test DT.isleftofline(pts, 1, DT.UpperBoundingIndex, 2) == 1
    @test DT.isleftofline(pts, DT.UpperBoundingIndex, 1, 2) == -1
end

@testset "Testing Boolean intriangle values" begin
    e = [
        1 1 1
        1 1 0
        1 1 -1
        1 0 1
        1 0 0
        1 0 -1
        1 -1 1
        1 -1 0
        1 -1 -1
        0 1 1
        0 1 0
        0 1 -1
        0 0 1
        0 0 0
        0 0 -1
        0 -1 1
        0 -1 0
        0 -1 -1
        -1 1 1
        -1 1 0
        -1 1 -1
        -1 0 1
        -1 0 0
        -1 0 -1
        -1 -1 1
        -1 -1 0
        -1 -1 -1
    ]
    ev = [1, 0, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, -1, 0, -1]
    @test all(DT.isintriangle(e[i, 1], e[i, 2], e[i, 3]) == ev[i] for i in eachindex(ev))
end

@testset "Testing that a point is in a triangle, outside, or on the edge" begin
    p1 = Float64[5, 5]
    p2 = Float64[4.5, 2.5]
    p3 = Float64[2.5, 1.5]
    p4 = Float64[3, 3.5]
    p5 = Float64[0, 2]
    p6 = Float64[1, 5]
    p7 = Float64[1, 3]
    p8 = Float64[4, -1]
    p9 = Float64[-1, 4]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    T1 = DT.construct_triangle(NTuple{3,Int64}, 4, 1, 6)
    T2 = DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    T3 = DT.construct_triangle(NTuple{3,Int64}, 3, 2, 4)
    T4 = DT.construct_triangle(NTuple{3,Int64}, 8, 1, 2)
    T5 = DT.construct_triangle(NTuple{3,Int64}, 8, 2, 3)
    T6 = DT.construct_triangle(NTuple{3,Int64}, 8, 3, 5)
    T7 = DT.construct_triangle(NTuple{3,Int64}, 5, 3, 7)
    T8 = DT.construct_triangle(NTuple{3,Int64}, 3, 4, 7)
    T9 = DT.construct_triangle(NTuple{3,Int64}, 5, 7, 9)
    T10 = DT.construct_triangle(NTuple{3,Int64}, 7, 6, 9)
    T11 = DT.construct_triangle(NTuple{3,Int64}, 7, 4, 6)
    T = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11]
    pt1 = (
        A=(3.9551298987095, 4.7489988803935),
        B=(3.2585303811361, 4.3003415639903),
        C=(2.69180534989, 4.4066025073489),
        D=(3.2231100666833, 4.6781582514877),
        E=(2.0424329182538, 4.8198395092992)
    )
    pt2 = (
        A=(3.5, 3.5),
        B=(3.8384340736083, 3.2777159548861),
        C=(4.173119090818, 3.0606229707501),
        D=(4.0736181397556, 3.8475850382431),
        E=(4.390212074954, 3.6938108411468),
        F=(4.390212074954, 4.2546343834981),
        G=(4.7520337151806, 4.2998620885264),
        H=(4.7520337151806, 4.661683728753)
    )
    pt3 = (
        A=(3.114790793155, 2.9520764786821),
        B=(3.3952025643307, 2.8163933635972),
        C=(3.3137926952797, 2.4545717233705),
        D=(2.951971055053, 2.3098430672799),
        E=(2.9157888910304, 2.4726628053818),
        F=(3.0786086291324, 2.6535736254952),
        G=(3.901752860648, 2.6626191665008),
        H=(3.0, 2.0),
        I=(2.7258325299114, 1.8213838529739)
    )
    pt4 = (
        A=(4.5781396673455, 2.7004728027825),
        B=(4.6264236360138, 2.9231155471972),
        C=(4.6693427192745, 3.1618529478347),
        D=(4.70153203172, 3.3871781349531),
        E=(4.5325381413811, 2.437593417811),
        F=(4.4708419591939, 2.1988560171735),
        G=(4.4520648602673, 2.0352270122422),
        H=(4.3501320375233, 1.3082850395147)
    )
    pt5 = (
        A=(3.433718968984, 1.1294534677817),
        B=(3.7382811110409, 1.4137114670348),
        C=(3.8499538964617, 0.9162599683419),
        D=(3.6875207540314, 0.5304812550699),
        E=(3.0377881843101, 1.2614303960064),
        F=(3.7484331824428, -0.0481868148381),
        G=(4.0, 2.0)
    )
    pt6 = (
        A=(2.8956591846836, 0.4086563982472),
        B=(2.4286639001964, 0.90610789694),
        C=(2.0936455439339, 1.2309741818007),
        D=(1.5149774740259, 1.4035593956329),
        E=(0.824636618697, 1.586296680867),
        F=(1.3322401887918, 1.2715824674082),
        G=(1.6063461166429, 1.0177806823609),
        H=(2.0225810441206, 0.713218540304),
        I=(2.3169911147756, 0.4391126124529),
        J=(2.6317053282343, 0.1954628988074),
        K=(3.1697651125348, -0.1395554574552)
    )
    pt7 = (
        A=(1.0581342609406, 2.2766375361959),
        B=(0.9363094041178, 2.550743464047),
        C=(1.2916319031842, 2.3070937504015),
        D=(1.7281709734657, 1.9923795369428),
        E=(1.0, 2.0),
        F=(0.5809869050515, 2.1142043937655)
    )
    pt8 = (
        A=(1.5454336882315, 2.845153534702),
        B=(1.9921248299149, 2.5405913926451),
        C=(2.1545579723453, 2.936522177319),
        D=(2.1951662579528, 2.4187665358224),
        E=(2.4185118287945, 2.8756097489077),
        F=(2.4185118287945, 2.5101351784394),
        G=(2.4489680430002, 2.0431398939523),
        H=(2.6317053282343, 2.9162180345153)
    )
    pt9 = (
        A=(-0.5458930205588, 3.5557985328346),
        B=(0.0733833349568, 3.2512363907778),
        C=(0.3170330486022, 2.9771304629266),
        D=(0.0835354063587, 2.6421121066641),
        E=(0.0, 2.4187665358224),
        F=(0.3576413342098, 2.4695268928319),
        G=(-0.2210267356982, 2.9872825343285),
        H=(-0.4849805921475, 3.2410843193759),
        I=(0.5099224052383, 2.9162180345153)
    )
    pt10 = (
        A=(0.3576413342098, 4.1649228169483),
        B=(0.0, 4.0),
        C=(0.4794661910326, 3.6573192468536),
        D=(0.7028117618743, 3.3527571047967),
        E=(0.6520514048648, 4.4796370304071),
        F=(-0.3530036639228, 4.1243145313408),
        G=(0.0, 3.7689920322744)
    )
    pt11 = (
        A=(1.3931526172031, 4.0735541743313),
        B=(2.0022769013168, 3.718231675265),
        C=(1.3931526172031, 3.5354943900309),
        D=(2.1444059009434, 3.4339736760119),
        E=(1.3017839745861, 4.3172038879768),
        F=(1.3017839745861, 4.5507015302204),
        G=(1.7992354732789, 3.9923376031161),
        H=(1.6875626878581, 3.5151902472271),
        I=(1.4337609028107, 3.809600317882)
    )
    test_pts = [pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11]
    for i in eachindex(T)
        @test DT.isoriented(T[i], pts) == 1
        for r in eachindex(test_pts[i])
            push!(pts, [test_pts[i][r]...])
            @test DT.isintriangle(T[i], pts, length(pts)) == 1
            pop!(pts)
        end
    end
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.BoundingTriangle, pts, length(pts)) == 1
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.0, 5.7])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.7, 29.5])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [2.422, 188.2])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [172.3, 178.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [49.1, 1720.0])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
    push!(pts, [-2.02, 13.4])
    @test DT.isintriangle(DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, length(pts)) == 1
end

@testset "Test if edge is on the bounding triangle" begin
    @test DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    @test DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)
    @test !DT.edge_on_bounding_triangle(2, 7)
    @test !DT.edge_on_bounding_triangle(3, 9)
    @test !DT.edge_on_bounding_triangle(DT.LowerLeftBoundingIndex, 2)
    @test !DT.edge_on_bounding_triangle(20, DT.LowerLeftBoundingIndex)
    @test !DT.edge_on_bounding_triangle(10, DT.LowerRightBoundingIndex)
    @test !DT.edge_on_bounding_triangle(DT.LowerRightBoundingIndex, 18)
    @test !DT.edge_on_bounding_triangle(DT.UpperBoundingIndex, 8)
    @test !DT.edge_on_bounding_triangle(91, DT.UpperBoundingIndex)
end

@testset "Test if a given edge is legal" begin
    p1 = [-1.0, 4.0]
    p2 = [4.0, 6.0]
    p3 = [2.0, 2.0]
    p4 = [6.0, -3.0]
    p5 = [7.0, 3.0]
    p6 = [-2.5519459826976, 13.6106262700637]
    p7 = [-21.0502221073507, -23.3458204355075]
    p8 = [23.7055438911111, 19.1906513812123]
    p9 = [31.7813088302274, -22.9759380718838]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 3, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 4, 3))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 4, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 5, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 7, 4))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 7, 9, 4))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 4, 9, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 5, 9, 8))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 2, 5, 8))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 2, 8, 6))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 6, 1, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 6, 7, 1))
    @test DT.islegal(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, adj, pts)
    @test DT.islegal(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, adj, pts)
    @test DT.islegal(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, adj, pts)
    @test DT.islegal(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, adj, pts)
    @test all(DT.islegal(i, j, adj, pts) for (i, j) in ((1, 3), (3, 2), (2, 1), (1, 4), (4, 3), (3, 1), (3, 4), (4, 5), (5, 3), (3, 5), (5, 2), (2, 3)))
    @test_broken !DT.islegal(2, DT.LowerRightBoundingIndex, 17, 3, pts)
    @test_broken !DT.islegal(DT.LowerRightBoundingIndex, 2, 17, 3, pts)
    @test_broken DT.islegal(17, 5, DT.LowerRightBoundingIndex, 5, pts)
    @test_broken DT.islegal(20, 21, 17, DT.LowerRightBoundingIndex, pts)
    @test_broken !DT.islegal(2, DT.LowerLeftBoundingIndex, 1117, 36, pts)
    @test_broken !DT.islegal(DT.LowerLeftBoundingIndex, 223, 5117, 613, pts)
    @test_broken DT.islegal(147, 235, DT.LowerLeftBoundingIndex, 15, pts)
    @test_broken DT.islegal(2610, 1221, 1517, DT.LowerLeftBoundingIndex, pts)
    @test_broken !DT.islegal(2, DT.UpperBoundingIndex, 517, 63, pts)
    @test_broken !DT.islegal(DT.UpperBoundingIndex, 2, 117, 36, pts)
    @test_broken DT.islegal(157, 51, DT.UpperBoundingIndex, 35, pts)
    @test_broken DT.islegal(205, 2121, 174, DT.UpperBoundingIndex, pts)
    @test_broken !DT.islegal(-2, 3, 1, -5, pts)
    @test_broken DT.islegal(-5, 53, 51, -3, pts)
    @test_broken !DT.islegal(4, -2, 51, -3, pts)
    @test_broken DT.islegal(4, -12, 51, -7, pts)
    @test_broken !DT.islegal(-2, 3, -5, 17, pts)
    @test_broken DT.islegal(-5, 53, -3, 18, pts)
    @test_broken !DT.islegal(4, -2, -3, 172, pts)
    @test_broken DT.islegal(4, -12, -7, 80, pts)
end

@testset "Can we find the correct edge that a point is on for a given triangle?" begin
    p1 = [2.0, 3.5]
    p2 = [0.0, 0.0]
    p3 = [3.0, 0.0]
    p4 = [17.2, -2.5]
    p5 = [0.0, 3.0]
    T = DT.construct_triangle(NTuple{3,Int64}, 2, 3, 5)
    pts = [p1, p2, p3, p4, p5]
    push!(pts, [1.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [2.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [1.5, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [1.0, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [0.5, 0.0])
    @test DT.find_edge(T, pts, length(pts)) == (2, 3)
    push!(pts, [2.5, 0.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [2.0, 1.0])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [1.5, 1.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [1.0, 2.0])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [0.5, 2.5])
    @test DT.find_edge(T, pts, length(pts)) == (3, 5)
    push!(pts, [0.0, 2.5])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 2.2])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 2.0])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 1.5])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 0.8])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
    push!(pts, [0.0, 0.2])
    @test DT.find_edge(T, pts, length(pts)) == (5, 2)
end

@testset "Testing if we can identify boundary and invalid edges" begin
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    tri = triangulate(pts; randomise=false)
    p7 = @SVector[2.0, 1.0]
    newtri = deepcopy(tri)
    push!(points(newtri), p7)
    DT.split_triangle!(triangles(newtri),
        pointlocation(newtri), adjacent(newtri),
        adjacent2vertex(newtri), graph(newtri),
        (1, 3, 5), 7)
    @test DT.is_boundary_edge((5, 2), adjacent(newtri))
    @test !DT.is_boundary_edge((2, 5), adjacent(newtri))
    @test DT.is_boundary_edge((2, 6), adjacent(newtri))
    @test DT.is_boundary_edge((6, 4), adjacent(newtri))
    @test DT.is_boundary_edge((4, 5), adjacent(newtri))
    @test !DT.is_boundary_edge((1, 7), adjacent(newtri))
    @test !DT.is_boundary_edge((7, 3), adjacent(newtri))
    @test !DT.is_boundary_edge((6, 3), adjacent(newtri))
    @test !DT.is_boundary_edge((3, 6), adjacent(newtri))
    @test all(DT.is_valid_edge(ij, adjacent(newtri)) for ij in edges((newtri)))
    @test !DT.is_valid_edge(7, 4, adjacent(newtri))
    @test !DT.is_valid_edge((2, 7), adjacent(newtri))
    @test !DT.is_valid_edge(7, 6, adjacent(newtri))
    @test !DT.is_valid_edge((7, 6), adjacent(newtri))
    @test !DT.is_valid_edge(2, 2, adjacent(newtri))
    @test !DT.is_valid_edge((4, 11), adjacent(newtri))
end

@testset "Testing if we can identify what triangles contain a point in their circumdisk" begin
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    tri = triangulate(pts; randomise=false)
    p11 = (6.0, 2.5)
    push!(pts, p11)
    bw_in = NTuple{3,Int64}[]
    for T in triangles(tri)
        DT.isincircle(T, pts, 11) == 1 && push!(bw_in, T)
    end
    @test bw_in == [
        (8, 1, 7),
        (10, 8, 9),
        (8, 7, 9),
        (10, 4, 8)
    ]
end

############################################
##
## UNCONSTRAINED INCREMENTAL TRIANGULATION
##
############################################
@testset "Initialisation" begin
    p1 = [2.7, 13.0]
    p2 = [-2.0, 3.0]
    p3 = [7.0, 2.0]
    p4 = [18.0, 3.5]
    pts = [p1, p2, p3, p4]
    Tri = UnconstrainedTriangulation(pts; method=:berg)
    @test adjacent(Tri) == Tri.adjacent
    @test adjacent2vertex(Tri) == Tri.adjacent2vertex
    @test graph(Tri) == Tri.graph
    @test pointlocation(Tri) == Tri.pointlocation
    @test triangles(Tri) == Tri.triangles
    @test points(Tri) == Tri.points
    @test triangles(Tri) == Set{NTuple{3,Int64}}([DT.BoundingTriangle])
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri, (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, (DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test neighbours(Tri, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test pointlocation(Tri).graph.N[DT.BoundingTriangle] == pointlocation(Tri).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == Set()
    @test points(Tri) == pts
    @test adjacent(Tri) isa DT.Adjacent{Int64,NTuple{2,Int64}}
    @test adjacent2vertex(Tri) isa DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}
    @test graph(Tri) isa DT.DelaunayGraph{Int64}
    @test pointlocation(Tri) isa DT.HistoryGraph{NTuple{3,Int64}}
    @test triangles(Tri) isa Set{NTuple{3,Int64}}
    @test points(Tri) isa Vector{Vector{Float64}}
    @test Tri isa UnconstrainedTriangulation{NTuple{3,Int64},NTuple{2,Int64},Set{NTuple{3,Int64}},Set{NTuple{2,Int64}},Vector{Vector{Float64}},Int64,DT.HistoryGraph{NTuple{3,Int64}}}
    @test typeof(Tri) <: DT.AbstractTriangulation
    @test typeof(Tri) <: DT.AbstractUnconstrainedTriangulation

    Tri3 = UnconstrainedTriangulation(pts; IntegerType=Int16, method=:berg)
    @test adjacent(Tri3) == Tri3.adjacent
    @test adjacent2vertex(Tri3) == Tri3.adjacent2vertex
    @test graph(Tri3) == Tri3.graph
    @test pointlocation(Tri3) == Tri3.pointlocation
    @test triangles(Tri3) == Tri3.triangles
    @test points(Tri3) == Tri3.points
    @test triangles(Tri3) == Set{NTuple{3,Int16}}([NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))])
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex))) == DT.UpperBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex))) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerRightBoundingIndex))) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex))) == DT.BoundaryIndex
    @test adjacent2vertex(Tri3, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))])
    @test adjacent2vertex(Tri3, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))])
    @test adjacent2vertex(Tri3, DT.UpperBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex))])
    @test neighbours(Tri3, DT.LowerLeftBoundingIndex) == Set{Int16}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.LowerRightBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.UpperBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test pointlocation(Tri3).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == pointlocation(Tri3).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == Set()
    @test points(Tri3) == pts
    @test adjacent(Tri3) isa DT.Adjacent{Int16,NTuple{2,Int16}}
    @test adjacent2vertex(Tri3) isa DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test graph(Tri3) isa DT.DelaunayGraph{Int16}
    @test pointlocation(Tri3) isa DT.HistoryGraph{NTuple{3,Int16}}
    @test triangles(Tri3) isa Set{NTuple{3,Int16}}
    @test points(Tri3) isa Vector{Vector{Float64}}
    @test Tri3 isa UnconstrainedTriangulation{NTuple{3,Int16},NTuple{2,Int16},Set{NTuple{3,Int16}},Set{NTuple{2,Int16}},Vector{Vector{Float64}},Int16,DT.HistoryGraph{NTuple{3,Int16}}}
    @test typeof(Tri) <: DT.AbstractTriangulation
    @test typeof(Tri) <: DT.AbstractUnconstrainedTriangulation
end

@testset "Testing that the point location data structure can correctly locate triangles" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))

    p = (33.373, 15.2287)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-31.0689, 52.90257)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (3.63, 1.679)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test flag == 1

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))

    p = (3.63, 1.679)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (27.706, 0.968)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (-13.6689, 1.3567)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = (-3.56804, 1.745279)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = (5.95, -2.91669)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (32.9507, -4.2764)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (41.4976, 46.81)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-33.48, 23.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))

    p = (-10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-36.59, 13.594)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (35.86, 34.379)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (15.66, 7.766)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (5.173, -3.305)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (-6.09, -3.305)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (14.57, 0.48686)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (9.89, 2.198)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (-5.735, 3.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-3.7957, 3.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-11.21, 2.54)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-3.68, 1.057)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3)
    @test flag == 1
    p = (0.0, 0.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test flag == 0
    p = (0.916, 1.408)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
    p = (0.0, 2.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
    p = (2.5057, 2.8986)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
end

@testset "Testing that we can correctly add a point into a triangulation" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)

    # Now do the actual test 
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == 4
    @test DT.get_edge(adj, 1, 4) == 3
    @test DT.get_edge(adj, 4, 3) == 1
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (3, 2),
        (4, 3),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 3),
        (DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (1, 4),
        (2, 1)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 3, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 1, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3])
    @test !DT.islegal(3, 1, adj, pts)
end

@testset "Testing that we can correctly flip an edge" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    p = (1.88, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    @test flag == 1
    p = (-0.9802, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    @test flag == 1
    p = (-3.642, 2.8615)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = (4.3036, -3.0977)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (4, 2),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 4),
        (DT.LowerRightBoundingIndex, 1),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 2),
        (2, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1, 4])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3, 2])
    @test DT.islegal(r, k, adj, pts)
end

@testset "Testing legalise edge" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    @test !DT.islegal(i, j, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    p = (1.88, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    @test flag == 1
    p = (-0.9802, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    @test flag == 1
    p = (-3.642, 2.8615)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = (4.3036, -3.0977)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (4, 2),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 4),
        (DT.LowerRightBoundingIndex, 1),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 2),
        (2, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1, 4])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3, 2])
    @test DT.islegal(r, k, adj, pts)
end

@testset "Testing that we can remove the bounding triangle" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    DT.remove_bounding_triangle!(T, adj, adj2v, DG)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    ])
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 2, 3) == DT.BoundaryIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, 2) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 1) == DT.BoundaryIndex
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 3, 4) == DT.BoundaryIndex
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test length(adj.adjacent) == 10
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (4, 2)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (1, 4),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (3, 2),
        (2, 1)
    ])
    @test DT.get_edge(adj2v, DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (4, 1),
        (3, 4),
        (2, 3)
    ])
    @test length(adj2v.adjacent2vertex) == 5
    @test DT.get_neighbour(DG, 1) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([1, 4, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([2, 3, 1])
end

@testset "Testing the steps for a small problem" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    pts = [p1, p2, p3, p4]

    DTri = DT.UnconstrainedTriangulation(pts; method=:berg)
    T, HG, adj, adj2v, DG, root = triangles(DTri), pointlocation(DTri),
    adjacent(DTri), adjacent2vertex(DTri), graph(DTri), DT.BoundingTriangle

    r = 1
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, r, root)
    @test Tᵢⱼₖ == DT.BoundingTriangle
    @test interior_flag == 1
    i, j, k = Tᵢⱼₖ
    @test i == DT.LowerRightBoundingIndex
    @test j == DT.UpperBoundingIndex
    @test k == DT.LowerLeftBoundingIndex
    DT.split_triangle!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, 1) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1),
        (1, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]

    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, r, root)
    i, j, k = Tᵢⱼₖ
    @test Tᵢⱼₖ == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test interior_flag == 1
    DT.split_triangle!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 2) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 2) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, 1) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 2, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 2, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, DT.LowerRightBoundingIndex) == 1
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1),
        (1, 2),
        (2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (2, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 1),
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, r, root)
    i, j, k = Tᵢⱼₖ
    DT.split_triangle!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test Tᵢⱼₖ == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    @test interior_flag == 1
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 2, 1, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 2) == 3
    @test DT.get_edge(adj, 2, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, 3) == 2
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1)
        (1, 2)
        (2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 3),
        (3, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (3, 2),
        (2, DT.LowerRightBoundingIndex),
        (DT.UpperBoundingIndex, 3)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (1, 3),
        (3, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (1, DT.UpperBoundingIndex),
        (DT.LowerLeftBoundingIndex, 2),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (2, 1)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, 3, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([1, 2, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    DTri = triangulate(pts; method=:berg)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)
end

@testset "Random triangulations" begin
    p1 = (-6.88, 3.61)
    p2 = (-6.08, -0.43)
    p3 = (-0.3, 2.01)
    p4 = (5.1, -1.27)
    p5 = (6.18, 1.87)
    p6 = (3.08, 4.43)
    p7 = (-1.34, 4.83)
    p8 = (-1.68, -0.77)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; method=:berg)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)
    pts = [[-6.88, 3.61], [-6.08, -0.43], [-0.3, 2.01], [5.1, -1.27], [6.18, 1.87], [3.08, 4.43], [-1.34, 4.83], [-1.68, -0.77]]
    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)

    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; randomise=false, trim=false, method=:berg)
    DTri2 = triangulate([p1, p2, p3, p4, p5, p6, p7]; randomise=false, trim=false, method=:berg)
    add_point!(DTri2, (-1.68, -0.77))
    @test triangles(DTri2) == triangles(DTri)

    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; randomise=false, trim=true, method=:berg)
    DTri2 = triangulate([p1, p2, p3, p4, p5, p6, p7]; randomise=false, trim=false, method=:berg)
    add_point!(DTri2, (-1.68, -0.77))
    DT.remove_bounding_triangle!(DTri2)
    @test triangles(DTri2) == triangles(DTri)

    for _ in 1:250
        x = rand(100)
        y = rand(100)
        pts = [(x, y) for (x, y) in zip(x, y)]
        DTri = triangulate(pts; method=:berg)
        @test DT.is_delaunay(DTri)
        @test all(DT.isoriented(T, pts) == 1 for T in triangles(DTri))
        DTri = triangulate(pts; method=:berg, trim=false)
        @test DT.is_delaunay(DTri)
        @test all(DT.isoriented(T, pts) == 1 for T in triangles(DTri))
        @test 1 == num_points(DTri) - num_edges(DTri) + num_triangles(DTri) # Euler's formula
    end
end

import CairoMakie: poly!, Figure, Axis, lines!, save
fig = Figure()
ax = Axis(fig[1, 1])
x = rand(100)
y = rand(100)
pts = [(x, y) for (x, y) in zip(x, y)]
DTri = triangulate(pts; method=:berg)
Tmat = zeros(Int64, num_triangles(DTri), 3)
for (i, T) in enumerate(triangles(DTri))
    Tmat[i, :] = [geti(T), getj(T), getk(T)]
end
pmat = zeros(num_points(DTri), 2)
for (i, p) in enumerate(points(DTri))
    pmat[i, :] = [getx(p), gety(p)]
end
poly!(ax, pmat, Tmat, strokewidth=2)
for (i, j) in adjacent2vertex(DTri, DelaunayTriangulation.BoundaryIndex)
    p = DT.get_point(DTri, i)
    q = DT.get_point(DTri, j)
    lines!(ax, [getx(p), getx(q)], [gety(p), gety(q)], color=:red, linewidth=5)
end
save("figures/test_triangulation.png", fig)

@testset "Custom integer types" begin
    p1 = (-6.88, 3.61)
    p2 = (-6.08, -0.43)
    p3 = (-0.3, 2.01)
    p4 = (5.1, -1.27)
    p5 = (6.18, 1.87)
    p6 = (3.08, 4.43)
    p7 = (-1.34, 4.83)
    p8 = (-1.68, -0.77)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; IntegerType=Int16, randomise=false, method=:berg)
    @test points(DTri) == pts
    @test adjacent(DTri) isa DT.Adjacent{Int16,NTuple{2,Int16}}
    @test adjacent2vertex(DTri) isa DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test graph(DTri) isa DT.DelaunayGraph{Int16}
    @test pointlocation(DTri) isa DT.HistoryGraph{NTuple{3,Int16}}
    @test triangles(DTri) isa Set{NTuple{3,Int16}}
    @test points(DTri) isa Vector{NTuple{2,Float64}}
    @test DTri isa UnconstrainedTriangulation{NTuple{3,Int16},NTuple{2,Int16},Set{NTuple{3,Int16}},
        Set{NTuple{2,Int16}},Vector{NTuple{2,Float64}},Int16}
    @test DT.is_delaunay(DTri)
    @test DT.get_point(DTri, 1) == p1
    @test DT.get_point(DTri, 2) == p2
    @test DT.get_point(DTri, 3) == p3
    @test DT.get_point(DTri, 4) == p4
    @test DT.get_point(DTri, 5) == p5
    @test DT.get_point(DTri, 6) == p6
    @test DT.get_point(DTri, 7) == p7
    @test DT.get_point(DTri, 8) == p8
end

############################################
##
## UTILITIES
##
############################################ 
@testset "Can we correctly find the root of a DAG?" begin
    D = DT.HistoryGraph{NTuple{3,Int64}}()
    add!(DT.graph(D), (1, 2, 3))
    add!(DT.graph(D), (4, 5, 6))
    add!(DT.graph(D), (7, 8, 9))
    add!(DT.graph(D), (10, 11, 12))
    add!(DT.graph(D), (13, 14, 15))
    add!(DT.graph(D), (1, 2, 3), (4, 5, 6))
    add!(DT.graph(D), (1, 2, 3), (7, 8, 9))
    add!(DT.graph(D), (7, 8, 9), (10, 11, 12))
    add!(DT.graph(D), (7, 8, 9), (4, 5, 6))
    add!(DT.graph(D), (4, 5, 6), (13, 14, 15))
    @test DT.find_root(D; method=:brute) == (1, 2, 3)
    @test all(DT.find_root(D; method=:rng) == (1, 2, 3) for _ in 1:10)
end

############################################
##
## DELETING A TRIANGLE 
##
############################################
@testset "Deleting a triangle" begin
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    tri = triangulate(pts; randomise=false)
    p7 = @SVector[2.0, 1.0]
    newtri = deepcopy(tri)
    push!(points(newtri), p7)
    DT.split_triangle!(triangles(newtri),
        pointlocation(newtri), adjacent(newtri),
        adjacent2vertex(newtri), graph(newtri),
        (1, 3, 5), 7)
    newtri2 = deepcopy(newtri)
    adj = adjacent(newtri2)
    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        1, 3, 7)
    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (6, 2, 3),
        (5, 3, 2),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 7
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2)
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == 5
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 3) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.is_boundary_edge(2, 6, adjacent(newtri2))
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 2, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
        (5, 3),
        (3, 6)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
        (2, 5),
        (6, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
        (3, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (2, 3),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (5, 2),
        (2, 6),
        (6, 4)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
        3, 6, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 2, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 2) == 5
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == 2
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        3, 2, 5)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (6, 2, 3),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 6
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2)
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.is_boundary_edge(2, 6, adjacent(newtri2))
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 2, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
        (3, 6)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
        (6, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (2, 3),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (2, 6),
        (6, 4),
        (5, 3),
        (3, 2)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
        3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        3, 6, 2)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 5
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2) &&
          (3, 6, 2) ∉ triangles(newtri2) &&
          (6, 2, 3) ∉ triangles(newtri2) &&
          (2, 3, 6) ∉ triangles(newtri2)
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test !DT.is_valid_edge(7, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 7, adjacent(newtri2))
    @test DT.is_valid_edge(5, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 5, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 2, adjacent(newtri2))
    @test DT.is_valid_edge(3, 6, adjacent(newtri2))
    @test !DT.is_valid_edge(6, 2, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 3, adjacent(newtri2))
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 6, adjacent(newtri2))
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (6, 4),
        (5, 3),
        (3, 6)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        5, 1, 7)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 4
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2) &&
          (3, 6, 2) ∉ triangles(newtri2) &&
          (6, 2, 3) ∉ triangles(newtri2) &&
          (2, 3, 6) ∉ triangles(newtri2) &&
          (5, 1, 7) ∉ triangles(newtri2) &&
          (1, 7, 5) ∉ triangles(newtri2) &&
          (7, 5, 1) ∉ triangles(newtri2)
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 1) == 6
    @test !DT.is_valid_edge(7, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 7, adjacent(newtri2))
    @test DT.is_valid_edge(5, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 5, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 2, adjacent(newtri2))
    @test DT.is_valid_edge(3, 6, adjacent(newtri2))
    @test !DT.is_valid_edge(6, 2, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(5, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 7, adjacent(newtri2))
    @test !DT.is_valid_edge(7, 5, adjacent(newtri2))
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 6, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (6, 4),
        (5, 3),
        (3, 6)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (5, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (1, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 5)
    @test (7, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 5) == DT.DefaultAdjacentValue
end