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

    e = Edge{Int16}(7, 5)
    i, j = e
    @test i == 7 && j == 5
    @test_throws BoundsError i, j, k = e
    @test typeof(i) == typeof(j) == Int16
    @test length(e) == 2
    @test eltype(e) == Int16
    e = Edge{Int64}(17, 13)
    i, j = e
    @test i == 17 && j == 13
    @test_throws BoundsError i, j, k = e
    @test typeof(i) == typeof(j) == Int64
    @test length(e) == 2
    @test eltype(e) == Int64
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

    T = Triangle{Int16}(2, 3, 4)
    @test length(T) == 3
    @test eltype(T) == Int16

    T = Triangle(17, 32, 90)
    i, j, k = T
    @test i == 17
    @test j == 32
    @test k == 90
    @test length(T) == 3
    @test typeof(i) == typeof(j) == typeof(k) == Int64
    @test eltype(T) == Int64
    T = Triangle{Int16}(81, 101, 172)
    i, j, k = T
    @test i == 81
    @test j == 101
    @test k == 172
    @test length(T) == 3
    @test typeof(i) == typeof(j) == typeof(k) == Int16
    @test eltype(T) == Int16

    T = Triangle(17, 32, 90)
    @test T == T
    @test T == Triangle(32, 90, 17)
    @test T == Triangle(90, 17, 32)
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

    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    @test p₁ ≠ p₂
    @test Point(2.0, 3.0) == Point(2.0, 3.0)

    p = Point(2.0, 27.0)
    x, y = p
    @test x == 2.0
    @test y == 27.0
    @test length(p) == 2
    @test eltype(p) == Float64
    @test typeof(x) == typeof(y) == Float64

    p = Point{Float32,NTuple{2,Float32}}(17.0, 3.0)
    @test p == p
    @test p.coords == (17.0f0, 3.0f0)
    x, y = p
    @test x == 17.0
    @test y == 3.0
    @test length(p) == 2
    @test eltype(p) == Float32
    @test typeof(x) == typeof(y) == Float32
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
    T_set_param = Triangles{Int64,Triangle{Int64}}(T_set)
    @test T_set_param == T_set_struct
    @test Triangles(T₁, T₂, T₃, T₄, T₅, T₆).triangles == T_set

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

    T₁ = Triangle(1, 2, 3)
    T₂ = Triangle(5, 8, 2)
    T₃ = Triangle(17, 29, 15)
    T₄ = Triangle(-2, 0, 5)
    T₅ = Triangle(17, 30, 72)
    T₆ = Triangle(5, 50, 101)
    T_set = Set{Triangle{Int}}([T₁, T₂, T₃, T₄, T₅, T₆])
    T_tris = Triangles(T_set)
    T_testset = Set{Triangle{Int}}()
    for T in T_tris
        push!(T_testset, T)
    end
    @test T_set == T_testset
    @test length(T_tris) == 6
    @test eltype(T_set) == Triangle{Int64}
    @test eltype(Triangles{Int,Triangle{Int}}) == Triangle{Int64}
end

@testset "Points" begin
    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts_vec_struct = Points(pts_vec)
    @test size(pts_vec_struct) == (4,)
    @test eachindex(pts_vec_struct) == Base.OneTo(4)
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
    qts = Points(pts)
    @test getx.(points(qts)) == first.(pts)
    @test gety.(points(qts)) == last.(pts)
    rts = Points(p0, p1, p2, p3, p4, p5, p6, p7, p8)
    @test getx.(points(rts)) == first.(pts)
    @test gety.(points(rts)) == last.(pts)
    @test typeof(rts) == Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}

    pt_test = Vector{Point{Float64,Vector{Float64}}}()
    for pt in rts
        push!(pt_test, pt)
    end
    @test pt_test == rts.points
    @test length(rts) == 9
    @test eltype(rts) == Point{Float64,Vector{Float64}}
    @test eltype(Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}) == Point{Float64,Vector{Float64}}

    p₁ = Point(2.0, 5.0)
    p₂ = Point(5.0, 1.7)
    p₃ = Point(2.2, 2.2)
    p₄ = Point(-17.0, 5.0)
    pts_vec = [p₁, p₂, p₃, p₄]
    pts = Points(pts_vec)
    @test pts.points == pts_vec
    @test DT.xcentroid(pts) == (-17.0 + 5.0) / 2
    @test DT.ycentroid(pts) == (5.0 + 1.7) / 2
    @test DT.xmin(pts) == -17.0
    @test DT.xmax(pts) == 5.0
    @test DT.ymin(pts) == 1.7
    @test DT.ymax(pts) == 5.0
    @test DT.width(pts) == 22.0
    @test DT.height(pts) == 3.3
    @test DT.max_width_height(pts) == 22.0
    @test DT.get_point(pts, 1) == p₁
    @test DT.get_point(pts, 2) == p₂
    @test DT.get_point(pts, 3) == p₃
    @test DT.get_point(pts, 4) == p₄
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == Point(-6.0 + DT.BoundingTriangleShift * 22.0, 3.35 - 22.0)
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == Point(-6.0 - DT.BoundingTriangleShift * 22.0, 3.35 - 22.0)
    @test DT.get_point(pts, DT.UpperBoundingIndex) == Point(-6.0, 3.35 + DT.BoundingTriangleShift * 22.0)
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == pts.lower_right_bounding_triangle_coords
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == pts.lower_left_bounding_triangle_coords
    @test DT.get_point(pts, DT.UpperBoundingIndex) == pts.upper_bounding_triangle_coords
    @test_throws BoundsError DT.get_point(pts, 0)
    @test_throws BoundsError DT.get_point(pts, -5)
    @test_throws BoundsError DT.get_point(pts, 17)

    pts_copy = deepcopy(pts)
    Random.seed!(292991)
    shuffle!(pts_copy)
    pts_copy2 = deepcopy(pts_vec)
    Random.seed!(292991)
    shuffle!(pts_copy2)
    @test pts_copy.points == pts_copy2
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
    @test collect(DT.edges(adj)) == [Edge{Int64}(2, 3), Edge{Int64}(5, 7)]
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
@testset "Orientation of a triangle points" begin
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
    @test orient(Triangle(1, 2, 3), Points([[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]])) == 0
    pts = Points([0, -3.0], [3.0, 0.0], [0.0, 3.0], [-3.0, 0.0])
    @test orient(Triangle(1, 2, 3), pts) == orient(Triangle(2, 3, 4), pts) ==
          orient(Triangle(4, 1, 2), pts) == 1
    p₁ = Point(5.7044025422189, 1.801603986463)
    p₂ = Point(8.3797127128527, 5.8924221838871)
    p₃ = Point(2.8875415689061, 6.2038339497809)
    @test orient(p₁, p₂, p₃) == 1
    @test orient(p₁, p₂, Point(10.0, 4.0)) == -1
    p₁ = Point(5.0, 1.0)
    p₂ = Point(5.0, 6.0)
    @test orient(p₁, p₂, Point(5.0, 5.0)) == 0
    @test orient(p₁, p₂, Point(5.0, 2.0)) == 0
    @test orient(p₂, p₁, Point(5.0, 2.0)) == 0
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
    pts = Points(p0, p1, p2, p3, p4, p5, p6, p7, p8)
    @test incircle(pts, 5, 7, 6, 9) == 1
    @test incircle(pts, 5, 7, 6, 3) == -1
    @test incircle(pts, 5, 7, 6, 3) == -1
    @test incircle(pts, 5, 7, 6, 6) == 0
    @test incircle(pts, 3, 2, 1, 4) == 1
    @test incircle(pts, 3, 2, 1, 6) == 1
    @test incircle(pts, 3, 2, 1, 7) == 1
    @test incircle(pts, 3, 2, 1, 5) == -1
    @test incircle(pts, 3, 2, 1, 8) == -1
    @test incircle(Point(p4), Point(p6), Point(p5), Point(p8)) == 1
end

@testset "Testing if a point is to the left of a line" begin
    A = Point(4.6, 3.2)
    B = Point(3.2, 2.2)
    C = Point(3.4, 3.2)
    @test DT.leftofline(C, B, A) == 1
    @test DT.leftofline(C, A, B) == -1
    @test DT.leftofline(Points([A, B]), C, 2, 1) == 1
    @test DT.leftofline(Points([A, B]), C, 1, 2) == -1
    C = Point(5.8, 3.6)
    @test DT.leftofline(C, B, A) == -1
    A = Point(1.0, 7.0)
    B = Point(1.0, 1.0)
    C = Point(1.0, 5.0)
    @test DT.leftofline(C, B, A) == 0
    @test DT.leftofline(C, A, B) == 0
    @test DT.leftofline(Points(A, B), C, 2, 1) == 0
    @test DT.leftofline(Points(A, B), C, 1, 2) == 0
    A = Point(2.123933267613, 7.1892809338214)
    B = Point(-1.5542939635314, 3.3935384556756)
    C = Point(2.8172732249214, 5.085758012496)
    @test DT.leftofline(C, B, A) == -1
    @test DT.leftofline(C, A, B) == 1
    C = Point(-2.8172732249214, 5.085758012496)
    @test DT.leftofline(C, B, A) == 1
    @test DT.leftofline(C, A, B) == -1
    pts = Points([A, B])
    @test DT.leftofline(pts, C, 2, 1) == 1
    @test DT.leftofline(pts, C, 1, 2) == -1
end

@testset "Testing if a point is to the left of a line defined by edges of the bounding triangle" begin
    p1 = Point(2.0, 3.0)
    p2 = Point(5.7, 2.3)
    p3 = Point(17.0, -2.0)
    p = Point(0.0, 0.0)
    pts = Points(p1, p2, p3)
    i, j = DT.LowerRightBoundingIndex, DT.UpperBoundingIndex
    @test DT.leftofline(pts, p, i, j) == 1
    i, j = DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.leftofline(pts, p, i, j) == -1
    i, j = DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex
    @test DT.leftofline(pts, p, i, j) == 1
    i, j = DT.UpperBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.leftofline(pts, p, i, j) == -1
    i, j = DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex
    @test DT.leftofline(pts, p, i, j) == 1
    i, j = DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex
    @test DT.leftofline(pts, p, i, j) == -1
end

@testset "Testing if we can find where a point lies to a segment connecting to the bounding triangle" begin
    pᵢ = Point(17.0, 2.0)
    p = Point(2.0, -2.5)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerRightBoundingIndex) == -1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerRightBoundingIndex, 1) == 1
    pᵢ = Point(-8.7999865, 21.89396)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerRightBoundingIndex) == -1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerRightBoundingIndex, 1) == 1
    pᵢ = Point(-12.49835, -10.738)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerRightBoundingIndex) == 1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerRightBoundingIndex, 1) == -1
    pᵢ = Point(20.0, 20.0)
    p = Point(0.0, 20.0)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerRightBoundingIndex) == -1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerRightBoundingIndex, 1) == 1
    pᵢ = Point(20.569333, 8.405)
    p = Point(8.743, -17.13)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerLeftBoundingIndex) == 1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerLeftBoundingIndex, 1) == -1
    p = Point(-98.23, 26.9)
    @test DT.leftofline(Points(pᵢ), p, 1, DT.LowerLeftBoundingIndex) == -1
    @test DT.leftofline(Points(pᵢ), p, DT.LowerLeftBoundingIndex, 1) == 1
    @test DT.leftofline(Points(pᵢ), p, 1, DT.UpperBoundingIndex) == 1
    @test DT.leftofline(Points(pᵢ), p, DT.UpperBoundingIndex, 1) == -1
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
    @test all(DT.intriangle(e[i, 1], e[i, 2], e[i, 3]) == ev[i] for i in eachindex(ev))
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
    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8, p9)
    T1 = Triangle((4, 1, 6))
    T2 = Triangle((4, 2, 1))
    T3 = Triangle((3, 2, 4))
    T4 = Triangle((8, 1, 2))
    T5 = Triangle((8, 2, 3))
    T6 = Triangle((8, 3, 5))
    T7 = Triangle((5, 3, 7))
    T8 = Triangle((3, 4, 7))
    T9 = Triangle((5, 7, 9))
    T10 = Triangle((7, 6, 9))
    T11 = Triangle((7, 4, 6))
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
        @test orient(T[i], pts) == 1
        for p in test_pts[i]
            @test DT.intriangle(T[i], pts, Point(p)) == 1
        end
    end
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(2.0, 5.7)) == 1
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(-2.7, 29.5)) == 1
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(2.422, 188.2)) == 1
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(172.3, 178.0)) == 1
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(49.1, 1720.0)) == 1
    @test DT.intriangle(DT.BoundingTriangle, pts, Point(-2.02, 13.4)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(2.0, 5.7)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(-2.7, 29.5)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(2.422, 188.2)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(172.3, 178.0)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(49.1, 1720.0)) == 1
    @test DT.intriangle(Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex), pts, Point(-2.02, 13.4)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(2.0, 5.7)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(-2.7, 29.5)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(2.422, 188.2)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(172.3, 178.0)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(49.1, 1720.0)) == 1
    @test DT.intriangle(Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex), pts, Point(-2.02, 13.4)) == 1
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
    p1 = Point(-1.0, 4.0)
    p2 = Point(4.0, 6.0)
    p3 = Point(2.0, 2.0)
    p4 = Point(6.0, -3.0)
    p5 = Point(7.0, 3.0)
    p6 = Point(-2.5519459826976, 13.6106262700637)
    p7 = Point(-21.0502221073507, -23.3458204355075)
    p8 = Point(23.7055438911111, 19.1906513812123)
    p9 = Point(31.7813088302274, -22.9759380718838)
    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8, p9)
    adj = DT.Adjacent()
    DT.add_triangle!(adj, Triangle(1, 3, 2))
    DT.add_triangle!(adj, Triangle(1, 4, 3))
    DT.add_triangle!(adj, Triangle(3, 4, 5))
    DT.add_triangle!(adj, Triangle(3, 5, 2))
    DT.add_triangle!(adj, Triangle(1, 7, 4))
    DT.add_triangle!(adj, Triangle(7, 9, 4))
    DT.add_triangle!(adj, Triangle(4, 9, 5))
    DT.add_triangle!(adj, Triangle(5, 9, 8))
    DT.add_triangle!(adj, Triangle(2, 5, 8))
    DT.add_triangle!(adj, Triangle(2, 8, 6))
    DT.add_triangle!(adj, Triangle(6, 1, 2))
    DT.add_triangle!(adj, Triangle(6, 7, 1))
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

############################################
##
## MAIN TRIANGULATION FUNCTIONS
##
############################################
@testset "Initialisation" begin
    p1 = Point(2.7, 13.0)
    p2 = Point(-2.0, 3.0)
    p3 = Point(17.0, 2.0)
    p4 = Point(18.0, 3.5)
    pts = Points(p1, p2, p3, p4)
    Tri = DT.initialise_triangulation(pts)
    @test adjacent(Tri) == Tri.adjacent
    @test adjacent2vertex(Tri) == Tri.adjacent2vertex
    @test graph(Tri) == Tri.graph
    @test history(Tri) == Tri.history
    @test triangles(Tri) == Tri.triangles
    @test points(Tri) == Tri.points
    @test DT.root(Tri) == Tri.root
    @test DT.root(Tri) == DT.BoundingTriangle
    @test triangles(Tri).triangles == Set{Triangle{Int64}}([DT.BoundingTriangle])
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri, DT.UpperBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test neighbours(Tri, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test history(Tri).graph.N[DT.BoundingTriangle] == history(Tri).graph.N[DT.root(Tri)] == Set()
    @test points(Tri) == pts
    @test adjacent(Tri) isa DT.Adjacent{Int64,Edge{Int64}}
    @test adjacent2vertex(Tri) isa DT.Adjacent2Vertex{Int64,Edge{Int64}}
    @test graph(Tri) isa DT.DelaunayGraph{Int64}
    @test history(Tri) isa DT.HistoryDAG{Int64,Triangle{Int64}}
    @test triangles(Tri) isa Triangles{Int64,Triangle{Int64}}
    @test points(Tri) isa Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}
    @test DT.root(Tri) isa Triangle{Int64}
    @test Tri isa Triangulation{typeof(adjacent(Tri)),typeof(adjacent2vertex(Tri)),
        typeof(graph(Tri)),typeof(history(Tri)),typeof(triangles(Tri)),
        typeof(points(Tri)),typeof(DT.root(Tri))}

    pts = [p1, p2, p3, p4]
    Tri2 = DT.initialise_triangulation(pts)
    @test Tri2.adjacent.adjacent == Tri.adjacent.adjacent
    @test Tri2.points.points == Tri.points.points
    @test Tri2.adjacent2vertex.adjacent2vertex == Tri.adjacent2vertex.adjacent2vertex
    @test Tri2.graph.graph == Tri.graph.graph
    @test Tri2.history.graph == Tri.history.graph
    @test Tri2.triangles.triangles == Tri.triangles.triangles
    @test Tri2.root == Tri.root

    Tri3 = DT.initialise_triangulation(pts; IntegerType=Int16)
    @test adjacent(Tri3) == Tri3.adjacent
    @test adjacent2vertex(Tri3) == Tri3.adjacent2vertex
    @test graph(Tri3) == Tri3.graph
    @test history(Tri3) == Tri3.history
    @test triangles(Tri3) == Tri3.triangles
    @test points(Tri3) == Tri3.points
    @test DT.root(Tri3) == Tri3.root
    @test DT.root(Tri3) == Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    @test triangles(Tri3).triangles == Set{Triangle{Int16}}([Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri3, DT.LowerLeftBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri3, DT.LowerRightBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri3, DT.UpperBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test neighbours(Tri3, DT.LowerLeftBoundingIndex) == Set{Int16}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.LowerRightBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.UpperBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test history(Tri3).graph.N[Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)] == history(Tri3).graph.N[DT.root(Tri3)] == Set()
    @test points(Tri3).points == Points(pts).points
    @test adjacent(Tri3) isa DT.Adjacent{Int16,Edge{Int16}}
    @test adjacent2vertex(Tri3) isa DT.Adjacent2Vertex{Int16,Edge{Int16}}
    @test graph(Tri3) isa DT.DelaunayGraph{Int16}
    @test history(Tri3) isa DT.HistoryDAG{Int16,Triangle{Int16}}
    @test triangles(Tri3) isa Triangles{Int16,Triangle{Int16}}
    @test points(Tri3) isa Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}
    @test DT.root(Tri3) isa Triangle{Int16}
    @test Tri3 isa Triangulation{typeof(adjacent(Tri3)),typeof(adjacent2vertex(Tri3)),
        typeof(graph(Tri3)),typeof(history(Tri3)),typeof(triangles(Tri3)),
        typeof(points(Tri3)),typeof(DT.root(Tri3))}
end

@testset "Testing that the point location data structure can correctly locate triangles" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))

    p = Point(33.373, 15.2287)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-31.0689, 52.90257)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(3.63, 1.679)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test flag == 1

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))

    p = Point(3.63, 1.679)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(27.706, 0.968)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(-13.6689, 1.3567)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = Point(-3.56804, 1.745279)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = Point(5.95, -2.91669)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(32.9507, -4.2764)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(41.4976, 46.81)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-33.48, 23.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))

    p = Point(-10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-36.59, 13.594)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(35.86, 34.379)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(15.66, 7.766)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(5.173, -3.305)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(-6.09, -3.305)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(14.57, 0.48686)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(9.89, 2.198)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(-5.735, 3.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-3.7957, 3.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-11.21, 2.54)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-3.68, 1.057)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 3)
    @test flag == 1
    p = Point(0.0, 0.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test flag == 0
    p = Point(0.916, 1.408)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
    p = Point(0.0, 2.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
    p = Point(2.5057, 2.8986)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
end

@testset "Testing that we can correctly add a point into a triangulation" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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

    # Now do the actual test 
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, 2),
        Edge{Int64}(4, 3),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 3),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∉ triangles(T)
    @test Triangle(4, 3, 2) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(4, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(4, 3, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(3, 1, 4)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    @test out_neighbors(HG, Triangle(3, 2, 1)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    p = Point(1.88, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 2, 1)
    @test flag == 1
    p = Point(-0.9802, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 3, 2)
    @test flag == 1
    p = Point(-3.642, 2.8615)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = Point(4.3036, -3.0977)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)]
    @test in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test_throws KeyError DT.get_edge(adj, 3, 1) == 4
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
    @test_throws KeyError DT.get_edge(adj, 1, 3)
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(4, 2),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    @test !DT.islegal(i, j, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∉ triangles(T)
    @test Triangle(4, 3, 2) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(4, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(4, 3, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(3, 1, 4)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    @test out_neighbors(HG, Triangle(3, 2, 1)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    p = Point(1.88, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 2, 1)
    @test flag == 1
    p = Point(-0.9802, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 3, 2)
    @test flag == 1
    p = Point(-3.642, 2.8615)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = Point(4.3036, -3.0977)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)]
    @test in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test_throws KeyError DT.get_edge(adj, 3, 1) == 4
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
    @test_throws KeyError DT.get_edge(adj, 1, 3)
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(4, 2),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    DT.remove_bounding_triangle!(T, adj, adj2v, DG)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(4, 3, 2),
        Triangle(4, 2, 1)
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
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(4, 2)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(1, 4),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_edge(adj2v, DT.BoundaryIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(4, 1),
        Edge{Int64}(3, 4),
        Edge{Int64}(2, 3)
    ])
    @test length(adj2v.adjacent2vertex) == 5
    @test DT.get_neighbour(DG, 1) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([1, 4, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([2, 3, 1])
end

@testset "Testing the steps for a small problem" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    pts = Points(p1, p2, p3, p4)

    DTri = DT.initialise_triangulation(pts)
    T, HG, adj, adj2v, DG, root = triangles(DTri), history(DTri),
    adjacent(DTri), adjacent2vertex(DTri), graph(DTri), DT.root(DTri)

    r = 1
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    @test Tᵢⱼₖ == DT.BoundingTriangle
    @test interior_flag == 1
    i, j, k = Tᵢⱼₖ
    @test i == DT.LowerRightBoundingIndex
    @test j == DT.UpperBoundingIndex
    @test k == DT.LowerLeftBoundingIndex
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]

    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    @test Tᵢⱼₖ == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test interior_flag == 1
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(1, DT.LowerLeftBoundingIndex, 2)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1),
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, 1, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test Tᵢⱼₖ == Triangle(1, DT.LowerLeftBoundingIndex, 2)
    @test interior_flag == 1
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(3, 1, DT.UpperBoundingIndex),
        Triangle(2, 1, 3),
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1)
        Edge{Int64}(1, 2)
        Edge{Int64}(2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 3)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(1, 3),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.UpperBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, 3, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([1, 2, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)
end

DTri = triangulate(pts)
@test DT.is_delaunay(DTri)
DTri = triangulate(pts; trim=false)
@test DT.is_delaunay(DTri)
p1 = Point(-6.88, 3.61)
p2 = Point(-6.08, -0.43)
p3 = Point(-0.3, 2.01)
p4 = Point(5.1, -1.27)
p5 = Point(6.18, 1.87)
p6 = Point(3.08, 4.43)
p7 = Point(-1.34, 4.83)
p8 = Point(-1.68, -0.77)
pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
DTri = triangulate(pts)
@test DT.is_delaunay(DTri)
DTri = triangulate(pts; trim=false)
@test DT.is_delaunay(DTri)

for _ in 1:10000
    x = rand(100)
    y = rand(100)
    pts = Points([Point(x, y) for (x, y) in zip(x, y)])
    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false)
    @test DT.is_delaunay(DTri)
end
############################################
##
## UTILITY FUNCTIONS 
##
############################################
@testset "Can we correctly find the root of a DAG?" begin
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
    v = Vector{Point{Float64,Vector{Float64}}}(undef, 100)
    for _ in 1:500
        v .= [Point(rand(2)) for _ in 1:100]
        DT.partial_highest_point_sort!(v, 1)
        @test all(DT.is_point_higher(v[1], v[j]) for j in 2:lastindex(v))
    end
end

@testset "Can we count the number of negative values?" begin
    v = [-1, -2, 3, 5, 4.0]
    @test DT.num_less(0, v) == 2
    @test DT.num_less(7, v) == 5
    v = [2, 5, 10, 29.0]
    @test DT.num_less(0, v) == 0
    v = [2, 1, 5, 4]
    @test DT.num_less(10, v) == 4
    @test DT.num_less(5, v) == 3
end

@testset "Can we correctly identify a given triangulation as Delaunay?" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    pts = Points(p1, p2, p3, p4)
    DTri = DT.initialise_triangulation(pts)
    T, HG, adj, adj2v, DG, root = triangles(DTri), history(DTri),
    adjacent(DTri), adjacent2vertex(DTri), graph(DTri), DT.root(DTri)
    r = 1
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.is_delaunay(DTri)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test DT.is_delaunay(DTri)
end