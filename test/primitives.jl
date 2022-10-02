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
    @test Tuple(p) === (17.0, 3.0)
    @test typeof(Tuple(p)) == NTuple{2, Float64}
end