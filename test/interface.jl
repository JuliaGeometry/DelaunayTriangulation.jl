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

    pts = [
        (-2.26, -4.31),
        (-6.36, -1.55),
        (-3.32, 3.53)
    ]
    T = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 1, 2, 3, pts)
    @test T == (3, 2, 1)
    T = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 3, 2, 1, pts)
    @test T == (3, 2, 1)
    for _ in 1:5000
        pts = rand(SVector{2,Float64}, 3)
        T1 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 1, 2, 3, pts)
        T2 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 2, 3, 1, pts)
        T3 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 3, 1, 2, pts)
        T4 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 3, 2, 1, pts)
        T5 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 2, 1, 3, pts)
        T6 = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 1, 3, 2, pts)
        @test all(DT.isoriented(T, pts) == 1 for T in (T1, T2, T3, T4, T5, T6))
    end
end

@testset "Point" begin
    p = [1.0, 5.0]
    @test DT.getx(p) == 1.0
    @test DT.gety(p) == 5.0
    @test DT.number_type(p) == Float64
    p = (1.388182, 5.0001)
    @test DT.getx(p) == 1.388182
    @test DT.gety(p) == 5.0001
    @test DT.number_type(p) == Float64
    p = (1.398219f0, 5.09291f0)
    @test DT.number_type(p) == Float32
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
    @test DT.number_type(p) == Float64
    @test DT.get_point(p, 1) == (1.0, 2.0)
    @test DT.get_point(p, 2) == (5.0, 2.0)
    @test DT.get_point(p, 3) == (5.0, 1.0)
    @test DT.get_point(p, 4) == (17.0, 2.0)
    @test DT.get_point(p, 1, 2, 3, 4) == ((1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0))
    DT.add_point!(p, (1.0, 5.0))
    @test p == [(1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0), (1.0, 5.0)]
    DT.add_point!(p, (1.0, 2.0), (5.0, 17.0))
    @test p == [(1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0), (1.0, 5.0), (1.0, 2.0), (5.0, 17.0)]
    @test DT.number_type(p) == Float64

    p = [[1.0, 2.0], [5.0, 2.0], [5.0, 1.0], [17.0, 2.0]]
    @test DT.number_type(p) == Float64
    @test DT.get_point(p, 1) == (1.0, 2.0)
    @test DT.get_point(p, 2) == (5.0, 2.0)
    @test DT.get_point(p, 3) == (5.0, 1.0)
    @test DT.get_point(p, 4) == (17.0, 2.0)
    @test DT.get_point(p, 1, 2, 3, 4) == ((1.0, 2.0), (5.0, 2.0), (5.0, 1.0), (17.0, 2.0))

    p₁ = [2.0, 5.0]
    p₂ = [5.0, 1.7]
    p₃ = [2.2, 2.2]
    p₄ = [-17.0, 5.0]
    pts = [p₁, p₂, p₃, p₄]
    @test DT.number_type(p) == Float64
    @test DT.get_point(pts, DT.LowerRightBoundingIndex)[1] ≈ -6.0 + DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerRightBoundingIndex)[2] ≈ 3.35 - 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex)[1] ≈ -6.0 - DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex)[2] ≈ 3.35 - 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.UpperBoundingIndex)[1] ≈ -6.0 rtol = 1e-1
    @test DT.get_point(pts, DT.UpperBoundingIndex)[2] ≈ 3.35 + DT.BoundingTriangleShift * 22.0 rtol = 1e-1
    @test DT.get_point(pts, DT.LowerRightBoundingIndex) == DT.lower_right_bounding_triangle_coords(pts)
    @test DT.get_point(pts, DT.LowerLeftBoundingIndex) == DT.lower_left_bounding_triangle_coords(pts)
    @test DT.get_point(pts, DT.UpperBoundingIndex) == DT.upper_bounding_triangle_coords(pts)
    # @test_throws BoundsError DT.get_point(pts, 0)
    @test_throws BoundsError DT.get_point(pts, -5)
    @test_throws BoundsError DT.get_point(pts, 17)

    pts = [(1.0f0, 2.0f0), (2.1f0, 2.9f0)]
    @test DT.number_type(pts) == Float32
end
