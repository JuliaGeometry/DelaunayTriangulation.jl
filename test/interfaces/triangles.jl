using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArraysCore

i, j, k = 1, 5, 10
T1 = (i, j, k)
T2 = [i, j, k]
T3 = SVector{3,Int32}((i, j, k))

@test_throws "The" DT.construct_triangle(String, i, j, k)
@test DT.construct_triangle(NTuple{3,Int64}, i, j, k) == T1
@test DT.construct_triangle(Vector{Int64}, i, j, k) == T2
@test DT.construct_triangle(SVector{3,Int32}, i, j, k) == T3

@test_throws "The" DT.geti(String)
@test DT.geti(T1) == DT.geti(T2) == DT.geti(T3) == i

@test_throws "The" DT.getj(String)
@test DT.getj(T1) == DT.getj(T2) == DT.getj(T3) == j

@test_throws "The" DT.getk(String)
@test DT.getk(T1) == DT.getk(T2) == DT.getk(T3) == k

@test indices(T1) == indices(T2) == indices(T3) == (i, j, k)

@test_throws "The" DT.integer_type(String)
@test DT.integer_type(typeof(T1)) == DT.integer_type(typeof(T2)) == Int64
@test DT.integer_type(typeof(T3)) == typeof(T3[1]) == Int32

@test DT.triangle_edges(T1) ==
      DT.triangle_edges(T2) ==
      DT.triangle_edges(T3) ==
      DT.triangle_edges(i, j, k) ==
      ((i, j), (j, k), (k, i))

for T in (T1, T2, T3)
    for r in 0:2
        @test DT.rotate_triangle(T, r) ==
              DT.construct_triangle(typeof(T), ((i, j, k), (j, k, i), (k, i, j))[r+1]...)
    end
end

for T in (T1, T2, T3)
    @test DT.compare_triangles(T, (i, j, k))
    @test DT.compare_triangles(T, (j, k, i))
    @test DT.compare_triangles(T, (k, i, j))
    @test DT.compare_triangles(T, [i, j, k])
    @test !DT.compare_triangles(T, (k, j, i))
    @test !DT.compare_triangles(T, (1, 50, 323))
end

Ts1 = Set{typeof(T1)}(((1, 2, 3), (4, 6, 1), (10, 6, 1), (3, 10, 9), (5, 10, 3)))
Ts2 = Set{typeof(T2)}(([1, 2, 3], [4, 6, 1], [10, 6, 1], [3, 10, 9], [5, 10, 3]))
Ts3 = Set{typeof(T3)}((SVector{3,Int32}((1, 2, 3)),
    SVector{3,Int32}((4, 6, 1)),
    SVector{3,Int32}((10, 6, 1)),
    SVector{3,Int32}((3, 10, 9)),
    SVector{3,Int32}((5, 10, 3))))

@test_throws "The" DT.initialise_triangles(String)
for T in (T1, T2, T3)
    V = typeof(T)
    @test DT.initialise_triangles(Set{V}) == Set{V}()
end

@test_throws "The" DT.triangle_type(String)
for (T, Ts) in zip((T1, T2, T3), (Ts1, Ts2, Ts3))
    @test DT.triangle_type(typeof(Ts)) == typeof(T)
end
@test DT.triangle_type(Vector{NTuple{3,Int64}}) == NTuple{3,Int64}

@test_throws "The" num_triangles(String)
for Ts in (Ts1, Ts2, Ts3)
    @test num_triangles(Ts) == length(Ts) == 5
end

for Ts in (Ts1, Ts2, Ts3)
    V = eltype(Ts)
    V1 = DT.construct_triangle(V, 1, 2, 3)
    V2 = DT.construct_triangle(V, 6, 1, 4)
    V3 = DT.construct_triangle(V, 10, 6, 1)
    V4 = DT.construct_triangle(V, 9, 3, 10)
    V5 = DT.construct_triangle(V, 10, 3, 5)
    V6 = DT.construct_triangle(V, 10, 12, 5)
    V7 = DT.construct_triangle(V, 7, 17, 5)
    e1, p1 = DT.contains_triangle(V1, Ts)
    e2, p2 = DT.contains_triangle(V2..., Ts)
    e3, p3 = DT.contains_triangle(V3, Ts)
    e4, p4 = DT.contains_triangle(V4, Ts)
    e5, p5 = DT.contains_triangle(V5..., Ts)
    e6, p6 = DT.contains_triangle(V6, Ts)
    e7, p7 = DT.contains_triangle(V7..., Ts)
    @test all((p1, p2, p3, p4, p5))
    @test all((!p6, !p7))
    @test e1 == V1
    @test e2 == DT.construct_triangle(V, 4, 6, 1)
    @test e3 == V3
    @test e4 == DT.construct_triangle(V, 3, 10, 9)
    @test e5 == DT.construct_triangle(V, 5, 10, 3)
    @test e6 == V6
    @test e7 == V7
end

@test_throws "The" DT.add_to_triangles!(String, (1, 2, 3))
for Ts in (Ts1, Ts2, Ts3)
    V = eltype(Ts)
    V1 = DT.construct_triangle(V, 10, 12, 5)
    V2 = DT.construct_triangle(V, 7, 17, 5)
    V3 = DT.construct_triangle(V, 6, 17, 29)
    V4 = DT.construct_triangle(V, 13, 50, 101)
    V5 = DT.construct_triangle(V, 20, 25, 50)
    V6 = [17, 23, 507]
    V7 = (1, 100, 500)
    V8 = SVector{3,Int32}(100, 50, 901)
    G6 = DT.construct_triangle(V, V6...)
    G7 = DT.construct_triangle(V, V7...)
    G8 = DT.construct_triangle(V, V8...)
    DT.add_to_triangles!(Ts, V1)
    DT.add_triangle!(Ts, V2)
    DT.add_triangle!(Ts, V3, V4, V5)
    DT.add_triangle!(Ts, V6)
    DT.add_triangle!(Ts, V7)
    DT.add_to_triangles!(Ts, V8)
    for F in (V1, V2, V3, V4, V5, G6, G7, G8)
        @test F ∈ Ts
    end
    @test length(Ts) == 13
end

@test_throws "The" DT.delete_from_triangles!(String, (1, 2, 3))
for Ts in (Ts1, Ts2, Ts3)
    V = eltype(Ts)
    V1 = DT.construct_triangle(V, 10, 12, 5)
    V2 = DT.construct_triangle(V, 7, 17, 5)
    V3 = DT.construct_triangle(V, 6, 17, 29)
    V4 = DT.construct_triangle(V, 13, 50, 101)
    V5 = DT.construct_triangle(V, 20, 25, 50)
    G1 = DT.construct_triangle(V, 12, 5, 10)
    G2 = DT.construct_triangle(V, 5, 7, 17)
    G5 = DT.construct_triangle(V, 20, 25, 50)
    DT.delete_triangle!(Ts, G1)
    DT.delete_triangle!(Ts, V2)
    DT.delete_triangle!(Ts, V3)
    DT.delete_triangle!(Ts, V4...)
    DT.delete_triangle!(Ts, V5)
    for F in (V1, V2, V3, V4, V5)
        @test F ∉ Ts
    end
    @test length(Ts) == 8
end

@test_throws "The" DT.each_triangle(String)
for Ts in (Ts1, Ts2, Ts3)
    @test each_triangle(Ts) == Ts
end
T = [1 2 3; 4 5 6; 7 8 9]'
@test each_triangle(T) == eachcol(T)

points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
@test DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 1, 2, 3, points) ==
      (1, 2, 3)
@test DT.construct_positively_oriented_triangle(NTuple{3,Int64}, 2, 1, 3, points) ==
      (1, 2, 3)

T = Set{NTuple{3,Int64}}(((1, 2, 3),
    (2, 3, 4),
    (4, 5, 6),
    (6, 9, 11)))
V = Set{NTuple{3,Int64}}(((1, 2, 3),
    (2, 3, 4),
    (4, 5, 6),
    (6, 9, 11)))
@test DT.compare_triangle_collections(T, V)
V = Set{NTuple{3,Int64}}(((1, 2, 3),
    (2, 3, 4),
    (4, 5, 6)))
@test !DT.compare_triangle_collections(T, V)
V = Set{NTuple{3,Int64}}(((3, 1, 2),
    (3, 4, 2),
    (6, 4, 5),
    (6, 9, 11)))
@test DT.compare_triangle_collections(T, V)
V = Set{NTuple{3,Int64}}(((3, 1, 2),
    (3, 4, 2),
    (6, 4, 5),
    (6, 11, 9)))
@test !DT.compare_triangle_collections(T, V)

T = (1, 7, 5)
sortT = (1, 7, 5)
@test DT.sort_triangle(T) == sortT
T = [3, 1, 5]
@test DT.sort_triangle(T) == [1, 5, 3]
T = (5, 7, 1)
@test DT.sort_triangle(T) == (1, 5, 7)

T = Set{NTuple{3,Int64}}([
    (1, 2, 3),
    (10, 9, 3),
    (11, 5, 1),
    (193, 12, 10),
    (5, 3, 1),
    (19, 18, 17),
    (17, 5, 23),
    (20, 50, 72),
    (30, 31, 32),
    (20, 13, 37)
])
V = DT.sort_triangles(T)
@test V == Set{NTuple{3,Int64}}([
    (1, 2, 3),
    (3, 10, 9),
    (1, 11, 5),
    (10, 193, 12),
    (1, 5, 3),
    (17, 19, 18),
    (5, 23, 17),
    (20, 50, 72),
    (30, 31, 32),
    (13, 37, 20)
])
@test DT.compare_triangle_collections(T, V)

T = Set{NTuple{3,Int64}}([
    (1, 2, 3),
    (2, 3, 1),
    (3, 1, 2),
    (4, 5, 6),
    (11, 8, 3),
    (3, 11, 8),
    (18, 13, 1)
])
V = DT.remove_duplicate_triangles(T)
@test V == Set{NTuple{3,Int64}}([
    (1, 2, 3),
    (4, 5, 6),
    (1, 18, 13),
    (3, 11, 8),
])
T = Set{NTuple{3,Int64}}([
    (11, 9, 1),
    (1, 11, 9)
])
V = DT.remove_duplicate_triangles(T)
@test V == Set{NTuple{3,Int64}}([
    (1, 11, 9)
])