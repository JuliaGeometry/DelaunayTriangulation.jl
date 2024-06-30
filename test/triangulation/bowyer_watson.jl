using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using StableRNGs
using CairoMakie
using DataStructures

Base.include(@__MODULE__, "../helper_functions.jl")

@testset "Getting the correct order" begin
    points = rand(2, 50)
    randomise = false
    skip_points = ()
    IntegerType = Int
    rng = Random.default_rng()
    pt_order = DT.get_insertion_order(points, randomise, skip_points, IntegerType, rng)
    @test pt_order == 1:50
    randomise = true
    pt_order = DT.get_insertion_order(points, randomise, skip_points, IntegerType, rng)
    @test pt_order ≠ 1:50
    @test all(∈(pt_order), DT.each_point_index(points))
    skip_points = (37, 48, 11)
    randomise = false
    pt_order = DT.get_insertion_order(points, randomise, skip_points, IntegerType, rng)
    @test pt_order == setdiff(1:50, skip_points)
    randomise = true
    pt_order = DT.get_insertion_order(points, randomise, skip_points, IntegerType, rng)
    @test pt_order ≠ setdiff(1:50, skip_points)
    @test all(∉(pt_order), skip_points)
    @test all(∈(pt_order), setdiff(DT.each_point_index(points), skip_points))
end

@testset "Getting the initial triangle" begin
    tri = triangulate_rectangle(0.0, 10.0, 0.0, 20.0, 11, 21)
    insertion_order = collect(DT.each_point_index(tri))
    insertion_order_orig = deepcopy(insertion_order)
    initial_triangle = DT.get_initial_triangle(tri, insertion_order)
    @test initial_triangle == (10, 11, 12)
    @test DT.is_positively_oriented(DT.triangle_orientation(tri, initial_triangle))
end

@testset "Initialising the Bowyer-Watson algorithm" begin
    tri = triangulate_rectangle(0.0, 10.0, 0.0, 20.0, 11, 21, delete_ghosts=false)
    _tri = DT.Triangulation(tri.points)
    DT.initialise_bowyer_watson!(_tri, DT.each_point_index(_tri) |> collect)
    S = filter(!DT.is_ghost_triangle, _tri.triangles)
    u, v, w = first(S)
    BI = DT.𝒢
    @test get_triangles(_tri) == Set(((u, v, w), (w, v, BI), (v, u, BI), (u, w, BI)))
    @test get_adjacent(get_adjacent(_tri)) == Dict((u, v) => w,
        (v, w) => u,
        (w, u) => v,
        (w, v) => BI,
        (v, BI) => w,
        (BI, w) => v,
        (v, u) => BI,
        (u, BI) => v,
        (BI, v) => u,
        (u, w) => BI,
        (w, BI) => u,
        (BI, u) => w)
    @test get_adjacent2vertex(get_adjacent2vertex(_tri)) ==
          Dict(BI => Set{NTuple{2,Int}}([(w, v), (v, u), (u, w)]),
        u => Set{NTuple{2,Int}}([(v, w), (BI, v), (w, BI)]),
        v => Set{NTuple{2,Int}}([(w, u), (BI, w), (u, BI)]),
        w => Set{NTuple{2,Int}}([(u, v), (v, BI), (BI, u)]))
    @test DT.get_graph(_tri).vertices == Set([u, v, w, BI])
end

@testset "A basic triangulation" begin
    a, b, c, d, e, f, g = [-1.78, 5.77], [-4.96, 0.31], [0.08, -3.73], [7.74, -3.03],
    [8.0, 3.0], [-0.6, 1.57], [3.58, 6.15]
    points = [a, b, c, d, e, f, g]
    rng = StableRNG(8881)
    tri = DT.triangulate(points; rng, delete_ghosts=false)
    @static if VERSION ≥ v"1.9"
        @inferred DT.triangulate(points; rng, delete_ghosts=false)
    end
    BI = DT.𝒢
    @test get_triangles(tri) == Set{NTuple{3,Int}}(((3, 2, BI),
        (4, BI, 5),
        (5, 7, 6),
        (3, 6, 2),
        (5, BI, 7),
        (4, 6, 3),
        (4, 3, BI),
        (4, 5, 6),
        (1, 2, 6),
        (1, 6, 7),
        (1, 7, BI),
        (1, BI, 2)))
    @test validate_triangulation(tri)
end

@testset "A random triangulation" begin
    pts = randn(2, 500)
    tri = DT.triangulate(pts; delete_ghosts=false)
    @test validate_triangulation(tri)
end

@testset "Lots of collinearity" begin
    _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true)
    @test validate_triangulation(_tri)
    for _ in 1:100
        tri = DT.triangulate(_tri.points)
        @test validate_triangulation(tri)
    end
end

@testset "A detailed example" begin
    @time for _ in 1:500
        T = Set{NTuple{3,Int}}(((2, 9, 8),
            (2, 8, 6),
            (2, 6, 1),
            (1, 6, 10),
            (6, 8, 10),
            (10, 8, 3),
            (10, 3, 4),
            (10, 4, 5),
            (1, 10, 5),
            (8, 9, 3),
            (9, 7, 3),
            (3, 7, 4),
            (1, 5, DT.𝒢),
            (5, 4, DT.𝒢),
            (4, 7, DT.𝒢),
            (7, 9, DT.𝒢),
            (9, 2, DT.𝒢),
            (2, 1, DT.𝒢)))
        T2 = Set{NTuple{3,Int}}(((2, 9, 8),
            (2, 8, 6),
            (2, 6, 1),
            (1, 6, 10),
            (6, 8, 10),
            (10, 8, 3),
            (10, 3, 4),
            (10, 4, 5),
            (1, 10, 5),
            (8, 9, 3),
            (9, 7, 3),
            (3, 7, 4)))
        adj = DT.Adjacent(
            Dict(@_adj(2, 9, 8)...,
                @_adj(2, 8, 6)...,
                @_adj(2, 6, 1)...,
                @_adj(1, 6, 10)...,
                @_adj(6, 8, 10)...,
                @_adj(10, 8, 3)...,
                @_adj(10, 3, 4)...,
                @_adj(10, 4, 5)...,
                @_adj(1, 10, 5)...,
                @_adj(8, 9, 3)...,
                @_adj(9, 7, 3)...,
                @_adj(3, 7, 4)...,
                @_adj(1, 5, DT.𝒢)...,
                @_adj(5, 4, DT.𝒢)...,
                @_adj(4, 7, DT.𝒢)...,
                @_adj(7, 9, DT.𝒢)...,
                @_adj(9, 2, DT.𝒢)...,
                @_adj(2, 1, DT.𝒢)...))
        adj2 = DT.Adjacent(
            Dict(@_adj(2, 9, 8)...,
                @_adj(2, 8, 6)...,
                @_adj(2, 6, 1)...,
                @_adj(1, 6, 10)...,
                @_adj(6, 8, 10)...,
                @_adj(10, 8, 3)...,
                @_adj(10, 3, 4)...,
                @_adj(10, 4, 5)...,
                @_adj(1, 10, 5)...,
                @_adj(8, 9, 3)...,
                @_adj(9, 7, 3)...,
                @_adj(3, 7, 4)...,
                (1, 5) => DT.𝒢,
                (5, 4) => DT.𝒢,
                (4, 7) => DT.𝒢,
                (7, 9) => DT.𝒢,
                (9, 2) => DT.𝒢,
                (2, 1) => DT.𝒢))
        adj2v = DT.Adjacent2Vertex(Dict(DT.𝒢 => Set{NTuple{2,Int}}(((1, 5),
                (5, 4),
                (4, 7),
                (7, 9),
                (9, 2),
                (2, 1))),
            1 => Set{NTuple{2,Int}}(((2, 6), (6, 10),
                (10, 5),
                (5, DT.𝒢),
                (DT.𝒢, 2))),
            2 => Set{NTuple{2,Int}}(((9, 8), (8, 6), (6, 1),
                (1, DT.𝒢),
                (DT.𝒢, 9))),
            3 => Set{NTuple{2,Int}}(((9, 7), (7, 4), (4, 10),
                (10, 8), (8, 9))),
            4 => Set{NTuple{2,Int}}(((5, 10), (10, 3),
                (3, 7),
                (7, DT.𝒢),
                (DT.𝒢, 5))),
            5 => Set{NTuple{2,Int}}(((1, 10), (10, 4),
                (4, DT.𝒢),
                (DT.𝒢, 1))),
            6 => Set{NTuple{2,Int}}(((2, 8), (8, 10),
                (10, 1), (1, 2))),
            7 => Set{NTuple{2,Int}}(((4, 3), (3, 9),
                (9, DT.𝒢),
                (DT.𝒢, 4))),
            8 => Set{NTuple{2,Int}}(((3, 10), (10, 6),
                (6, 2), (2, 9), (9, 3))),
            9 => Set{NTuple{2,Int}}(((7, 3), (3, 8), (8, 2),
                (2, DT.𝒢),
                (DT.𝒢, 7))),
            10 => Set{NTuple{2,Int}}(((8, 3), (3, 4), (4, 5),
                (5, 1), (1, 6),
                (6, 8)))))
        adj2v2 = DT.Adjacent2Vertex(Dict(DT.𝒢 => Set{NTuple{2,Int}}(((1, 5),
                (5, 4),
                (4, 7),
                (7, 9),
                (9, 2),
                (2, 1))),
            1 => Set{NTuple{2,Int}}(((2, 6), (6, 10),
                (10, 5))),
            2 => Set{NTuple{2,Int}}(((9, 8), (8, 6),
                (6, 1))),
            3 => Set{NTuple{2,Int}}(((9, 7), (7, 4),
                (4, 10), (10, 8),
                (8, 9))),
            4 => Set{NTuple{2,Int}}(((5, 10), (10, 3),
                (3, 7))),
            5 => Set{NTuple{2,Int}}(((1, 10), (10, 4))),
            6 => Set{NTuple{2,Int}}(((2, 8), (8, 10),
                (10, 1), (1, 2))),
            7 => Set{NTuple{2,Int}}(((4, 3), (3, 9))),
            8 => Set{NTuple{2,Int}}(((3, 10), (10, 6),
                (6, 2), (2, 9),
                (9, 3))),
            9 => Set{NTuple{2,Int}}(((7, 3), (3, 8),
                (8, 2))),
            10 => Set{NTuple{2,Int}}(((8, 3), (3, 4),
                (4, 5), (5, 1),
                (1, 6), (6, 8)))))
        A = zeros(Int, 10, 10)
        A[1, [2, 6, 10, 5]] .= 1
        A[2, [1, 6, 8, 9]] .= 1
        A[3, [8, 9, 7, 4, 10]] .= 1
        A[4, [5, 10, 3, 7]] .= 1
        A[5, [1, 10, 4]] .= 1
        A[6, [2, 8, 10, 1]] .= 1
        A[7, [4, 3, 9]] .= 1
        A[8, [6, 2, 9, 3, 10]] .= 1
        A[9, [2, 8, 3, 7]] .= 1
        A[10, [1, 6, 8, 3, 4, 5]] .= 1
        B = zeros(Int, 11, 11)
        B[2:end, 2:end] .= A
        B[1, [1, 5, 4, 7, 9, 2].+1] .= 1
        B[[1, 5, 4, 7, 9, 2].+1, 1] .= 1
        graph = _make_graph_from_adjacency(B, Dict(1 => DT.𝒢, (2:11 .=> 1:10)...))
        a = [1.5, 4.0]
        b = [0.0, 3.5]
        c = [2.0, 1.5]
        d = [3.0, 2.5]
        e = [2.5, 3.5]
        f = [0.5, 3.0]
        g = [2.5, -2.0]
        h = [0.5, 1.5]
        i = [0.0, 0.5]
        j = [1.5, 3.0]
        pts = [a, b, c, d, e, f, g, h, i, j]
        ch = DT.ConvexHull(pts, [1, 2, 9, 7, 4, 5, 1])
        tri = DT.triangulate(pts; delete_ghosts=false)
        @test DT.compare_triangle_collections(T, get_triangles(tri)) &&
              get_adjacent(tri) == adj &&
              get_adjacent2vertex(tri) == adj2v &&
              DT.get_graph(tri) == graph &&
              get_convex_hull(tri) == ch
        @test validate_triangulation(tri)
        delete_ghost_triangles!(tri)
        @test DT.compare_triangle_collections(T2, get_triangles(tri)) &&
              get_adjacent(tri) == adj2 &&
              get_adjacent2vertex(tri) == adj2v2 &&
              get_graph(tri) == graph &&
              get_convex_hull(tri) == ch
        @test validate_triangulation(tri)
    end
    a = [1.5, 4.0]
    b = [0.0, 3.5]
    c = [2.0, 1.5]
    d = [3.0, 2.5]
    e = [2.5, 3.5]
    f = [0.5, 3.0]
    g = [2.5, -2.0]
    h = [0.5, 1.5]
    i = [0.0, 0.5]
    j = [1.5, 3.0]
    pts = [a, b, c, d, e, f, g, h, i, j]
    tri = DT.triangulate(pts; delete_ghosts=false)
    @test validate_triangulation(tri)
end

@testset "Issue #94" begin
    i = 265
    rng = StableRNG(i)
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5), (0.2, 0.8), (0.1, 0.785)]
    tri = triangulate(points; rng)
    @test validate_triangulation(tri)
end