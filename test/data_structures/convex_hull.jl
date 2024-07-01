using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

@testset "Basic getters" begin
    ch = DT.ConvexHull(rand(2, 500), [1, 2, 4, 5, 6, 10, 29, 45, 71, 1])
    @test DT.get_points(ch) == ch.points
    @test DT.get_vertices(ch) == ch.vertices == [1, 2, 4, 5, 6, 10, 29, 45, 71, 1]
end

@testset "Specific example" begin
    a = [-0.7, 4.99]
    b = [-4.36, 2.39]
    c = [-4.5, -1.43]
    d = [-0.76, -4.75]
    e = [5.42, -4.85]
    f = [2.88, -0.79]
    g = [-0.2, -1.55]
    h = [-2.24, 0.59]
    i = [0.0, 2.21]
    j = [4.08, 2.19]
    k = [6.02, -1.51]
    ℓ = [5.0, 4.0]
    m = [2.02, 4.85]
    points = [a, b, c, d, e, f, g, h, i, j, k, ℓ, m]
    ch = convex_hull(points)
    @test ch == ch
    @test ch == DelaunayTriangulation.ConvexHull(points, [13, 1, 2, 3, 4, 5, 11, 12, 13])
    @test DT.get_vertices(ch) == reverse([3, 2, 1, 13, 12, 11, 5, 4, 3])
    @test get_points(ch) == points
end

@testset "Broken example" begin
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (1.0, 2.0), (1.0, 3.0), (1.0, 4.0), (1.0, 5.0), (1.0, 6.0)]
    ch = convex_hull(points)
    @test DT.circular_equality(ch.vertices, [1, 2, 3, 4, 5, 6, 7, 8, 1])
end

@testset "Corner cases" begin
    ch = convex_hull([(0.0, 0.0)])
    @test ch.points == [(0.0, 0.0)] && ch.vertices == [1]

    ch = convex_hull([(0.0, 0.0), (1.0, 0.0)])
    @test ch.points == [(0.0, 0.0), (1.0, 0.0)] && ch.vertices == [1, 2, 1]

    ch = convex_hull([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
    @test ch.points == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)] && ch.vertices == [1, 2, 3, 1]
    ch = convex_hull([(0.0, 0.0), (1.0, 1.0), (1.0, 0.0)])
    @test ch.points == [(0.0, 0.0), (1.0, 1.0), (1.0, 0.0)] && ch.vertices == [1, 3, 2, 1]
end

@testset "Random tests" begin
    for _ in 1:500
        pts = rand(2, 6)
        tri = triangulate(pts)
        ch = convex_hull(pts)
        @test ch == get_convex_hull(tri)

        tri = triangulate_rectangle(0, 1, 0, 1, 11, 11)
        ch = convex_hull(tri.points)
        @test ch == get_convex_hull(tri)

        A = DT.polygon_features(ch.points, ch.vertices)[1]
        @test A ≥ 0
    end
end

@testset "empty!" begin
    tri = triangulate(rand(2,50))
    ch = get_convex_hull(tri)
    @test !isempty(DT.get_vertices(ch))
    empty!(ch)
    @test isempty(DT.get_vertices(ch))
end

@testset "Issue #109" begin
    points = rand(2, 50)
    @test convex_hull(points).vertices == convex_hull(vcat(points, points)).vertices
end