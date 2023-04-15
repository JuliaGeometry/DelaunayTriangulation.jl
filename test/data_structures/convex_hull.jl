using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

@testset "Basic getters" begin
    ch = ConvexHull(rand(2, 500), [1, 2, 4, 5, 6, 10, 29, 45, 71, 1])
    @test DT.get_points(ch) == ch.points
    @test DT.get_indices(ch) == ch.indices == [1, 2, 4, 5, 6, 10, 29, 45, 71, 1]
    @test num_points(ch) == 500
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
    @test ch == ConvexHull(points, [13, 1, 2, 3, 4, 5, 11, 12, 13])
    @test DT.get_indices(ch) == reverse([3, 2, 1, 13, 12, 11, 5, 4, 3])
    @test get_points(ch) == points
end