using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie

_functions.jl")

@testset "Random tests" begin
    for _ in 1:100
        pts = rand(2, 38)
        tri = triangulate(pts)
        @test validate_triangulation(tri)
        _tri = DT.triangulate(pts)
        @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri)) &&
              get_adjacent(tri) == get_adjacent(_tri) &&
              get_adjacent2vertex(tri) == get_adjacent2vertex(_tri) &&
              get_graph(tri) == get_graph(_tri) &&
              get_convex_hull(tri) == get_convex_hull(_tri)
        __tri = retriangulate(_tri)
        @inferred retriangulate(_tri)
        DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri)) &&
            get_adjacent(tri) == get_adjacent(_tri) &&
            get_adjacent2vertex(tri) == get_adjacent2vertex(_tri) &&
            get_graph(tri) == get_graph(_tri) &&
            get_convex_hull(tri) == get_convex_hull(_tri)
    end
end

@testset "Retriangulate should ignore deleted points" begin 
    points = [(0.0, 0.0), (0.87, 0.0), (1.0006, 0.7766), (0.0, 1.0), (0.5, 0.5)]
    tri = triangulate(points; skip_points = 5)
    _tri = retriangulate(tri)
    @test tri == _tri && !DelaunayTriangulation.has_vertex(_tri, 5) && validate_triangulation(_tri)
end

@testset "Lots of collinearity" begin
    _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true)
    @test validate_triangulation(_tri)
    for _ in 1:10
        tri = triangulate(_tri.points)
        @test validate_triangulation(tri)
    end
end