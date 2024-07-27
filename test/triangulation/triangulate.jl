using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie

@testset "Random tests" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:100
            pts = rand(2, 38)
            tri = triangulate(pts; predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
            _tri = DT.triangulate(pts; predicates=PT())
            @test tri == _tri
            __tri = retriangulate(_tri; predicates=PT())
            @inferred retriangulate(_tri; predicates=PT())
            @test __tri == _tri
        end
    end
end

@testset "Retriangulate should ignore deleted points" begin
    points = [(0.0, 0.0), (0.87, 0.0), (1.0006, 0.7766), (0.0, 1.0), (0.5, 0.5)]
    tri = triangulate(points; skip_points=5)
    _tri = retriangulate(tri)
    @test tri == _tri && !DelaunayTriangulation.has_vertex(_tri, 5) && validate_triangulation(_tri)
end

@testset "Lots of collinearity" begin
    for PT in (DT.Adaptive, DT.Exact)
        _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true, predicates=PT())
        @test validate_triangulation(_tri; predicates=PT())
        for _ in 1:10
            tri = triangulate(_tri.points; predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end