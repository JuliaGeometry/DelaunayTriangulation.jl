using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie
using StableRNGs

@testset "Random tests" begin
    rng = StableRNG(1)
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:100
            pts = rand(rng, 2, 38)
            tri = triangulate(pts; predicates=PT(), rng)
            @test validate_triangulation(tri; predicates=PT())
            _tri = DT.triangulate(pts; predicates=PT(), rng)
            @test tri == _tri
            __tri = retriangulate(_tri; predicates=PT(), rng)
            @inferred retriangulate(_tri; predicates=PT(), rng)
            @test __tri == _tri
        end
    end
end

@testset "Retriangulate should ignore deleted points" begin
    rng = StableRNG(2)
    points = [(0.0, 0.0), (0.87, 0.0), (1.0006, 0.7766), (0.0, 1.0), (0.5, 0.5)]
    tri = triangulate(points; skip_points=5, rng)
    _tri = retriangulate(tri; rng)
    @test tri == _tri && !DelaunayTriangulation.has_vertex(_tri, 5) && validate_triangulation(_tri)
end

@testset "Lots of collinearity" begin
    rng = StableRNG(3)
    for PT in (DT.AdaptiveKernel, DT.ExactKernel)
        _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true, predicates=PT())
        @test validate_triangulation(_tri; predicates=PT())
        for _ in 1:10
            tri = triangulate(_tri.points; predicates=PT(), rng)
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Signed zeros are duplicate points" begin
    # (0.0, 0.0) and (-0.0, 0.0) are == but hash differently
    points = [(0.0, 0.0), (-0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    @test !DT.points_are_unique(points)
    @test DT.find_duplicate_points(points) == Dict((0.0, 0.0) => [1, 2])
    tri = @test_logs (:warn, r"duplicate points") triangulate(points; rng = StableRNG(1))
    @test validate_triangulation(tri; print_result = false)
    @test DT.has_vertex(tri, 1) && !DT.has_vertex(tri, 2)
    @test DT.num_solid_vertices(tri) == 3
end
