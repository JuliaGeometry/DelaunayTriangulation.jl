using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StatsBase
using StructEquality

include("../helper_functions.jl")

DT.empty_representative_points!()

@struct_equal DT.RepresentativeCoordinates

@testset "Initialise" begin
    c = DT.RepresentativeCoordinates{Int64,Float64}()
    @test c == DT.RepresentativeCoordinates(0.0, 0.0, 0)
end

@testset "Given coordinates" begin
    c = DT.RepresentativeCoordinates(0.39182, 0.5919198, 5)
    @test getx(c) == c.x == 0.39182
    @test gety(c) == c.y == 0.5919198
    @test DT.getn(c) == c.n == 5
end

@testset "Conversion" begin
    c = DT.RepresentativeCoordinates(0.39182, 0.5919198, 5)
    newc = convert(DT.RepresentativeCoordinates{Int32,Float32}, c)
    @test newc.x isa Float32
    @test newc.y isa Float32
    @test newc.n isa Int32
    @test getx(newc) == 0.39182f0
    @test gety(newc) == 0.5919198f0
    @test DT.getn(newc) == 5
end

@testset "Resetting" begin
    c = DT.RepresentativeCoordinates(0.39182, 0.5919198, 5)
    DT.reset!(c)
    @test getx(c) == gety(c) == DT.getn(c) == 0.0
end

@testset "Adding and deleting points to representative coordinates" begin
    c = DT.RepresentativeCoordinates{Int64,Float64}()
    pts = [15randn(2) for _ in 1:200]
    foreach(p -> DT.add_point!(c, p), pts)
    mx, my = mean(pts)
    @test getx(c) ≈ mx
    @test gety(c) ≈ my
    @test DT.getn(c) == 200

    foreach(p -> DT.delete_point!(c, p), pts[151:200])
    mx, my = mean(pts[1:150])
    @test getx(c) ≈ mx
    @test gety(c) ≈ my
    @test DT.getn(c) == 150

    DT.compute_centroid!(c, pts)
    mx, my = mean(pts)
    @test getx(c) ≈ mx
    @test gety(c) ≈ my
    @test DT.getn(c) == 200
end

@testset "Getting a new representative point" begin
    DT.new_representative_point(2, Float64)
    @test 2 ∈ keys(DT.RepresentativePointList)
    @test DT.RepresentativeCoordinates(0.0, 0.0, 0) == DT.RepresentativePointList[2]
    DT.new_representative_point(3, Float32)
    @test 3 ∈ keys(DT.RepresentativePointList)
    @test DT.RepresentativeCoordinates(0.0, 0.0, 0) == DT.RepresentativePointList[3]

    DT.update_centroid_after_addition!(2, [2.0, 5.0])
    @test DT.RepresentativePointList[2] == DT.RepresentativeCoordinates(2.0, 5.0, 1)
    DT.update_centroid_after_addition!(3, [-17.0, 0.0])
    @test DT.RepresentativePointList[3] == DT.RepresentativeCoordinates(-17.0, 0.0, 1)
end

@testset "Using get! with a new representative point" begin
    centroid = get!(DT.RepresentativeCoordinates{Int64,Float64}, DT.RepresentativePointList,
        13)
    centroid.x = 0.292
    centroid.y = 0.991
    @test DT.RepresentativePointList[13] == DT.RepresentativeCoordinates(0.292, 0.991, 0)

    DT.update_centroid_after_addition!(7, [13.3, -5.0])
    @test DT.RepresentativePointList[7] == DT.RepresentativeCoordinates(13.3, -5.0, 1)
end

@testset "Resetting all representative points" begin
    DT.reset_representative_points!()
    @test DT.RepresentativeCoordinates(0.0, 0.0, 0) == DT.RepresentativePointList[2]
    @test DT.RepresentativeCoordinates(0.0, 0.0, 0) == DT.RepresentativePointList[3]

    DT.empty_representative_points!()
    @test_throws KeyError DT.RepresentativePointList[2]
    @test_throws KeyError DT.RepresentativePointList[3]
end

@testset "Representative points from a triangulation" begin
    tri = generate_mesh(-2.0, 2.0, -2.0, 2.0, 0.05; single_boundary=true,
        convert_result=true)
    DT.compute_representative_points!(tri)
    @test DT.RepresentativePointList[1].x ≈ 0.0 atol = 1e-12
    @test DT.RepresentativePointList[1].y ≈ 0.0 atol = 1e-12
    @test DT.RepresentativePointList[1].n == 0

    tri = generate_mesh(-2.0, 2.0, -2.0, 2.0, 0.05; single_boundary=false,
        convert_result=true)
    DT.compute_representative_points!(tri)
    @test DT.RepresentativePointList[1].x ≈ 0.0 atol = 1e-12
    @test DT.RepresentativePointList[1].y ≈ 0.0 atol = 1e-12
    @test DT.RepresentativePointList[1].n == 0
end

@testset "Polylabel, representative points, and distance" begin
    x, y, x1, x2, x3, x4, x5, y1, y2, y3, y4, y5 = complicated_geometry()
    tri = generate_mesh(x, y, 0.1; add_ghost_triangles=true, convert_result=true)
    for i in (1, 2)
        local x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
        if i == 1
            DT.compute_representative_points!(tri)
        else
            DT.compute_representative_points!(get_points(tri), get_boundary_nodes(tri))
        end
        x1, y1 = polylabel(tri.points, tri.boundary_nodes)
        x2, y2 = polylabel(tri.points, tri.boundary_nodes[2])
        x3, y3 = polylabel(tri.points, tri.boundary_nodes[3])
        x4, y4 = polylabel(tri.points, tri.boundary_nodes[4])
        x5, y5 = polylabel(tri.points, tri.boundary_nodes[5])
        d1 = DT.distance_to_polygon((x1, y1), tri.points, tri.boundary_nodes)
        d2 = DT.distance_to_polygon((x2, y2), tri.points, tri.boundary_nodes[2])
        d3 = DT.distance_to_polygon((x3, y3), tri.points, tri.boundary_nodes[3])
        d4 = DT.distance_to_polygon((x4, y4), tri.points, tri.boundary_nodes[4])
        d5 = DT.distance_to_polygon((x5, y5), tri.points, tri.boundary_nodes[5])
        @test all(≥(0), (d1, d2, d3, d4, d5)) # all points are inside their respective regions
        @test DT.RepresentativePointList[1].x ≈ x1 atol = 1e-14
        @test DT.RepresentativePointList[1].y ≈ y1 atol = 1e-14
        @test DT.RepresentativePointList[2].x ≈ x2 atol = 1e-14
        @test DT.RepresentativePointList[2].y ≈ y2 atol = 1e-14
        @test DT.RepresentativePointList[3].x ≈ x3 atol = 1e-14
        @test DT.RepresentativePointList[3].y ≈ y3 atol = 1e-14
        @test DT.RepresentativePointList[4].x ≈ x4 atol = 1e-14
        @test DT.RepresentativePointList[4].y ≈ y4 atol = 1e-14
        @test DT.RepresentativePointList[5].x ≈ x5 atol = 1e-14
        @test DT.RepresentativePointList[5].y ≈ y5 atol = 1e-14
    end
end