using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StatsBase

Base.include(@__MODULE__, "../helper_functions.jl")

@testset "Initialise" begin
    c = DT.RepresentativeCoordinates{Int,Float64}()
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
    c = DT.RepresentativeCoordinates{Int,Float64}()
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

@testset "Representative points from a triangulation" begin
    x = [0.0, 1.0, 1.0, 0.0, 0.0]
    y = [0.0, 0.0, 1.0, 1.0, 0.0]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes)
    DT.compute_representative_points!(tri)
    rep = DT.get_representative_point_list(tri)
    @test rep[1].x ≈ 0.5 atol = 1e-12
    @test rep[1].y ≈ 0.5 atol = 1e-12
    @test rep[1].n == 0

    x = [[0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]]
    y = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes)
    rep = DT.get_representative_point_list(tri)
    DT.compute_representative_points!(tri)
    @test rep[1].x ≈ 0.5 atol = 1e-12
    @test rep[1].y ≈ 0.5 atol = 1e-12
    @test rep[1].n == 0
end

@testset "Polylabel, representative points, and distance" begin
    x, y, x1, x2, x3, x4, x5, y1, y2, y3, y4, y5 = complicated_geometry()
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes, delete_ghosts=false)
    for i in (1, 2, 3)
        local x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
        if i == 3
            DT.compute_representative_points!(tri)
            rep = DT.get_representative_point_list(tri)
        elseif i == 2
            empty!(tri.representative_point_list)
            DT.compute_representative_points!(tri)
            rep = DT.get_representative_point_list(tri)
        else
            rep = DT.get_representative_point_list(tri)
        end
        x1, y1 = DT.pole_of_inaccessibility(tri.points, tri.boundary_nodes)
        x2, y2 = DT.pole_of_inaccessibility(tri.points, tri.boundary_nodes[2])
        x3, y3 = DT.pole_of_inaccessibility(tri.points, tri.boundary_nodes[3])
        x4, y4 = DT.pole_of_inaccessibility(tri.points, tri.boundary_nodes[4])
        x5, y5 = DT.pole_of_inaccessibility(tri.points, tri.boundary_nodes[5])
        d1 = DT.distance_to_polygon((x1, y1), tri.points, tri.boundary_nodes)
        d2 = DT.distance_to_polygon((x2, y2), tri.points, tri.boundary_nodes[2])
        d3 = DT.distance_to_polygon((x3, y3), tri.points, tri.boundary_nodes[3])
        d4 = DT.distance_to_polygon((x4, y4), tri.points, tri.boundary_nodes[4])
        d5 = DT.distance_to_polygon((x5, y5), tri.points, tri.boundary_nodes[5])
        @test all(≥(0), (d1, d2, d3, d4, d5)) # all points are inside their respective regions
        @test rep[1].x ≈ x1 atol = 1e-14
        @test rep[1].y ≈ y1 atol = 1e-14
        @test rep[2].x ≈ x2 atol = 1e-14
        @test rep[2].y ≈ y2 atol = 1e-14
        @test rep[3].x ≈ x3 atol = 1e-14
        @test rep[3].y ≈ y3 atol = 1e-14
        @test rep[4].x ≈ x4 atol = 1e-14
        @test rep[4].y ≈ y4 atol = 1e-14
        @test rep[5].x ≈ x5 atol = 1e-14
        @test rep[5].y ≈ y5 atol = 1e-14
        @test collect(get_point(tri, -1)) ≈ [x1, y1]
        @test collect(get_point(tri, -2)) ≈ [x1, y1]
        @test collect(get_point(tri, -3)) ≈ [x1, y1]
        @test collect(get_point(tri, -4)) ≈ [x1, y1]
        @test collect(get_point(tri, -5)) ≈ [x2, y2]
        @test collect(get_point(tri, -6)) ≈ [x3, y3]
        @test collect(get_point(tri, -7)) ≈ [x4, y4]
        @test collect(get_point(tri, -8)) ≈ [x4, y4]
        @test collect(get_point(tri, -9)) ≈ [x4, y4]
        @test collect(get_point(tri, -10)) ≈ [x4, y4]
        @test collect(get_point(tri, -11)) ≈ [x5, y5]
    end
end
