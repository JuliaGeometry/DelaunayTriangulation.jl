using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArraysCore
using StatsBase

global p1 = [1.3, 2.5]
global p2 = (1.3, 2.5)
global p3 = SVector{2,Float32}((1.3, 2.5))

@testset "Individual points" begin
    @testset "Getting coordinates" begin
        @test_throws "The" getx(String)
        for p in (p1, p2, p3)
            @test getx(p) == 1.3 || getx(p) == 1.3f0
            @inferred getx(p)
        end

        @test_throws "The" gety(String)
        for p in (p1, p2, p3)
            @test gety(p) == 2.5 || gety(p) == 2.5f0
            @inferred gety(p)
        end

        for p in (p1, p2, p3)
            @test getxy(p) == (p[1], p[2])
            @inferred getxy(p)
        end
    end
end

global pts1 = [[2.0, 3.5], [1.7, 23.3], [-1.0, 0.0]]
global pts2 = [(2.0, 3.5), (1.7, 23.3), (-1.0, 0.0)]
global pts3 = [2.0 1.7 -1.0; 3.5 23.3 0.0]

@testset "Collection of points" begin
    @testset "Getting points" begin
        @test_throws "The" DT.getpoint(String, 5)
        for pts in (pts1, pts2, pts3)
            @test DT.getpoint(pts, 1) == (2.0, 3.5)
            @test DT.getpoint(pts, 2) == (1.7, 23.3)
            @test DT.getpoint(pts, 3) == (-1.0, 0.0)
            @test get_point(pts, 1) == (2.0, 3.5)
            @test get_point(pts, 2, 3) == ((1.7, 23.3), (-1.0, 0.0))
            @inferred get_point(pts, 2, 3)
            @inferred get_point(pts, 1)
            @test DT.getpoint(pts, DT.getpoint(pts, 1)) == (2.0, 3.5)
            @test DT.getpoint(pts, DT.getpoint(pts, 2)) == (1.7, 23.3)
            @test DT.get_point(pts, DT.getpoint(pts, 3)) == (-1.0, 0.0)
            @inferred DT.get_point(pts, DT.getpoint(pts, 3))
            @inferred DT.getpoint(pts, DT.getpoint(pts, 3)) == (-1.0, 0.0)
            @test DT.getpoint(pts, (17.0, 23.0)) == (17.0, 23.0)
            @test DT.get_point(pts, (2.0, 5.3), (17.0, 5.3)) == ((2.0, 5.3), (17.0, 5.3))
            @inferred DT.get_point(pts, (2.0, 5.3), (17.0, 5.3)) == ((2.0, 5.3), (17.0, 5.3))
            @test DT.get_point(pts, 1, 2, (17.0, -2.0), (57.0, 23.0)) ==
                  ((2.0, 3.5), (1.7, 23.3), (17.0, -2.0), (57.0, 23.0))
            @inferred DT.get_point(pts, 1, 2, (17.0, -2.0), (57.0, 23.0)) ==
                      ((2.0, 3.5), (1.7, 23.3), (17.0, -2.0), (57.0, 23.0))
        end
    end

    @testset "Each point index" begin
        @test_throws "The" each_point_index(String)
        for pts in (pts1, pts2, pts3)
            @test each_point_index(pts) == 1:3
            @inferred each_point_index(pts)
        end
    end

    @testset "Each point" begin
        @test_throws "The" each_point(String)
        @test each_point(pts1) == pts1
        @test each_point(pts2) == pts2
        @test each_point(pts3) == eachcol(pts3)
        @inferred each_point(pts3)
    end

    @testset "Number of points" begin
        @test_throws "The" num_points(String)
        for pts in (pts1, pts2, pts3)
            @test num_points(pts) == 3
            @inferred num_points(pts)
        end
    end

    @testset "Getting points corresponding to ghost vertices" begin
        bn1 = [[[1, 2], [3, 4], [5, 6], [10, 12]],
            [[13, 25, 50], [75, 17, 5, 10]],
            [[17, 293, 101], [29, 23]]]
        bn2 = [[13, 25, 50], [75, 17, 5, 10]]
        bn3 = [17, 293, 101, 29, 23]
        map1 = DT.construct_boundary_map(bn1)
        map2 = DT.construct_boundary_map(bn2)
        map3 = DT.construct_boundary_map(bn3)
        rep = DT.get_empty_representative_points()
        rep[1] = DT.RepresentativeCoordinates(0.5, 0.3, 2)
        rep[2] = DT.RepresentativeCoordinates(2.5, 7.3, 7)
        rep[3] = DT.RepresentativeCoordinates(6.5, -0.6, 13)
        pts = rand(2, 500)
        @test DT.get_point(pts, rep, map1, 2) == Tuple(pts[:, 2])
        @inferred DT.get_point(pts, rep, map1, 2)
        @test DT.get_point(pts, rep, map1, 3, 4, 5, -1, -2, -7) ==
              (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]), (0.5, 0.3), (0.5, 0.3),
            (6.5, -0.6))
        @test DT.get_point(pts, rep, map1, -1) == (0.5, 0.3)
        @inferred DT.get_point(pts, rep, map1, -1)
        @test DT.get_point(pts, rep, map1, -2) == (0.5, 0.3)
        @test DT.get_point(pts, rep, map1, -5) == (2.5, 7.3)
        @test DT.get_point(pts, rep, map1, -6) == (2.5, 7.3)
        @test DT.get_point(pts, rep, map1, -7) == (6.5, -0.6)
        @test DT.get_point(pts, rep, map1, -8) == (6.5, -0.6)
        @test DT.get_point(pts, rep, map2, 2) == Tuple(pts[:, 2])
        @test DT.get_point(pts, rep, map2, 3, 4, 5) ==
              (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]))
        @test DT.get_point(pts, rep, map2, -1) == (0.5, 0.3)
        @test DT.get_point(pts, rep, map2, -2) == (0.5, 0.3)
        @test_throws KeyError DT.get_point(pts, rep, map2, -3)
        @test_throws KeyError DT.get_point(pts, rep, map2, -6)
        @test_throws KeyError DT.get_point(pts, rep, map2, -7)
        @test_throws KeyError DT.get_point(pts, rep, map2, -8)
        @test DT.get_point(pts, rep, map3, 2) == Tuple(pts[:, 2])
        @test DT.get_point(pts, rep, map3, 3, 4, 5) ==
              (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]))
        @test DT.get_point(pts, rep, map3, -1) == (0.5, 0.3)
        @test_throws KeyError DT.get_point(pts, rep, map3, -2)
        @test_throws KeyError DT.get_point(pts, rep, map3, -3)
        @test_throws KeyError DT.get_point(pts, rep, map3, -6)
        @test_throws KeyError DT.get_point(pts, rep, map3, -7)
        @test_throws KeyError DT.get_point(pts, rep, map3, -8)
        @test DT.get_point(pts, rep, map3, (2.0, 5.0)) == (2.0, 5.0)
        @test DT.get_point(pts, rep, map2, (2.0, 5.0), (7.0, 13.3), 5) ==
              ((2.0, 5.0), (7.0, 13.3), (pts[1, 5], pts[2, 5]))
    end

    @testset "Testing if a set of points contains no duplicates" begin
        pts = [(2.0, 3.0), (2.5, 3.7), (2.0, 17.0)]
        @test DT.points_are_unique(pts)
        pts = [2.0 2.5 2.0; 3.0 3.7 17.0]
        @test DT.points_are_unique(pts)
        pts = [(2.0, 3.0), (17.3, -5.0), (2.0, 3.0)]
        @test !DT.points_are_unique(pts)
        pts = [2.0 17.3 2.0; 3.0 -5.0 3.0]
        @test !DT.points_are_unique(pts)
        @inferred DT.points_are_unique(pts)
    end

    @testset "Sorting points lexicographically" begin
        A = [2.0 5.0 1.0 10.0 17.0 23.0 5.0 5.0
            7.0 2.0 0.0 7.5 2.0 -2.5 3.5 3.0]
        idx = DT.lexicographic_order(A)
        @test idx == [3, 1, 2, 8, 7, 4, 5, 6]
        A = [(2.0, 7.0),
            (5.0, 2.0),
            (1.0, 0.0),
            (10.0, 7.5),
            (17.0, 2.0),
            (23.0, -2.5),
            (5.0, 3.5),
            (5.0, 3.0)]
        idx = DT.lexicographic_order(A)
        @test idx == [3, 1, 2, 8, 7, 4, 5, 6]
    end

    @testset "push_point!" begin
        p1 = [(1.0, 2.0), (3.0, 4.0)]
        DT.push_point!(p1, 0.5, 3.0)
        @test p1 == [(1.0, 2.0), (3.0, 4.0), (0.5, 3.0)]
        p1 = [[1.0, 2.0], [5.0, 13.7]]
        DT.push_point!(p1, 13.9, 25.0)
        @test p1 == [[1.0, 2.0], [5.0, 13.7], [13.9, 25.0]]
        p1 = [(1.0, 5.0), (17.0, 3.0)]
        DT.push_point!(p1, (17.9, 15.0))
        @test p1 == [(1.0, 5.0), (17.0, 3.0), (17.9, 15.0)]
    end

    @testset "pop_point!" begin
        p1 = [(1.0, 2.0), (3.0, 4.0), (10.0, 10.0)]
        DT.pop_point!(p1)
        @test p1 == [(1.0, 2.0), (3.0, 4.0)]
    end

    @testset "mean_points" begin
        for _ in 1:500
            points = [(rand(), rand()) for _ in 1:500]
            px = first.(points)
            py = last.(points)
            cx = mean(px)
            cy = mean(py)
            _cx, _cy = DT.mean_points(points)
            @test cx ≈ _cx && cy ≈ _cy
            verts = rand(1:500, 50)
            unique!(verts)
            cx = mean(px[verts])
            cy = mean(py[verts])
            _cx, _cy = DT.mean_points(points, verts)
            @test cx ≈ _cx && cy ≈ _cy
        end
    end

    @testset "set_point!" begin
        points = [(1.0, 2.0), (5.0, 9.0)]
        DT.set_point!(points, 2, 3.0, 4.0)
        @test points == [(1.0, 2.0), (3.0, 4.0)]
        DT.set_point!(points, 1, (6.0, 7.7))
        @test points == [(6.0, 7.7), (3.0, 4.0)]
        points = [1.0 10.0; 7.8 5.5]
        DT.set_point!(points, 2, 3.0, 4.0)
        @test points == [1.0 3.0; 7.8 4.0]
        DT.set_point!(points, 1, (6.0, 7.7))
        @test points == [6.0 3.0; 7.7 4.0]
    end
end