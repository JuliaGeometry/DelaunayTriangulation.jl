using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArrays
using StatsBase
import GeometryBasics: Point2f

global p1 = [1.3, 2.5]
global p2 = (1.3, 2.5)
global p3 = SVector{2,Float32}((1.3, 2.5))

@testset "Individual points" begin
    @testset "Getting coordinates" begin
        for p in (p1, p2, p3)
            @test getx(p) == 1.3 || getx(p) == 1.3f0
            @inferred getx(p)
        end

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
global pts4 = ((2.0, 3.5), (1.7, 23.3), (-1.0, 0.0))

@testset "Collection of points" begin
    @testset "Getting points" begin
        for pts in (pts1, pts2, pts3, pts4)
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
        for pts in (pts1, pts2, pts3, pts4)
            @test DT.each_point_index(pts) == 1:3
            @inferred DT.each_point_index(pts)
        end
    end

    @testset "Each point" begin
        @test DT.each_point(pts1) == pts1
        @test DT.each_point(pts2) == pts2
        @test DT.each_point(pts3) == eachcol(pts3)
        @test DT.each_point(pts4) == pts4
        @inferred DT.each_point(pts4)
        @inferred DT.each_point(pts3)
    end

    @testset "Number of points" begin
        for pts in (pts1, pts2, pts3, pts4)
            @test DT.num_points(pts) == 3
            @inferred DT.num_points(pts)
        end
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
        A = [
            2.0 5.0 1.0 10.0 17.0 23.0 5.0 5.0
            7.0 2.0 0.0 7.5 2.0 -2.5 3.5 3.0
        ]
        idx = DT.lexicographic_order(A)
        @test idx == [3, 1, 2, 8, 7, 4, 5, 6]
        A = [
            (2.0, 7.0),
            (5.0, 2.0),
            (1.0, 0.0),
            (10.0, 7.5),
            (17.0, 2.0),
            (23.0, -2.5),
            (5.0, 3.5),
            (5.0, 3.0),
        ]
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

        points = [[1.0, 5.0], [6.5, 2.3]]
        DT.set_point!(points, 2, 3.4, 6.7)
        @test points == [[1.0, 5.0], [3.4, 6.7]]

        points = [SVector{2,Float64}(1.0, 5.4), SVector{2,Float64}(6.5, 2.3)]
        DT.set_point!(points, 2, 3.4, 6.7)
        @test points == [SVector{2,Float64}(1.0, 5.4), SVector{2,Float64}(3.4, 6.7)]

        points = [Float32[1.0, 5.0], Float32[2.3, -6.7]]
        DT.set_point!(points, 1, 2.3, 6.9)
        @test points == [Float32[2.3, 6.9], Float32[2.3, -6.7]]

        points = [Point2f(2.0, 10.0), Point2f(3.0, 4.0)]
        DT.set_point!(points, 2, 3.0, 4.0)
        @test points == [Point2f(2.0, 10.0), Point2f(3.0, 4.0)]

        points = rand(3, 50)
        _points = copy(points)
        DT.set_point!(points, 7, 3.0, 4.0)
        @test points[:, 7] == [3.0, 4.0, _points[3, 7]]
    end
end

@testset "find_point_index" begin
    points1 = [(1.0, 2.0), (5.0, 9.0), (3.0, 4.0)]
    points2 = [1.0 5.0 3.0; 2.0 9.0 4.0]
    points3 = rand(2, 75)
    points4 = [SVector{2,Float32}((1.0, 2.0)), SVector{2,Float32}((5.0, 9.0)), SVector{2,Float32}((3.0, 4.0))]
    points5 = [Point2f(1.0, 2.0), Point2f(5.0, 9.0), Point2f(3.0, 4.0)]
    points6 = [Float32[1.0, 2.0], Float32[5.0, 9.0], Float32[3.0, 4.0]]
    @test DT.find_point_index(points1, 1.0, 2.0) == 1
    @test DT.find_point_index(points1, 5.0, 9.0) == 2
    @test DT.find_point_index(points1, 3.0, 4.0) == 3
    @test DT.find_point_index(points2, 1.0, 2.0) == 1
    @test DT.find_point_index(points2, 5.0, 9.0) == 2
    @test DT.find_point_index(points2, 3.0, 4.0) == 3
    for i in 1:75
        @test DT.find_point_index(points3, points3[1, i], points3[2, i]) == i
    end
    @test DT.find_point_index(points4, 1.0, 2.0) == 1
    @test DT.find_point_index(points4, 5.0, 9.0) == 2
    @test DT.find_point_index(points4, 3.0, 4.0) == DT.find_point_index(points4, (3.0, 4.0)) == 3
    @test DT.find_point_index(points5, 1.0, 2.0) == 1
    @test DT.find_point_index(points5, 5.0, 9.0) == 2
    @test DT.find_point_index(points5, 3.0, 4.0) == 3
    @test DT.find_point_index(points6, 1.0, 2.0) == 1
    @test DT.find_point_index(points6, 5.0, 9.0) == 2
    @test DT.find_point_index(points6, 3.0, 4.0) == 3
    for _ in 1:1000
        @test DT.find_point_index(points1, rand(2)...) == DT.∅
        @test DT.find_point_index(points2, rand(2)...) == DT.∅
        @test DT.find_point_index(points3, rand(2)...) == DT.∅
        @test DT.find_point_index(points4, rand(2)...) == DT.∅
        @test DT.find_point_index(points5, rand(2)...) == DT.∅
        @test DT.find_point_index(points6, rand(2)...) == DT.∅
        @test DT.find_point_index(points1, rand(2)) == DT.∅
    end
end

@testset "find_duplicate_points" begin
    points = rand(2, 50)
    dict = DT.find_duplicate_points(points)
    @test isempty(dict)
    points = [
        (1.0, 2.0),
        (2.5, 2.5),
        (17.5, 17.5),
        (1.0, 2.0),
        (17.5, 2.0),
        (17.5, 17.5),
        (0.0, 0.0),
        (0.0, 0.0),
        (20.0, 20.0),
    ]
    dict = DT.find_duplicate_points(points)
    @test dict == Dict(
        (1.0, 2.0) => [1, 4],
        (17.5, 17.5) => [3, 6],
        (0.0, 0.0) => [7, 8],
    )
end

@testset "getz/_getz/getxyz/_getxyz" begin
    p = (1.0, 2.0, 3.0)
    q = (1.0f0, 2.0f0, 3.0f0)
    r = [1.0f0, 4.0f0, 5.0f0]
    s = @SVector [3, 7, 5]
    t = [5.0, 13.0, -5.0]

    @test DT.getz(p) === 3.0
    @test DT._getz(p) === 3.0
    @test DT.getxyz(p) === (1.0, 2.0, 3.0)
    @test DT._getxyz(p) === (1.0, 2.0, 3.0)

    @test DT.getz(q) === 3.0f0
    @test DT._getz(q) === 3.0
    @test DT.getxyz(q) === (1.0f0, 2.0f0, 3.0f0)
    @test DT._getxyz(q) === (1.0, 2.0, 3.0)

    @test DT.getz(r) === 5.0f0
    @test DT._getz(r) === 5.0
    @test DT.getxyz(r) === (1.0f0, 4.0f0, 5.0f0)
    @test DT._getxyz(r) === (1.0, 4.0, 5.0)

    @test DT.getz(s) === 5
    @test DT._getz(s) === 5.0
    @test DT.getxyz(s) === (3, 7, 5)
    @test DT._getxyz(s) === (3.0, 7.0, 5.0)

    @test DT.getz(t) === -5.0
    @test DT._getz(t) === -5.0
    @test DT.getxyz(t) === (5.0, 13.0, -5.0)
    @test DT._getxyz(t) === (5.0, 13.0, -5.0)
end

@testset "is_point2/is_point3" begin
    @test DT.is_point2((1.0, 2.0))
    @test !DT.is_point2((1.0, 2.0, 3.0))
    @test !DT.is_point2((1.0,))
    @test DT.is_point2([1.0, 2.0])
    @test !DT.is_point2([1.0, 2.0, 3.0])
    @test DT.is_point2(SVector{2,Float64}(1.0, 2.0))
    @test DT.is_point2(Point2f(1.0, 2.0))
    @test !DT.is_point2([[1.0, 1.0]])
    @test !DT.is_point2([[1.0, 2.0], [1.0, 2.0]])

    @test DT.is_point3((1.0, 2.0, 3.0))
    @test DT.is_point3([1.0, 2.0, 3.0])
    @test !DT.is_point3([[1.0, 2.0, 3.0], [3.0, 4.0, 5.0], [1.0, 2.0, 3.0]])
    @test DT.is_point3(SVector{3,Float64}(1.0, 2.0, 3.0))
end

@testset "is_planar" begin
    @test DT.is_planar(rand(2, 50))
    @test DT.is_planar([(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)])
    @test DT.is_planar([SVector{2,Float64}((1.0, 2.0)), SVector{2,Float64}((2.0, 3.0)), SVector{2,Float64}((3.0, 4.0))])
    @test !DT.is_planar([(1.0, 2.0, 3.0), (2.0, 3.0, 4.0), (3.0, 4.0, 5.0)])
    @test !DT.is_planar(rand(3, 50))
    @test !DT.is_planar([SVector{3,Float64}((1.0, 2.0, 3.0)), SVector{3,Float64}((2.0, 3.0, 4.0)), SVector{3,Float64}((3.0, 4.0, 5.0))])
end