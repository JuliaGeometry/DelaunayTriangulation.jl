using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using StaticArraysCore

p1 = [1.3, 2.5]
p2 = (1.3, 2.5)
p3 = SVector{2,Float32}((1.3, 2.5))

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

pts1 = [[2.0, 3.5], [1.7, 23.3], [-1.0, 0.0]]
pts2 = [(2.0, 3.5), (1.7, 23.3), (-1.0, 0.0)]
pts3 = [2.0 1.7 -1.0; 3.5 23.3 0.0]

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

@test_throws "The" each_point_index(String)
for pts in (pts1, pts2, pts3)
    @test each_point_index(pts) == 1:3
    @inferred each_point_index(pts)
end

@test_throws "The" each_point(String)
@test each_point(pts1) == pts1
@test each_point(pts2) == pts2
@test each_point(pts3) == eachcol(pts3)
@inferred each_point(pts3)

@test_throws "The" num_points(String)
for pts in (pts1, pts2, pts3)
    @test num_points(pts) == 3
    @inferred num_points(pts)
end

bn1 = [[[1, 2], [3, 4], [5, 6], [10, 12]],
       [[13, 25, 50], [75, 17, 5, 10]],
       [[17, 293, 101], [29, 23]]]
bn2 = [[13, 25, 50], [75, 17, 5, 10]]
bn3 = [17, 293, 101, 29, 23]
map1 = DT.construct_boundary_map(bn1)
map2 = DT.construct_boundary_map(bn2)
map3 = DT.construct_boundary_map(bn3)
DT.RepresentativePointList[1] = DT.RepresentativeCoordinates(0.5, 0.3, 2)
DT.RepresentativePointList[2] = DT.RepresentativeCoordinates(2.5, 7.3, 7)
DT.RepresentativePointList[3] = DT.RepresentativeCoordinates(6.5, -0.6, 13)
pts = rand(2, 500)
@test DT.get_point(pts, map1, 2) == Tuple(pts[:, 2])
@inferred DT.get_point(pts, map1, 2)
@test DT.get_point(pts, map1, 3, 4, 5, -1, -2, -7) ==
      (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]), (0.5, 0.3), (0.5, 0.3),
       (6.5, -0.6))
@test DT.get_point(pts, map1, -1) == (0.5, 0.3)
@inferred DT.get_point(pts, map1, -1)
@test DT.get_point(pts, map1, -2) == (0.5, 0.3)
@test DT.get_point(pts, map1, -5) == (2.5, 7.3)
@test DT.get_point(pts, map1, -6) == (2.5, 7.3)
@test DT.get_point(pts, map1, -7) == (6.5, -0.6)
@test DT.get_point(pts, map1, -8) == (6.5, -0.6)
@test DT.get_point(pts, map2, 2) == Tuple(pts[:, 2])
@test DT.get_point(pts, map2, 3, 4, 5) ==
      (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]))
@test DT.get_point(pts, map2, -1) == (0.5, 0.3)
@test DT.get_point(pts, map2, -2) == (0.5, 0.3)
@test_throws KeyError DT.get_point(pts, map2, -3)
@test_throws KeyError DT.get_point(pts, map2, -6)
@test_throws KeyError DT.get_point(pts, map2, -7)
@test_throws KeyError DT.get_point(pts, map2, -8)
@test DT.get_point(pts, map3, 2) == Tuple(pts[:, 2])
@test DT.get_point(pts, map3, 3, 4, 5) ==
      (Tuple(pts[:, 3]), Tuple(pts[:, 4]), Tuple(pts[:, 5]))
@test DT.get_point(pts, map3, -1) == (0.5, 0.3)
@test_throws KeyError DT.get_point(pts, map3, -2)
@test_throws KeyError DT.get_point(pts, map3, -3)
@test_throws KeyError DT.get_point(pts, map3, -6)
@test_throws KeyError DT.get_point(pts, map3, -7)
@test_throws KeyError DT.get_point(pts, map3, -8)
@test DT.get_point(pts, map3, (2.0, 5.0)) == (2.0, 5.0)
@test DT.get_point(pts, map2, (2.0, 5.0), (7.0, 13.3), 5) ==
      ((2.0, 5.0), (7.0, 13.3), (pts[1, 5], pts[2, 5]))

pts = [(2.0, 3.0), (2.5, 3.7), (2.0, 17.0)]
@test DT.points_are_unique(pts)
pts = [2.0 2.5 2.0; 3.0 3.7 17.0]
@test DT.points_are_unique(pts)
pts = [(2.0, 3.0), (17.3, -5.0), (2.0, 3.0)]
@test !DT.points_are_unique(pts)
pts = [2.0 17.3 2.0; 3.0 -5.0 3.0]
@test !DT.points_are_unique(pts)
@inferred DT.points_are_unique(pts)

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