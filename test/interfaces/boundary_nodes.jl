using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using DataStructures
using StaticArraysCore

include("../helper_functions.jl")

global bn1 = [[[1, 2], [3, 4], [5, 6], [10, 12]],
       [[13, 25, 50], [75, 17, 5, 10]],
       [[17, 293, 101], [29, 23]]]
global bn2 = [[13, 25, 50], [75, 17, 5, 10]]
global bn3 = [17, 293, 101, 29, 23]
global map1 = DT.construct_boundary_map(bn1)
global map2 = DT.construct_boundary_map(bn2)
global map3 = DT.construct_boundary_map(bn3)
global idx = DT.BoundaryIndex

@testset "Testing number of segments/curves" begin
       @test_throws "The" DT.has_multiple_curves(String)
       @test DT.has_multiple_curves(bn1)
       @test !DT.has_multiple_curves(bn2)
       @test !DT.has_multiple_curves(bn3)
       @test_throws "The" DT.has_multiple_segments(String)
       @test DT.has_multiple_segments(bn1)
       @test DT.has_multiple_segments(bn2)
       @test !DT.has_multiple_segments(bn3)
end

@testset "Getting number of segments/curves" begin
       @test_throws "The" DT.num_curves(String)
       @test DT.num_curves(bn1) == 3
       @test_throws "The" DT.num_segments(String)
       @test DT.num_segments(bn2) == 2
end

@testset "Number of boundary edges" begin
       @test_throws "The" DT.num_boundary_edges(String)
       @test DT.num_boundary_edges(bn3) == 4
       @test DT.num_boundary_edges(bn1[1][1]) == 1
       @test DT.num_boundary_edges(bn2[2]) == 3
       @test DT.num_boundary_edges(Int64[]) == 0
end

@testset "Getting boundary nodes" begin
       @test_throws "The" DT.getboundarynodes(String, [1, 2])
       @test get_boundary_nodes(bn1, 1) == bn1[1]
       @test get_boundary_nodes(bn1, 2) == bn1[2]
       @test get_boundary_nodes(bn2, 2) == bn2[2]
       @test get_boundary_nodes(bn3, 4) == bn3[4]
       @test get_boundary_nodes(bn1, (1, 2)) == bn1[1][2]
       @test get_boundary_nodes(bn3, bn3) == bn3
end

@testset "Getting each boundary node" begin
       @test_throws "The" DT.each_boundary_node(String)
       @test DT.each_boundary_node(bn3) == bn3
end

@testset "Constructing the boundary map" begin
       map1 = DT.construct_boundary_map(bn1)
       map2 = DT.construct_boundary_map(bn2)
       map3 = DT.construct_boundary_map(bn3)
       idx = DT.BoundaryIndex
       @test map1 ==
             OrderedDict(idx => (1, 1), idx - 1 => (1, 2), idx - 2 => (1, 3), idx - 3 => (1, 4),
              idx - 4 => (2, 1), idx - 5 => (2, 2),
              idx - 6 => (3, 1), idx - 7 => (3, 2))
       @test map2 == OrderedDict(idx => 1, idx - 1 => 2)
       @test map3 == OrderedDict(idx => bn3)
end

@testset "Mapping a boundary index" begin
       @test DT.map_boundary_index(map1, idx - 4) == (2, 1)
       @test DT.map_boundary_index(map2, idx - 1) == 2
       @test DT.map_boundary_index(map3, idx) == bn3
end

@testset "Getting a curve index" begin
       @test DT.get_curve_index(map1, idx - 4) == 2
       @test DT.get_curve_index(map2, idx - 1) == 1
       @test DT.get_curve_index(3) == 1
       @test DT.get_curve_index((5, 7)) == 5
       @test DT.get_curve_index(map3, idx) == 1
end

@testset "Getting a segment index" begin
       @test DT.get_segment_index(map1, idx - 4) == 1
       @test DT.get_segment_index(map2, idx - 1) == 2
       @test DT.get_segment_index(3) == 3
       @test DT.get_segment_index((5, 7)) == 7
       @test DT.get_segment_index(map3, idx) == 1
end

@testset "Number of boundary segments" begin
       @test DT.num_outer_boundary_segments(bn1) == 4
       @test DT.num_outer_boundary_segments(bn2) == 2
       @test DT.num_outer_boundary_segments(bn3) == 1
end

@testset "Testing if a boundary has multiple segments from a boundary map" begin
       @test DT.has_multiple_segments(map1)
       @test DT.has_multiple_segments(map2)
       @test !DT.has_multiple_segments(map3)
end

@testset "Getting boundary index ranges" begin
       d1 = DT.construct_boundary_index_ranges(bn1)
       d2 = DT.construct_boundary_index_ranges(bn2)
       d3 = DT.construct_boundary_index_ranges(bn3)
       boundary_nodes = [[[1, 2, 3, 4], [4, 5, 6, 1]],
              [[18, 19, 20, 25, 26, 30]],
              [[50, 51, 52, 53, 54, 55], [55, 56, 57, 58],
                     [58, 101, 103, 105, 107, 120],
                     [120, 121, 122, 50]]]
       d4 = DT.construct_boundary_index_ranges(boundary_nodes)
       @test d4 == OrderedDict(-1 => -2:-1,
              -2 => -2:-1,
              -3 => -3:-3,
              -4 => -7:-4,
              -5 => -7:-4,
              -6 => -7:-4,
              -7 => -7:-4)
       @test d1 == OrderedDict(-1 => -4:-1,
              -2 => -4:-1,
              -3 => -4:-1,
              -4 => -4:-1,
              -5 => -6:-5,
              -6 => -6:-5,
              -7 => -8:-7,
              -8 => -8:-7)
       @test d2 == OrderedDict(-1 => -2:-1, -2 => -2:-1)
       @test d3 == OrderedDict(-1 => -1:-1)

       x, y = complicated_geometry()
       tri = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)
       @test tri.boundary_index_ranges == OrderedDict(-1 => -4:-1,
              -2 => -4:-1,
              -3 => -4:-1,
              -4 => -4:-1,
              -5 => -5:-5,
              -6 => -6:-6,
              -7 => -10:-7,
              -8 => -10:-7,
              -9 => -10:-7,
              -10 => -10:-7,
              -11 => -11:-11)
end
