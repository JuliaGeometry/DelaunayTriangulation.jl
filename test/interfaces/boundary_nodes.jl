using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Test
using DataStructures
using StaticArraysCore

Base.include(@__MODULE__, "../helper_functions.jl")

global bn1 = [[[1, 2], [3, 4], [5, 6], [10, 12]],
       [[13, 25, 50], [75, 17, 5, 10]],
       [[17, 293, 101], [29, 23]]]
global bn2 = [[13, 25, 50], [75, 17, 5, 10]]
global bn3 = [17, 293, 101, 29, 23]
global map1 = DT.construct_ghost_vertex_map(bn1)
global map2 = DT.construct_ghost_vertex_map(bn2)
global map3 = DT.construct_ghost_vertex_map(bn3)
global idx = DT.𝒢

@testset "Testing number of sections/curves" begin
       @test DT.has_multiple_curves(bn1)
       @test !DT.has_multiple_curves(bn2)
       @test !DT.has_multiple_curves(bn3)
       @test !DT.has_multiple_curves(Int[])
       @test !DT.has_multiple_sections(Int[])
       @test DT.has_multiple_sections(bn1)
       @test DT.has_multiple_sections(bn2)
       @test !DT.has_multiple_sections(bn3)
end

@testset "Getting number of sections/curves" begin
       @test DT.num_curves(bn1) == 3
       @test DT.num_curves(bn2) == 1
       @test DT.num_curves(bn3) == 1
       @test DT.num_sections(bn2) == 2
end

@testset "Number of boundary edges" begin
       @test DT.num_boundary_edges(bn3) == 4
       @test DT.num_boundary_edges(bn1[1][1]) == 1
       @test DT.num_boundary_edges(bn2[2]) == 3
       @test DT.num_boundary_edges(Int[]) == 0
end

@testset "Getting boundary nodes" begin
       @test get_boundary_nodes(bn1, 1) == bn1[1]
       @test get_boundary_nodes(bn1, 2) == bn1[2]
       @test get_boundary_nodes(bn2, 2) == bn2[2]
       @test get_boundary_nodes(bn3, 4) == bn3[4]
       @test get_boundary_nodes(bn1, (1, 2)) == bn1[1][2]
       @test get_boundary_nodes(bn3, bn3) == bn3
end

@testset "Getting each boundary node" begin
       @test DT.each_boundary_node(bn3) == bn3
end

@testset "Constructing the ghost vertex map" begin
       map1 = DT.construct_ghost_vertex_map(bn1)
       map2 = DT.construct_ghost_vertex_map(bn2)
       map3 = DT.construct_ghost_vertex_map(bn3)
       idx = DT.𝒢
       @test map1 ==
             Dict(idx => (1, 1), idx - 1 => (1, 2), idx - 2 => (1, 3), idx - 3 => (1, 4),
              idx - 4 => (2, 1), idx - 5 => (2, 2),
              idx - 6 => (3, 1), idx - 7 => (3, 2))
       @test map2 == Dict(idx => 1, idx - 1 => 2)
       @test map3 == Dict(idx => bn3)
end

@testset "Getting a curve index" begin
       @test DT.get_curve_index(map1, idx - 4) == 2
       @test DT.get_curve_index(map2, idx - 1) == 1
       @test DT.get_curve_index(3) == 1
       @test DT.get_curve_index((5, 7)) == 5
       @test DT.get_curve_index(map3, idx) == 1
end

@testset "Getting a section index" begin
       @test DT.get_section_index(map1, idx - 4) == 1
       @test DT.get_section_index(map2, idx - 1) == 2
       @test DT.get_section_index(3) == 3
       @test DT.get_section_index((5, 7)) == 7
       @test DT.get_section_index(map3, idx) == 1
end

@testset "Getting ghost vertex ranges" begin
       for r in 1:500 # this example used to StackOverflow randomly, so let's just be sure it doesn't come back
              d1 = DT.construct_ghost_vertex_ranges(bn1)
              d2 = DT.construct_ghost_vertex_ranges(bn2)
              d3 = DT.construct_ghost_vertex_ranges(bn3)
              boundary_nodes = [[[1, 2, 3, 4], [4, 5, 6, 1]],
                     [[18, 19, 20, 25, 26, 30]],
                     [[50, 51, 52, 53, 54, 55], [55, 56, 57, 58],
                            [58, 101, 103, 105, 107, 120],
                            [120, 121, 122, 50]]]
              d4 = DT.construct_ghost_vertex_ranges(boundary_nodes)
              @test d4 == Dict(-1 => -2:-1,
                     -2 => -2:-1,
                     -3 => -3:-3,
                     -4 => -7:-4,
                     -5 => -7:-4,
                     -6 => -7:-4,
                     -7 => -7:-4)
              @test d1 == Dict(-1 => -4:-1,
                     -2 => -4:-1,
                     -3 => -4:-1,
                     -4 => -4:-1,
                     -5 => -6:-5,
                     -6 => -6:-5,
                     -7 => -8:-7,
                     -8 => -8:-7)
              @test d2 == Dict(-1 => -2:-1, -2 => -2:-1)
              @test d3 == Dict(-1 => -1:-1)

              x, y = complicated_geometry()
              boundary_nodes, points = convert_boundary_points_to_indices(x, y)
              tri = triangulate(points; boundary_nodes)
              @test tri.ghost_vertex_ranges == Dict(-1 => -4:-1,
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
end

@testset "construct_boundary_edge_map" begin
       bn = [1, 2, 3, 4, 5, 6, 7, 1]
       bn_map = DT.construct_boundary_edge_map(bn)
       for (ij, (index, k)) in bn_map
              S = get_boundary_nodes(bn, index)
              @test get_boundary_nodes(S, k) == ij[1]
              @test get_boundary_nodes(S, k + 1) == ij[2]
       end
       bn = [[1, 2, 3, 4], [4, 5, 6, 7, 8], [8, 9, 10, 1]]
       bn_map = DT.construct_boundary_edge_map(bn)
       for (ij, (index, k)) in bn_map
              S = get_boundary_nodes(bn, index)
              @test get_boundary_nodes(S, k) == ij[1]
              @test get_boundary_nodes(S, k + 1) == ij[2]
       end
       bn = [
              [[1, 2, 3, 4, 5], [5, 6, 7], [7, 8], [8, 9, 10, 1]],
              [[13, 14, 15, 16, 17], [17, 18, 19, 20], [20, 13]]
       ]
       bn_map = DT.construct_boundary_edge_map(bn)
       for (ij, (index, k)) in bn_map
              S = get_boundary_nodes(bn, index)
              @test get_boundary_nodes(S, k) == ij[1]
              @test get_boundary_nodes(S, k + 1) == ij[2]
       end
       bn = Int[]
       bn_map = DT.construct_boundary_edge_map(bn)
       @test bn_map == Dict{Tuple{Int32,Int32},Tuple{Vector{Int},Int}}()
end

@testset "insert_boundary_node!" begin
       bn = [1, 2, 3, 4, 5, 6, 7, 1]
       DT.insert_boundary_node!(bn, (bn, 5), 17)
       DT.insert_boundary_node!(bn, (bn, 1), 13)
       @test bn == [13, 1, 2, 3, 4, 17, 5, 6, 7, 1]
       bn = [[1, 2, 3, 4], [4, 5, 6, 7, 8], [8, 9, 10, 1]]
       DT.insert_boundary_node!(bn, (1, 2), 9)
       DT.insert_boundary_node!(bn, (1, 4), 18)
       DT.insert_boundary_node!(bn, (2, 4), 23)
       DT.insert_boundary_node!(bn, (3, 1), 5)
       @test bn == [[1, 9, 2, 18, 3, 4], [4, 5, 6, 23, 7, 8], [5, 8, 9, 10, 1]]
       bn = [
              [[1, 2, 3, 4, 5], [5, 6, 7], [7, 8], [8, 9, 10, 1]],
              [[13, 14, 15, 16, 17], [17, 18, 19, 20], [20, 13]]
       ]
       DT.insert_boundary_node!(bn, ((1, 1), 1), 17)
       DT.insert_boundary_node!(bn, ((1, 2), 3), 38)
       DT.insert_boundary_node!(bn, ((1, 3), 2), 50)
       DT.insert_boundary_node!(bn, ((1, 4), 3), 67)
       DT.insert_boundary_node!(bn, ((2, 1), 3), 500)
       DT.insert_boundary_node!(bn, ((2, 2), 3), 87)
       DT.insert_boundary_node!(bn, ((2, 3), 2), 671)
       @test bn == [
              [[17, 1, 2, 3, 4, 5], [5, 6, 38, 7], [7, 50, 8], [8, 9, 67, 10, 1]],
              [[13, 14, 500, 15, 16, 17], [17, 18, 87, 19, 20], [20, 671, 13]]
       ]
end

@testset "delete_boundary_node!" begin
       bn = [1, 2, 3, 4, 5, 6, 7, 1]
       DT.delete_boundary_node!(bn, (bn, 5))
       DT.delete_boundary_node!(bn, (bn, 1))
       @test bn == [2, 3, 4, 6, 7, 1]
       bn = [[1, 2, 3, 4], [4, 5, 6, 7, 8], [8, 9, 10, 1]]
       DT.delete_boundary_node!(bn, (1, 2))
       DT.delete_boundary_node!(bn, (1, 2))
       DT.delete_boundary_node!(bn, (2, 4))
       DT.delete_boundary_node!(bn, (3, 1))
       @test bn == [[1, 4], [4, 5, 6, 8], [9, 10, 1]]
       bn = [
              [[1, 2, 3, 4, 5], [5, 6, 7], [7, 8], [8, 9, 10, 1]],
              [[13, 14, 15, 16, 17], [17, 18, 19, 20], [20, 13]]
       ]
       DT.delete_boundary_node!(bn, ((1, 1), 1))
       DT.delete_boundary_node!(bn, ((1, 2), 3))
       DT.delete_boundary_node!(bn, ((1, 3), 2))
       DT.delete_boundary_node!(bn, ((1, 4), 3))
       DT.delete_boundary_node!(bn, ((2, 1), 3))
       DT.delete_boundary_node!(bn, ((2, 2), 3))
       DT.delete_boundary_node!(bn, ((2, 3), 2))
       @test bn == [
              [[2, 3, 4, 5], [5, 6], [7], [8, 9, 1]],
              [[13, 14, 16, 17], [17, 18, 20], [20]]
       ]
end