using ..DelaunayTriangulation
const DT = DelaunayTriangulation
import SimpleGraphs: relabel, UndirectedGraph
using DataStructures

include("../helper_functions.jl")



@testset "Small example" begin
    tri = example_triangulation()
    DT.add_triangle!(tri, 6, 2, 3)
    DT.split_triangle!(tri, 4, 6, 1, 8)
    DT.split_triangle!(tri, 1, 5, 4, 7)
    @testset "First flip" begin
        true_T = Set{NTuple{3,Int}}([
            (5, 6, 3),
            (3, 2, 5),
            (4, 1, 7),
            (5, 4, 7),
            (5, 1, 6),
            (1, 5, 7),
            (6, 2, 3),
            (6, 1, 8),
            (4, 6, 8),
            (1, 4, 8)
        ])
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict(
                (5, 6) => 3, (6, 3) => 5, (3, 5) => 6,
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 7, (1, 7) => 4, (7, 4) => 1,
                (5, 4) => 7, (4, 7) => 5, (7, 5) => 4,
                (5, 1) => 6, (1, 6) => 5, (6, 5) => 1,
                (1, 5) => 7, (5, 7) => 1, (7, 1) => 5,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (6, 1) => 8, (1, 8) => 6, (8, 6) => 1,
                (4, 6) => 8, (6, 8) => 4, (8, 4) => 6,
                (1, 4) => 8, (4, 8) => 1, (8, 1) => 4,
                (4, 5) => DT.BoundaryIndex,
                (5, 2) => DT.BoundaryIndex,
                (2, 6) => DT.BoundaryIndex,
                (6, 4) => DT.BoundaryIndex,
            )
        )
        true_adj2v = Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(6, 5), (5, 7), (7, 4), (4, 8), (8, 6)]),
            2 => Set{NTuple{2,Int}}([(5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(5, 6), (6, 2), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(6, 8), (8, 1), (1, 7), (7, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 7), (7, 1), (1, 6), (6, 3), (3, 2)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 5), (5, 1), (1, 8), (8, 4)]),
            7 => Set{NTuple{2,Int}}([(1, 5), (5, 4), (4, 1)]),
            8 => Set{NTuple{2,Int}}([(1, 4), (4, 6), (6, 1)])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 0 1 1 1 0 0
                    0 0 0 0 1 1 1 1 1
                    1 0 0 1 0 1 1 0 0
                    0 0 1 0 0 1 1 0 0
                    1 1 0 0 0 1 1 1 1
                    1 1 1 1 1 0 1 1 0
                    1 1 1 1 1 1 0 0 1
                    0 1 0 0 1 1 0 0 0
                    0 1 0 0 1 0 1 0 0
                ]
            ), Dict(1:9 .=> [-1, 1, 2, 3, 4, 5, 6, 7, 8]))
        DT.flip_edge!(tri, 1, 3)
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test (get_graph ∘ get_graph)(tri) == true_DG
    end

    @testset "Second flip" begin
        true_T = Set{NTuple{3,Int}}([
            (5, 6, 3),
            (3, 2, 5),
            (5, 4, 7),
            (5, 1, 6),
            (1, 5, 7),
            (6, 2, 3),
            (6, 1, 8),
            (4, 6, 8),
            (8, 1, 7),
            (8, 7, 4)
        ])
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict(
                (5, 6) => 3, (6, 3) => 5, (3, 5) => 6,
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (5, 4) => 7, (4, 7) => 5, (7, 5) => 4,
                (5, 1) => 6, (1, 6) => 5, (6, 5) => 1,
                (1, 5) => 7, (5, 7) => 1, (7, 1) => 5,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (6, 1) => 8, (1, 8) => 6, (8, 6) => 1,
                (4, 6) => 8, (6, 8) => 4, (8, 4) => 6,
                (8, 1) => 7, (1, 7) => 8, (7, 8) => 1,
                (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
                (4, 5) => DT.BoundaryIndex,
                (5, 2) => DT.BoundaryIndex,
                (2, 6) => DT.BoundaryIndex,
                (6, 4) => DT.BoundaryIndex,
            )
        )
        true_adj2v = Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(6, 5), (5, 7), (7, 8), (8, 6)]),
            2 => Set{NTuple{2,Int}}([(5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(5, 6), (6, 2), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(6, 8), (8, 7), (7, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 7), (7, 1), (1, 6), (6, 3), (3, 2)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 5), (5, 1), (1, 8), (8, 4)]),
            7 => Set{NTuple{2,Int}}([(1, 5), (5, 4), (4, 8), (8, 1)]),
            8 => Set{NTuple{2,Int}}([(1, 7), (7, 4), (4, 6), (6, 1)])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 0 1 1 1 0 0
                    0 0 0 0 0 1 1 1 1
                    1 0 0 1 0 1 1 0 0
                    0 0 1 0 0 1 1 0 0
                    1 0 0 0 0 1 1 1 1
                    1 1 1 1 1 0 1 1 0
                    1 1 1 1 1 1 0 0 1
                    0 1 0 0 1 1 0 0 1
                    0 1 0 0 1 0 1 1 0
                ]
            ), Dict(1:9 .=> [-1, 1, 2, 3, 4, 5, 6, 7, 8]))
        DT.flip_edge!(tri, 1, 4)
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test (get_graph ∘ get_graph)(tri) == true_DG
    end
end

@testset "Another example" begin
    tri = example_triangulation()
    DT.add_triangle!(tri, 6, 2, 3)
    DT.split_triangle!(tri, 1, 3, 5, 7)
    i, j = 5, 1
    r = 7
    e = DT.get_adjacent(tri, j, i)
    DT.flip_edge!(tri, i, j)
    true_T = Set{NTuple{3,Int}}([
        (3, 2, 5),
        (1, 3, 7),
        (3, 5, 7),
        (6, 3, 1),
        (4, 6, 1),
        (6, 2, 3),
        (7, 5, 4),
        (7, 4, 1)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (7, 5) => 4, (5, 4) => 7, (4, 7) => 5,
            (7, 4) => 1, (4, 1) => 7, (1, 7) => 4,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int}}([(6, 3), (3, 7), (7, 4), (4, 6)]),
        2 => Set{NTuple{2,Int}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int}}([(2, 5), (5, 7), (7, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int}}([(6, 1), (1, 7), (7, 5)]),
        5 => Set{NTuple{2,Int}}([(4, 7), (7, 3), (3, 2)]),
        6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int}}([(3, 5), (5, 4), (4, 1), (1, 3)])
    )
    true_DG = relabel(UndirectedGraph([
            0 0 1 0 1 1 1 0
            0 0 0 1 1 0 1 1
            1 0 0 1 0 1 1 0
            0 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 1
            1 0 1 1 1 0 0 1
            1 1 1 1 1 0 0 0
            0 1 0 1 1 1 0 0
        ]), Dict(1:8 .=> [-1, (1:7)...]))
    DT.clear_empty_features!(tri)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
    @test (get_graph ∘ get_graph)(tri) == true_DG
    @test all(DT.is_positively_oriented(DT.triangle_orientation(tri, T...)) for T in each_triangle(tri))
end