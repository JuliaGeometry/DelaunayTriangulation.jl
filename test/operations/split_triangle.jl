using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StaticArrays


@testset "Specific example" begin
    tri = example_triangulation()
    p9 = @SVector[-1.0, 1.0]
    p10 = @SVector[4.0, 1.0]
    push!(get_points(tri), p9, p10)
    DT.add_triangle!(tri, 5, 2, 8)
    DT.add_triangle!(tri, 6, 2, 3)

    @testset "Interior splitting" begin
        true_T = Set{NTuple{3,Int}}([
            (3, 2, 5),
            (4, 1, 5),
            (5, 2, 8),
            (6, 3, 1),
            (4, 6, 1),
            (1, 3, 7), (3, 5, 7), (5, 1, 7),
            (6, 2, 3),
        ])
        true_adj = DefaultDict(DT.∅,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
                (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (4, 5) => DT.𝒢,
                (5, 8) => DT.𝒢,
                (8, 2) => DT.𝒢,
                (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢
            )
        )
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(5, 4), (4, 6), (6, 3), (3, 7), (7, 5)]),
            2 => Set{NTuple{2,Int}}([(8, 5), (5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(6, 1), (1, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 4)]),
            7 => Set{NTuple{2,Int}}([(1, 3), (3, 5), (5, 1)]),
            8 => Set{NTuple{2,Int}}([(5, 2)])
        )
        true_DG = _make_graph_from_adjacency([
                0 0 1 0 1 1 1 0 1
                0 0 0 1 1 1 1 1 0
                1 0 0 1 0 1 1 0 1
                0 1 1 0 0 1 1 1 0
                1 1 0 0 0 1 1 0 0
                1 1 1 1 1 0 0 1 1
                1 1 1 1 1 0 0 0 0
                0 1 0 1 0 1 0 0 0
                1 0 1 0 0 1 0 0 0
            ], Dict(1:9 .=> [-1, (1:8)...]))
        DT.split_triangle!(tri, 1, 3, 5, 7)
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG
    end

    @testset "Splitting a triangle with one boundary edge" begin
        true_T = Set{NTuple{3,Int}}([
            (3, 2, 5),
            (4, 1, 5),
            (5, 2, 8),
            (6, 3, 1),
            (4, 6, 9), (6, 1, 9), (1, 4, 9),
            (1, 3, 7), (3, 5, 7), (5, 1, 7),
            (6, 2, 3),
        ])
        true_adj = DefaultDict(DT.∅,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
                (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
                (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
                (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
                (4, 5) => DT.𝒢,
                (5, 8) => DT.𝒢,
                (8, 2) => DT.𝒢,
                (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢
            )
        )
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
            2 => Set{NTuple{2,Int}}([(8, 5), (5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(6, 9), (9, 1), (1, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
            7 => Set{NTuple{2,Int}}([(1, 3), (3, 5), (5, 1)]),
            8 => Set{NTuple{2,Int}}([(5, 2)]),
            9 => Set{NTuple{2,Int}}([(6, 1), (1, 4), (4, 6)])
        )
        true_DG = _make_graph_from_adjacency([
                0 0 1 0 1 1 1 0 1 0
                0 0 0 1 1 1 1 1 0 1
                1 0 0 1 0 1 1 0 1 0
                0 1 1 0 0 1 1 1 0 0
                1 1 0 0 0 1 1 0 0 1
                1 1 1 1 1 0 0 1 1 0
                1 1 1 1 1 0 0 0 0 1
                0 1 0 1 0 1 0 0 0 0
                1 0 1 0 0 1 0 0 0 0
                0 1 0 0 1 0 1 0 0 0
            ], Dict(1:10 .=> [-1, (1:9)...]))
        DT.split_triangle!(tri, 4, 6, 1, 9)
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG
    end

    @testset "Splitting two boundary edges" begin
        true_T = Set{NTuple{3,Int}}([
            (3, 2, 5),
            (4, 1, 5),
            (5, 2, 10), (2, 8, 10), (8, 5, 10),
            (6, 3, 1),
            (4, 6, 9), (6, 1, 9), (1, 4, 9),
            (1, 3, 7), (3, 5, 7), (5, 1, 7),
            (6, 2, 3),
        ])
        true_adj = DefaultDict(DT.∅,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (5, 2) => 10, (2, 10) => 5, (10, 5) => 2,
                (2, 8) => 10, (8, 10) => 2, (10, 2) => 8,
                (8, 5) => 10, (5, 10) => 8, (10, 8) => 5,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
                (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
                (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
                (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
                (4, 5) => DT.𝒢,
                (5, 8) => DT.𝒢,
                (8, 2) => DT.𝒢,
                (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢
            )
        )
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
            2 => Set{NTuple{2,Int}}([(8, 10), (10, 5), (5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(6, 9), (9, 1), (1, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 10), (10, 8)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
            7 => Set{NTuple{2,Int}}([(1, 3), (3, 5), (5, 1)]),
            8 => Set{NTuple{2,Int}}([(5, 10), (10, 2)]),
            9 => Set{NTuple{2,Int}}([(6, 1), (1, 4), (4, 6)]),
            10 => Set{NTuple{2,Int}}([(5, 2), (2, 8), (8, 5)])
        )
        true_DG = _make_graph_from_adjacency([
                0 0 1 0 1 1 1 0 1 0 0
                0 0 0 1 1 1 1 1 0 1 0
                1 0 0 1 0 1 1 0 1 0 1
                0 1 1 0 0 1 1 1 0 0 0
                1 1 0 0 0 1 1 0 0 1 0
                1 1 1 1 1 0 0 1 1 0 1
                1 1 1 1 1 0 0 0 0 1 0
                0 1 0 1 0 1 0 0 0 0 0
                1 0 1 0 0 1 0 0 0 0 1
                0 1 0 0 1 0 1 0 0 0 0
                0 0 1 0 0 1 0 0 1 0 0
            ], Dict(1:11 .=> [-1, (1:10)...]))
        DT.split_triangle!(tri, 5, 2, 8, 10)
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG
    end
end

@testset "complete_split_triangle_and_legalise!, also with constrained edges" begin
    for _ in 1:100
        tri = fixed_shewchuk_example_constrained()
        function test_fnc(tri, i, j, k)
            p, q, r = get_point(tri, i, j, k)
            c = DT.triangle_centroid(p, q, r)
            DT.push_point!(tri, c)
            DT.complete_split_triangle_and_legalise!(tri, i, j, k, DT.num_points(tri))
            validate_statistics(tri)
            return validate_triangulation(tri)
        end
        @test test_fnc(tri, 10, 5, 8)
        @test test_fnc(tri, 4, 5, 9)
        @test test_fnc(tri, 9, 13, 12)
        @test test_fnc(tri, 9, 13, 14)
        @test test_fnc(tri, 10, 12, 11)
        @test test_fnc(tri, 3, 4, 10)
        @test test_fnc(tri, 17, 4, 9)
        @test test_fnc(tri, 10, 18, 9)
        @test test_fnc(tri, 10, 19, 9)
        @test test_fnc(tri, 16, 14, 12)
        @test test_fnc(tri, 16, 12, 8)
        @test test_fnc(tri, 10, 17, 19)
        @test test_fnc(tri, 3, 17, 10)
        @test test_fnc(tri, 3, 2, 17)
        add_segment!(tri, 4, 21)
        i, j, k = 21, 9, 4
        p, q, r = get_point(tri, i, j, k)
        c = DT.triangle_centroid(p, q, r)
        DT.push_point!(tri, c)
        DT.complete_split_triangle_and_legalise!(tri, i, j, k, DT.num_points(tri))
        validate_statistics(tri)
        @test validate_triangulation(tri)
        @test DT.contains_segment(tri, 4, 21)
        @test DT.edge_exists(tri, 4, 21) && DT.edge_exists(tri, 21, 4)
        add_segment!(tri, 1, 11)
        i, j, k = 11, 1, 9
        p, q, r = get_point(tri, i, j, k)
        c = DT.triangle_centroid(p, q, r)
        DT.push_point!(tri, c)
        DT.complete_split_triangle_and_legalise!(tri, i, j, k, DT.num_points(tri))
        validate_statistics(tri)
        @test validate_triangulation(tri)
        @test DT.contains_segment(tri, 1, 11)
        @test DT.edge_exists(tri, 1, 11) && DT.edge_exists(tri, 11, 1)
    end
end