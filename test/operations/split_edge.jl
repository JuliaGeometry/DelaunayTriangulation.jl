using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StableRNGs

include("../helper_functions.jl")

@testset "Specific example" begin
    tri = example_triangulation()
    DT.add_triangle!(tri, 6, 2, 3)
    DT.split_triangle!(tri, 1, 3, 5, 7)
    DT.flip_edge!(tri, 1, 5)
    p9 = @SVector[2.5, 0.5]
    p10 = @SVector[1.5, 2.0]
    push!(get_points(tri), p9, p10)

    @testset "Interior edge" begin
        i, j, r = 1, 7, 8
        DT.split_edge!(tri, i, j, r)
        DT.split_edge!(tri, j, i, r)
        true_T = Set{NTuple{3,Int}}([
            (3, 2, 5),
            (1, 8, 4),
            (7, 8, 3),
            (3, 5, 7),
            (6, 3, 1),
            (4, 6, 1),
            (4, 7, 5),
            (6, 2, 3),
            (8, 7, 4),
            (8, 1, 3),
        ])
        true_adj = DefaultDict(DT.∅,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
                (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
                (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
                (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
                (4, 5) => DT.𝒢,
                (5, 2) => DT.𝒢,
                (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢,
            ))
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
            2 => Set{NTuple{2,Int}}([(5, 3), (3, 6)]),
            3 => Set{NTuple{2,Int}}([(6, 2), (2, 5), (5, 7), (7, 8), (8, 1), (1, 6)]),
            4 => Set{NTuple{2,Int}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
            5 => Set{NTuple{2,Int}}([(4, 7), (7, 3), (3, 2)]),
            6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 4)]),
            7 => Set{NTuple{2,Int}}([(3, 5), (5, 4), (4, 8), (8, 3)]),
            8 => Set{NTuple{2,Int}}([(1, 3), (3, 7), (7, 4), (4, 1)])
        )
        true_DG = _make_graph_from_adjacency(
            [
                0 0 1 0 1 1 1 0 0
                0 0 0 1 1 0 1 0 1
                1 0 0 1 0 1 1 0 0
                0 1 1 0 0 1 1 1 1
                1 1 0 0 0 1 1 1 1
                1 0 1 1 1 0 0 1 0
                1 1 1 1 1 0 0 0 0
                0 0 0 1 1 1 0 0 1
                0 1 0 1 1 0 0 1 0
            ], Dict(1:9 .=> [-1, (1:8)...]))
        DT.clear_empty_features!(tri)
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG

        @testset "Interior edge for a triangle that is on the boundary" begin
            i, j, r = 5, 3, 9
            DT.split_edge!(tri, i, j, r)
            DT.split_edge!(tri, j, i, r)
            true_T = Set{NTuple{3,Int}}([
                (1, 8, 4),
                (7, 8, 3),
                (6, 3, 1),
                (4, 6, 1),
                (4, 7, 5),
                (6, 2, 3),
                (8, 7, 4),
                (8, 1, 3),
                (9, 5, 7),
                (3, 9, 7),
                (9, 3, 2),
                (5, 9, 2)
            ])
            true_adj = DefaultDict(DT.∅,
                Dict(
                    (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
                    (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
                    (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                    (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                    (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
                    (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                    (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
                    (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
                    (7, 9) => 5, (9, 5) => 7, (5, 7) => 9,
                    (3, 9) => 7, (9, 7) => 3, (7, 3) => 9,
                    (9, 3) => 2, (3, 2) => 9, (2, 9) => 3,
                    (5, 9) => 2, (9, 2) => 5, (2, 5) => 9,
                    (4, 5) => DT.𝒢,
                    (5, 2) => DT.𝒢,
                    (2, 6) => DT.𝒢,
                    (6, 4) => DT.𝒢,
                ))
            true_adj2v = Dict(
                DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
                1 => Set{NTuple{2,Int}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
                2 => Set{NTuple{2,Int}}([(5, 9), (9, 3), (3, 6)]),
                3 => Set{NTuple{2,Int}}([(2, 9), (9, 7), (7, 8), (8, 1), (1, 6), (6, 2)]),
                4 => Set{NTuple{2,Int}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
                5 => Set{NTuple{2,Int}}([(4, 7), (7, 9), (9, 2)]),
                6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 4)]),
                7 => Set{NTuple{2,Int}}([(9, 5), (5, 4), (4, 8), (8, 3), (3, 9)]),
                8 => Set{NTuple{2,Int}}([(1, 3), (3, 7), (7, 4), (4, 1)]),
                9 => Set{NTuple{2,Int}}([(2, 5), (5, 7), (7, 3), (3, 2)])
            )
            true_DG = _make_graph_from_adjacency(
                [
                    0 0 1 0 1 1 1 0 0 0
                    0 0 0 1 1 0 1 0 1 0
                    1 0 0 1 0 1 1 0 0 1
                    0 1 1 0 0 0 1 1 1 1
                    1 1 0 0 0 1 1 1 1 0
                    1 0 1 0 1 0 0 1 0 1
                    1 1 1 1 1 0 0 0 0 0
                    0 0 0 1 1 1 0 0 1 1
                    0 1 0 1 1 0 0 1 0 0
                    0 0 1 1 0 1 0 1 0 0
                ], Dict(1:10 .=> [-1, (1:9)...]))
            DT.clear_empty_features!(tri)
            @test get_triangles(tri) == true_T
            @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
            @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
            @test get_graph(tri) == true_DG
        end
    end
end

@testset "complete_split_edge_and_legalise!" begin
    @testset "Free edges and constrained segments" begin
        for i in 1:100
            tri = fixed_shewchuk_example_constrained()
            DT.push_point!(tri, 4.0, 1.5)
            DT.complete_split_edge_and_legalise!(tri, 9, 10, DT.num_points(tri))
            validate_triangulation(tri)
            DT.push_point!(tri, 4.0, 1.3)
            DT.complete_split_edge_and_legalise!(tri, 12, 9, DT.num_points(tri))
            validate_triangulation(tri)
            DT.push_point!(tri, 4.0 + 1e-13, 2.2 - 1e-9)
            DT.complete_split_edge_and_legalise!(tri, 12, 10, DT.num_points(tri))
            validate_triangulation(tri)
            DT.push_point!(tri, 4.0, 0.0)
            DT.complete_split_edge_and_legalise!(tri, 4, 5, DT.num_points(tri))
            @test DT.is_ghost_vertex(get_adjacent(tri, 15, 4))
            @test DT.is_ghost_vertex(get_adjacent(tri, 5, 15))
            add_segment!(tri, 2, 9)
            DT.push_point!(tri, 2.0, 1.0)
            DT.complete_split_edge_and_legalise!(tri, 2, 9, DT.num_points(tri))
            @test compare_edge_vectors(collect(get_interior_segments(tri)), [(16, 9), (2, 16)])
            @test compare_edge_vectors(collect(get_all_segments(tri)), [(16, 9), (2, 16)])
            @test !DT.has_boundary_nodes(tri)
            @test validate_triangulation(tri)
            validate_statistics(tri)
            DT.push_point!(tri, 0.5, 1.000000000000000691982)
            DT.complete_split_edge_and_legalise!(tri, 2, 16, DT.num_points(tri))
            @test validate_triangulation(tri)
            validate_statistics(tri)
            @test compare_edge_vectors(collect(get_interior_segments(tri)), [(16, 9), (2, 17), (17, 16)])
            @test compare_edge_vectors(collect(get_all_segments(tri)), [(16, 9), (2, 17), (17, 16)])
            DT.push_point!(tri, 0.5, -0.0000000000000000000001293)
            DT.complete_split_edge_and_legalise!(tri, 1, 4, DT.num_points(tri))
            @test validate_triangulation(tri)
            validate_statistics(tri)
        end
    end

    @testset "Boundary segments" begin
        for i in 1:100
            curve_1 = [[
                (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
                (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
                (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
                (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
                (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
                (-4.0, 0.0), (0.0, 0.0),
            ]]
            curve_2 = [[
                (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
                (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
                (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
            ]]
            curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
            curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
            curves = [curve_1, curve_2, curve_3, curve_4]
            points = [
                (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
                (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
                (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
                (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
                (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
                (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
                (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
                (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
                (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
                (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
                (7.0, 7.8), (7.0, 7.4)]
            boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
            uncons_tri = triangulate(points)
            cons_tri = triangulate(points; boundary_nodes=boundary_nodes)
            tri = cons_tri
            @test validate_triangulation(tri)
            validate_statistics(tri)
            add_ghost_triangles!(tri)

            function test_fnc(tri, i, j)
                p, q = get_point(tri, i, j)
                px, py = DT.getxy(p)
                qx, qy = DT.getxy(q)
                mx, my = (px + qx) / 2, (py + qy) / 2
                DT.push_point!(tri, mx, my)
                DT.complete_split_edge_and_legalise!(tri, i, j, DT.num_points(tri))
                flag = validate_triangulation(tri)
                validate_statistics(tri)
                return flag
            end

            @test test_fnc(tri, 17, 57)
            @test test_fnc(tri, 98, 20)
            @test test_fnc(tri, 99, 20)
            @test test_fnc(tri, 55, 17)

            i, j = 97, 96
            p, q = get_point(tri, i, j)
            px, py = DT.getxy(p)
            qx, qy = DT.getxy(q)
            mx, my = (px + qx) / 2, (py + qy) / 2
            DT.push_point!(tri, mx, my)
            DT.complete_split_edge_and_legalise!(tri, i, j, DT.num_points(tri))
            @test get_adjacent(tri, 97, 102) == DT.𝒢 - 3
            @test get_adjacent(tri, 102, 96) == DT.𝒢 - 3
            flag = validate_triangulation(tri)
            @test flag
            @test get_adjacent2vertex(tri, DT.𝒢) == Set((
                (70, 69), (69, 68), (68, 67), (67, 66), (66, 65), (65, 64), (64, 63),
                (63, 62), (62, 61), (61, 60), (60, 59), (59, 58), (58, 57),
                (57, 56), (56, 55), (55, 78), (78, 77), (77, 76), (76, 75),
                (75, 74), (74, 73), (73, 72), (72, 71), (71, 70),
            ))
            @test get_adjacent2vertex(tri, DT.𝒢 - 1) == Set((
                (79, 89), (89, 88), (88, 87), (87, 86), (86, 85), (85, 84),
                (84, 83), (83, 82), (82, 81), (81, 80), (80, 79)
            ))
            @test get_adjacent2vertex(tri, DT.𝒢 - 2) == Set((
                (90, 93), (93, 92), (92, 91), (91, 90)
            ))
            @test get_adjacent2vertex(tri, DT.𝒢 - 3) == Set((
                (94, 97), (97, 102), (102, 96), (96, 95), (95, 94)
            ))
            validate_statistics(tri)

            @test test_fnc(tri, 87, 86)
            @test test_fnc(tri, 93, 92)
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 93, DT.num_points(tri))
            @test test_fnc(tri, 78, 55)
            @test test_fnc(tri, 55, 56)
            @test test_fnc(tri, 58, 59)
            @test test_fnc(tri, 62, 63)
            @test test_fnc(tri, 65, 66)
            @test test_fnc(tri, 94, 95)
            add_segment!(tri, 101, 1)
            @test test_fnc(tri, 15, 16)
            @test test_fnc(tri, 15, 23)
            @test test_fnc(tri, 78, 77)
        end
    end
end