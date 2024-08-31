using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures
using ..DelaunayTriangulation: Triangulation
using StaticArrays


global tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)

@testset "Deleting a triangle" begin
    all_i = Set{NTuple{3, Int}}()
    for (a, b, c) in [
            ("a1", "f", "g"),
            ("a1", "g", "i"),
            ("a1", "i", "f"),
            ("f", "i", "ℓ"),
            ("f", "v", "e"),
            ("f", "w", "v"),
            ("e", "v", "z"),
        ]
        i, j, k = index_map[a], index_map[b], index_map[c]
        push!(all_i, (i, j, k))
        DT.delete_triangle!(tri, i, j, k; protect_boundary = true, update_ghost_edges = false)
        @test !DT.contains_triangle(tri, i, j, k)[2]
        @test !DT.contains_triangle(tri, j, k, i)[2]
        @test !DT.contains_triangle(tri, k, i, j)[2]
        @test DT.get_adjacent(tri, i, j) == DT.∅
        @test DT.get_adjacent(tri, j, k) == DT.∅
        @test DT.get_adjacent(tri, k, i) == DT.∅
        @test !DT.contains_edge(i, j, DT.get_adjacent2vertex(tri, k))
        @test !DT.contains_edge(j, k, DT.get_adjacent2vertex(tri, i))
        @test !DT.contains_edge(k, i, DT.get_adjacent2vertex(tri, j))
    end
    for T in each_triangle(tri)
        i, j, k = T
        if !DT.contains_triangle(i, j, k, all_i)[2]
            @test DT.contains_triangle(tri, i, j, k)[2]
            @test DT.contains_triangle(tri, j, k, i)[2]
            @test DT.contains_triangle(tri, k, i, j)[2]
            @test DT.get_adjacent(tri, i, j) == k
            @test DT.get_adjacent(tri, j, k) == i
            @test DT.get_adjacent(tri, k, i) == j
            @test DT.contains_edge(i, j, DT.get_adjacent2vertex(tri, k))
            @test DT.contains_edge(j, k, DT.get_adjacent2vertex(tri, i))
            @test DT.contains_edge(k, i, DT.get_adjacent2vertex(tri, j))
        end
    end
end

@testset "Simpler test" begin
    tri = example_triangulation()
    DT.add_triangle!(tri, 1, 3, 7)
    DT.add_triangle!(tri, 5, 2, 8)
    DT.add_triangle!(tri, 6, 2, 3)

    @testset "Deleting an interior triangle" begin
        i, j, k = 1, 3, 7
        DT.delete_triangle!(tri, i, j, k)
        true_T = Set{NTuple{3, Int}}(
            [
                (3, 2, 5),
                (4, 1, 5),
                (6, 3, 1),
                (4, 6, 1),
                (5, 1, 3),
                (5, 2, 8),
                (6, 2, 3),
            ],
        )
        true_adj =
            Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
            (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
            (4, 5) => DT.𝒢, (5, 8) => DT.𝒢,
            (8, 2) => DT.𝒢, (2, 6) => DT.𝒢,
            (6, 4) => DT.𝒢,
        )
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2, Int}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
            2 => Set{NTuple{2, Int}}([(5, 3), (8, 5), (3, 6)]),
            3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
            4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
            5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
            6 => Set{NTuple{2, Int}}([(3, 1), (1, 4), (2, 3)]),
            7 => Set{NTuple{2, Int}}([]),
            8 => Set{NTuple{2, Int}}([(5, 2)]),
        )
        true_DG = _make_graph_from_adjacency(
            [
                0 0 1 0 1 1 1 0 1
                0 0 0 1 1 1 1 0 0
                1 0 0 1 0 1 1 0 1
                0 1 1 0 0 1 1 0 0
                1 1 0 0 0 1 1 0 0
                1 1 1 1 1 0 0 0 1
                1 1 1 1 1 0 0 0 0
                0 0 0 0 0 0 0 0 0
                1 0 1 0 0 1 0 0 0
            ], Dict(1:9 .=> [-1, (1:8)...]),
        )
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG
    end

    @testset "Deleting a triangle with two boundary edges" begin
        for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
            local true_T, true_adj, true_adj2v, true_DG
            _tri = deepcopy(tri)
            DT.delete_triangle!(_tri, i, j, k)
            true_T = Set{NTuple{3, Int}}(
                [
                    (3, 2, 5),
                    (4, 1, 5),
                    (6, 3, 1),
                    (4, 6, 1),
                    (5, 1, 3),
                    (6, 2, 3),
                ],
            )
            true_adj =
                Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
                (4, 5) => DT.𝒢,
                (5, 2) => DT.𝒢, (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢,
            )
            true_adj2v = Dict(
                DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
                1 => Set{NTuple{2, Int}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
                2 => Set{NTuple{2, Int}}([(5, 3), (3, 6)]),
                3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
                4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
                5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2)]),
                6 => Set{NTuple{2, Int}}([(3, 1), (1, 4), (2, 3)]),
            )
            true_DG = _make_graph_from_adjacency(
                [
                    0 0 1 0 1 1 1
                    0 0 0 1 1 1 1
                    1 0 0 1 0 1 1
                    0 1 1 0 0 1 1
                    1 1 0 0 0 1 1
                    1 1 1 1 1 0 0
                    1 1 1 1 1 0 0
                ], Dict(1:7 .=> [-1, (1:6)...]),
            )
            DT.clear_empty_features!(_tri)
            @test get_triangles(_tri) == true_T
            @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
            @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
            @test get_graph(_tri) == true_DG
        end
        i, j, k = 5, 2, 8
        DT.delete_triangle!(tri, i, j, k)
    end

    @testset "Deleting a triangle with a single boundary edge" begin
        for (i, j, k) in ((6, 2, 3), (2, 3, 6), (3, 6, 2))
            local true_T, true_adj, true_adj2v, true_DG
            _tri = deepcopy(tri)
            DT.delete_triangle!(_tri, i, j, k)
            true_T = Set{NTuple{3, Int}}(
                [
                    (3, 2, 5),
                    (4, 1, 5),
                    (6, 3, 1),
                    (4, 6, 1),
                    (5, 1, 3),
                ],
            )
            true_adj =
                Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (4, 5) => DT.𝒢,
                (5, 2) => DT.𝒢,
                (6, 4) => DT.𝒢,
                (2, 3) => DT.𝒢, (3, 6) => DT.𝒢,
            )
            true_adj2v = Dict(
                DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
                1 => Set{NTuple{2, Int}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
                2 => Set{NTuple{2, Int}}([(5, 3)]),
                3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (1, 6)]),
                4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
                5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2)]),
                6 => Set{NTuple{2, Int}}([(3, 1), (1, 4)]),
            )
            true_DG = _make_graph_from_adjacency(
                [
                    0 0 1 1 1 1 1
                    0 0 0 1 1 1 1
                    1 0 0 1 0 1 0
                    1 1 1 0 0 1 1
                    1 1 0 0 0 1 1
                    1 1 1 1 1 0 0
                    1 1 0 1 1 0 0
                ], Dict(1:7 .=> [-1, (1:6)...]),
            )
            DT.clear_empty_features!(_tri)
            @test get_triangles(_tri) == true_T
            @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
            @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
            @test get_graph(_tri) == true_DG
        end
    end
end

@testset "Deleting the only triangle of a triangulation" begin
    tri = example_empty_triangulation()
    DT.add_triangle!(tri, 1, 2, 3)
    true_T = Set{NTuple{3, Int}}([])
    true_adj = DefaultDict(DT.∅, Dict())
    true_adj2v = Dict(
        DT.𝒢 => Set{NTuple{2, Int}}(),
        1 => Set{NTuple{2, Int}}(),
        2 => Set{NTuple{2, Int}}(),
        3 => Set{NTuple{2, Int}}(),
    )
    true_DG = _make_graph_from_adjacency([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], Dict(1:4 .=> [-1, 1, 2, 3]))
    DT.delete_triangle!(tri, 1, 2, 3)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
    @test get_graph(tri) == true_DG
end

@testset "Testing all the boundary deletion cases" begin
    p1 = @SVector[0.0, 0.0]
    p2 = @SVector[1.0, 0.0]
    p3 = @SVector[0.0, 1.0]
    pts = [p1, p2, p3]
    tri = Triangulation(pts)
    DT.add_triangle!(tri, 1, 2, 3; update_ghost_edges = true)
    p4 = @SVector[1.7, 1.7]
    push!(pts, p4)
    DT.add_triangle!(tri, 3, 2, 4; update_ghost_edges = true)
    p5 = @SVector[1.0, 3.0]
    p6 = @SVector[3.0, 1.0]
    push!(pts, p5, p6)
    DT.add_triangle!(tri, 3, 4, 5; update_ghost_edges = true)
    DT.add_triangle!(tri, 4, 2, 6; update_ghost_edges = true)
    DT.add_triangle!(tri, 5, 4, 6; update_ghost_edges = true)
    DT.delete_triangle!(tri, 5, 4, 6; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (3, 2, 4),
            (3, 4, 5),
            (4, 2, 6),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 5, DT.𝒢),
            (5, 4, DT.𝒢),
            (4, 6, DT.𝒢),
            (6, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
            (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
            (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 5) => DT.𝒢, (5, DT.𝒢) => 3, (DT.𝒢, 3) => 5,
            (5, 4) => DT.𝒢, (4, DT.𝒢) => 5, (DT.𝒢, 5) => 4,
            (4, 6) => DT.𝒢, (6, DT.𝒢) => 4, (DT.𝒢, 4) => 6,
            (6, 2) => DT.𝒢, (2, DT.𝒢) => 6, (DT.𝒢, 6) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(1, 3), (3, 5), (5, 4), (4, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (3, DT.𝒢), (DT.𝒢, 2)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (4, 3), (6, 4), (1, DT.𝒢), (DT.𝒢, 6)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (2, 4), (4, 5), (DT.𝒢, 1), (5, DT.𝒢)]),
            4 => Set{NTuple{2, Int}}([(3, 2), (5, 3), (2, 6), (DT.𝒢, 5), (6, DT.𝒢)]),
            5 => Set{NTuple{2, Int}}([(3, 4), (DT.𝒢, 3), (4, DT.𝒢)]),
            6 => Set{NTuple{2, Int}}([(4, 2), (DT.𝒢, 4), (2, DT.𝒢)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 3, 5, 4, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4]...)
    DT.add_neighbour!(true_DG, 6, [2, 4]...)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.delete_triangle!(tri, 4, 2, 6; update_ghost_edges = true)
    DT.delete_triangle!(tri, 3, 4, 5; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 4, DT.𝒢),
            (4, 2, DT.𝒢),
            (3, 2, 4),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 4) => DT.𝒢, (4, DT.𝒢) => 3, (DT.𝒢, 3) => 4,
            (4, 2) => DT.𝒢, (2, DT.𝒢) => 4, (DT.𝒢, 4) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(2, 1), (1, 3), (3, 4), (4, 2)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (DT.𝒢, 2), (3, DT.𝒢)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (1, DT.𝒢), (DT.𝒢, 4), (4, 3)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (DT.𝒢, 1), (4, DT.𝒢), (2, 4)]),
            4 => Set{NTuple{2, Int}}([(DT.𝒢, 3), (2, DT.𝒢), (3, 2)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 3, 4, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4]...)
    DT.add_neighbour!(true_DG, 4, [2, 3]...)
    DT.clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.delete_triangle!(tri, 3, 2, 4; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 2) => DT.𝒢, (2, DT.𝒢) => 3, (DT.𝒢, 3) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(2, 1), (1, 3), (3, 2)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (DT.𝒢, 2), (3, DT.𝒢)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (1, DT.𝒢), (DT.𝒢, 3)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (DT.𝒢, 1), (2, DT.𝒢)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 2, 3]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3]...)
    DT.add_neighbour!(true_DG, 3, [1, 2]...)
    DT.clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.delete_triangle!(tri, 1, 2, 3; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}()
    true_adj = DT.Adjacent{Int, NTuple{2, Int}}()
    true_adj2v = DT.Adjacent2Vertex{Int, Set{NTuple{2, Int}}}()
    true_DG = DT.Graph{Int}()
    DT.clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
end

@testset "A larger example" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    tri = Triangulation(pts)
    for (i, j, k) in (
            (1, 2, 6),
            (1, 6, 8),
            (9, 1, 8),
            (9, 8, 10),
            (10, 8, 11),
            (8, 7, 11),
            (8, 6, 7),
            (6, 2, 3),
            (6, 3, 4),
            (6, 4, 7),
            (7, 4, 5),
            (11, 7, 5),
            (10, 11, 5),
        )
        DT.add_triangle!(tri, i, j, k; update_ghost_edges = true)
    end
    [
    DT.delete_triangle!(tri, i, j, k; update_ghost_edges = true) for (i, j, k) in (
        (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1),
    )
    ]
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 6, 8),
            (6, 7, 8),
            (6, 4, 7),
            (8, 7, 11),
            (11, 7, 5),
            (10, 8, 11),
            (1, 8, DT.𝒢),
            (8, 10, DT.𝒢),
            (10, 11, DT.𝒢),
            (11, 5, DT.𝒢),
            (5, 7, DT.𝒢),
            (7, 4, DT.𝒢),
            (4, 6, DT.𝒢),
            (6, 1, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 6) => 8, (6, 8) => 1, (8, 1) => 6,
            (6, 7) => 8, (7, 8) => 6, (8, 6) => 7,
            (6, 4) => 7, (4, 7) => 6, (7, 6) => 4,
            (8, 7) => 11, (7, 11) => 8, (11, 8) => 7,
            (11, 7) => 5, (7, 5) => 11, (5, 11) => 7,
            (10, 8) => 11, (8, 11) => 10, (11, 10) => 8,
            (1, 8) => DT.𝒢, (8, DT.𝒢) => 1, (DT.𝒢, 1) => 8,
            (8, 10) => DT.𝒢, (10, DT.𝒢) => 8, (DT.𝒢, 8) => 10,
            (10, 11) => DT.𝒢, (11, DT.𝒢) => 10, (DT.𝒢, 10) => 11,
            (11, 5) => DT.𝒢, (5, DT.𝒢) => 11, (DT.𝒢, 11) => 5,
            (5, 7) => DT.𝒢, (7, DT.𝒢) => 5, (DT.𝒢, 5) => 7,
            (7, 4) => DT.𝒢, (4, DT.𝒢) => 7, (DT.𝒢, 7) => 4,
            (4, 6) => DT.𝒢, (6, DT.𝒢) => 4, (DT.𝒢, 4) => 6,
            (6, 1) => DT.𝒢, (1, DT.𝒢) => 6, (DT.𝒢, 6) => 1,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex{Int, Set{NTuple{2, Int}}}()
    for (ij, k) in true_adj.adjacent
        DT.add_adjacent2vertex!(true_adj2v, k, ij)
    end
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, 1, 8, 10, 11, 5, 7, 4, 6, 1)
    DT.add_neighbour!(true_DG, 1, 6, 8)
    DT.add_neighbour!(true_DG, 4, 6, 7)
    DT.add_neighbour!(true_DG, 5, 7, 11)
    DT.add_neighbour!(true_DG, 6, 1, 8, 7, 4)
    DT.add_neighbour!(true_DG, 7, 4, 6, 8, 11, 5)
    DT.add_neighbour!(true_DG, 8, 1, 6, 7, 11, 10)
    DT.add_neighbour!(true_DG, 10, 8, 11)
    DT.add_neighbour!(true_DG, 11, 10, 8, 7, 5)
    DT.clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
end
