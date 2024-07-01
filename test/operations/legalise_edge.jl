using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StaticArrays


@testset "Specific example" begin
    tri = example_triangulation()
    DT.add_triangle!(tri, 6, 2, 3)
    DT.split_triangle!(tri, 1, 3, 5, 7)
    @test DT.is_legal(DT.is_legal(tri, 1, 3))
    @test DT.is_legal(DT.is_legal(tri, 3, 5))
    i, j, r = 5, 1, 7
    e = DT.get_adjacent(tri, j, i)
    DT.legalise_edge!(tri, i, j, r)
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
    true_adj = DefaultDict(DT.∅,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (7, 5) => 4, (5, 4) => 7, (4, 7) => 5,
            (7, 4) => 1, (4, 1) => 7, (1, 7) => 4,
            (4, 5) => DT.𝒢,
            (5, 2) => DT.𝒢,
            (2, 6) => DT.𝒢,
            (6, 4) => DT.𝒢,
        )
    )
    true_adj2v = Dict(
        DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int}}([(6, 3), (3, 7), (7, 4), (4, 6)]),
        2 => Set{NTuple{2,Int}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int}}([(2, 5), (5, 7), (7, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int}}([(6, 1), (1, 7), (7, 5)]),
        5 => Set{NTuple{2,Int}}([(4, 7), (7, 3), (3, 2)]),
        6 => Set{NTuple{2,Int}}([(2, 3), (3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int}}([(3, 5), (5, 4), (4, 1), (1, 3)])
    )
    true_DG = _make_graph_from_adjacency([
            0 0 1 0 1 1 1 0
            0 0 0 1 1 0 1 1
            1 0 0 1 0 1 1 0
            0 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 1
            1 0 1 1 1 0 0 1
            1 1 1 1 1 0 0 0
            0 1 0 1 1 1 0 0
        ], Dict(1:8 .=> [-1, (1:7)...]))
    DT.clear_empty_features!(tri)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
    @test get_graph(tri) == true_DG
    @test all(DT.is_positively_oriented(DT.triangle_orientation(tri, T...)) for T in each_triangle(tri))
    for ((i, j), v) in get_adjacent(tri).adjacent
        @test DT.is_legal(tri, i, j) |> DT.is_legal
    end
end