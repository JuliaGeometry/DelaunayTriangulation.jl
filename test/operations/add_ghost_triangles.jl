using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StaticArrays


@testset "Adding ghost triangles" begin
    tri, label_map, index_map = simple_geometry()

    DT.add_ghost_triangles!(tri)
    outer_edges = [
        ("a", "b") => DT.𝒢,
        ("b", "c") => DT.𝒢,
        ("c", "d") => DT.𝒢,
        ("d", "e") => DT.𝒢,
        ("e", "f") => DT.𝒢,
        ("f", "g") => DT.𝒢,
        ("g", "h") => DT.𝒢,
        ("h", "a") => DT.𝒢,
    ]
    inner_edges_1 = [
        ("k", "j") => DT.𝒢 - 1,
        ("j", "i") => DT.𝒢 - 1,
        ("i", "ℓ") => DT.𝒢 - 1,
        ("ℓ", "k") => DT.𝒢 - 1,
    ]
    inner_edges_2 = [
        ("r", "q") => DT.𝒢 - 2,
        ("q", "p") => DT.𝒢 - 2,
        ("p", "o") => DT.𝒢 - 3,
        ("o", "n") => DT.𝒢 - 3,
        ("n", "m") => DT.𝒢 - 3,
        ("m", "r") => DT.𝒢 - 3,
    ]

    for ((a, b), k) in outer_edges
        i = index_map[a]
        j = index_map[b]
        @test DT.get_adjacent(tri, j, i) == k
        @test DT.get_adjacent(tri, i, k) == j
        @test DT.get_adjacent(tri, k, j) == i
        @test DT.contains_triangle(tri, j, i, k)[2]
        @test DT.contains_triangle(tri, j, i, k)[2]
        @test DT.contains_triangle(tri, i, k, j)[2]
        @test DT.contains_edge(j, i, DT.get_adjacent2vertex(tri, k))
        @test DT.contains_edge(i, k, DT.get_adjacent2vertex(tri, j))
        @test DT.contains_edge(k, j, DT.get_adjacent2vertex(tri, i))
    end
end
