using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

include("../test_setup.jl")
include("../helper_functions.jl")

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

@testset "Adding ghost triangles" begin
    tri, label_map, index_map = simple_geometry()

    fig = Figure()
    ax = Axis(fig[1, 1])
    xlims!(ax, -2, 22)
    ylims!(ax, -2, 22)
    triplot!(ax, tri)

    DT.add_ghost_triangles!(tri)
    outer_edges = [("a", "b") => DT.BoundaryIndex,
        ("b", "c") => DT.BoundaryIndex,
        ("c", "d") => DT.BoundaryIndex,
        ("d", "e") => DT.BoundaryIndex,
        ("e", "f") => DT.BoundaryIndex,
        ("f", "g") => DT.BoundaryIndex,
        ("g", "h") => DT.BoundaryIndex,
        ("h", "a") => DT.BoundaryIndex]
    inner_edges_1 = [("k", "j") => DT.BoundaryIndex - 1,
        ("j", "i") => DT.BoundaryIndex - 1,
        ("i", "ℓ") => DT.BoundaryIndex - 1,
        ("ℓ", "k") => DT.BoundaryIndex - 1]
    inner_edges_2 = [("r", "q") => DT.BoundaryIndex - 2,
        ("q", "p") => DT.BoundaryIndex - 2,
        ("p", "o") => DT.BoundaryIndex - 3,
        ("o", "n") => DT.BoundaryIndex - 3,
        ("n", "m") => DT.BoundaryIndex - 3,
        ("m", "r") => DT.BoundaryIndex - 3]

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

    ax = Axis(fig[1, 2])
    triplot!(ax, tri; show_ghost_edges=true, recompute_centers=true, ghost_edge_extension_factor=2.0)
    xlims!(ax, -2, 22)
    ylims!(ax, -2, 22)
    SAVE_FIGURE && save("$save_path/test_add_ghost_triangles.png", fig)
end