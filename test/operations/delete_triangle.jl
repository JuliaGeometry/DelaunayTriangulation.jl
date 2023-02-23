using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)

fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, tri; show_ghost_edges=true, ghost_edge_extension_factor=2.0)
xlims!(ax, -2, 22)
ylims!(ax, -2, 22)
all_i = Set{NTuple{3,Int64}}()
for (a, b, c) in [("a1", "f", "g"),
                  ("a1", "g", "i"),
                  ("a1", "i", "f"),
                  ("f", "i", "ℓ"),
                  ("f", "v", "e"),
                  ("f", "w", "v"),
                  ("e", "v", "z")]
    i, j, k = index_map[a], index_map[b], index_map[c]
    push!(all_i, (i, j, k))
    DT.delete_triangle!(tri, i, j, k)
    @test !DT.contains_triangle(tri, i, j, k)[2]
    @test !DT.contains_triangle(tri, j, k, i)[2]
    @test !DT.contains_triangle(tri, k, i, j)[2]
    @test DT.get_adjacent(tri, i, j) == DT.DefaultAdjacentValue
    @test DT.get_adjacent(tri, j, k) == DT.DefaultAdjacentValue
    @test DT.get_adjacent(tri, k, i) == DT.DefaultAdjacentValue
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

ax = Axis(fig[1, 2])
triplot!(ax, tri; show_ghost_edges=true, ghost_edge_extension_factor=2.0)
xlims!(ax, -2, 22)
ylims!(ax, -2, 22)

save("$save_path/test_delete_triangles.png", fig)
