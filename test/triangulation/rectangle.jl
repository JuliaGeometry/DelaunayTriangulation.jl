using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

include("../test_setup.jl")

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

@testset "Multiple boundaries" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny)
    fig, ax, sc = triplot(tri; recompute_centers=true, show_ghost_edges=true)
    xlims!(ax, a - 0.5, b + 0.5)
    ylims!(ax, c - 0.5, d + 0.5)
    lines!(ax, tri.points[:, tri.boundary_nodes[1]]; linewidth=4)
    lines!(ax, tri.points[:, tri.boundary_nodes[2]]; linewidth=4)
    lines!(ax, tri.points[:, tri.boundary_nodes[3]]; linewidth=4)
    lines!(ax, tri.points[:, tri.boundary_nodes[4]]; linewidth=4)
    SAVE_FIGURE && save("$save_path/rectangular_triangulation_1.png", fig)
    @test tri.boundary_nodes[1] == 1:20
    @test tri.boundary_nodes[2] == 20:20:200
    @test tri.boundary_nodes[3] == 200:-1:181
    @test tri.boundary_nodes[4] == 181:-20:1
end

@testset "Single boundary" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true)
    fig, ax, sc = triplot(tri; recompute_centers=true, show_ghost_edges=true)
    xlims!(ax, a - 0.5, b + 0.5)
    ylims!(ax, c - 0.5, d + 0.5)
    SAVE_FIGURE && save("$save_path/rectangular_triangulation_2.png", fig)
    bn = reduce(vcat, [1:20, 20:20:200, 200:-1:181, 181:-20:1])
    unique!(bn)
    push!(bn, 1)
    @test tri.boundary_nodes == bn
end
