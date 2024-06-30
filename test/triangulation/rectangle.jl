using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie

Base.include(@__MODULE__, "../helper_functions.jl")

@testset "Multiple boundaries" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny)
    @test validate_triangulation(tri)
end

@testset "Single boundary" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true)
    @test validate_triangulation(tri)
    bn = reduce(vcat, [1:20, 20:20:200, 200:-1:181, 181:-20:1])
    unique!(bn)
    push!(bn, 1)
    @test tri.boundary_nodes == bn
end
