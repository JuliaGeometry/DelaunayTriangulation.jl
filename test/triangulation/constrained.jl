using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
include("../helper_functions.jl")

@testset "Test random constrained Delaunay triangulations" begin
    rng = StableRNG(191919)
    for i in 1:250
        @show i
        points, edges, mat_edges = get_random_vertices_and_constrained_edges(40, 200, 20, rng)
        tri = triangulate(points; edges, rng)
        @test validate_triangulation(tri)
        empty!(get_all_constrained_edges(tri))
        @test !validate_triangulation(tri)
    end
    for i in 1:25
        @show i
        points, edges, mat_edges = get_random_vertices_and_constrained_edges(200, 1000, 200, rng)
        tri = triangulate(points; edges, rng)
        @test validate_triangulation(tri)
        empty!(get_all_constrained_edges(tri))
        @test !validate_triangulation(tri)
    end
end

@testset "Testing Shewchuk's PSLG example" begin
    pts, C = second_shewchuk_example_constrained()
    for _ in 1:500
        tri = triangulate(pts; edges=C, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Random parabolas" begin
    for i in 1:50
        @show i
        pts = [(2rand(rng) - 1, rand(rng)) for _ in 1:500]
        x = LinRange(-1, 1, 250)
        a = LinRange(0.0, 1.0, 8)
        C = Set{NTuple{2,Int64}}()
        for i in eachindex(a)
            y = a[i] * x .^ 2
            append!(pts, zip(x, y))
            push!(C, [(i, i + 1) for i in (500+250*(i-1)+1):(500+250*(i-1)+249)]...)
        end
        tri = triangulate(pts; edges=C, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Random collection of straight lines" begin
    for i in 1:10
        @show i
        pts = NTuple{2,Float64}[]
        C = Set{NTuple{2,Int64}}()
        j = 1
        for i in 1:100
            push!(pts, (2i / 101 - 1, 2rand() - 1))
            push!(pts, (2i / 101 - 1, 2rand() - 1))
            push!(C, (j, j + 1))
            j += 2
        end
        x1 = LinRange(-1, 1 - 1e-12, 100)
        y1 = LinRange(-1, -1, 100)
        x2 = LinRange(1, 1, 100)
        y2 = LinRange(-1, 1 - 1e-12, 100)
        x3 = LinRange(1, -1 + 1e-12, 100)
        y3 = LinRange(1, 1, 100)
        x4 = LinRange(-1, -1, 100)
        y4 = LinRange(1, -1 + 1e-12, 100)
        append!(pts, zip(x1, y1), zip(x2, y2), zip(x3, y3), zip(x4, y4))
        push!(C, [(j, j + 1) for j in 201:299]...)
        push!(C, [(j, j + 1) for j in 301:399]...)
        push!(C, [(j, j + 1) for j in 401:499]...)
        push!(C, [(j, j + 1) for j in 501:599]...)
        tri = triangulate(pts; edges=C)
        @test validate_triangulation(tri)
    end
end