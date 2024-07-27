using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using Test

#

@testset "Multiple boundaries" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny)
    @test validate_triangulation(tri)
    for PT in subtypes(DT.AbstractPredicateKernel)
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; predicates=PT())
        if PT() == DT.Fast()
            @test_broken validate_triangulation(tri; predicates=PT())
        else
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Single boundary" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        a, b, c, d = 2.0, 10.0, -5.0, 7.5
        nx = 20
        ny = 10
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true, predicates=PT())
        PT() != DT.Fast() && @test validate_triangulation(tri; predicates=PT())
        bn = reduce(vcat, [1:20, 20:20:200, 200:-1:181, 181:-20:1])
        unique!(bn)
        push!(bn, 1)
        @test tri.boundary_nodes == bn
    end
end
