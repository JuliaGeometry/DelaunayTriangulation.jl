using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using Test
using .HelperFunctions

@testset "Multiple boundaries" begin
    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri = DT.triangulate_rectangle(a, b, c, d, nx, ny)
    @test num_solid_vertices(tri) == nx * ny 
    @test num_ghost_vertices(tri) == 4
    @test DT.num_vertices(tri) == nx * ny + 4
    @test validate_triangulation(tri)
    for PT in subtypes(DT.AbstractPredicateKernel)
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; predicates = PT())
        if PT() == DT.FastKernel()
            @test_broken validate_triangulation(tri; predicates = PT())
        else
            @test validate_triangulation(tri; predicates = PT())
        end
    end
end

@testset "Single boundary" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        a, b, c, d = 2.0, 10.0, -5.0, 7.5
        nx = 20
        ny = 10
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; single_boundary = true, predicates = PT())
        @test num_solid_vertices(tri) == nx * ny
        @test num_ghost_vertices(tri) == 1
        @test DT.num_vertices(tri) == nx * ny + 1
        PT() != DT.FastKernel() && @test validate_triangulation(tri; predicates = PT())
        bn = reduce(vcat, [1:20, 20:20:200, 200:-1:181, 181:-20:1])
        unique!(bn)
        push!(bn, 1)
        @test tri.boundary_nodes == bn
    end
end

points = [(-1.0,0.0),(1.0,0.0),(0.0,1.0)]
triangles = Set(((1,2,3),))
boundary_nodes = [1, 2, 3, 1]
delete_ghosts = false
predicates = AdaptiveKernel()

_bn = copy(boundary_nodes)
tri = Triangulation(points; boundary_nodes = _bn)

boundary_nodes = _bn
Es = DT.edges_type(tri)
points = DT.get_points(tri)
polygon_hierarchy = DT.get_polygon_hierarchy(tri)
DT.construct_polygon_hierarchy!(polygon_hierarchy, points, boundary_nodes)
for τ in DT.each_triangle(triangles)
    DT.is_ghost_triangle(τ) && continue 
    DT.add_triangle!(tri, τ)
end
