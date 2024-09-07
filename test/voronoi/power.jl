using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes
using DataStructures
using StableRNGs
import GeometryBasics: Point2f
using Random
using StaticArrays
using LinearAlgebra
using StructEquality
using GeometryBasics
@struct_equal DT.Queue

for PT in subtypes(DT.AbstractPredicateKernel)
    for _ in 1:10
        A = (-1.0, 7.0)
        B = (4.0, 4.0)
        C = (-2.0, -1.0)
        D = (-1.0, 3.0)
        E = (3.0, -1.0)
        F = (1.0, 4.0)
        G = (-3.0, 5.0)
        points = [A, B, C, D, E, F, G]
        weights = [0.0, 0.1, 0.2, -2.0, 17.0, 0.9, 1.5]
        tri = triangulate(points; weights, randomise=false, predicates=PT())
        vorn = voronoi(tri)
        @test DT.is_weighted(vorn)
        for (i, p) in DT.get_generators(vorn)
            @test get_point(tri, i) == get_generator(vorn, i) == p
        end
        @test validate_tessellation(vorn; predicates=PT())
        @test DT.get_triangulation(vorn) == tri
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
            V = DT.sort_triangle(V)
            c = DT.get_triangle_to_circumcenter(vorn, V)
            c = get_polygon_point(vorn, c)
            i, j, k = triangle_vertices(V)
            p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
            cx, cy = DT.triangle_orthocenter(p, q, r, get_weight(tri, i), get_weight(tri, j), get_weight(tri, k))
            @test cx == c[1] && cy == c[2]
        end
        for (c, V) in circumcenter_to_triangle
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        for (V, c) in triangle_to_circumcenter
            @test DT.get_circumcenter_to_triangle(vorn, c) == V
            @test DT.get_triangle_to_circumcenter(vorn, V) == c
        end
        @test isempty(DT.get_boundary_polygons(vorn))
        poly1 = [7, 3, -3, -5, 7]
        poly2 = [3, 6, -1, -3, 3]
        poly3 = [1, 5, -4, -2, 1]
        poly4 = [4, 2, 5, 1, 4]
        poly5 = [1, -2, -1, 6, 4, 1]
        poly6 = [6, 3, 7, 2, 4, 6]
        poly7 = [2, 7, -5, -4, 5, 2]
        polys = [poly1, poly2, poly3, poly4, poly5, poly6, poly7]
        for (i, poly) in enumerate(polys)
            @test DT.circular_equality(get_polygon(vorn, i), poly)
        end
        polypoints = sort(reinterpret(reshape, Float64, DT.get_polygon_points(vorn))', dims=1)
        _polypoints = sort([
                2.6333333333333 3.3633333333333
                -0.765 5.14
                2.6333333333333 7.4055555555556
                -3.38 1.745
                -1.025 4.1
                -1.18 1.195
                -0.10833333333333 2.2666666666667
            ], dims=1)
        @test polypoints ≈ _polypoints
    end
end