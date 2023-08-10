using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes
using DataStructures
using StableRNGs
import GeometryBasics: Point2f
using StaticArrays
using LinearAlgebra

include("../helper_functions.jl")

@testset "Unconstrained test" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    pts = [A, B, C, D, E, F, G]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)

    vorn = voronoi(tri)
    for (i, p) in DT.get_generators(vorn)
        @test get_point(tri, i) == get_generator(vorn, i) == p
    end
    @test validate_tessellation(vorn)
    @test DT.get_triangulation(vorn) == tri
    circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
    for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
        V = DT.rotate_triangle_to_standard_form(V)
        c = DT.get_triangle_to_circumcenter(vorn, V)
        c = get_polygon_point(vorn, c)
        i, j, k = indices(V)
        p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
        cx, cy = DT.triangle_circumcenter(p, q, r)
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
    @test DT.circular_equality(
        get_polygon(vorn, 1),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (7, 4, 1), (4, 6, 1), (6, 2, 1), (1, 2, -1), (7, 1, -1), (7, 4, 1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 2),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (2, 5, -1), (1, 2, -1), (6, 2, 1), (6, 5, 2), (2, 5, -1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 3),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (5, 3, -1), (5, 4, 3), (4, 7, 3), (3, 7, -1), (5, 3, -1)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 4),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (5, 6, 4), (4, 6, 1), (7, 4, 1), (4, 7, 3), (5, 4, 3), (5, 6, 4)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 5),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (6, 5, 2), (5, 6, 4), (5, 4, 3), (5, 3, -1), (2, 5, -1), (6, 5, 2)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 6),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (6, 5, 2), (6, 2, 1), (4, 6, 1), (5, 6, 4), (6, 5, 2)
        ])
    )
    @test DT.circular_equality(
        get_polygon(vorn, 7),
        DT.get_triangle_to_circumcenter.(Ref(vorn), [
            (7, 4, 1), (7, 1, -1), (3, 7, -1), (4, 7, 3), (7, 4, 1)
        ])
    )

    xmin, xmax, ymin, ymax = DT.polygon_bounds(vorn, 0.5)
    fig, ax, sc = voronoiplot(vorn)
    triplot!(ax, tri)
    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)
    fig
end

@testset "Smaller example, checking ray coordinates" begin
    tri = example_triangulation()
    tri = triangulate(get_points(tri))
    vorn = voronoi(tri)
    @test validate_tessellation(vorn)
    bbox = DT.polygon_bounds(vorn, 0.1)
    xmin, xmax, ymin, ymax = bbox
    c1 = DT.get_polygon_coordinates(vorn, 1, bbox)
    c2 = DT.get_polygon_coordinates(vorn, 2, bbox)
    c3 = DT.get_polygon_coordinates(vorn, 3, bbox)
    c4 = DT.get_polygon_coordinates(vorn, 4, bbox)
    c5 = DT.get_polygon_coordinates(vorn, 5, bbox)
    c6 = DT.get_polygon_coordinates(vorn, 6, bbox)
    c7 = DT.get_polygon_coordinates(vorn, 7, bbox)
    @test DT.circular_equality(collect.(c1), collect.([
        (-1.5, 0.5)
        (0.166666666666, -1.1666666666666665)
        (1.0, 0.5)
        (1.0, 3.0)
        (-1.5, 0.5)
    ]), ≈)
    @test DT.circular_equality(collect.(c2), collect.([
        (0.5, -3.2)
        (5.7, -3.2)
        (5.7, -1.700000000000001)
        (3.5, 0.5)
        (0.5, -2.5)
        (0.5, -3.2)
    ]), ≈)
    @test DT.circular_equality(collect.(c3), collect.([
        (3.5, 0.5)
        (1.0, 0.5)
        (0.16666666666666666, -1.1666666666666665)
        (0.5, -2.5)
        (3.5, 0.5)
    ]), ≈)
    @test DT.circular_equality(collect.(c4), collect.([
        (1.5, 5.2)
        (-2.7, 5.2)
        (-2.7, 0.9000000000000001)
        (-1.5, 0.5)
        (1.0, 3.0)
        (1.5, 4.5)
        (1.5, 5.2)
    ]), ≈)
    @test DT.circular_equality(collect.(c5), collect.([
        (1.5, 4.5)
        (3.5, 0.5)
        (5.7, 2.6999999999999997)
        (5.7, 5.2)
        (1.5, 5.2)
        (1.5, 4.5)
    ]), ≈)
    @test DT.circular_equality(collect.(c6), collect.([
        (-2.7, 0.9000000000000001)
        (-2.7, -3.2)
        (0.5, -3.2)
        (0.5, -2.5)
        (0.16666666666666666, -1.1666666666666665)
        (-1.5, 0.5)
        (-2.7, 0.9000000000000001)
    ]), ≈)
    @test DT.circular_equality(collect.(c7), collect.([
        (1.5, 4.5)
        (1.0, 3.0)
        (1.0, 0.5)
        (3.5, 0.5)
        (1.5, 4.5)
    ]), ≈)
end

@testset "delete/add_polygon_adjacencies" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    pts = [A, B, C, D, E, F, G]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)
    vorn = voronoi(tri)
    @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
          get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
          5
    DT.delete_polygon_adjacent!(vorn, 5)
    @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
          get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
          DT.DefaultAdjacentValue
    DT.add_polygon_adjacent!(vorn, 5)
    @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
          get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
          5
end

@testset "Voronoi point location" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    pts = [A, B, C, D, E, F, G]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)
    vor = voronoi(tri)
    @test validate_tessellation(vor)
    xmin, xmax, ymin, ymax = DT.polygon_bounds(get_points(tri), get_convex_hull_indices(tri))
    p = NTuple{2,Float64}[]
    n = 100000
    while length(p) ≤ n # only going to test points that are inside the polygon
        pt = (xmin + rand() * (xmax - xmin), ymin + rand() * (ymax - ymin))
        if DT.distance_to_polygon(pt, get_points(tri), get_convex_hull_indices(tri)) ≥ 0
            push!(p, pt)
        end
    end
    for p in p
        u = jump_and_march(vor, p)
        all_dists = [norm(p .- get_generator(vor, i)) for i in sort(collect(each_generator(vor)))]
        @test findmin(all_dists)[2] == u
    end

    points = [
        (0.0, 0.0), (-1.0, 1.0), (-0.5, 1.0), (0.0, 1.0), (0.5, 1.0), (1.0, 1.0),
        (1.0, 0.8), (1.0, 0.0), (1.0, -0.5), (1.0, -1.0),
        (0.1, -1.0), (-0.8, -1.0), (-1.0, -1.0),
        (-1.0, -0.7), (-1.0, -0.1), (-1.0, 0.6),
        (-0.1, -0.8), (0.2, -0.8),
        (-0.6, -0.4), (0.9, 0.0), (-0.5, 0.5), (-0.4, 0.6), (-0.1, 0.8)
    ]
    tri = triangulate(points, delete_ghosts=false)
    vorn = voronoi(tri)
    @test validate_tessellation(vorn)
    xg = LinRange(-1, 1, 250)
    yg = LinRange(-1, 1, 250)
    x = vec([x for x in xg, _ in yg])
    y = vec([y for _ in xg, y in yg])
    for (ξ, η) in zip(x, y)
        p = (ξ, η)
        u = jump_and_march(vorn, p)
        @inferred jump_and_march(vorn, p)
        all_dists = [norm(p .- get_generator(vorn, i)) for i in sort(collect(each_generator(vorn)))]
        k = findmin(all_dists)[2]
        @test k == u
        for m in each_point_index(tri)
            u = jump_and_march(vorn, p, try_points=m)
            @test u == k
        end
    end
end

@testset "Clipping a simple VoronoiTessellation" begin
    for _ in 1:1000
        A = (-1.0, 7.0)
        B = (4.0, 4.0)
        C = (-2.0, -1.0)
        D = (-1.0, 3.0)
        E = (3.0, -1.0)
        F = (1.0, 4.0)
        G = (-3.0, 5.0)
        pts = [A, B, C, D, E, F, G]
        tri = triangulate(pts; delete_ghosts=false, randomise=true)
        triplot(tri)
        #lock_convex_hull!(tri)
        vorn = voronoi(tri, true)
        for (i, p) in DT.get_generators(vorn)
            @test get_point(tri, i) == get_generator(vorn, i) == p
        end
        @test DT.get_triangulation(vorn) == tri
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(vorn, V)
            c = get_polygon_point(vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
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
        voronoiplot(vorn)
        orig_pt = collect.([(0.5, 0.5)
            (2.5, 7.166666666666667)
            (1.1666666666666665, 1.1666666666666667)
            (-4.3, 1.7000000000000002)
            (-1.0, 5.0)
            (-0.75, 5.0)
            (2.5, 1.7000000000000002)
            (0.5, -1.0)
            (3.5, 1.5)
            (3.0, -1.0)
            (-2.7142857142857144, 3.2857142857142856)
            (-2.369565217391304, 1.2173913043478262)
            (2.5000000000000004, 4.8999999999999995)
            (0.7105263157894739, 5.973684210526315)
            (-2.0, 6.0)
            (-3.0, 5.0)
            (4.0, 4.0)
            (-2.0, -1.0)
            (-1.0, 7.0)])
        @test validate_tessellation(vorn)
        @test isempty(DT.get_unbounded_polygons(vorn))
        @test all([0 ≤ get_area(vorn, i) < Inf for i in each_polygon_index(vorn)])
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 5))), getindex.(Ref(orig_pt), [8, 10, 9, 7, 3, 1, 8]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 4))), getindex.(Ref(orig_pt), [12, 1, 3, 6, 5, 11, 12]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 6))), getindex.(Ref(orig_pt), [3, 7, 13, 14, 6, 3]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 7))), getindex.(Ref(orig_pt), [11, 5, 15, 16, 11]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 2))), getindex.(Ref(orig_pt), [7, 9, 17, 13, 7]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 3))), getindex.(Ref(orig_pt), [18, 8, 1, 12, 18]), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 1))), getindex.(Ref(orig_pt), [5, 6, 14, 19, 15, 5]), ≈)
        @test allunique(DT.get_polygon_points(vorn))
        for i in each_polygon_index(vorn)
            C = get_polygon(vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(vorn, v), get_points(tri), get_convex_hull_indices(tri))
                @test δ ≥ -1e-15
            end
        end
    end
end

@testset "Another example" begin
    for _ in 1:100
        tri = fixed_shewchuk_example_constrained()
        vorn = voronoi(tri, false)
        @test validate_tessellation(vorn)
        for (i, p) in DT.get_generators(vorn)
            @test get_point(tri, i) == get_generator(vorn, i) == p
        end
        @test DT.get_triangulation(vorn) == tri
        @test DT.get_unbounded_polygons(vorn) == Set((3, 10, 11, 7, 6, 5, 4, 1, 2))
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(vorn, V)
            c = get_polygon_point(vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
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
        @test DT.circular_equality(
            get_polygon(vorn, 1),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 2, 1), (1, 2, -1), (4, 1, -1), (4, 2, 1)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 2),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 2, 1), (4, 3, 2), (2, 3, -1), (1, 2, -1), (4, 2, 1)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 3),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (4, 3, 2), (4, 10, 3), (3, 10, -1), (2, 3, -1), (4, 3, 2)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 4),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (5, 9, 4), (9, 10, 4), (4, 10, 3), (4, 3, 2), (4, 2, 1), (4, 1, -1), (5, 4, -1), (5, 9, 4)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 5),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (6, 8, 5), (8, 10, 5), (10, 9, 5), (5, 9, 4), (5, 4, -1), (6, 5, -1), (6, 8, 5)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 6),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (7, 8, 6), (6, 8, 5), (6, 5, -1), (7, 6, -1), (7, 8, 6)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 7),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (11, 8, 7), (7, 8, 6), (7, 6, -1), (11, 7, -1), (11, 8, 7)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 8),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (7, 8, 6), (11, 8, 7), (11, 10, 8), (8, 10, 5), (6, 8, 5), (7, 8, 6)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 9),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (10, 9, 5), (9, 10, 4), (5, 9, 4), (10, 9, 5)
            ])
        )
        @test DT.circular_equality(
            get_polygon(vorn, 10),
            DT.get_triangle_to_circumcenter.(Ref(vorn), [
                (9, 10, 4), (10, 9, 5), (8, 10, 5), (11, 10, 8), (10, 11, -1), (3, 10, -1), (4, 10, 3), (9, 10, 4)
            ])
        )

        _vorn = voronoi(tri, true)
        @test validate_tessellation(_vorn)
        for (i, p) in DT.get_generators(_vorn)
            @test get_point(tri, i) == get_generator(_vorn, i) == p
        end
        @test DT.get_triangulation(_vorn) == tri
        circumcenter_to_triangle = DT.get_circumcenter_to_triangle(_vorn)
        triangle_to_circumcenter = DT.get_triangle_to_circumcenter(_vorn)
        for V in DT.each_solid_triangle(DT.get_triangulation(_vorn))
            V = DT.rotate_triangle_to_standard_form(V)
            c = DT.get_triangle_to_circumcenter(_vorn, V)
            c = get_polygon_point(_vorn, c)
            i, j, k = indices(V)
            p, q, r = get_point(DT.get_triangulation(_vorn), i, j, k)
            cx, cy = DT.triangle_circumcenter(p, q, r)
            @test cx == c[1] && cy == c[2]
        end
        for (c, V) in circumcenter_to_triangle
            @test DT.get_circumcenter_to_triangle(_vorn, c) == V
            @test DT.get_triangle_to_circumcenter(_vorn, V) == c
        end
        for (V, c) in triangle_to_circumcenter
            @test DT.get_circumcenter_to_triangle(_vorn, c) == V
            @test DT.get_triangle_to_circumcenter(_vorn, V) == c
        end
        orig_pt = [
            (7.25, 0.25)
            (2.375, 1.75)
            (5.625, 1.75)
            (2.0, 2.05)
            (5.815217391304348, 1.9021739130434783)
            (6.0, 2.333333333333333)
            (1.625, 1.75)
            (4.0, -1.5)
            (7.0, 0.125)
            (8.5, 1.5)
            (1.0, 0.5)
            (2.0, 2.5)
            (0.0, 1.75)
            (0.0, 2.5)
            (6.0, 2.5)
            (4.0, 2.5)
            (0.0, 0.5)
            (0.0, 1.0)
            (8.0, 1.6666666666666665)
            (8.0, 2.5)
            (1.0, 0.0)
            (0.0, 0.0)
            (8.0, 1.0)
            (8.0, 0.25)
            (8.0, 0.5)
            (3.25, 0.0)
            (2.0, 0.0)
            (7.0, 0.0)
            (8.0, 0.0)
            (4.75, 0.0)
            (6.0, 0.0)
        ]
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 1))), collect.(getindex.(Ref(orig_pt), [22, 21, 11, 17, 22])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 2))), collect.(getindex.(Ref(orig_pt), [18, 17, 11, 7, 13, 18])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 3))), collect.(getindex.(Ref(orig_pt), [13, 7, 4, 12, 14, 13])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 4))), collect.(getindex.(Ref(orig_pt), [11, 21, 27, 26, 2, 4, 7, 11])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 5))), collect.(getindex.(Ref(orig_pt), [30, 31, 28, 9, 5, 3, 30])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 6))), collect.(getindex.(Ref(orig_pt), [28, 29, 24, 1, 9, 28])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 7))), collect.(getindex.(Ref(orig_pt), [1, 24, 25, 23, 1])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 8))), collect.(getindex.(Ref(orig_pt), [9, 1, 23, 19, 6, 5, 9])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 9))), collect.(getindex.(Ref(orig_pt), [26, 30, 3, 2, 26])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 10))), collect.(getindex.(Ref(orig_pt), [4, 2, 3, 5, 6, 15, 16, 12, 4])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 11))), collect.(getindex.(Ref(orig_pt), [19, 20, 15, 6, 19])), ≈)
        @test isempty(DT.get_unbounded_polygons(_vorn))
        @test all([0 ≤ get_area(_vorn, i) < Inf for i in each_polygon_index(_vorn)])
        @test allunique(DT.get_polygon_points(_vorn))
        for i in each_polygon_index(_vorn)
            C = get_polygon(_vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(_vorn, v), get_points(tri), get_convex_hull_indices(tri))
                @test δ ≥ 0
            end
        end
    end
end

@testset "initialise_clipping_arrays" begin
    a = (4.0, 3.0)
    b = (0.0, 3.0)
    c = (0.0, 0.0)
    d = (4.0, 0.0)
    e = (2.0, 1.5)
    pts = [a, b, c, d, e]
    tri = triangulate(pts, delete_ghosts=false, randomise=false)
    vorn = voronoi(tri)
    lock_convex_hull!(tri)
    edges_to_process,
    polygon_edge_queue,
    boundary_sites,
    segment_intersections,
    processed_pairs,
    intersected_edge_cache,
    exterior_circumcenters,
    left_edge_intersectors,
    right_edge_intersectors,
    current_edge_intersectors,
    equal_circumcenter_mapping = DT.initialise_clipping_arrays(vorn)
    @test edges_to_process == Set(((1, 2), (4, 1), (3, 4), (2, 3)))
    @test polygon_edge_queue == Queue{Tuple{NTuple{2,Int},Int}}()
    @test boundary_sites == Dict{Int,Set{Int}}()
    @test segment_intersections == NTuple{2,Int}[]
    @test processed_pairs == Set{Tuple{NTuple{2,Int},Int}}()
    @test intersected_edge_cache == Pair{NTuple{2,Int},NTuple{2,Int}}[]
    @test left_edge_intersectors == Set{NTuple{2,Int}}()
    @test right_edge_intersectors == Set{NTuple{2,Int}}()
    @test current_edge_intersectors == Set{NTuple{2,Int}}()
    @test equal_circumcenter_mapping == Dict{Int,Int}()
end

@testset "enqueue_new_edge" begin
    for _ in 1:500
        a = (3.0, 3.0)
        b = (0.0, 3.0)
        c = (0.0, 0.0)
        d = (4.0, 0.0)
        e = (1.0, 1.5)
        pts = [a, b, c, d, e]
        tri = triangulate(pts, delete_ghosts=false)
        vorn = voronoi(tri)
        lock_convex_hull!(tri)
        edges_to_process, polygon_edge_queue = DT.initialise_clipping_arrays(vorn)
        e = (1, 2)
        DT.enqueue_new_edge!(polygon_edge_queue, vorn, e)
        _e, _polygon = dequeue!(polygon_edge_queue)
        @test _e == e && _polygon == 2
        e = (2, 3)
        DT.enqueue_new_edge!(polygon_edge_queue, vorn, e)
        _e, _polygon = dequeue!(polygon_edge_queue)
        @test _e == e && _polygon == 5
        e = (3, 4)
        DT.enqueue_new_edge!(polygon_edge_queue, vorn, e)
        _e, _polygon = dequeue!(polygon_edge_queue)
        @test _e == e && _polygon == 5
        e = (4, 1)
        DT.enqueue_new_edge!(polygon_edge_queue, vorn, e)
        _e, _polygon = dequeue!(polygon_edge_queue)
        @test _e == e && _polygon == 4
    end
end

@testset "Segment classification" begin
    @test DT.is_segment_between_two_ghosts(-1, -2)
    @test !DT.is_segment_between_two_ghosts(1, 2)
    @test DT.is_ray_going_in(-1, 2)
    @test !DT.is_ray_going_in(1, 2)
    @test !DT.is_ray_going_in(1, -2)
    @test DT.is_ray_going_out(1, -2)
    @test !DT.is_ray_going_out(1, 2)
    @test !DT.is_ray_going_out(-1, 2)
    @test DT.is_finite_segment(1, 2)
    @test !DT.is_finite_segment(-1, 2)
    @test !DT.is_finite_segment(1, -2)
    @test !DT.is_finite_segment(-1, -2)
end

@testset "add_to_intersected_edge_cache" begin
    u, v, a, b = 1, 7, 5, 9
    intersected_edge_cache = Pair{NTuple{2,Int},NTuple{2,Int}}[]
    DT.add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)
    @test intersected_edge_cache == [(u, v) => (a, b)]
    u, v, a, b = -2, 5, 10, 17
    DT.add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)
    @test intersected_edge_cache == [(1, 7) => (5, 9), (u, v) => (a, b)]
end

@testset "More detailed test with a tessellation that has a finite edge going completely through the interior" begin
    a = (3.0, 3.0)
    b = (0.0, 3.0)
    c = (0.0, 0.0)
    d = (4.0, 0.0)
    e = (1.0, 1.5)
    pts = [a, b, c, d, e]
    tri = triangulate(pts, delete_ghosts=false, randomise=false)
    vorn = voronoi(tri)
    lock_convex_hull!(tri)
    edges_to_process,
    polygon_edge_queue,
    boundary_sites,
    segment_intersections,
    processed_pairs,
    intersected_edge_cache,
    exterior_circumcenters,
    left_edge_intersectors,
    right_edge_intersectors,
    current_edge_intersectors,
    equal_circumcenter_mapping = DT.initialise_clipping_arrays(vorn)
    e = DT.convert_to_boundary_edge(vorn, first(edges_to_process))
    DT.enqueue_new_edge!(polygon_edge_queue, vorn, e)

    empty!(intersected_edge_cache)
    e, incident_polygon = dequeue!(polygon_edge_queue)
    push!(processed_pairs, (e, incident_polygon))
    left_edge, right_edge = DT.get_neighbouring_boundary_edges(vorn, e)
    @test left_edge == (1, 4)
    @test right_edge == (3, 2)
    polygon_vertices = get_polygon(vorn, incident_polygon)
    nedges = num_boundary_edges(polygon_vertices)

    ℓ = 1
    u = get_boundary_nodes(polygon_vertices, ℓ)
    v = get_boundary_nodes(polygon_vertices, ℓ + 1)
    DT.get_circumcenter_to_triangle.(Ref(vorn), (u, v))
    @test DT.is_ray_going_out(u, v)
    @test !any(isnan, DT.process_ray_intersection!(vorn, v, u, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping))
    @test intersected_edge_cache == [(v, u) => e]
    @test segment_intersections == [(1.5, 3.0)]
    @test boundary_sites == Dict(incident_polygon => Set(1))

    ℓ = 2
    u = get_boundary_nodes(polygon_vertices, ℓ)
    v = get_boundary_nodes(polygon_vertices, ℓ + 1)
    DT.get_circumcenter_to_triangle.(Ref(vorn), (u, v))
    @test DT.is_segment_between_two_ghosts(u, v)

    ℓ = 3
    u = get_boundary_nodes(polygon_vertices, ℓ)
    v = get_boundary_nodes(polygon_vertices, ℓ + 1)
    DT.get_circumcenter_to_triangle.(Ref(vorn), (u, v))
    @test DT.is_ray_going_in(u, v)
    @test all(isnan, DT.process_ray_intersection!(vorn, u, v, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping))

    ℓ = 4
    u = get_boundary_nodes(polygon_vertices, ℓ)
    v = get_boundary_nodes(polygon_vertices, ℓ + 1)
    DT.get_circumcenter_to_triangle.(Ref(vorn), (u, v))
    @test DT.is_finite_segment(u, v)
    @test all(isnan, DT.process_segment_intersection!(vorn, u, v, e, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping))
    @test all(isnan, DT.process_segment_intersection!(vorn, u, v, left_edge, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping))
    @test !any(isnan, DT.process_segment_intersection!(vorn, u, v, right_edge, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping))
    @test collect.(segment_intersections) ≈ collect.([(1.5, 3.0), (0.0, 1.916666666666666666666666)])
    @test boundary_sites == Dict(incident_polygon => Set((1, 2)))
    @test intersected_edge_cache == [(-3, 3) => e, (u, v) => right_edge]

    DT.classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, e)
    @test left_edge_intersectors == Set{NTuple{2,Int}}()
    @test right_edge_intersectors == Set([(u, v)])
    @test current_edge_intersectors == Set([(-3, 3)])

    DT.process_intersection_points!(polygon_edge_queue, vorn, incident_polygon,
        left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
        left_edge, right_edge, e, processed_pairs, segment_intersections, boundary_sites)

    _queue = Queue{Tuple{NTuple{2,Int},Int}}()
    enqueue!(_queue, ((3, 2), 3))
    enqueue!(_queue, ((3, 2), 2))
    enqueue!(_queue, ((3, 2), 5))
    enqueue!(_queue, ((2, 1), 1))
    @test _queue == polygon_edge_queue

    a = (3.0, 3.0)
    b = (0.0, 3.0)
    c = (0.0, 0.0)
    d = (4.0, 0.0)
    e = (1.0, 1.5)
    pts = [a, b, c, d, e]
    tri = triangulate(pts, delete_ghosts=false, randomise=false)
    vorn = voronoi(tri)
    lock_convex_hull!(tri)
    boundary_sites, segment_intersections, exterior_circumcenters = DT.find_all_intersections(vorn)
    @test collect.(segment_intersections) ≈ collect.([
        (1.5, 3.0)
        (0.0, 1.9166666666666665)
        (0.0, 3.0)
        (0.0, 1.0833333333333333)
        (1.625, 0.0)
        (0.0, 0.0)
        (2.125, 0.0)
        (3.5, 1.5)
        (3.0, 3.0)
        (4.0, 0.0)
    ])
    @test boundary_sites[4] == Set((7, 10, 8))
    @test boundary_sites[2] == Set((2, 3, 1))
    @test boundary_sites[1] == Set((9, 8, 1))
    @test boundary_sites[3] == Set((5, 4, 6))
    @test boundary_sites[5] == Set((5, 4, 7, 2))
    @test exterior_circumcenters == Set((2, 1))

    DT.clip_voronoi_tessellation!(vorn)
    @test isempty(vorn.unbounded_polygons)
    @test validate_tessellation(vorn)
    @test DT.points_are_unique(vorn.polygon_points)
    @test collect.(DT.get_polygon_points(vorn)) ≈ collect.([
        (2.0, -0.25)
        (-0.625, 1.5)
        (1.5, 2.9166666666666665)
        (2.75, 1.25)
        (1.5, 3.0)
        (0.0, 1.9166666666666665)
        (0.0, 3.0)
        (0.0, 1.0833333333333333)
        (1.625, 0.0)
        (0.0, 0.0)
        (2.125, 0.0)
        (3.5, 1.5)
        (3.0, 3.0)
        (4.0, 0.0)
    ])
    @test vorn.polygons == Dict(
        5 => [8, 9, 11, 4, 3, 6, 8],
        4 => [11, 14, 12, 4, 11],
        2 => [6, 3, 5, 7, 6],
        3 => [10, 9, 8, 10],
        1 => [4, 12, 13, 5, 3, 4]
    )
    @test DT.get_boundary_polygons(vorn) == Set((4, 2, 1, 3, 5))
end

@testset "Single triangle" begin
    for _ in 1:1000
        points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        tri = triangulate(points)
        vorn = voronoi(tri)
        @test validate_tessellation(vorn)
        vorn = voronoi(tri, true)
        @test validate_tessellation(vorn)
        @test vorn.polygon_points == [
            (0.5, 0.5),
            (0.5, 0.5),
            (0.5, 0.0),
            (0.0, 0.5),
            (1.0, 0.0),
            (0.0, 1.0),
            (0.0, 0.0)
        ]
    end

    for _ in 1:100
        points = [
            0.290978 0.830755 0.0139574
            0.386411 0.630008 0.803881
        ]
        tri = triangulate(points, delete_ghosts=false)
        vorn = voronoi(tri, true)
        @test validate_tessellation(vorn)
        @test collect.(sort(vorn.polygon_points)) ≈ collect.(sort([
            (0.43655799581398663, 0.7836598194374879)
            (0.5608665, 0.5082095)
            (0.4713751999728193, 0.7065097477648392)
            (0.1524677, 0.595146)
            (0.35698797851173, 0.7308595369379512)
            (0.830755, 0.630008)
            (0.0139574, 0.803881)
            (0.290978, 0.386411)
        ]))
    end

    for _ in 1:1000
        pts = rand(2, 3)
        tri = triangulate(pts)
        vorn = voronoi(tri, true)
        @test validate_tessellation(vorn)
    end
end

@testset "Another previously broken example with a non-boundary generator's intersecting edges not being previously detected" begin
    for _ in 1:100
        a = (0.0, 0.0)
        b = (6.0, 0.0)
        c = (6.0, 6.0)
        d = (0.0, 6.0)
        e = (5.0, 3.0)
        f = (0.2, 5.8)
        g = (0.1, 0.2)
        h = (0.2, 0.1)
        pts = [a, b, c, d, e, f, g, h]
        tri = triangulate(pts)
        vorn = voronoi(tri)
        _vorn = voronoi(tri, true)
        @test validate_tessellation(vorn)
        @test validate_tessellation(_vorn)
        orig_pt = [
            (3.1112716763005785, 0.703757225433526)
            (3.000000000000019, -5.750000000000036)
            (3.12093023255814, 5.293023255813954)
            (2.999999999999998, 8.799999999999999)
            (2.2045454545454386, 2.2045454545454577)
            (-2.7482456140350915, 3.0517543859649106)
            (0.08333333333333333, 0.08333333333333333)
            (10.0, 3.0)
            (1.7664948453608247, 2.9711340206185564)
            (-5.750000000000036, 3.000000000000019)
            (0.12499999999999911, 0.0)
            (3.0991379310344853, 0.0)
            (0.0, 0.125)
            (0.0, 0.0)
            (6.0, 1.666666666666666)
            (6.0, 0.0)
            (0.0, 5.800000000000001)
            (0.19999999999999837, 6.0)
            (0.0, 6.0)
            (0.0, 3.0026785714285706)
            (3.096551724137931, 6.0)
            (6.0, 4.333333333333334)
            (6.0, 6.0)
        ]
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 1))), collect.(getindex.(Ref(orig_pt), [14, 11, 7, 13, 14])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 2))), collect.(getindex.(Ref(orig_pt), [12, 16, 15, 1, 12])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 3))), collect.(getindex.(Ref(orig_pt), [3, 22, 23, 21, 3])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 4))), collect.(getindex.(Ref(orig_pt), [17, 18, 19, 17])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 5))), collect.(getindex.(Ref(orig_pt), [5, 1, 15, 22, 3, 9, 5])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 6))), collect.(getindex.(Ref(orig_pt), [20, 9, 3, 21, 18, 17, 20])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 7))), collect.(getindex.(Ref(orig_pt), [13, 7, 5, 9, 20, 13])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(_vorn), get_polygon(_vorn, 8))), collect.(getindex.(Ref(orig_pt), [7, 11, 12, 1, 5, 7])), ≈)
        @test isempty(DT.get_unbounded_polygons(_vorn))
        @test all([0 ≤ get_area(_vorn, i) < Inf for i in each_polygon_index(_vorn)])
        @test allunique(DT.get_polygon_points(_vorn))
        for i in each_polygon_index(_vorn)
            C = get_polygon(_vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(_vorn, v), get_points(tri), get_convex_hull_indices(tri))
                @test δ ≥ -1e-14
            end
        end
        @test DT.get_boundary_polygons(_vorn) == Set((4, 6, 3, 5, 2, 8, 1, 7))
    end
end

@testset "Varying size" begin
    for n in 3:200
        for j in 1:250
            rng = StableRNG(n + 4j)
            pts = rand(rng, 2, n)
            tri = triangulate(pts; rng)
            vorn = voronoi(tri)
            flag1 = validate_tessellation(vorn)
            vorn = voronoi(tri, true)
            flag2 = validate_tessellation(vorn)
            @test flag1
            @test flag2
            (!flag1 || !flag2) && @show n, j
        end
    end
end

@testset "Centroidal tessellation" begin
    flag = 0
    tot = 0
    for i in 1:50
        p1 = randn(2, 50)
        p2 = rand(SVector{2,Float64}, 30)
        p3 = rand(Point2f, 250)
        p4 = randn(Float32, 2, 70)
        p5 = randn(2, 4)
        p6 = rand(SVector{2,Float64}, 7)
        p7 = rand(Point2f, 5)
        p8 = randn(Float32, 2, 15)
        _pts = (p1, p2, p3, p4, p5, p6, p7, p8)
        for jj in eachindex(_pts)
            points = _pts[jj]
            println("Starting centroidal test at $((i, jj)).")
            tri = triangulate(points)
            vorn = voronoi(tri, true)
            @test validate_tessellation(vorn, check_convex=!(jj ∈ (3, 4, 7, 8)))
            smooth_vorn = centroidal_smooth(vorn, maxiters=5000)
            @test validate_tessellation(smooth_vorn, check_convex=!(jj ∈ (3, 4, 7, 8)))
            for i in each_polygon_index(smooth_vorn)
                p = get_generator(smooth_vorn, i)
                c = DT.get_centroid(smooth_vorn, i)
                px, py = getxy(p)
                cx, cy = getxy(c)
                dx, dy = px - cx, py - cy
                ℓ = sqrt(dx^2 + dy^2)
                _flag = ℓ ≤ 1e-1
                flag += _flag
                tot += 1
            end
        end
    end
    @test flag / tot > 0.9
end

@testset "Lattice" begin
    for _ in 1:100
        tri = triangulate_rectangle(0, 1, 0, 1, 11, 11, add_ghost_triangles=true)
        tri = triangulate(tri.points)
        vorn = voronoi(tri)
        @test validate_tessellation(vorn; check_convex=false)
        vorn = voronoi(tri, true)
        @test validate_tessellation(vorn; check_convex=false, check_adjacent=false)
    end
end

@testset "Example that used to previously break: Plotting a Voronoi tile with unbounded edges intersecting over non-neighbouring sides" begin
    pts = [0.508812 0.656662 0.785124 0.63427 0.444969 0.609253 0.0826304 0.265388 0.830807 0.658346
        0.647732 0.482994 0.809909 0.482046 0.0170022 0.821742 0.835057 0.591724 0.881006 0.97652]
    tri = triangulate(pts)
    vorn = voronoi(tri)
    clip = DT.get_polygon_coordinates(vorn, 9, DT.polygon_bounds(vorn, 0.1))
    @test DT.circular_equality(collect.(clip), collect.([
        (0.8374509236290323, 1.0964579554357115)
        (0.7271857522962732, 0.8973620981454822)
        (1.7423743304675243, 0.24505806539308744)
        (2.0053766736553027, 0.1299847935961671)
        (2.0053766736553027, 1.0964579554357115)
        (0.8374509236290323, 1.0964579554357115)
    ]), ≈)
end

@testset "Example that used to previously break: Plotting a Voronoi tile with unbounded edges intersecting over non-neighbouring sides, needing THREE corners" begin
    pts = [0.279402 0.874842 0.163028
        0.274178 0.831658 0.223709]
    tri = triangulate(pts)
    vorn = voronoi(tri)
    clip = DT.get_polygon_coordinates(vorn, 2, DT.polygon_bounds(vorn, 2.0))
    @test DT.circular_equality(collect.(clip), collect.([
        (-2.551580488601045, 4.122780970352073)
        (-0.3314903102646037, 1.5233996567840244)
        (3.2875066205292076, -2.3420224793856605)
        (3.2875066205292076, 4.122780970352073)
        (-2.551580488601045, 4.122780970352073)
    ]), ≈)
end

@testset "Issue #72" begin
    points = [
        Float32[0.32965052, 0.7966664],
        Float32[0.015137732, 0.31555605],
        Float32[0.54775107, 0.7222514],
        Float32[0.687552, 0.6982844],
        Float32[0.65762305, 0.5177773],
        Float32[0.9649102, 0.8300816],
        Float32[0.12174326, 0.82220316],
        Float32[0.007668495, 0.7747718],
        Float32[0.3144052, 0.5493178],
        Float32[0.6848137, 0.12493092],
        Float32[0.39197737, 0.6912688],
        Float32[0.41400427, 0.025964081],
        Float32[0.35919905, 0.7255166],
        Float32[0.80712754, 0.3415957],
        Float32[0.542679, 0.51094216],
        Float32[0.092720866, 0.90151125],
        Float32[0.90992355, 0.8814645],
        Float32[0.02194357, 0.00064593554],
        Float32[0.9616154, 0.10633117],
        Float32[0.0044495463, 0.97074896],
        Float32[0.4309939, 0.5323847],
        Float32[0.90867966, 0.55974066],
        Float32[0.580766, 0.7668439],
        Float32[0.8563475, 0.88143903],
        Float32[0.18311942, 0.8367877]
    ]
    tri = triangulate(points)
    vorn = voronoi(tri)
    @test validate_tessellation(vorn)

    points = [
        [0.01877582f0, 0.33188105f0],
        [0.57645035f0, 0.58131695f0],
        [0.14916456f0, 0.37567925f0],
        [0.87054604f0, 0.29108626f0],
        [0.15384161f0, 0.80313444f0],
        [0.66470474f0, 0.2547925f0],
        [0.69431657f0, 0.42567456f0],
        [0.43695337f0, 0.42649788f0],
        [0.3316936f0, 0.18936294f0],
        [0.98043495f0, 0.8360868f0],
        [0.5788496f0, 0.103449225f0],
        [0.5252029f0, 0.96790665f0],
        [0.33206534f0, 0.88216203f0],
        [0.07115775f0, 0.5983915f0],
        [0.29895544f0, 0.103566706f0],
        [0.52547264f0, 0.57929194f0],
        [0.19257814f0, 0.30570602f0],
        [0.12954468f0, 0.11141288f0],
        [0.28790158f0, 0.39447558f0],
        [0.6525599f0, 0.6425986f0]
    ]
    tri = triangulate(points)
    vorn = voronoi(tri)
    @test validate_tessellation(vorn)
end

@testset "Polygon bounds with generators only" begin
    points = rand(2, 50)
    tri = triangulate(points)
    vorn = voronoi(tri)
    xmin, xmax, ymin, ymax = DT.polygon_bounds(vorn, 0.1; include_polygon_vertices=false)
    _xmin, _xmax = extrema(points[1, :])
    _ymin, _ymax = extrema(points[2, :])
    @test xmin == _xmin - 0.1(_xmax - _xmin)
    @test xmax == _xmax + 0.1(_xmax - _xmin)
    @test ymin == _ymin - 0.1(_ymax - _ymin)
    @test ymax == _ymax + 0.1(_ymax - _ymin)
end

@testset "Position of Voronoi polygons relative to box" begin
    a = (-3.0, 7.0)
    b = (1.0, 6.0)
    c = (-1.0, 3.0)
    d = (-2.0, 4.0)
    e = (3.0, -2.0)
    f = (5.0, 5.0)
    g = (-4.0, -3.0)
    h = (3.0, 8.0)
    points = [a, b, c, d, e, f, g, h]
    tri = triangulate(points)
    vorn = voronoi(tri)
    a, b, c, d = -8.0, 6.0, -2.0, 10.0
    bounding_box = (a, b, c, d)
    results = Dict(
        1 => DT.has_multiple_intersections,
        2 => DT.is_inside,
        3 => DT.is_inside,
        4 => DT.has_multiple_intersections,
        5 => DT.has_multiple_intersections,
        6 => DT.has_multiple_intersections,
        7 => DT.has_multiple_intersections,
        8 => DT.has_multiple_intersections,
    )
    for (i, cert_f) in results
        @test cert_f(DT.polygon_position_relative_to_box(vorn, bounding_box, i))
    end
end

A = (-3.0, 7.0)
B = (1.0, 6.0)
C = (-1.0, 3.0)
D = (-2.0, 4.0)
E = (3.0, -2.0)
F = (5.0, 5.0)
G = (-4.0, -3.0)
H = (3.0, 8.0)
points = [A, B, C, D, E, F, G, H]
tri = triangulate(points)
vorn = voronoi(tri)
a, b, c, d = -4.0, 6.0, -2.0, 4.0
bounding_box = (a, b, c, d)
pa = NTuple{2,Float64}[]
pb = [(0.75, 4.0), (2.357142857142857, 2.9285714285714284), (2.625, 4.0), (0.75, 4.0)]
pc = [(2.710526315789474, 1.868421052631579), (2.357142857142857, 2.9285714285714284), (0.75, 4.0), (-1.0, 4.0), (-4.0, 1.0), (-4.0, 0.75), (-0.7307692307692308, -0.8846153846153846), (2.710526315789474, 1.868421052631579)]
pd = [(-4.0, 4.0), (-4.0, 1.0), (-1.0, 4.0), (-4.0, 4.0)]
collect.(pd) ≈ collect.(pa)
pe = get_polygon_coordinates(vorn, 7; bounding_box)

fig, ax, sc = voronoiplot(vorn)
lines!(ax, [(a, c), (b, c), (b, d), (a, d), (a, c)], color=:red, linewidth=6)
lines!(ax, pd, color=:white, linewidth=4)
fig

i = 7
poly, clip_poly = DT.get_clipping_poly_structs(vorn, i, bounding_box)
vertices, clip_vertices, clip_points = poly.vertices, clip_poly.vertices, clip_poly.points
new_vertices, new_point_list = DT.get_new_polygon_indices(vorn, vertices)
output_vertices = new_vertices
output_points = deepcopy(new_point_list)
q = clip_points[1]
p = clip_points[2]
input_vertices = output_vertices
T = typeof(q)
I = eltype(input_vertices)
output_vertices = I[]
output_points = T[]
s_vertex = input_vertices[end]
vertex = input_vertices[1]
DT._clip_unbounded_polygon_edge_to_ray_going_in!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[2]
DT._clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[3]
DT._clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[4]
DT._clip_unbounded_polygon_edge_to_ray_going_out!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)


i = 7
poly, clip_poly = DT.get_clipping_poly_structs(vorn, i, bounding_box)
vertices, clip_vertices, clip_points = poly.vertices, clip_poly.vertices, clip_poly.points
new_vertices, new_point_list = DT.get_new_polygon_indices(vorn, vertices)
output_vertices = new_vertices
output_points = deepcopy(new_point_list)
q = clip_points[end] 
p = clip_points[1] 
input_vertices = output_vertices
T = typeof(q)
I = eltype(input_vertices)
output_vertices = I[]
output_points = T[]
s_vertex = input_vertices[end]
vertex = input_vertices[1]
DT._clip_unbounded_polygon_edge_to_ray_going_in!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[2]
DT._clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[3]
DT._clip_unbounded_polygon_edge_to_finite_segment!(s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
s_vertex = vertex
vertex = input_vertices[4]
DT._clip_unbounded_polygon_edge_to_ray_going_out!(vorn, i, s_vertex, vertex, new_point_list, q, p, output_vertices, output_points)
