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

@testset "Unconstrained test" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:10
            A = (-1.0, 7.0)
            B = (4.0, 4.0)
            C = (-2.0, -1.0)
            D = (-1.0, 3.0)
            E = (3.0, -1.0)
            F = (1.0, 4.0)
            G = (-3.0, 5.0)
            pts = [A, B, C, D, E, F, G]
            tri = triangulate(pts; delete_ghosts=false, randomise=false)

            vorn = voronoi(tri; predicates=PT())
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
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (7, 4, 1), (4, 6, 1), (6, 2, 1), (1, 2, -1), (7, 1, -1), (7, 4, 1),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 2),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (2, 5, -1), (1, 2, -1), (6, 2, 1), (6, 5, 2), (2, 5, -1),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 3),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (5, 3, -1), (5, 4, 3), (4, 7, 3), (3, 7, -1), (5, 3, -1),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 4),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (5, 6, 4), (4, 6, 1), (7, 4, 1), (4, 7, 3), (5, 4, 3), (5, 6, 4),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 5),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (6, 5, 2), (5, 6, 4), (5, 4, 3), (5, 3, -1), (2, 5, -1), (6, 5, 2),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 6),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (6, 5, 2), (6, 2, 1), (4, 6, 1), (5, 6, 4), (6, 5, 2),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 7),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (7, 4, 1), (7, 1, -1), (3, 7, -1), (4, 7, 3), (7, 4, 1),
                    ],
                ),
            )
        end
    end
end

@testset "Smaller example, checking ray coordinates" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:10
            tri = example_triangulation()
            tri = triangulate(get_points(tri))
            vorn = voronoi(tri; predicates=PT())
            @test validate_tessellation(vorn; predicates=PT())
            bbox = DT.polygon_bounds(vorn, 0.1)
            xmin, xmax, ymin, ymax = bbox
            c1 = DT.get_polygon_coordinates(vorn, 1, bbox, predicates=PT())
            c2 = DT.get_polygon_coordinates(vorn, 2, bbox, predicates=PT())
            c3 = DT.get_polygon_coordinates(vorn, 3, bbox, predicates=PT())
            c4 = DT.get_polygon_coordinates(vorn, 4, bbox, predicates=PT())
            c5 = DT.get_polygon_coordinates(vorn, 5, bbox, predicates=PT())
            c6 = DT.get_polygon_coordinates(vorn, 6, bbox, predicates=PT())
            c7 = DT.get_polygon_coordinates(vorn, 7, bbox, predicates=PT())
            @test all(DT.is_circular, (c1, c2, c3, c4, c5, c6, c7))
            @test DT.circular_equality(
                collect.(c1), collect.(
                    [
                        (-1.5, 0.5)
                        (0.166666666666, -1.1666666666666665)
                        (1.0, 0.5)
                        (1.0, 3.0)
                        (-1.5, 0.5)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c2), collect.(
                    [
                        (0.5, -3.2)
                        (5.7, -3.2)
                        (5.7, -1.700000000000001)
                        (3.5, 0.5)
                        (0.5, -2.5)
                        (0.5, -3.2)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c3), collect.(
                    [
                        (3.5, 0.5)
                        (1.0, 0.5)
                        (0.16666666666666666, -1.1666666666666665)
                        (0.5, -2.5)
                        (3.5, 0.5)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c4), collect.(
                    [
                        (1.5, 5.2)
                        (-2.7, 5.2)
                        (-2.7, 0.9000000000000001)
                        (-1.5, 0.5)
                        (1.0, 3.0)
                        (1.5, 4.5)
                        (1.5, 5.2)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c5), collect.(
                    [
                        (1.5, 4.5)
                        (3.5, 0.5)
                        (5.7, 2.6999999999999997)
                        (5.7, 5.2)
                        (1.5, 5.2)
                        (1.5, 4.5)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c6), collect.(
                    [
                        (-2.7, 0.9000000000000001)
                        (-2.7, -3.2)
                        (0.5, -3.2)
                        (0.5, -2.5)
                        (0.16666666666666666, -1.1666666666666665)
                        (-1.5, 0.5)
                        (-2.7, 0.9000000000000001)
                    ],
                ), ≈,
            )
            @test DT.circular_equality(
                collect.(c7), collect.(
                    [
                        (1.5, 4.5)
                        (1.0, 3.0)
                        (1.0, 0.5)
                        (3.5, 0.5)
                        (1.5, 4.5)
                    ],
                ), ≈,
            )
        end
    end
end

@testset "delete/add_polygon_adjacencies" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:10
            A = (-1.0, 7.0)
            B = (4.0, 4.0)
            C = (-2.0, -1.0)
            D = (-1.0, 3.0)
            E = (3.0, -1.0)
            F = (1.0, 4.0)
            G = (-3.0, 5.0)
            pts = [A, B, C, D, E, F, G]
            tri = triangulate(pts; delete_ghosts=false, randomise=false)
            vorn = voronoi(tri; predicates=PT())
            @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
                  get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
                  5
            DT.delete_polygon_adjacent!(vorn, 5)
            @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
                  get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
                  DT.∅
            DT.add_polygon_adjacent!(vorn, 5)
            @test get_adjacent(vorn, 1, -2) == get_adjacent(vorn, -2, -1) ==
                  get_adjacent(vorn, -1, 7) == get_adjacent(vorn, 7, 3) == get_adjacent(vorn, 3, 1) ==
                  5
        end
    end
end

@testset "Voronoi point location" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        A = (-1.0, 7.0) .+ (1.0e-6rand(), 1.0e-6rand()) # perturb to allow the tests to work even without ExactPredicates
        B = (4.0, 4.0) .+ (1.0e-6rand(), 1.0e-6rand())
        C = (-2.0, -1.0) .+ (1.0e-6rand(), 1.0e-6rand())
        D = (-1.0, 3.0) .+ (1.0e-6rand(), 1.0e-6rand())
        E = (3.0, -1.0) .+ (1.0e-6rand(), 1.0e-6rand())
        F = (1.0, 4.0) .+ (1.0e-6rand(), 1.0e-6rand())
        G = (-3.0, 5.0) .+ (1.0e-6rand(), 1.0e-6rand())
        pts = [A, B, C, D, E, F, G]
        tri = triangulate(pts; delete_ghosts=false, randomise=false)
        vor = voronoi(tri, predicates=PT())
        @test validate_tessellation(vor, predicates=PT())
        xmin, xmax, ymin, ymax = DT.polygon_bounds(get_points(tri), get_convex_hull_vertices(tri))
        p = NTuple{2,Float64}[]
        n = 10000
        while length(p) ≤ n # only going to test points that are inside the polygon
            pt = (xmin + rand() * (xmax - xmin), ymin + rand() * (ymax - ymin))
            if DT.distance_to_polygon(pt, get_points(tri), get_convex_hull_vertices(tri)) ≥ 0
                push!(p, pt)
            end
        end
        for p in p
            u = get_nearest_neighbour(vor, p; predicates=PT())
            all_dists = [norm(p .- get_generator(vor, i)) for i in sort(collect(each_generator(vor)))]
            @test findmin(all_dists)[2] == u
        end

        if PT() != DT.FastKernel()
            points = [
                (0.0, 0.0), (-1.0, 1.0), (-0.5, 1.0), (0.0, 1.0), (0.5, 1.0), (1.0, 1.0),
                (1.0, 0.8), (1.0, 0.0), (1.0, -0.5), (1.0, -1.0),
                (0.1, -1.0), (-0.8, -1.0), (-1.0, -1.0),
                (-1.0, -0.7), (-1.0, -0.1), (-1.0, 0.6),
                (-0.1, -0.8), (0.2, -0.8),
                (-0.6, -0.4), (0.9, 0.0), (-0.5, 0.5), (-0.4, 0.6), (-0.1, 0.8),
            ]
            tri = triangulate(points, delete_ghosts=false)
            vorn = voronoi(tri, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            xg = LinRange(-1, 1, 50)
            yg = LinRange(-1, 1, 50)
            x = vec([x for x in xg, _ in yg])
            y = vec([y for _ in xg, y in yg])
            for (ξ, η) in zip(x, y)
                p = (ξ, η)
                u = get_nearest_neighbour(vorn, p; predicates=PT())
                @inferred get_nearest_neighbour(vorn, p; predicates=PT())
                all_dists = [norm(p .- get_generator(vorn, i)) for i in sort(collect(each_generator(vorn)))]
                k = findmin(all_dists)[2]
                @test k == u
                for m in DT.each_point_index(tri)
                    u = get_nearest_neighbour(vorn, p, try_points=m; predicates=PT())
                    @test u == k
                end
            end
        end
    end
end

@testset "Clipping a simple VoronoiTessellation" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:333
            A = (-1.0, 7.0)
            B = (4.0, 4.0)
            C = (-2.0, -1.0)
            D = (-1.0, 3.0)
            E = (3.0, -1.0)
            F = (1.0, 4.0)
            G = (-3.0, 5.0)
            pts = [A, B, C, D, E, F, G]
            tri = triangulate(pts; delete_ghosts=false, randomise=true)
            #lock_convex_hull!(tri)
            vorn = voronoi(tri; clip=true, predicates=PT())
            for (i, p) in DT.get_generators(vorn)
                @test get_point(tri, i) == get_generator(vorn, i) == p
            end
            @test DT.get_triangulation(vorn) == tri
            circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
            triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
            for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
                V = DT.sort_triangle(V)
                c = DT.get_triangle_to_circumcenter(vorn, V)
                c = get_polygon_point(vorn, c)
                i, j, k = triangle_vertices(V)
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
            orig_pt = collect.(
                [
                    (0.5, 0.5)
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
                    (-1.0, 7.0)
                ],
            )
            @test validate_tessellation(vorn, predicates=PT())
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
            for i in DT.each_polygon_index(vorn)
                @test DT.has_polygon(vorn, i)
                C = get_polygon(vorn, i)
                for (j, v) in pairs(C)
                    δ = DT.distance_to_polygon(get_polygon_point(vorn, v), get_points(tri), get_convex_hull_vertices(tri))
                    @test δ ≥ -1.0e-15
                end
            end
            @test !DT.has_polygon(vorn, 200)
        end
    end
end

@testset "Another example" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:100
            tri = fixed_shewchuk_example_constrained()
            vorn = voronoi(tri; clip=false, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            for (i, p) in DT.get_generators(vorn)
                @test get_point(tri, i) == get_generator(vorn, i) == p
            end
            @test DT.get_triangulation(vorn) == tri
            @test DT.get_unbounded_polygons(vorn) == Set((3, 10, 11, 7, 6, 5, 4, 1, 2))
            circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
            triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
            for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
                V = DT.sort_triangle(V)
                c = DT.get_triangle_to_circumcenter(vorn, V)
                c = get_polygon_point(vorn, c)
                i, j, k = triangle_vertices(V)
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
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (4, 2, 1), (1, 2, -1), (4, 1, -1), (4, 2, 1),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 2),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (4, 2, 1), (4, 3, 2), (2, 3, -1), (1, 2, -1), (4, 2, 1),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 3),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (4, 3, 2), (4, 10, 3), (3, 10, -1), (2, 3, -1), (4, 3, 2),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 4),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (5, 9, 4), (9, 10, 4), (4, 10, 3), (4, 3, 2), (4, 2, 1), (4, 1, -1), (5, 4, -1), (5, 9, 4),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 5),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (6, 8, 5), (8, 10, 5), (10, 9, 5), (5, 9, 4), (5, 4, -1), (6, 5, -1), (6, 8, 5),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 6),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (7, 8, 6), (6, 8, 5), (6, 5, -1), (7, 6, -1), (7, 8, 6),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 7),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (11, 8, 7), (7, 8, 6), (7, 6, -1), (11, 7, -1), (11, 8, 7),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 8),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (7, 8, 6), (11, 8, 7), (11, 10, 8), (8, 10, 5), (6, 8, 5), (7, 8, 6),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 9),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (10, 9, 5), (9, 10, 4), (5, 9, 4), (10, 9, 5),
                    ],
                ),
            )
            @test DT.circular_equality(
                get_polygon(vorn, 10),
                DT.get_triangle_to_circumcenter.(
                    Ref(vorn), [
                        (9, 10, 4), (10, 9, 5), (8, 10, 5), (11, 10, 8), (10, 11, -1), (3, 10, -1), (4, 10, 3), (9, 10, 4),
                    ],
                ),
            )

            _vorn = voronoi(tri; clip=true, predicates=PT())
            @test validate_tessellation(_vorn, predicates=PT())
            for (i, p) in DT.get_generators(_vorn)
                @test get_point(tri, i) == get_generator(_vorn, i) == p
            end
            @test DT.get_triangulation(_vorn) == tri
            circumcenter_to_triangle = DT.get_circumcenter_to_triangle(_vorn)
            triangle_to_circumcenter = DT.get_triangle_to_circumcenter(_vorn)
            for V in DT.each_solid_triangle(DT.get_triangulation(_vorn))
                V = DT.sort_triangle(V)
                c = DT.get_triangle_to_circumcenter(_vorn, V)
                c = get_polygon_point(_vorn, c)
                i, j, k = triangle_vertices(V)
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
                    δ = DT.distance_to_polygon(get_polygon_point(_vorn, v), get_points(tri), get_convex_hull_vertices(tri))
                    @test δ ≥ 0
                end
            end
        end
    end
end

@testset "initialise_clipping_arrays" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        a = (4.0, 3.0)
        b = (0.0, 3.0)
        c = (0.0, 0.0)
        d = (4.0, 0.0)
        e = (2.0, 1.5)
        pts = [a, b, c, d, e]
        tri = triangulate(pts, delete_ghosts=false, randomise=false)
        vorn = voronoi(tri, predicates=PT())
        lock_convex_hull!(tri; predicates=PT())
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
        @test polygon_edge_queue == DT.Queue{Tuple{NTuple{2,Int},Int}}()
        @test boundary_sites == Dict{Int,Set{Int}}()
        @test segment_intersections == NTuple{2,Int}[]
        @test processed_pairs == Set{Tuple{NTuple{2,Int},Int}}()
        @test intersected_edge_cache == Pair{NTuple{2,Int},NTuple{2,Int}}[]
        @test left_edge_intersectors == Set{NTuple{2,Int}}()
        @test right_edge_intersectors == Set{NTuple{2,Int}}()
        @test current_edge_intersectors == Set{NTuple{2,Int}}()
        @test equal_circumcenter_mapping == Dict{Int,Int}()
    end
end

@testset "enqueue_new_edge" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:180
            a = (3.0, 3.0)
            b = (0.0, 3.0)
            c = (0.0, 0.0)
            d = (4.0, 0.0)
            e = (1.0, 1.5)
            pts = [a, b, c, d, e]
            tri = triangulate(pts, delete_ghosts=false)
            vorn = voronoi(tri; predicates=PT())
            lock_convex_hull!(tri; predicates=PT())
            edges_to_process, polygon_edge_queue = DT.initialise_clipping_arrays(vorn)
            e = (1, 2)
            DT.enqueue_new_edge!(polygon_edge_queue, vorn, e, Random.default_rng(), PT())
            _e, _polygon = popfirst!(polygon_edge_queue)
            @test _e == e && _polygon == 2
            e = (2, 3)
            DT.enqueue_new_edge!(polygon_edge_queue, vorn, e, Random.default_rng(), PT())
            _e, _polygon = popfirst!(polygon_edge_queue)
            @test _e == e && _polygon == 5
            e = (3, 4)
            DT.enqueue_new_edge!(polygon_edge_queue, vorn, e, Random.default_rng(), PT())
            _e, _polygon = popfirst!(polygon_edge_queue)
            @test _e == e && _polygon == 5
            e = (4, 1)
            DT.enqueue_new_edge!(polygon_edge_queue, vorn, e, Random.default_rng(), PT())
            _e, _polygon = popfirst!(polygon_edge_queue)
            @test _e == e && _polygon == 4
        end
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
    for PT in subtypes(DT.AbstractPredicateKernel)
        a = (3.0, 3.0)
        b = (0.0, 3.0)
        c = (0.0, 0.0)
        d = (4.0, 0.0)
        e = (1.0, 1.5)
        pts = [a, b, c, d, e]
        tri = triangulate(pts, delete_ghosts=false, randomise=false)
        vorn = voronoi(tri, predicates=PT())
        lock_convex_hull!(tri, predicates=PT())
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
        e = DT.convert_to_edge_adjoining_ghost_vertex(vorn, first(edges_to_process))
        DT.enqueue_new_edge!(polygon_edge_queue, vorn, e, Random.default_rng(), PT())

        empty!(intersected_edge_cache)
        e, incident_polygon = popfirst!(polygon_edge_queue)
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
        @test !any(isnan, DT.process_ray_intersection!(vorn, v, u, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, PT()))
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
        @test all(isnan, DT.process_ray_intersection!(vorn, u, v, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, PT()))

        ℓ = 4
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        DT.get_circumcenter_to_triangle.(Ref(vorn), (u, v))
        @test DT.is_finite_segment(u, v)
        @test all(isnan, DT.process_segment_intersection!(vorn, u, v, e, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, PT()))
        @test all(isnan, DT.process_segment_intersection!(vorn, u, v, left_edge, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, PT()))
        @test !any(isnan, DT.process_segment_intersection!(vorn, u, v, right_edge, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, PT()))
        @test collect.(segment_intersections) ≈ collect.([(1.5, 3.0), (0.0, 1.916666666666666666666666)])
        @test boundary_sites == Dict(incident_polygon => Set((1, 2)))
        @test intersected_edge_cache == [(-3, 3) => e, (u, v) => right_edge]

        DT.classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, e)
        @test left_edge_intersectors == Set{NTuple{2,Int}}()
        @test right_edge_intersectors == Set([(u, v)])
        @test current_edge_intersectors == Set([(-3, 3)])

        DT.process_intersection_points!(
            polygon_edge_queue, vorn, incident_polygon,
            left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
            left_edge, right_edge, e, processed_pairs, segment_intersections, boundary_sites,
        )

        _queue = DT.Queue{Tuple{NTuple{2,Int},Int}}()
        push!(_queue, ((3, 2), 3))
        push!(_queue, ((3, 2), 2))
        push!(_queue, ((3, 2), 5))
        push!(_queue, ((2, 1), 1))
        @test _queue == polygon_edge_queue

        a = (3.0, 3.0)
        b = (0.0, 3.0)
        c = (0.0, 0.0)
        d = (4.0, 0.0)
        e = (1.0, 1.5)
        pts = [a, b, c, d, e]
        tri = triangulate(pts, delete_ghosts=false, randomise=false)
        vorn = voronoi(tri; predicates=PT())
        lock_convex_hull!(tri; predicates=PT())
        boundary_sites, segment_intersections, exterior_circumcenters = DT.find_all_intersections(vorn, predicates=PT())
        @test collect.(segment_intersections) ≈ collect.(
            [
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
            ],
        )
        @test boundary_sites[4] == Set((7, 10, 8))
        @test boundary_sites[2] == Set((2, 3, 1))
        @test boundary_sites[1] == Set((9, 8, 1))
        @test boundary_sites[3] == Set((5, 4, 6))
        @test boundary_sites[5] == Set((5, 4, 7, 2))
        @test exterior_circumcenters == Set((2, 1))

        DT.clip_voronoi_tessellation!(vorn, predicates=PT())
        @test isempty(vorn.unbounded_polygons)
        @test validate_tessellation(vorn, predicates=PT())
        @test DT.points_are_unique(vorn.polygon_points)
        @test collect.(DT.get_polygon_points(vorn)) ≈ collect.(
            [
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
            ],
        )
        @test vorn.polygons == Dict(
            5 => [8, 9, 11, 4, 3, 6, 8],
            4 => [11, 14, 12, 4, 11],
            2 => [6, 3, 5, 7, 6],
            3 => [10, 9, 8, 10],
            1 => [4, 12, 13, 5, 3, 4],
        )
        @test DT.get_boundary_polygons(vorn) == Set((4, 2, 1, 3, 5))
    end
end

@testset "Single triangle" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:333
            points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
            tri = triangulate(points)
            vorn = voronoi(tri, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            vorn = voronoi(tri, clip=true, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            @test sort(unique(vorn.polygon_points)) == sort(unique([
                (0.5, 0.5),
                (0.5, 0.5),
                (0.5, 0.0),
                (0.0, 0.5),
                (1.0, 0.0),
                (0.0, 1.0),
                (0.0, 0.0),
            ]))
            @test DT.circular_equality(collect.(vorn.polygon_points[vorn.polygons[2]]), [[1.0, 0.0], [0.5, 0.5], [0.5, 0.0], [1.0, 0.0]])
            @test DT.circular_equality(collect.(vorn.polygon_points[vorn.polygons[1]]), [[0.0, 0.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5], [0.0, 0.0]])
            @test DT.circular_equality(collect.(vorn.polygon_points[vorn.polygons[3]]), [[0.0, 1.0], [0.0, 0.5], [0.5, 0.5], [0.0, 1.0]])
        end

        for _ in 1:100
            points = [
                0.290978 0.830755 0.0139574
                0.386411 0.630008 0.803881
            ]
            tri = triangulate(points, delete_ghosts=false)
            vorn = voronoi(tri, clip=true, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            @test unique(collect.(sort(vorn.polygon_points))) ≈ collect.(
                sort(
                    [
                        (0.4365579958139866, 0.7836598194374877)
                        (0.013957399999999976, 0.803881)
                        (0.15246769999999996, 0.5951460000000001)
                        (0.3569879785117299, 0.7308595369379514)
                        (0.29097799999999996, 0.386411)
                        (0.5608665, 0.5082095)
                        (0.47137519997281924, 0.7065097477648392)
                        (0.15246769999999998, 0.595146)
                        (0.8307550000000001, 0.6300080000000001)
                    ],
                ),
            )
        end

        for _ in 1:1000
            pts = rand(2, 3)
            tri = triangulate(pts)
            vorn = voronoi(tri, clip=true, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
        end
    end
end

@testset "Another previously broken example with a non-boundary generator's intersecting edges not being previously detected" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:43
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
            vorn = voronoi(tri, predicates=PT())
            _vorn = voronoi(tri, clip=true, predicates=PT())
            @test validate_tessellation(vorn, predicates=PT())
            @test validate_tessellation(_vorn, predicates=PT())
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
                    δ = DT.distance_to_polygon(get_polygon_point(_vorn, v), get_points(tri), get_convex_hull_vertices(tri))
                    @test δ ≥ -1.0e-14
                end
            end
            @test DT.get_boundary_polygons(_vorn) == Set((4, 6, 3, 5, 2, 8, 1, 7))
        end
    end
end

@testset "Varying size" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for n in 3:4:200
            for j in 1:15:250
                rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
                pts = rand(rng, 2, n)
                tri = triangulate(pts; rng)
                vorn = voronoi(tri, predicates=PT())
                flag1 = validate_tessellation(vorn, predicates=PT())
                vorn = voronoi(tri, clip=true, predicates=PT())
                flag2 = validate_tessellation(vorn, predicates=PT())
                @test flag1
                @test flag2
            end
        end
    end
end

@testset "Centroidal tessellation" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        flag = 0
        tot = 0
        for i in 1:25
            @info "Testing centroidal tessellation: Run: $i"
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
                tri = triangulate(points, predicates=PT())
                vorn = voronoi(tri, clip=true, predicates=PT())
                @test validate_tessellation(vorn, check_convex=!(jj ∈ (3, 4, 7, 8)), predicates=PT())
                for smooth_vorn in (centroidal_smooth(vorn, maxiters=5000, predicates=PT()), voronoi(tri, clip=true, smooth=true, maxiters=5000, predicates=PT()))
                    @test validate_tessellation(smooth_vorn, check_convex=!(jj ∈ (3, 4, 7, 8)), predicates=PT())
                    for i in each_polygon_index(smooth_vorn)
                        p = get_generator(smooth_vorn, i)
                        c = DT.get_centroid(smooth_vorn, i)
                        px, py = DT.getxy(p)
                        cx, cy = DT.getxy(c)
                        dx, dy = px - cx, py - cy
                        ℓ = sqrt(dx^2 + dy^2)
                        _flag = ℓ ≤ 1.0e-1
                        flag += _flag
                        tot += 1
                    end
                end
            end
        end
        @test flag / tot > 0.9
    end
end

@testset "Lattice" begin
    for PT in (DT.ExactKernel, DT.AdaptiveKernel)
        for _ in 1:50
            tri = triangulate_rectangle(0, 1, 0, 1, 11, 11, delete_ghosts=false, predicates=PT())
            tri = triangulate(tri.points, predicates=PT())
            vorn = voronoi(tri, predicates=PT())
            @test validate_tessellation(vorn; check_convex=false, predicates=PT())
            vorn = voronoi(tri, clip=true, predicates=PT())
            @test validate_tessellation(vorn; check_convex=false, check_adjacent=false, predicates=PT())
        end
    end
end

@testset "Example that used to previously break: Plotting a Voronoi tile with unbounded edges intersecting over non-neighbouring sides" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        pts = [
            0.508812 0.656662 0.785124 0.63427 0.444969 0.609253 0.0826304 0.265388 0.830807 0.658346
            0.647732 0.482994 0.809909 0.482046 0.0170022 0.821742 0.835057 0.591724 0.881006 0.97652
        ]
        tri = triangulate(pts)
        vorn = voronoi(tri, predicates=PT())
        clip = DT.get_polygon_coordinates(vorn, 9, DT.polygon_bounds(vorn, 0.1), predicates=PT())
        @test DT.circular_equality(
            collect.(clip), collect.(
                [
                    (0.8374509236290323, 1.0964579554357115)
                    (0.7271857522962732, 0.8973620981454822)
                    (1.7423743304675243, 0.24505806539308744)
                    (2.0053766736553027, 0.1299847935961671)
                    (2.0053766736553027, 1.0964579554357115)
                    (0.8374509236290323, 1.0964579554357115)
                ],
            ), ≈,
        )
    end
end

@testset "Example that used to previously break: Plotting a Voronoi tile with unbounded edges intersecting over non-neighbouring sides, needing THREE corners" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        pts = [
            0.279402 0.874842 0.163028
            0.274178 0.831658 0.223709
        ]
        tri = triangulate(pts, predicates=PT())
        vorn = voronoi(tri, predicates=PT())
        clip = DT.get_polygon_coordinates(vorn, 2, DT.polygon_bounds(vorn, 2.0), predicates=PT())
        @test DT.circular_equality(
            collect.(clip), collect.(
                [
                    (-2.551580488601045, 4.122780970352073)
                    (-0.3314903102646037, 1.5233996567840244)
                    (3.2875066205292076, -2.3420224793856605)
                    (3.2875066205292076, 4.122780970352073)
                    (-2.551580488601045, 4.122780970352073)
                ],
            ), ≈,
        )
    end
end

@testset "Issue #72" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
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
            Float32[0.18311942, 0.8367877],
        ]
        tri = triangulate(points)
        vorn = voronoi(tri, predicates=PT())
        @test validate_tessellation(vorn, predicates=PT())

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
            [0.6525599f0, 0.6425986f0],
        ]
        tri = triangulate(points)
        vorn = voronoi(tri, predicates=PT())
        @test validate_tessellation(vorn, predicates=PT())
    end
end

@testset "Polygon bounds with generators only" begin
    points = rand(2, 50)
    tri = triangulate(points)
    vorn = voronoi(tri, predicates=rt())
    xmin, xmax, ymin, ymax = DT.polygon_bounds(vorn, 0.1; include_polygon_vertices=false)
    _xmin, _xmax = extrema(points[1, :])
    _ymin, _ymax = extrema(points[2, :])
    @test xmin == _xmin - 0.1(_xmax - _xmin)
    @test xmax == _xmax + 0.1(_xmax - _xmin)
    @test ymin == _ymin - 0.1(_ymax - _ymin)
    @test ymax == _ymax + 0.1(_ymax - _ymin)
end

@testset "grow_polygon_outside_of_box" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
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
        vorn = voronoi(tri, predicates=PT())
        a, b, c, d = -4.0, 6.0, -2.0, 4.0
        bounding_box = (a, b, c, d)
        new_vertices, new_points = DT.grow_polygon_outside_of_box(vorn, 1, bounding_box, PT())
        @test collect.(new_points) ≈ collect.(
            [
                (-1.1363636363636365, 5.954545454545455)
                (-0.2999999999999998, 9.3)
                (-2.3999999999999986, 21.900000000000006)
                (-61.96153846153845, 7.846153846153847)
                (-10.807692307692307, 2.730769230769231)
            ],
        )
        @test new_vertices == [1, 2, 3, 4, 5]
        new_vertices, new_points = DT.grow_polygon_outside_of_box(vorn, 5, bounding_box, PT())
        @test collect.(new_points) ≈ collect.(
            [
                (2.710526315789474, 1.868421052631579)
                (-0.7307692307692308, -0.8846153846153846)
                (1.1153846153846159, -13.807692307692308)
                (13.026315789473683, -1.0789473684210529)
            ],
        )
        @test new_vertices == [1, 2, 3, 4]
        new_vertices, new_points = DT.grow_polygon_outside_of_box(vorn, 6, bounding_box, PT())
        @test collect.(new_points) ≈ collect.(
            [
                (23.34210526315789, -4.026315789473685)
                (17.5, 15.499999999999995)
                (3.1, 5.9)
                (2.357142857142857, 2.9285714285714284)
                (2.710526315789474, 1.868421052631579)
            ],
        )
        @test new_vertices == [1, 2, 3, 4, 5]
        new_vertices, new_points = DT.grow_polygon_outside_of_box(vorn, 7, bounding_box, PT())
        @test collect.(new_points) ≈ collect.(
            [
                (-0.7307692307692308, -0.8846153846153846)
                (-4.166666666666666, 0.8333333333333335)
                (-10.807692307692307, 2.730769230769231)
                (-61.96153846153845, 7.846153846153847)
                (1.1153846153846159, -13.807692307692308)
            ],
        )
        @test new_vertices == [1, 2, 3, 4, 5]
        new_vertices, new_points = DT.grow_polygon_outside_of_box(vorn, 8, bounding_box, PT())
        @test collect.(new_points) ≈ collect.(
            [
                (17.5, 15.499999999999995)
                (-4.799999999999997, 36.30000000000001)
                (-0.2999999999999998, 9.3)
                (3.1, 5.9)
            ],
        )
        @test new_vertices == [1, 2, 3, 4]
        @inferred DT.grow_polygon_outside_of_box(vorn, 8, bounding_box, PT())
    end
end

@testset "Clipping polygons to arbitrary bounding box" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
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
        vorn = voronoi(tri, predicates=PT())
        a, b, c, d = -4.0, 6.0, -2.0, 4.0
        bounding_box = (a, b, c, d)
        _pa = get_polygon_coordinates(vorn, 1, bounding_box; predicates=PT())
        _pb = get_polygon_coordinates(vorn, 2, bounding_box; predicates=PT())
        _pc = get_polygon_coordinates(vorn, 3, bounding_box; predicates=PT())
        _pd = get_polygon_coordinates(vorn, 4, bounding_box; predicates=PT())
        _pe = get_polygon_coordinates(vorn, 5, bounding_box; predicates=PT())
        _pf = get_polygon_coordinates(vorn, 6, bounding_box; predicates=PT())
        _pg = get_polygon_coordinates(vorn, 7, bounding_box; predicates=PT())
        _ph = get_polygon_coordinates(vorn, 8, bounding_box; predicates=PT())
        pa = NTuple{2,Float64}[]
        pb = [(0.75, 4.0), (2.357142857142857, 2.9285714285714284), (2.625, 4.0), (0.75, 4.0)]
        pc = [(2.710526315789474, 1.868421052631579), (2.357142857142857, 2.9285714285714284), (0.75, 4.0), (-1.0, 4.0), (-4.0, 1.0), (-4.0, 0.75), (-0.7307692307692308, -0.8846153846153846), (2.710526315789474, 1.868421052631579)]
        pd = [(-4.0, 4.0), (-4.0, 1.0), (-1.0, 4.0), (-4.0, 4.0)]
        pe = [(6.0, 0.9285714285714279), (2.710526315789474, 1.868421052631579), (-0.7307692307692308, -0.8846153846153846), (-0.5714285714285712, -2.0), (6.0, -2.0), (6.0, 0.9285714285714279)]
        pf = [(6.0, 0.9285714285714284), (6.0, 4.0), (2.625, 4.0), (2.357142857142857, 2.9285714285714284), (2.710526315789474, 1.868421052631579), (6.0, 0.9285714285714284)]
        pg = [(-0.5714285714285721, -2.0), (-0.7307692307692308, -0.8846153846153846), (-4.0, 0.75), (-4.0, -2.0), (-0.5714285714285721, -2.0)]
        ph = NTuple{2,Float64}[]
        @test DT.circular_equality(collect.(pa), collect.(_pa), ≈)
        @test DT.circular_equality(collect.(pb), collect.(_pb), ≈)
        @test DT.circular_equality(collect.(pc), collect.(_pc), ≈)
        @test DT.circular_equality(collect.(pd), collect.(_pd), ≈)
        @test DT.circular_equality(collect.(pe), collect.(_pe), ≈)
        @test DT.circular_equality(collect.(pf), collect.(_pf), ≈)
        @test DT.circular_equality(collect.(pg), collect.(_pg), ≈)
        @test DT.circular_equality(collect.(ph), collect.(_ph), ≈)

        # test a small example 
        points = [(0.0, 1.0), (-1.0, 2.0), (-2.0, -1.0)]
        tri = triangulate(points)
        vorn = voronoi(tri, predicates=PT())
        bb = (-1.0, 0.0, -1.0, 2.0)
        coord1 = get_polygon_coordinates(vorn, 1, bb; predicates=PT())
        coord2 = get_polygon_coordinates(vorn, 2, bb; predicates=PT())
        coord3 = get_polygon_coordinates(vorn, 3, bb; predicates=PT())
        _coord1 = [(0.0, 2.0), (0.0, 2.0), (-1.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (0.0, 2.0)]
        _coord2 = [(-1.0, 2.0), (-1.0, 1.0), (0.0, 2.0), (-1.0, 2.0)]
        _coord3 = [(-1.0, -1.0), (0.0, -1.0), (0.0, -1.0), (-1.0, 0.0), (-1.0, -1.0)]
        @test DT.circular_equality(collect.(coord1), collect.(_coord1), ≈)
        @test DT.circular_equality(collect.(coord2), collect.(_coord2), ≈)
        @test DT.circular_equality(collect.(coord3), collect.(_coord3), ≈)

        # example that used to break for G because its unbounded edge 
        # didn't initially start inside the bounding box
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
        vorn = voronoi(tri, predicates=PT())
        bounding_box = (0.0, 5.0, -15.0, 15.0)
        _pa = get_polygon_coordinates(vorn, 1, bounding_box; predicates=PT())
        _pb = get_polygon_coordinates(vorn, 2, bounding_box; predicates=PT())
        _pc = get_polygon_coordinates(vorn, 3, bounding_box; predicates=PT())
        _pd = get_polygon_coordinates(vorn, 4, bounding_box; predicates=PT())
        _pe = get_polygon_coordinates(vorn, 5, bounding_box; predicates=PT())
        _pf = get_polygon_coordinates(vorn, 6, bounding_box; predicates=PT())
        _pg = get_polygon_coordinates(vorn, 7, bounding_box; predicates=PT())
        # For _pg: This case used to be broken because the initial unbounded ray did not touch the bounding box, 
        # but then later it does! So just using the Liang-Barsky algorithm by itself is not sufficient.
        # To fix this, I added the stuff about maximum_distance_to_box inside grow_polygon_outside_of_box
        _ph = get_polygon_coordinates(vorn, 8, bounding_box; predicates=PT())
        pa = NTuple{2,Float64}[]
        pb = [(0.0, 4.499999999999998), (2.357142857142857, 2.9285714285714284), (3.1, 5.9), (0.0, 9.000000000000002), (0.0, 4.499999999999998)]
        pc = [(0.0, -0.3000000000000007), (2.710526315789474, 1.868421052631579), (2.357142857142857, 2.9285714285714284), (0.0, 4.499999999999998), (0.0, -0.3000000000000007)]
        pd = NTuple{2,Float64}[]
        pe = [(5.0, 1.2142857142857117), (2.710526315789474, 1.868421052631579), (0.0, -0.3000000000000007), (0.0, -6.0), (1.2857142857142865, -15.0), (5.0, -15.0), (5.0, 1.2142857142857117)]
        pf = [(5.0, 1.2142857142857153), (5.0, 7.166666666666664), (3.1, 5.9), (2.357142857142857, 2.9285714285714284), (2.710526315789474, 1.868421052631579), (5.0, 1.2142857142857153)]
        pg = [(1.2857142857142865, -15.0), (0.0, -6.0), (0.0, -15.0), (1.2857142857142865, -15.0)]
        ph = [(5.0, 7.166666666666664), (5.0, 15.0), (0.0, 15.0), (0.0, 9.000000000000002), (3.1, 5.9), (5.0, 7.166666666666664)]
        @test DT.circular_equality(collect.(pa), collect.(_pa), ≈)
        @test DT.circular_equality(collect.(pb), collect.(_pb), ≈)
        @test DT.circular_equality(collect.(pc), collect.(_pc), ≈)
        @test DT.circular_equality(collect.(pd), collect.(_pd), ≈)
        @test DT.circular_equality(collect.(pe), collect.(_pe), ≈)
        @test DT.circular_equality(collect.(pf), collect.(_pf), ≈)
        @test DT.circular_equality(collect.(pg), collect.(_pg), ≈)
        @test DT.circular_equality(collect.(ph), collect.(_ph), ≈)
    end
end

@testset "maximum_distance_to_box" begin
    a, b, c, d = -8.0, -2.0, 0.0, 7.0
    A = (-5.34, 4.51)
    δ = DT.maximum_distance_to_box(a, b, c, d, A)
    @test sqrt(δ) ≈ 5.6121029926401
    A = (-8.0, 5.0)
    δ = DT.maximum_distance_to_box(a, b, c, d, A)
    @test sqrt(δ) ≈ 7.8102496759067
    A = (-12.0, 10.0)
    δ = DT.maximum_distance_to_box(a, b, c, d, A)
    @test sqrt(δ) ≈ 14.142135623731
end

@testset "Previously broken example" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        # The points were duplicated from refine and not included in tri, but 
        # retriangulate kept trying to use them
        points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], predicates=PT())
        refine!(tri; max_area=1.0e-2, min_angle=29.871, predicates=PT())
        vorn = voronoi(tri, predicates=PT())
        smooth_vorn = centroidal_smooth(vorn; maxiters=250, predicates=PT())
        @test validate_tessellation(smooth_vorn, predicates=PT())
    end
end

@testset "toggle_inf_warn" begin
    @test DT.INF_WARN[]
    DT.toggle_inf_warn!()
    @test !DT.INF_WARN[]
    DT.toggle_inf_warn!()
    @test DT.INF_WARN[]
end

@testset "Issue #72" begin
    function is_point_in_polygon(point::Point, polygon::Polygon)
        # Check if the point is on the boundary of the polygon
        point in polygon && return false

        n = length(coordinates(polygon))
        inside = false
        p1 = polygon[1]
        for i in 2:(n+1)
            p2 = polygon[mod1(i, n)]
            if ((p1[2] > point[2]) != (p2[2] > point[2])) && (point[1] < (p2[1] - p1[1]) * (point[2] - p1[2]) / (p2[2] - p1[2]) + p1[1])
                inside = !inside
            end
            p1 = copy(p2)
        end
        return inside
    end
    b_vec = reverse([Point(0.0, 0.0), Point(0.0, 1.0), Point(1.0, 1.0), Point(1.0, 0.0), Point(0.5, 0.5)])
    p = Polygon(b_vec)
    p2 = coordinates(p)
    push!(p2, p2[1])
    boundary_nodes, points2 = convert_boundary_points_to_indices(getxy.(p2))
    triangulation = triangulate(points2; boundary_nodes)
    minx, miny = extrema(p2)[1]
    maxx, maxy = extrema(p2)[2]
    for i in range(minx + 0.05, maxx - 0.05, step=0.05)
        for j in range(maxy - 0.05, miny + 0.05, step=-0.05)
            if is_point_in_polygon(Point(i, j), p)
                try
                    add_point!(triangulation, (i, j))
                catch
                    @warn "There was a problem adding this point."
                end
            end
        end
    end
    vorn = voronoi(triangulation, clip=false)
    @test validate_tessellation(vorn, check_convex=false)
end

@testset "Clipping to a generic convex polygon" begin
    a, b, c, d, e, f, g, h, i = (-5.0, 7.0),
    (-8.0, 6.0), (-7.0, 4.0), (-4.0, 3.0), (-6.0, 6.0),
    (-1.0, 6.0), (-4.0, 6.0), (-3.0, 5.0), (-2.0, 8.0)
    points = [a, b, c, d, e, f, g, h, i]
    J, K, L, M = (-10.0, 2.0), (2.0, 2.0), (2.0, 12.0), (-10.0, 12.0)
    clip_points = (J, K, L, M)
    clip_vertices = (1, 2, 3, 4, 1)
    clip_polygon = (clip_points, clip_vertices)
    rng = StableRNG(123)
    tri = triangulate(points; rng)
    @test_logs (:warn, "Clip polygon provided but clip = false. Ignoring clip polygon.") voronoi(tri, clip_polygon=clip_polygon)
    vorn = voronoi(tri, clip=false, clip_polygon=clip_polygon, rng=rng)
    @test DT.find_all_exterior_circumcenters(vorn, clip_polygon...) == Set{Int}()
    int, ints = DT.find_all_intersecting_polygons(vorn, clip_polygon..., DT.find_all_exterior_circumcenters(vorn, clip_polygon...))
    @test int == Set((5, 7, 8))
    @test ints == Set((4, 6, 2, 9, 3, 1))
    vorn = voronoi(tri, clip=true, clip_polygon=clip_polygon, rng=rng)
    pt1 = [
        (-5.0, 6.0)
        (-3.5, 7.5)
        (-5.0, 12.0)
        (-8.333333333333334, 12.0)
        (-7.0, 8.0)
        (-5.0, 6.0)
    ]
    pt2 = [
        (-7.0, 8.0)
        (-8.333333333333334, 12.0)
        (-10.0, 12.0)
        (-10.0, 3.75)
        (-7.0, 5.25)
        (-7.0, 8.0)
    ]
    pt3 = [
        (-7.0, 5.25)
        (-10.0, 3.75)
        (-10.0, 2.0)
        (-6.0, 2.0)
        (-5.214285714285714, 4.357142857142857)
        (-7.0, 5.25)
    ]
    pt4 = [
        (-5.214285714285714, 4.357142857142857)
        (-6.0, 2.0)
        (0.0, 2.0)
        (-0.5, 2.5)
        (-4.5, 4.5)
        (-5.0, 4.5)
        (-5.214285714285714, 4.357142857142857)
    ]
    pt5 = [
        (-5.0, 4.5)
        (-5.0, 6.0)
        (-7.0, 8.0)
        (-7.0, 5.25)
        (-5.214285714285714, 4.357142857142857)
        (-5.0, 4.5)
    ]
    pt6 = [
        (0.0, 2.0)
        (2.0, 2.0)
        (2.0, 8.75)
        (-2.5, 6.5)
        (-0.5, 2.5)
        (0.0, 2.0)
    ]
    pt7 = [
        (-5.0, 4.5)
        (-4.5, 4.5)
        (-2.5, 6.5)
        (-3.5, 7.5)
        (-5.0, 6.0)
        (-5.0, 4.5)
    ]
    pt8 = [
        (-0.5, 2.5)
        (-2.5, 6.5)
        (-4.5, 4.5)
        (-0.5, 2.5)
    ]
    pt9 = [
        (2.0, 8.75)
        (2.0, 12.0)
        (-5.0, 12.0)
        (-3.5, 7.5)
        (-2.5, 6.5)
        (2.0, 8.75)
    ]
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[1]...), pt1, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[2]...), pt2, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[3]...), pt3, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[4]...), pt4, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[5]...), pt5, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[6]...), pt6, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[7]...), pt7, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[8]...), pt8, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[9]...), pt9, ⪧)
    @test validate_tessellation(vorn; check_area=false)
    @test isempty(vorn.unbounded_polygons)
    @test vorn.boundary_polygons == Set((4, 6, 2, 9, 3, 1))

    vorn = voronoi(tri, clip=true, clip_polygon=clip_polygon, rng=rng, smooth=true)
    pt1 = [
        (-5.1590315116169645, 8.028359015649276)
        (-2.0329830757675342, 9.155154224395774)
        (-1.5930583264171334, 12.0)
        (-6.067362675362922, 12.0)
        (-6.059020473726599, 8.97155876785828)
        (-5.1590315116169645, 8.028359015649276)
    ]
    pt2 = [
        (-10.0, 12.0)
        (-10.0, 8.234478977169747)
        (-6.059020473726599, 8.97155876785828)
        (-6.067362675362922, 12.0)
        (-10.0, 12.0)
    ]
    pt3 = [
        (-10.0, 5.764961735169299)
        (-10.0, 2.0)
        (-6.066696497738481, 2.0)
        (-6.058578287005368, 5.028158005937482)
        (-10.0, 5.764961735169299)
    ]
    pt4 = [
        (-6.058578287005368, 5.028158005937482)
        (-6.066696497738481, 2.0)
        (-1.592374540035335, 2.0)
        (-2.032629690319553, 4.8454051776002)
        (-5.158841519410821, 5.971443006049375)
        (-6.058578287005368, 5.028158005937482)
    ]
    pt5 = [
        (-5.158841519410821, 5.971443006049375)
        (-5.1590315116169645, 8.028359015649276)
        (-6.059020473726599, 8.97155876785828)
        (-10.0, 8.234478977169747)
        (-10.0, 5.764961735169299)
        (-6.058578287005368, 5.028158005937482)
        (-5.158841519410821, 5.971443006049375)
    ]
    pt6 = [
        (-1.592374540035335, 2.0)
        (2.0, 2.0)
        (2.0, 7.000229455337406)
        (-0.15539885081025462, 7.000135534420414)
        (-2.032629690319553, 4.8454051776002)
        (-1.592374540035335, 2.0)
    ]
    pt7 = [
        (-5.158841519410821, 5.971443006049375)
        (-2.032629690319553, 4.8454051776002)
        (-0.15539885081025462, 7.000135534420414)
        (-2.0329830757675342, 9.155154224395774)
        (-5.1590315116169645, 8.028359015649276)
        (-5.158841519410821, 5.971443006049375)
    ]
    @test !DT.has_polygon(vorn, 8)
    pt9 = [
        (2.0, 7.000229455337406)
        (2.0, 12.0)
        (-1.5930583264171334, 12.0)
        (-2.0329830757675342, 9.155154224395774)
        (-0.15539885081025462, 7.000135534420414)
        (2.0, 7.000229455337406)
    ]
    @test validate_tessellation(vorn; check_area=false)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[1]...), pt1, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[2]...), pt2, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[3]...), pt3, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[4]...), pt4, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[5]...), pt5, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[6]...), pt6, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[7]...), pt7, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[9]...), pt9, ⪧)
    for index in each_polygon_index(vorn)
        c = DT.get_centroid(vorn, index)
        @test c ⪧ get_generator(vorn, index) rtol = 1e-2 atol = 1e-4
    end

    J, K, L, M = (-2.0, 12.0), (-8.0, 4.0), (-4.0, 2.0), (-2.0, 10.0)
    clip_points = (J, K, L, M)
    clip_vertices = (1, 2, 3, 4, 1)
    clip_polygon = (clip_points, clip_vertices)
    rng = StableRNG(123)
    tri = triangulate(points; rng)
    vorn = voronoi(tri, clip=false, clip_polygon=clip_polygon, rng=rng)
    ext = DT.find_all_exterior_circumcenters(vorn, clip_polygon...)
    @test ext == Set((7, 2, 1))
    int, ints = DT.find_all_intersecting_polygons(vorn, clip_polygon..., ext)
    @test int == Set()
    @test ints == Set((5, 4, 7, 2, 9, 8, 3, 6, 1))
    vorn = voronoi(tri, clip=true, clip_polygon=clip_polygon, rng=rng)
    @test validate_tessellation(vorn; check_area=false)
    pt1 = [
        (-5.857142857142858, 6.857142857142857)
        (-5.0, 6.0)
        (-3.5, 7.5)
        (-4.076923076923077, 9.23076923076923)
        (-5.857142857142858, 6.857142857142857)
    ]
    pt2 = [
        (-7.0, 5.333333333333333)
        (-7.1, 5.2)
        (-7.0, 5.25)
        (-7.0, 5.333333333333333)
    ]
    pt3 = [
        (-7.0, 5.25)
        (-7.1, 5.2)
        (-8.0, 4.0)
        (-5.7142857142857135, 2.8571428571428568)
        (-5.214285714285714, 4.357142857142857)
        (-7.0, 5.25)
    ]
    pt4 = [
        (-5.214285714285714, 4.357142857142857)
        (-5.7142857142857135, 2.8571428571428568)
        (-4.0, 2.0)
        (-3.5, 4.0)
        (-4.5, 4.5)
        (-5.0, 4.5)
        (-5.214285714285714, 4.357142857142857)
    ]
    pt5 = [
        (-5.0, 4.5)
        (-5.0, 6.0)
        (-5.857142857142858, 6.857142857142857)
        (-7.0, 5.333333333333333)
        (-7.0, 5.25)
        (-5.214285714285714, 4.357142857142857)
        (-5.0, 4.5)
    ]
    pt7 = [
        (-5.0, 4.5)
        (-4.5, 4.5)
        (-3.0, 6.0)
        (-2.8, 6.8)
        (-3.5, 7.5)
        (-5.0, 6.0)
        (-5.0, 4.5)
    ]
    pt8 = [
        (-3.5, 4.0)
        (-3.0, 6.0)
        (-4.5, 4.5)
        (-3.5, 4.0)
    ]
    pt9 = [
        (-2.0, 10.0)
        (-2.0, 12.0)
        (-4.076923076923077, 9.23076923076923)
        (-3.5, 7.5)
        (-2.8, 6.8)
        (-2.0, 10.0)
    ]
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[1]...), pt1, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[2]...), pt2, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[3]...), pt3, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[4]...), pt4, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[5]...), pt5, ⪧)
    @test !DT.has_polygon(vorn, 6)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[7]...), pt7, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[8]...), pt8, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[9]...), pt9, ⪧)
    @test vorn.boundary_polygons == Set((5, 4, 7, 2, 9, 8, 3, 1))
    @test isempty(vorn.unbounded_polygons)

    vorn = voronoi(tri, clip=true, clip_polygon=clip_polygon, rng=rng, smooth=true)
    pt1 = [
        (-2.364836958330981, 8.540652166676075)
        (-4.029198089942065, 9.294402546743912)
        (-5.173208375837535, 7.7690554988832865)
        (-4.17127917698101, 6.722638795624013)
        (-2.855277812679444, 6.578888749282224)
        (-2.364836958330981, 8.540652166676075)
    ]
    pt2 = [
        (-6.659997738183602, 5.786669682421865)
        (-5.185697692423832, 4.903140482793648)
        (-4.17127917698101, 6.722638795624013)
        (-5.173208375837535, 7.7690554988832865)
        (-6.659997738183602, 5.786669682421865)
    ]
    pt3 = [
        (-8.0, 4.0)
        (-6.026707081932514, 3.013353540966257)
        (-5.120548573308532, 4.754263677601372)
        (-5.185697692423832, 4.903140482793648)
        (-6.659997738183602, 5.786669682421865)
        (-8.0, 4.0)
    ]
    pt4 = [
        (-4.0, 2.0)
        (-3.4550765743281895, 4.179693702687242)
        (-5.120548573308532, 4.754263677601372)
        (-6.026707081932514, 3.013353540966257)
        (-4.0, 2.0)
    ]
    @test !DT.has_polygon(vorn, 5)
    @test !DT.has_polygon(vorn, 6)
    @test !DT.has_polygon(vorn, 7)
    pt8 = [
        (-3.4550765743281895, 4.179693702687242)
        (-2.855277812679444, 6.578888749282224)
        (-4.17127917698101, 6.722638795624013)
        (-5.185697692423832, 4.903140482793648)
        (-5.120548573308532, 4.754263677601372)
        (-3.4550765743281895, 4.179693702687242)
    ]
    pt9 = [
        (-2.0, 10.0)
        (-2.0, 12.0)
        (-4.029198089942065, 9.294402546743912)
        (-2.364836958330981, 8.540652166676075)
        (-2.0, 10.0)
    ]
    @test validate_tessellation(vorn; check_area=false)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[1]...), pt1, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[2]...), pt2, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[3]...), pt3, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[4]...), pt4, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[8]...), pt8, ⪧)
    @test DT.circular_equality(DT.get_polygon_point(vorn, vorn.polygons[9]...), pt9, ⪧)
    for index in each_polygon_index(vorn)
        c = DT.get_centroid(vorn, index)
        @test c ⪧ get_generator(vorn, index) rtol = 1e-2 atol = 1e-4
    end
end

@testset "==" begin
    tri = triangulate(rand(2, 100))
    vorn = voronoi(tri, clip=true)
    vorn2 = voronoi(tri, clip=true)
    @test vorn == vorn
    @test vorn == vorn2
    g1 = vorn.generators[1]
    vorn.generators[1] = vorn.generators[2]
    @test vorn != vorn2
    vorn.generators[1] = g1
    @test vorn == vorn2
    p1 = vorn.polygon_points[1]
    vorn.polygon_points[1] = vorn.polygon_points[2]
    @test vorn != vorn2
    vorn.polygon_points[1] = p1
    @test vorn == vorn2
    poly1 = vorn.polygons[1]
    vorn.polygons[1] = vorn.polygons[2]
    @test vorn != vorn2
    vorn.polygons[1] = poly1
    @test vorn == vorn2
    T1 = vorn.circumcenter_to_triangle[1]
    vorn.circumcenter_to_triangle[1] = vorn.circumcenter_to_triangle[2]
    @test vorn != vorn2
    vorn.circumcenter_to_triangle[1] = T1
    @test vorn == vorn2
    T2 = vorn.circumcenter_to_triangle[2]
    p1 = vorn.triangle_to_circumcenter[T1]
    vorn.triangle_to_circumcenter[T1] = vorn.triangle_to_circumcenter[T2]
    @test vorn != vorn2
    vorn.triangle_to_circumcenter[T1] = p1
    @test vorn == vorn2
    push!(vorn.unbounded_polygons, 5)
    @test vorn != vorn2
    delete!(vorn.unbounded_polygons, 5)
    (u, v), w = first(vorn.adjacent.adjacent)
    vorn.adjacent.adjacent[(u, v)] = w + 1
    @test vorn != vorn2
    vorn.adjacent.adjacent[(u, v)] = w
    @test vorn == vorn2
    push!(vorn.boundary_polygons, 5000)
    @test vorn != vorn2
    delete!(vorn.boundary_polygons, 5000)
    @test vorn == vorn2
end

@testset "copy/deepcopy" begin
    vorn = voronoi(triangulate(rand(2, 100)), clip=true, smooth=true)
    vorn1 = deepcopy(vorn)
    vorn2 = copy(vorn)
    @test typeof(vorn1) == typeof(vorn) == typeof(vorn2)
    for _vorn in (vorn1, vorn2)
        @test vorn == _vorn && !(vorn === _vorn)
        for f in fieldnames(typeof(vorn))
            @test getfield(vorn, f) == getfield(_vorn, f)
            @test !(getfield(vorn, f) === getfield(_vorn, f))
        end
    end
end

@testset "Issue #206" begin
    Random.seed!(123)
    pts = [(0.1, 0.2), (1.1, 0.2), (0.1, 1.7)]
    tri = triangulate(pts)
    vorn = voronoi(tri, clip=true, clip_polygon=(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], [1, 2, 3, 4, 1]
    ))
    @test validate_tessellation(vorn; check_area=false)
    A = 0.0
    for i in each_polygon_index(vorn)
        A += get_area(vorn, i)
    end
    @test A ≈ 1.0
    vorn = voronoi(tri)
    poly = get_polygon_coordinates(vorn, 3, (0.0, 1.0, 0.0, 1.0))
    @test poly ⪧ [(0.675, 1.0), (0.0, 1.0), (0.0, 0.95), (0.6, 0.95), (0.675, 1.0)]
end