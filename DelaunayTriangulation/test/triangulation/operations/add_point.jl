using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using StableRNGs
using StaticArrays


rng = StableRNG(8888)

@testset "Adding a point" begin
    rng = StableRNG(1)
    for _ in 1:150
        pts = tuple.(rand(rng, 50), rand(rng, 50))
        tri = DT.triangulate(pts, delete_ghosts = false, skip_points = 4:50; rng)
        for i in setdiff(DT.each_point_index(pts), DT.get_vertices(tri))
            add_point!(tri, i; predicates = rt(rng), rng)
            @test validate_triangulation(tri)
        end
        @test num_solid_vertices(tri) == DT.num_points(tri)
        DT.convex_hull!(tri; reconstruct = false)
        DT.clear_empty_features!(tri)
        _tri = triangulate(pts; rng)
        DT.compute_representative_points!(tri)
        DT.clear_empty_features!(_tri)
        @test tri == _tri

        n = DT.num_points(tri)
        ctr = 1
        for _ in 1:50
            add_point!(tri, rand(rng), rand(rng); rng)
            @test num_solid_vertices(tri) == DT.num_points(tri) == n + ctr
            ctr += 1
        end
    end
end

@testset "Returning the triangle" begin
    rng = StableRNG(2)
    tri = triangulate([rand(rng, 2) for _ in 1:50]; rng)
    V = add_point!(tri, 0.5, 0.3; rng)
    @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V..., (0.5, 0.3)))
end

@testset "Peeking" begin
    rng = StableRNG(3)
    for _ in 1:250
        pts = [(rand(rng), rand(rng)) for _ in 1:50]
        tri = triangulate(pts; rng)
        _tri = deepcopy(tri)
        history = DT.InsertionEventHistory(tri)
        add_point!(tri, 0.4, 0.4, store_event_history = Val(true), event_history = history, peek = Val(true); rng)
        @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
        DT.clear_empty_features!(tri)
        @test get_adjacent(tri) == get_adjacent(_tri)
        @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
        @test get_graph(tri) == get_graph(_tri)
        @test get_convex_hull(tri) == get_convex_hull(_tri)
        @test length(history.added_triangles) > 0
        @test length(history.deleted_triangles) > 0
        for T in DT.each_added_triangle(history)
            @test !DT.contains_triangle(T, get_triangles(tri))[2]
        end
        for T in history.deleted_triangles
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(tri, T, (0.4, 0.4)))
        end
        @test get_points(tri) == get_points(_tri)

        tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 11, 6, delete_ghosts = false)
        _tri = deepcopy(tri)
        history = DT.InsertionEventHistory(tri)
        for (x, y) in [
                (0.0, 0.23), (0.0, 0.8), (1.0, 1.0), (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (0.39881, 0.0), (0.5, 1.0), (0.0, 0.4), (1.0, 0.881),
                [(0.0, rand(rng)) for _ in 1:50]..., [(rand(rng), 0.0) for _ in 1:50]..., [(rand(rng), 1.0) for _ in 1:50]..., [(1.0, rand(rng)) for _ in 1:50]...,
                [(rand(rng), rand(rng)) for _ in 1:500]...,
            ]
            empty!(history)
            if rand(rng) < 1 / 2
                add_point!(tri, x, y, store_event_history = Val(true), event_history = history, peek = Val(true); rng)
            else
                add_point!(tri, (x, y), store_event_history = Val(true), event_history = history, peek = Val(true); rng)
            end
            DT.clear_empty_features!(tri)
            @test tri == _tri
            @test length(history.added_triangles) > 0
            @test length(history.deleted_triangles) > 0
            for T in DT.each_added_triangle(history)
                @test !DT.contains_triangle(T, get_triangles(tri))[2]
            end
            for T in history.deleted_triangles
                @test !DT.is_outside(DT.point_position_relative_to_circumcircle(tri, T, (x, y)))
            end
            @test get_points(tri) == get_points(_tri)
        end
    end
end

@testset "Peeking concurrency" begin
    rng = StableRNG(4)
    pts = [(rand(rng), rand(rng)) for _ in 1:50]
    tri = triangulate(pts, delete_ghosts = false; rng)
    Base.Threads.@threads for _ in 1:500
        history = DT.InsertionEventHistory(tri)
        add_point!(tri, rand(rng, 2), store_event_history = Val(true), event_history = history, peek = Val(true); rng)
    end
    @test DT.num_points(tri) == 50
end
