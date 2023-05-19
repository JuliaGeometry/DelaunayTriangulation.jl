using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using StableRNGs

include("../helper_functions.jl")

rng = StableRNG(8888)

@testset "Adding a point" begin
    pts = rand(rng, 2, 500)
    tri = DT.initialise_bowyer_watson(pts)
    for i in setdiff(each_point_index(pts), get_vertices(tri))
        add_point!(tri, i)
        @test validate_triangulation(tri)
    end
    convex_hull!(tri; reconstruct=false)
    delete_ghost_triangles!(tri)
    clear_empty_features!(tri)
    _tri = triangulate(pts; rng)
    clear_empty_features!(_tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    @test get_convex_hull(tri) == get_convex_hull(_tri)
end

@testset "Returning the triangle" begin
    tri = triangulate([rand(2) for _ in 1:50])
    V = add_point!(tri, 0.5, 0.3)
    @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V..., (0.5, 0.3)))
end

@testset "Peeking" begin
    for _ in 1:250
        pts = [(rand(), rand()) for _ in 1:50]
        tri = triangulate(pts)
        _tri = deepcopy(tri)
        history = DT.initialise_event_history(tri)
        add_point!(tri, 0.4, 0.4, store_event_history=Val(true), event_history=history, peek=Val(true))
        @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
        clear_empty_features!(tri)
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

        tri = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 11, 6, add_ghost_triangles=true)
        _tri = deepcopy(tri)
        history = DT.initialise_event_history(tri)
        for (x, y) in [(0.0, 0.23), (0.0, 0.8), (1.0, 1.0), (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (0.39881, 0.0), (0.5, 1.0), (0.0, 0.4), (1.0, 0.881),
            [(0.0, rand()) for _ in 1:50]..., [(rand(), 0.0) for _ in 1:50]..., [(rand(), 1.0) for _ in 1:50]..., [(1.0, rand()) for _ in 1:50]...,
            [(rand(), rand()) for _ in 1:500]...]
            empty!(history)
            if rand() < 1 / 2
                add_point!(tri, x, y, store_event_history=Val(true), event_history=history, peek=Val(true))
            else
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history, peek=Val(true))
            end
            clear_empty_features!(tri)
            @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
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
                @test !DT.is_outside(DT.point_position_relative_to_circumcircle(tri, T, (x, y)))
            end
            @test get_points(tri) == get_points(_tri)
        end
    end
end

@testset "Peeking concurrency" begin
    pts = [(rand(), rand()) for _ in 1:50]
    tri = triangulate(pts, delete_ghosts = false)
    Base.Threads.@threads for _ in 1:500
        history = DT.initialise_event_history(tri)
        add_point!(tri, rand(2), store_event_history=Val(true), event_history=history, peek=Val(true))
    end
    @test num_points(tri) == 50
end