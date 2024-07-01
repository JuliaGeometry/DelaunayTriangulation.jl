using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using CairoMakie
using StableRNGs
using ReferenceTests
using StatsBase



@testset "Triangulating random convex polygons" begin
    for n in Iterators.flatten([3:20, 25:50:1000])
        points = rand(2, n)
        S = get_random_convex_polygon(points)
        skip_points = setdiff(axes(points, 2), S)
        for delete_ghosts in (false, true)
            tri_bowyer = triangulate(points; skip_points, delete_ghosts)
            tri_chew = triangulate_convex(points, S; delete_ghosts)
            @test validate_triangulation(tri_bowyer)
            @test validate_triangulation(tri_chew)
            @test DT.compare_triangle_collections(get_triangles(tri_bowyer), get_triangles(tri_chew))
            for (ij, k) in get_adjacent(tri_chew).adjacent
                if DT.edge_exists(k)
                    @test get_adjacent(tri_bowyer, ij) == k
                end
            end
            for (w, E) in get_adjacent2vertex(tri_chew).adjacent2vertex
                for ij in each_edge(E)
                    @test DT.contains_edge(ij, get_adjacent2vertex(tri_bowyer, w))
                end
            end
            for i in each_vertex(tri_chew)
                @test get_neighbours(tri_chew, i) == get_neighbours(tri_bowyer, i)
            end
            @test DT.get_convex_hull_vertices(tri_chew) == [S; S[begin]]
            @test tri_chew.representative_point_list[1].x ≈ mean(points[1, s] for s in S)
            @test tri_chew.representative_point_list[1].y ≈ mean(points[2, s] for s in S)
            @test tri_chew.representative_point_list[1].n ≈ length(S)
            @test compare_trees(DT.get_polygon_hierarchy(tri_bowyer), DT.construct_polygon_hierarchy(points))
            @test compare_trees(DT.get_polygon_hierarchy(tri_chew), DT.construct_polygon_hierarchy(points))
        end
    end
end

@testset "Triangulating a small polygon with some collinearities" begin
    for T in (Float64, Float32)
        p1 = T[8.0, 4.0]
        p2 = T[10.0, 4.0]
        p3 = T[12.0, 4.0]
        p4 = T[14.0, 4.0]
        p5 = T[14.0, 6.0]
        p6 = T[14.0, 8.0]
        p7 = T[14.0, 10.0]
        p8 = T[12.0, 10.0]
        p9 = T[10.0, 10.0]
        p10 = T[8.0, 10.0]
        p11 = T[8.0, 8.0]
        p12 = T[8.0, 6.0]
        pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12]
        for _ in 1:100
            tri_chew = triangulate_convex(pts, 1:12)
            @test validate_triangulation(tri_chew)
        end
    end
end

@testset "Triangulating more random polygons, smaller size" begin
    for _ in 1:1500
        pts = rand(2, 50)
        p1 = [0.0, 0.0]
        p2 = [1.0, 0.0]
        p3 = [0.0, 1.0]
        pts[:, 11] .= p1
        pts[:, 27] .= p2
        pts[:, 5] .= p3
        S = [11, 27, 5, 11]
        @test_throws AssertionError("S must not be circular.") triangulate_convex(pts, S)
        pop!(S)
        tri_chew = triangulate_convex(pts, S)
        tri_bowyer = triangulate(pts; skip_points=setdiff(1:50, [11, 27, 5]), delete_ghosts=false)
        @test get_convex_hull(tri_chew) == get_convex_hull(tri_bowyer)
        @test DT.compare_triangle_collections(get_triangles(tri_chew), get_triangles(tri_bowyer))
        @test (get_adjacent ∘ get_adjacent)(tri_chew) == (get_adjacent ∘ get_adjacent)(tri_bowyer)
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_chew) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bowyer)
        @test get_graph(tri_chew) == get_graph(tri_bowyer)

        pts[:, 28] .= [1.01, 1.01]
        S = [11, 27, 28, 5]
        tri_chew = triangulate_convex(pts, S)
        tri_bowyer = triangulate(pts; skip_points=setdiff(1:50, [11, 27, 5, 28]), delete_ghosts=false)
        @test get_convex_hull(tri_chew) == get_convex_hull(tri_bowyer)
        @test DT.compare_triangle_collections(get_triangles(tri_chew), get_triangles(tri_bowyer))
        @test (get_adjacent ∘ get_adjacent)(tri_chew) == (get_adjacent ∘ get_adjacent)(tri_bowyer)
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_chew) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bowyer)
        @test get_graph(tri_chew) == get_graph(tri_bowyer)
    end
end