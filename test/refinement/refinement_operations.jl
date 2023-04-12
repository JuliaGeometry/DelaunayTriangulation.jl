using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using StableRNGs

include("../helper_functions.jl")

@testset "Tracking the Bowyer-Watson history" begin
    for i in 1:1000
        curve_1 = [[
            (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
            (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
            (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
            (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
            (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
            (-4.0, 0.0), (0.0, 0.0),
        ]]
        curve_2 = [[
            (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
            (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
            (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
        ]]
        curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
        curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
        curves = [curve_1, curve_2, curve_3, curve_4]
        points = [
            (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
            (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
            (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
            (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
            (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
            (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
            (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
            (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
            (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
            (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
            (7.0, 7.8), (7.0, 7.4)]
        boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
        uncons_tri = triangulate(points)
        cons_tri = triangulate(points; boundary_nodes=boundary_nodes)
        add_ghost_triangles!(uncons_tri)
        add_ghost_triangles!(cons_tri)
        add_edge!(cons_tri, 16, 1)
        add_edge!(cons_tri, 16, 20)
        add_edge!(uncons_tri, 16, 1)
        add_edge!(uncons_tri, 33, 1)
        tri = cons_tri

        storage = DT.initialise_event_history(tri)
        @test storage.added_triangles == Set{DT.triangle_type(tri)}()
        @test storage.deleted_triangles == Set{DT.triangle_type(tri)}()
        @test storage.added_segments == Set{DT.edge_type(tri)}()
        @test storage.deleted_segments == Set{DT.edge_type(tri)}()
        @test storage.added_boundary_segments == Set{DT.edge_type(tri)}()
        @test storage.deleted_boundary_segments == Set{DT.edge_type(tri)}()
        add_triangle!(storage, (1, 2, 3))
        add_triangle!(storage, (2, 5, 7))
        delete_triangle!(storage, (17, 3, 5))
        delete_triangle!(storage, (5, 7, 15))
        add_edge!(storage, (1, 2))
        add_edge!(storage, (2, 5))
        DT.delete_edge!(storage, (5, 7))
        DT.split_boundary_edge!(storage, 2, 17, 10)
        DT.split_boundary_edge!(storage, 18, 23, 50)
        @test storage.added_triangles == Set{DT.triangle_type(tri)}([(1, 2, 3), (2, 5, 7)])
        @test storage.deleted_triangles == Set{DT.triangle_type(tri)}([(17, 3, 5), (5, 7, 15)])
        @test storage.added_segments == Set{DT.edge_type(tri)}([(1, 2), (2, 5)])
        @test storage.deleted_segments == Set{DT.edge_type(tri)}([(5, 7)])
        @test storage.added_boundary_segments == Set{DT.edge_type(tri)}([(2, 10), (10, 17), (18, 50), (50, 23)])
        @test storage.deleted_boundary_segments == Set{DT.edge_type(tri)}([(2, 17), (18, 23)])
        @test DT.has_segment_changes(storage)
        @test DT.each_added_triangle(storage) == each_triangle(storage.added_triangles)
        @test DT.each_added_segment(storage) == each_edge(storage.added_segments)
        @test DT.each_added_boundary_segment(storage) == each_edge(storage.added_boundary_segments)
        empty!(storage)
        @test storage.added_triangles == Set{DT.triangle_type(tri)}()
        @test storage.deleted_triangles == Set{DT.triangle_type(tri)}()
        @test storage.added_segments == Set{DT.edge_type(tri)}()
        @test storage.deleted_segments == Set{DT.edge_type(tri)}()
        @test storage.added_boundary_segments == Set{DT.edge_type(tri)}()
        @test storage.deleted_boundary_segments == Set{DT.edge_type(tri)}()
        @test !DT.has_segment_changes(storage)

        cur_triangles = deepcopy(get_triangles(tri))
        cur_constrained_edges = deepcopy(get_all_constrained_edges(tri))
        cur_boundary = deepcopy(keys(get_boundary_edge_map(tri)))
        add_point!(tri, 2.0, 2.0; store_event_history=Val(true), event_history=storage)
        new_triangles = get_triangles(tri)
        new_constrained_edges = get_all_constrained_edges(tri)
        new_boundary = keys(get_boundary_edge_map(tri))
        added_triangles = setdiff(new_triangles, cur_triangles)
        deleted_triangles = setdiff(cur_triangles, new_triangles)
        @test isempty(storage.added_segments)
        @test isempty(storage.deleted_segments)
        @test isempty(storage.added_boundary_segments)
        @test isempty(storage.deleted_boundary_segments)
        @test DT.compare_triangle_collections(storage.added_triangles, Set{DT.triangle_type(tri)}(added_triangles))
        @test DT.compare_triangle_collections(storage.deleted_triangles, Set{DT.triangle_type(tri)}(deleted_triangles))
        @test !DT.has_segment_changes(storage)
        empty!(storage)

        cur_triangles = deepcopy(get_triangles(tri))
        cur_constrained_edges = deepcopy(get_all_constrained_edges(tri))
        cur_boundary = deepcopy(keys(get_boundary_edge_map(tri)))
        add_point!(tri, 2.0, 14.0; store_event_history=Val(true), event_history=storage)
        new_triangles = get_triangles(tri)
        new_constrained_edges = get_all_constrained_edges(tri)
        new_boundary = keys(get_boundary_edge_map(tri))
        added_triangles = setdiff(new_triangles, cur_triangles)
        deleted_triangles = setdiff(cur_triangles, new_triangles)
        added_segments = setdiff(new_constrained_edges, cur_constrained_edges)
        deleted_segments = setdiff(cur_constrained_edges, new_constrained_edges)
        @test isempty(storage.added_boundary_segments)
        @test isempty(storage.deleted_boundary_segments)
        @test DT.compare_triangle_collections(storage.added_triangles, Set{DT.triangle_type(tri)}(added_triangles))
        @test DT.compare_triangle_collections(storage.deleted_triangles, Set{DT.triangle_type(tri)}(deleted_triangles))
        @test compare_edge_vectors(storage.added_segments, added_segments)
        @test compare_edge_vectors(storage.deleted_segments, deleted_segments)

        empty!(storage)

        cur_triangles = deepcopy(get_triangles(tri))
        cur_constrained_edges = deepcopy(get_all_constrained_edges(tri))
        cur_boundary = deepcopy(keys(get_boundary_edge_map(tri)))
        add_point!(tri, 5.0, 0.0; store_event_history=Val(true), event_history=storage)
        new_triangles = get_triangles(tri)
        new_constrained_edges = get_all_constrained_edges(tri)
        new_boundary = keys(get_boundary_edge_map(tri))
        added_triangles = setdiff(new_triangles, cur_triangles)
        deleted_triangles = setdiff(cur_triangles, new_triangles)
        added_segments = setdiff(new_constrained_edges, cur_constrained_edges)
        deleted_segments = setdiff(cur_constrained_edges, new_constrained_edges)
        added_boundary_segments = setdiff(new_boundary, cur_boundary)
        deleted_boundary_segments = setdiff(cur_boundary, new_boundary)
        @test DT.compare_triangle_collections(storage.added_triangles, Set{DT.triangle_type(tri)}(added_triangles))
        @test DT.compare_triangle_collections(storage.deleted_triangles, Set{DT.triangle_type(tri)}(deleted_triangles))
        @test compare_edge_vectors(storage.added_segments, added_segments)
        @test compare_edge_vectors(storage.deleted_segments, deleted_segments)
        @test compare_edge_vectors(storage.added_boundary_segments, added_boundary_segments)
        @test compare_edge_vectors(storage.deleted_boundary_segments, deleted_boundary_segments)
    end
end

@testset "Circumcenter insertion" begin
    curve_1 = [[
        (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
        (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
        (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
        (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
        (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
        (-4.0, 0.0), (0.0, 0.0),
    ]]
    curve_2 = [[
        (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
        (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
        (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
    ]]
    curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
    curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
    curves = [curve_1, curve_2, curve_3, curve_4]
    points = [
        (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
        (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
        (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
        (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
        (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
        (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
        (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
        (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
        (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
        (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
        (7.0, 7.8), (7.0, 7.4)]
    boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
    uncons_tri = triangulate(points)
    cons_tri = triangulate(points; boundary_nodes=boundary_nodes)
    add_ghost_triangles!(uncons_tri)
    add_ghost_triangles!(cons_tri)
    add_edge!(cons_tri, 16, 1)
    add_edge!(cons_tri, 16, 20)
    add_edge!(uncons_tri, 16, 1)
    add_edge!(uncons_tri, 33, 1)
    tri = cons_tri
    targets = DT.RefinementTargets(; min_angle=30.0)
    queue = DT.initialise_refinement_queue(tri, targets)
    storage = DT.initialise_event_history(tri)
    T = (16, 55, 17)
    stats1 = deepcopy(statistics(tri))
    @test !DT.try_circumcenter_insertion!(tri, T, storage, queue)
    stats2 = deepcopy(statistics(tri))
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    @test stats1 == stats2

    T = (63, 92, 62)
    stats1 = deepcopy(statistics(tri))
    @test !DT.try_circumcenter_insertion!(tri, T, storage, queue)
    stats2 = deepcopy(statistics(tri))
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    @test stats1 == stats2

    T = (19, 16, 97)
    stats1 = deepcopy(statistics(tri))
    @test !DT.try_circumcenter_insertion!(tri, T, storage, queue)
    stats2 = deepcopy(statistics(tri))
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    @test stats1 == stats2
end