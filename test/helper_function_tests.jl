using ..DelaunayTriangulation
using Test
const DT = DelaunayTriangulation

@testset "Check that validate_triangulation breaks on triangulations that it should break on" begin
    for _ in 1:20
        a = [0.0, 0.0]
        b = [3.0, 0.0]
        c = [5.0, 0.0]
        d = [6.0, 0.0]
        e = [2.0, 4.0]
        f = [2.0, 2.0]
        g = [2.0, -2.0]
        h = [6.0, 4.0]
        pts = [a, b, c, d, e, f, g, h]
        tri = triangulate(pts; delete_ghosts=false)
        @test validate_triangulation(tri; print_result=false)
        DT.delete_triangle!(tri, 6, 2, 3)
        @test !validate_triangulation(tri; print_result=false)
        DT.add_triangle!(tri, 6, 3, 2)
        @test !test_state(test_triangle_orientation(tri))
        @test !validate_triangulation(tri; print_result=false)
        tri.points[7] = [2, -0.5]
        @test !validate_triangulation(tri; print_result=false)
        tri.points[7] = [2.0, -2.0]
        tri.adjacent.adjacent[(6, 3)] = DT.∅
        @test !test_state(test_each_edge_has_two_incident_triangles(tri))
        @test !validate_triangulation(tri; print_result=false)
        tri.adjacent.adjacent[(6, 3)] = 2
        tri.adjacent.adjacent[(6, 17)] = 10
        @test !test_state(test_adjacent2vertex_map_matches_adjacent_map(tri))
        @test !test_state(test_adjacent_map_matches_adjacent2vertex_map(tri))
        tri = triangulate(pts; delete_ghosts=false)
        DT.delete_vertex!(tri, 3)
        @test !validate_triangulation(tri; print_result=false)
        tri = triangulate(pts; delete_ghosts=false)
        push!(tri.triangles, (11, 17, 20))
        @test any(!test_state, test_iterators(tri))
        tri = triangulate(pts; delete_ghosts=true)
        push!(tri.triangles, (11, 17, -20))
        @test any(!test_state, test_iterators(tri))
        tri = triangulate(pts; delete_ghosts=false)
        push!(tri.triangles, (11, 17, 20000))
        @test any(!test_state, test_iterators(tri))
    end
end

@testset "Validating constrained triangulations" begin
    _x, _y = complicated_geometry()
    x = _x
    y = _y
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_1 = triangulate(points; boundary_nodes)
    A = get_area(tri_1)
    refine!(tri_1; max_area=1e-3A, use_circumcenter=true)
    @test validate_triangulation(tri_1)

    boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
    tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_area(tri_2)
    refine!(tri_2; max_area=1e-3A, use_circumcenter=true)
    @test validate_triangulation(tri_2)

    x = [0.0, 2.0, 2.0, 0.0, 0.0]
    y = [0.0, 0.0, 2.0, 2.0, 0.0]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_area(tri_3)
    refine!(tri_3; max_area=1e-3A, use_circumcenter=true)
    @test validate_triangulation(tri_3)

    x = reverse(reverse.(_x[2]))
    y = reverse(reverse.(_y[2]))
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_area(tri_4)
    refine!(tri_4; max_area=1e-3A, use_circumcenter=true)
    @test validate_triangulation(tri_4)
end

@testset "is_conformal" begin
    tri = triangulate(rand(2, 50))
    @test is_conformal(tri)

    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.48, 0.5)]
    tri = triangulate(points)
    @test is_conformal(tri)
    add_segment!(tri, 1, 3)
    @test !is_conformal(tri)
end

@testset "test_visibility" begin
    points = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (1.0, 1.0), (3.5, 5.5), (4.5, 4.0)]
    boundary_nodes = [1, 2, 3, 4, 1]
    tri = triangulate(points; boundary_nodes, segments=Set([(2, 4)]))
    @test DT.is_visible(DT.test_visibility(tri, 4, 2, 3))
    @test DT.is_invisible(DT.test_visibility(tri, 5, 1, 7))
end


