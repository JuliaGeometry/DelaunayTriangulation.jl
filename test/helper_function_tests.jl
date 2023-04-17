using ..DelaunayTriangulation
using Test
const DT = DelaunayTriangulation
include("./helper_functions.jl")

@testset "Check that validate_triangulation breaks on triangulations that it should break on" begin
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
    @test validate_triangulation(tri)
    @test test_planarity(tri)
    DT.delete_triangle!(tri, 6, 2, 3)
    @test !test_planarity(tri)
    @test !validate_triangulation(tri)
    DT.add_triangle!(tri, 6, 3, 2)
    @test test_planarity(tri)
    @test !test_triangle_orientation(tri)
    @test !validate_triangulation(tri)
    DT.delete_triangle!(tri, 6, 3, 2)
    DT.add_triangle!(tri, 6, 2, 3)
    tri.points[7] = [2, -0.5]
    @test !test_delaunay_criterion(tri)
    @test !validate_triangulation(tri)
    tri.points[7] = [2.0, -2.0]
    tri.adjacent.adjacent[(6, 3)] = DT.DefaultAdjacentValue
    @test !test_each_edge_has_two_incident_triangles(tri)
    @test !validate_triangulation(tri)
    tri.adjacent.adjacent[(6, 3)] = 2
    tri.adjacent.adjacent[(6, 17)] = 10
    @test !test_adjacent2vertex_map_matches_adjacent_map(tri)
    @test !test_adjacent_map_matches_adjacent2vertex_map(tri)
    tri = triangulate(pts; delete_ghosts=false)
    DT.delete_vertex!(tri, 3)
    @test !validate_triangulation(tri)
    tri = triangulate(pts; delete_ghosts=false)
    push!(tri.triangles, (11, 17, 20))
    @test !test_iterators(tri)
    @test !validate_triangulation(tri)
    tri = triangulate(pts; delete_ghosts=true)
    push!(tri.triangles, (11, 17, -20))
    @test !test_iterators(tri)
    @test !validate_triangulation(tri)
    tri = triangulate(pts; delete_ghosts=false)
    push!(tri.triangles, (11, 17, 20000))
    @test !test_iterators(tri)
    @test !validate_triangulation(tri)
end

@testset "Validating constrained triangulations" begin
    _x, _y = complicated_geometry()
    x = _x
    y = _y
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_1 = triangulate(points; boundary_nodes)
    A = get_total_area(tri_1)
    refine!(tri_1; max_area=1e-3A)
    @test validate_triangulation(tri_1)

    boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
    tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_2)
    refine!(tri_2; max_area=1e-3A)
    @test validate_triangulation(tri_2)

    x = [0.0, 2.0, 2.0, 0.0, 0.0]
    y = [0.0, 0.0, 2.0, 2.0, 0.0]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_3)
    refine!(tri_3; max_area=1e-3A)
    @test validate_triangulation(tri_3)

    x = reverse(reverse.(_x[2]))
    y = reverse(reverse.(_y[2]))
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_4)
    refine!(tri_4; max_area=1e-3A)
    @test validate_triangulation(tri_4)

    a, b = 0.0, 5.0
    c, d = 3.0, 7.0
    nx = 3
    ny = 3
    tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
    @test validate_triangulation(tri_5)
    tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
    @test validate_triangulation(tri_6)
    push!(tri_6.triangles, (11, 1191919, 2000))
    @test !validate_triangulation(tri_6)
    push!(tri_5.triangles, (-1, -2, -3))
    @test !validate_triangulation(tri_5)
end
