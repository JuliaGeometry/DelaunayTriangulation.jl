using ..DelaunayTriangulation
const DT = DelaunayTriangulation

include("../helper_functions.jl")

@testset "Computing statistics" begin
    _x, _y = complicated_geometry()
    x = _x
    y = _y
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_1 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_1)
    refine!(tri_1; max_area=1e-3A)
    boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
    tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_2)
    refine!(tri_2; max_area=1e-3A)
    boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
    tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_3)
    refine!(tri_3; max_area=1e-3A)
    boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
    tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri_4)
    refine!(tri_4; max_area=1e-3A)
    a, b = 0.0, 5.0
    c, d = 3.0, 7.0
    nx = 3
    ny = 3
    tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
    tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
    for (i, tri) in enumerate((tri_1, tri_2, tri_3, tri_4, tri_5, tri_6))
        @show i
        validate_statistics(tri)
        @inferred statistics(tri)
        stats = statistics(tri)
        validate_statistics(tri, stats)
        delete_ghost_triangles!(tri)
        validate_statistics(tri)
        @inferred statistics(tri)
        stats = statistics(tri)
        validate_statistics(tri, stats)
    end

    pts = [
        (-7.36, 12.55), (-9.32, 8.59), (-9.0, 3.0), (-6.32, -0.27),
        (-4.78, -1.53), (2.78, -1.41), (-5.42, 1.45), (7.86, 0.67),
        (10.92, 0.23), (9.9, 7.39), (8.14, 4.77), (13.4, 8.61),
        (7.4, 12.27), (2.2, 13.85), (-3.48, 10.21), (-4.56, 7.35),
        (3.44, 8.99), (3.74, 5.87), (-2.0, 8.0), (-2.52, 4.81),
        (1.34, 6.77), (1.24, 4.15)
    ]
    boundary_points = [
        (0.0, 0.0), (2.0, 1.0), (3.98, 2.85), (6.0, 5.0),
        (7.0, 7.0), (7.0, 9.0), (6.0, 11.0), (4.0, 12.0),
        (2.0, 12.0), (1.0, 11.0), (0.0, 9.13), (-1.0, 11.0),
        (-2.0, 12.0), (-4.0, 12.0), (-6.0, 11.0), (-7.0, 9.0),
        (-6.94, 7.13), (-6.0, 5.0), (-4.0, 3.0), (-2.0, 1.0), (0.0, 0.0)
    ]
    boundary_nodes, pts = convert_boundary_points_to_indices(boundary_points; existing_points=pts)
    uncons_tri = triangulate(pts, delete_ghosts=false)
    cons_tri = triangulate(pts; boundary_nodes, delete_ghosts=false)
    add_point!(cons_tri, 0.0, 5.0)
    add_edge!(cons_tri, 40, 26)
    add_edge!(cons_tri, 39, 27)
    add_edge!(cons_tri, 38, 28)
    add_point!(cons_tri, -3.0, 12.0)
    validate_statistics(uncons_tri)
    validate_statistics(cons_tri)
    delete_ghost_triangles!(uncons_tri)
    delete_ghost_triangles!(cons_tri)
    validate_statistics(uncons_tri)
    validate_statistics(cons_tri)
end
