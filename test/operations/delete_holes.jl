using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs

include("../helper_functions.jl")

@testset "Simple domain" begin
    p1 = (0.2, 0.2)
    p2 = (0.8, 0.2)
    p3 = (0.8, 0.8)
    p4 = (0.2, 0.8)
    rng = StableRNG(1919)
    existing_pts = [(rand(rng), rand(rng)) for _ in 1:15]
    push!(existing_pts, (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0))
    nodes, points = convert_boundary_points_to_indices([p1, p2, p3, p4, p1], existing_points=existing_pts)
    tri = triangulate(points; boundary_nodes=nodes, delete_holes=false, delete_ghosts=false, rng)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((17, 1, 19, 12, 3, 8, 16, 15, 5, 11, 18, 9, 4))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (1, 22, 4),
        (1, 4, 19),
        (1, 19, 17),
        (1, 17, 12),
        (1, 12, 22),
        (3, 8, 23),
        (3, 23, 17),
        (3, 17, 8),
        (4, 11, 19),
        (4, 22, 21),
        (4, 21, 11),
        (5, 18, 11),
        (5, 11, 21),
        (5, 21, 9),
        (5, 9, 15),
        (5, 15, 16),
        (8, 16, 20),
        (8, 17, 16),
        (8, 20, 23),
        (9, 20, 15),
        (9, 21, 20),
        (5, 16, 18),
        (11, 18, 19),
        (12, 23, 22),
        (12, 17, 23),
        (15, 20, 16),
        (17, 19, DT.BoundaryIndex),
        (19, 18, DT.BoundaryIndex),
        (18, 16, DT.BoundaryIndex),
        (16, 17, DT.BoundaryIndex)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    p1 = (0.2, 0.2)
    p2 = (0.8, 0.2)
    p3 = (0.8, 0.8)
    p4 = (0.2, 0.8)
    rng = StableRNG(1919)
    existing_pts = [(rand(rng), rand(rng)) for _ in 1:15]
    push!(existing_pts, (0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0))
    nodes, points = convert_boundary_points_to_indices([p1, p2, p3, p4, p1], existing_points=existing_pts)
    tri = triangulate(points; boundary_nodes=nodes, rng)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end

@testset "Multiply-connected" begin
    a = (0.0, 0.0)
    b = (5.0, 0.0)
    c = (5.0, 5.0)
    d = (0.0, 5.0)
    e = (1.0, 4.0)
    f = (0.0, 3.0)
    g = (1.5, 2.0)
    h = (4.0, 2.0)
    i = (4.0, 1.0)
    j = (1.5, 1.0)
    k = (2.0, 1.5)
    ℓ = (2.4, 1.0)
    m = (2.8, 3.6)
    n = (2.2, 2.8)
    o = (7.2, 2.4)
    p = (7.3, 3.9)
    q = (2.0, 1.4)
    r = (3.0, 1.6)
    s = (3.4, 1.4)
    t = (3.6, 2.2)
    u = (1.0, 2.8)
    v = (0.8, 1.8)
    w = (0.8, 4.0)
    z = (0.4, 4.4)
    c1 = (5.0, 4.6)
    d1 = (4.9, 4.8)
    a1 = (5.6, 1.8)
    b1 = (6.0, 0.6)
    e1 = (5.1, 4.4)
    f1 = (5.0, 4.3)
    pts = [[[d, e, f, a, b, f1, e1, c1, d1, c, d]], [[g, h, i, ℓ, k, j, g]]]
    existing_points = [m, n, o, p, q, r, s, t, u, v, w, z, a1, b1]
    nodes, points = convert_boundary_points_to_indices(pts; existing_points=existing_points)
    tri = triangulate(points; boundary_nodes=nodes, delete_holes=false, delete_ghosts=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((12, 11, 14, 3, 13, 4, 6, 7, DT.BoundaryIndex))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (15, 17, 12),
        (12, 17, 11),
        (12, 11, 16),
        (11, 17, 16),
        (25, 30, 29),
        (29, 28, 6),
        (25, 29, 6),
        (25, 6, 26),
        (6, 7, 26),
        (6, 28, 7),
        (7, 28, 27),
        (26, 7, 27),
        (20, 19, 13),
        (13, 19, 14),
        (13, 14, 3),
        (20, 13, 3),
        (20, 3, 4),
        (21, 20, 4),
        (24, 23, 22),
        (24, 22, 21),
        (24, 21, 4),
        (15, 24, DT.BoundaryIndex),
        (24, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex),
        (3, 14, DT.BoundaryIndex),
        (14, 19, DT.BoundaryIndex),
        (19, 18, DT.BoundaryIndex),
        (18, 17, DT.BoundaryIndex),
        (17, 15, DT.BoundaryIndex),
        (12, 16, 15)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    tri = triangulate(points; boundary_nodes=nodes)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)
end

@testset "Interior holes that were already triangles" begin
    p1 = (0.0, 0.0)
    p2 = (1.0, 0.0)
    p3 = (1.0, 1.0)
    p4 = (0.0, 1.0)
    p5 = (0.2, 0.2)
    p6 = (0.8, 0.2)
    p7 = (0.8, 0.8)
    pts = [p1, p2, p3, p4, p5, p6, p7]
    tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]], delete_holes=false, delete_ghosts=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((DT.BoundaryIndex,))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    @test DT.compare_triangle_collections(triangles_to_delete, Set((
        (1, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex),
        (3, 2, DT.BoundaryIndex),
        (2, 1, DT.BoundaryIndex),
        (5, 6, 7)
    )))
    tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]])
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)

    tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7], [7, 6, 5]]], delete_holes=false, delete_ghosts=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((DT.BoundaryIndex,))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    @test DT.compare_triangle_collections(triangles_to_delete, Set((
        (1, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex),
        (3, 2, DT.BoundaryIndex),
        (2, 1, DT.BoundaryIndex),
        (5, 6, 7)
    )))
    tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7], [7, 6, 5]]])
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)

    tri = triangulate(pts, boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]], delete_holes=false, delete_ghosts=false)
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((DT.BoundaryIndex,))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    @test DT.compare_triangle_collections(triangles_to_delete, Set((
        (1, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex),
        (3, 2, DT.BoundaryIndex),
        (2, 1, DT.BoundaryIndex),
        (5, 6, 7)
    )))
    tri = triangulate(pts; boundary_nodes=[[[1, 2, 3, 4, 1]], [[5, 7, 6, 5]]])
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)
end

@testset "A previously broken example" begin
    a = 4 / 5
    t = LinRange(0, 2π, 6)
    x = @. a * (2cos(t) + cos(2t))
    y = @. a * (2sin(t) - sin(2t))
    tri = generate_mesh(x, y, 5.0)
    points = get_points(tri)
    bn_nodes = get_boundary_nodes(tri)
    _tri = triangulate(points; boundary_nodes=bn_nodes, delete_ghosts=false, delete_holes=false)
    tri = _tri
    points_to_process = DT.find_all_points_to_delete(tri)
    @test points_to_process == Set((DT.BoundaryIndex,))
    triangles_to_delete = DT.find_all_triangles_to_delete(tri, points_to_process)
    true_triangles_to_delete = Set((
        (3, 2, 1),
        (4, 1, 5),
        (3, 1, DT.BoundaryIndex),
        (1, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex)
    ))
    @test DT.compare_triangle_collections(triangles_to_delete, true_triangles_to_delete)

    tri = triangulate(points; boundary_nodes=bn_nodes)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)
end