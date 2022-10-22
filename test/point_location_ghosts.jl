@testset "Can we correctly identify intersections and triangles for boundary points and ghost triangles?" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    for (i, j, k) in (
        (1, 2, 6),
        (1, 6, 8),
        (9, 1, 8),
        (9, 8, 10),
        (10, 8, 11),
        (8, 7, 11),
        (8, 6, 7),
        (6, 2, 3),
        (6, 3, 4),
        (6, 4, 7),
        (7, 4, 5),
        (11, 7, 5),
        (10, 11, 5)
    )
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end
    p12 = @SVector[3.382, 4.3599]
    push!(pts, p12)

    # Test our interior checks
    DT.compute_centroid!(pts[1:11])
    k = 9
    pts[12] = @SVector[3.0, 1.0]
    q = DT.get_point(pts, 12)
    p = DT.get_point(pts, k)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (8, 10, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (8, 11, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (10, 11, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (5, 7, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (8, 9, true, false)
    pts[12] = @SVector[4.14, -1.95]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (8, 10, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (11, 5, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (11, 7, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (5, 7, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (6, 8, true, false)
    pts[12] = @SVector[5.83, -5.2]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (8, 10, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (11, 5, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (3, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (6, 8, true, false)
    pts[12] = @SVector[-2.84, 2.25]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (1, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (9, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (11, 7, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (7, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (6, 2, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (6, 8, false, true)
    pts[12] = @SVector[-4.19, 2.89]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (1, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (9, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (11, 7, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (6, 3, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (6, 2, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (2, 6, false, true)
    pts[12] = @SVector[-4.375, -2.44]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (1, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (8, 11, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (7, 4, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (6, 3, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (3, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (2, 6, true, false)
    pts[12] = @SVector[1.79, 2.34]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (8, 10, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (9, 8, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (10, 11, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (7, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (8, 9, true, false)
    pts[12] = @SVector[0.34, 3.25]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (1, 8, false, true)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (9, 8, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (10, 11, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (7, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (8, 9, false, true)
    pts[12] = @SVector[4.6, 5.49]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (7, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (4, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (0, 0, false, false)
    pts[12] = @SVector[0.0, 5.0]
    q = DT.get_point(pts, 12)
    @test DT.check_interior_edge_intersections(q, adj, DG, 9, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 10, pts) == (0, 0, false, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 5, pts) == (10, 11, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 4, pts) == (7, 6, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 3, pts) == (6, 2, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 2, pts) == (6, 1, true, false)
    @test DT.check_interior_edge_intersections(q, adj, DG, 1, pts) == (0, 0, false, false)

    # Test the straight line search
    for k in [9, 10, 5, 4, 3, 2, 1]
        pts[12] = @SVector[3.382, 4.3599]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (9, 10)
        pts[12] = @SVector[15.27, 9.77]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (9, 10)
        pts[12] = @SVector[18.0, 0.0]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (10, 5)
        pts[12] = @SVector[20.257, 7.27]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (10, 5)
        pts[12] = @SVector[12.0, -8.0]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (5, 4)
        pts[12] = @SVector[0.466, -5.18]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (5, 4)
        pts[12] = @SVector[-2.746, -4.614]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (4, 3)
        pts[12] = @SVector[-6.94, -3.37]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (4, 3)
        pts[12] = @SVector[-4.999, -2.749]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (4, 3)
        pts[12] = @SVector[-7.046, -1.45]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (3, 2)
        pts[12] = @SVector[-73.36, 1.30258]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (3, 2)
        pts[12] = @SVector[-8.0, 2.0]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (3, 2)
        pts[12] = @SVector[-7.17, 2.96]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (2, 1)
        pts[12] = @SVector[-6.0, 4.0]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (2, 1)
        pts[12] = @SVector[-9.39, 7.46]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (2, 1)
        pts[12] = @SVector[-3.0, 4.17]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (1, 9)
        pts[12] = @SVector[0.0, 6.0]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (1, 9)
        pts[12] = @SVector[2.47, 7.027]
        q = DT.get_point(pts, 12)
        @test DT.straight_line_search_ghost_triangles(q, adj, k, pts) == (9, 10)
    end
end