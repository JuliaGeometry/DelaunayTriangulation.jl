@testset "Can we pick the correct ghost triangle that a point resides in?" begin
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
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    [DT.delete_triangle!(i, j, k, T, adj, adj2v, DG) for (i, j, k) in (
        (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1)
    )]
    [DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true) for (i, j, k) in (
        (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1)
    )]
    DT.compute_centroid!(pts)
    p12 = @SVector[3.538447, 3.99844]
    push!(pts, p12)
    @test DT.isleftofline(pts, 9, 10, 12) == 1
    @test DT.isleftofline(pts, 10, DT.BoundaryIndex, 12) == 1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 9, 12) == 1
    @test DT.isintriangle((9, 10, DT.BoundaryIndex), pts, 12) == 1
    @test DT.isintriangle((9, 8, 10), pts, 12) == -1
    @test DT.locate_triangle(T, pts, 12) == (9, 10, DT.BoundaryIndex)
    p13 = @SVector[-1.182399, 4.5127]
    push!(pts, p13)
    @test DT.isleftofline(pts, 1, 9, 13) == 1
    @test DT.isleftofline(pts, 9, 1, 13) == -1
    @test DT.isleftofline(pts, 9, DT.BoundaryIndex, 13) == 1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 1, 13) == 1
    @test DT.isintriangle((1, 9, DT.BoundaryIndex), pts, 13) == 1
    @test DT.isintriangle((9, DT.BoundaryIndex, 1), pts, 13) == 1
    @test DT.isintriangle((DT.BoundaryIndex, 1, 9), pts, 13) == 1
    @test DT.isintriangle((1, 6, 8), pts, 13) == -1
    @test DT.locate_triangle(T, pts, 13) == (1, 9, DT.BoundaryIndex)
    p14 = @SVector[-4.85877, 3.086]
    push!(pts, p14)
    @test DT.isleftofline(pts, 2, 1, 14) == 1
    @test DT.isleftofline(pts, 1, DT.BoundaryIndex, 14) == 1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 2, 14) == 1
    @test DT.isleftofline(pts, 10, 5, 14) == -1
    @test DT.isintriangle((2, 1, DT.BoundaryIndex), pts, 14) == 1
    @test DT.isintriangle((1, DT.BoundaryIndex, 2), pts, 14) == 1
    @test DT.isintriangle((DT.BoundaryIndex, 2, 1), pts, 14) == 1
    @test DT.isintriangle((4, 3, DT.BoundaryIndex), pts, 14) == -1
    @test DT.locate_triangle(T, pts, 14) == (2, 1, DT.BoundaryIndex)
    p15 = @SVector[-2.0, -5.0]
    push!(pts, p15)
    @test DT.isleftofline(pts, 5, 4, 15) == 1
    @test DT.isleftofline(pts, 4, DT.BoundaryIndex, 15) == 1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 4, 15) == -1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 5, 15) == 1
    @test DT.isleftofline(pts, 5, DT.BoundaryIndex, 15) == -1
    @test DT.isintriangle((5, 4, DT.BoundaryIndex), pts, 15) == 1
    @test DT.isintriangle((4, DT.BoundaryIndex, 5), pts, 15) == 1
    @test DT.isintriangle((DT.BoundaryIndex, 5, 4), pts, 15) == 1
    @test DT.isintriangle((6, 7, 8), pts, 15) == -1
    @test DT.isintriangle((10, 5, DT.BoundaryIndex), pts, 15) == -1
    @test DT.isintriangle((1, 9, DT.BoundaryIndex), pts, 15) == -1
    @test DT.isintriangle((9, DT.BoundaryIndex, 1), pts, 15) == -1
    @test DT.isintriangle((DT.BoundaryIndex, 1, 9), pts, 15) == -1
    @test DT.locate_triangle(T, pts, 15) == (5, 4, DT.BoundaryIndex)
    p16 = @SVector[16.27, 0.92]
    push!(pts, p16)
    @test DT.isleftofline(pts, 10, 5, 16) == 1
    @test DT.isleftofline(pts, 5, DT.BoundaryIndex, 16) == 1
    @test DT.isleftofline(pts, DT.BoundaryIndex, 10, 16) == 1
    @test DT.isintriangle((10, 5, DT.BoundaryIndex), pts, 16) == 1
    @test DT.isintriangle((6, 7, 8), pts, 16) == -1
    @test DT.isintriangle((1, 9, DT.BoundaryIndex), pts, 16) == -1
    @test DT.isintriangle((2, 6, 3), pts, 16) == -1
    @test DT.isintriangle((3, 2, DT.BoundaryIndex), pts, 16) == -1
    @test DT.isleftofline(pts, 10, DT.BoundaryIndex, 16) == -1
    @test DT.locate_triangle(T, pts, 16) == (10, 5, DT.BoundaryIndex)
end