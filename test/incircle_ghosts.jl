@testset "Can we detect if a point is in a ghost triangle's circumdisk?" begin
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
    p12 = @SVector[-1.86, 5.99]
    push!(pts, p12)
    @test DT.isincircle(pts, 1, 9, 0, 12) == 1
    @test DT.isincircle(pts, 9, 0, 1, 12) == 1
    @test DT.isincircle(pts, 0, 1, 9, 12) == 1
    @test DT.isincircle(pts, 1, 0, 9, 12) == -1
    p13 = @SVector[-1.86, 102.9]
    push!(pts, p13)
    @test DT.isincircle(pts, 1, 9, 0, 13) == 1
    @test DT.isincircle(pts, 9, 0, 1, 12) == 1
    @test DT.isincircle(pts, 0, 1, 9, 12) == 1
    @test DT.isincircle(pts, 1, 0, 9, 13) == -1
    p14 = @SVector[3.54, 4.684]
    push!(pts, p14)
    @test DT.isincircle(pts, 9, 10, 0, 14) == 1
    @test DT.isincircle(pts, 10, 0, 9, 14) == 1
    @test DT.isincircle(pts, 0, 9, 10, 14) == 1
    @test DT.isincircle(pts, 9, 10, 0, 9) == 0
    @test DT.isincircle(pts, 10, 0, 9, 9) == 0
    @test DT.isincircle(pts, 0, 9, 10, 9) == 0
    @test DT.isincircle(pts, 10, 9, 8, 14) == -1
    p15 = @SVector[1.57, 2.514]
    push!(pts, p15)
    @test DT.isincircle(pts, 9, 10, 0, 15) == -1
    @test DT.isincircle(pts, 10, 0, 9, 15) == -1
    @test DT.isincircle(pts, 0, 9, 10, 15) == -1
    p16 = @SVector[6.77, 0.269]
    push!(pts, p16)
    @test DT.isincircle(pts, 10, 5, 0, 16) == 1
    @test DT.isincircle(pts, 5, 0, 10, 16) == 1
    @test DT.isincircle(pts, 0, 10, 5, 16) == 1
    p17 = @SVector[4.21754, 3.00067]
    push!(pts, p17)
    @test DT.isincircle(pts, 10, 5, 0, 17) == -1
    @test DT.isincircle(pts, 5, 0, 10, 17) == -1
    @test DT.isincircle(pts, 0, 10, 5, 17) == -1
    p18 = @SVector[4.816, -3.696112]
    push!(pts, p18)
    @test DT.isincircle(pts, 5, 4, 0, 18) == 1
    @test DT.isincircle(pts, 4, 0, 5, 18) == 1
    @test DT.isincircle(pts, 0, 5, 4, 18) == 1
    p19 = @SVector[6.499685, -2.935]
    push!(pts, p19)
    @test DT.isincircle(pts, 5, 4, 0, 19) == -1
    @test DT.isincircle(pts, 4, 0, 5, 19) == -1
    @test DT.isincircle(pts, 0, 5, 4, 19) == -1
    p20 = @SVector[2.79587, -4.020351]
    push!(pts, p20)
    @test DT.isincircle(pts, 5, 4, 0, 20) == 1
    @test DT.isincircle(pts, 4, 0, 5, 20) == 1
    @test DT.isincircle(pts, 0, 5, 4, 20) == 1
    p21 = @SVector[0.0, -4.0]
    push!(pts, p21)
    @test DT.isincircle(pts, 5, 4, 0, 21) == -1
    @test DT.isincircle(pts, 4, 0, 5, 21) == -1
    @test DT.isincircle(pts, 0, 5, 4, 21) == -1
    p22 = @SVector[-2.815953, -4.25729]
    push!(pts, p22)
    @test DT.isincircle(pts, 4, 3, 0, 22) == 1
    @test DT.isincircle(pts, 3, 0, 4, 22) == 1
    @test DT.isincircle(pts, 0, 4, 3, 22) == 1
    p23 = @SVector[-1.4317, -4.3196]
    push!(pts, p23)
    @test DT.isincircle(pts, 4, 3, 0, 23) == -1
    @test DT.isincircle(pts, 3, 0, 4, 23) == -1
    @test DT.isincircle(pts, 0, 4, 3, 23) == -1
    p24 = @SVector[-5.04821, -2.54880]
    push!(pts, p24)
    @test DT.isincircle(pts, 4, 3, 0, 24) == 1
    @test DT.isincircle(pts, 3, 0, 4, 24) == 1
    @test DT.isincircle(pts, 0, 4, 3, 24) == 1
    @test DT.isincircle(pts, 4, 3, 0, 4) == 0
    @test DT.isincircle(pts, 3, 0, 4, 4) == 0
    @test DT.isincircle(pts, 0, 4, 3, 4) == 0
    @test DT.isincircle(pts, 4, 3, 0, 3) == 0
    @test DT.isincircle(pts, 3, 0, 4, 3) == 0
    @test DT.isincircle(pts, 0, 4, 3, 3) == 0
    p25 = @SVector[-6.3327007, -1.7257]
    push!(pts, p25)
    @test DT.isincircle(pts, 4, 3, 0, 25) == 1
    @test DT.isincircle(pts, 3, 0, 4, 25) == 1
    @test DT.isincircle(pts, 0, 4, 3, 25) == 1
    p26 = @SVector[-6.444937, -0.54101]
    push!(pts, p26)
    @test DT.isincircle(pts, 3, 2, 0, 26) == 1
    @test DT.isincircle(pts, 2, 0, 3, 26) == 1
    @test DT.isincircle(pts, 0, 3, 2, 26) == 1
    p27 = @SVector[-5.310, 2.87596]
    push!(pts, p27)
    @test DT.isincircle(pts, 2, 1, 0, 27) == 1
    @test DT.isincircle(pts, 1, 0, 2, 27) == 1
    @test DT.isincircle(pts, 0, 2, 1, 27) == 1
    p28 = @SVector[-5.247746, 0.905588]
    push!(pts, p28)
    @test DT.isincircle(pts, 2, 1, 0, 28) == -1
    @test DT.isincircle(pts, 1, 0, 2, 28) == -1
    @test DT.isincircle(pts, 0, 2, 1, 28) == -1
    @test DT.isincircle(pts, 8, 7, 11, 28) == -1
    @test DT.isincircle(pts, 8, 11, 10, 28) == -1
    @test DT.isincircle(pts, 6, 4, 7, 28) == -1
end

@testset "Circumdisk of a ghost triangle on the edge" begin
    p8 = (6.0, 6.0)
    p9 = (10.0, -2.0)
    p13 = (8.0, 2.0)
    pts = [p8, p9, p13]
    @test DT.isincircle(pts, 1, 2, 0, 3) == 1
    push!(pts, (2.0, 14.0))
    @test DT.isincircle(pts, 1, 2, 0, 4) == -1
    push!(pts, (12.0, -6.0))
    @test DT.isincircle(pts, 1, 2, 0, 5) == -1
    push!(pts, (34.0, -6.0))
    @test DT.isincircle(pts, 1, 2, 0, 6) == 1
end